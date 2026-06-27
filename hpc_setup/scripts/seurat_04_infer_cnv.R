#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 04: Inferred copy-number variation
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-23
# Last update:  2026-05-23
# Input:        --rds_in        Seurat .rds from seurat_03_annotate_reference.R
#               --gene_order    inferCNV gene_order_file.txt (hg38)
#               --out_dir       output directory for inferCNV run + summaries
#               --ref_groups    comma-separated cell_type labels to use as the
#                               diploid reference (default: "B cells,T cells")
#               --query_groups  comma-separated cell_type labels to include as
#                               malignant/ambiguous query (default: everything
#                               except ref_groups and HS5 stromal / MSC)
# Output:       <out_dir>/infercnv/                  (full inferCNV run dir)
#                 ├── infercnv.png / .pdf            (chromosome heatmap)
#                 ├── HMM_CNV_predictions.HMMi6...   (HMM calls)
#                 └── 21_denoise.txt etc.            (intermediate)
#               <out_dir>/tables/per_cell_cnv_score.tsv
#               <out_dir>/tables/cluster_cnv_summary.tsv
#               <out_dir>/tables/cd16_classification.tsv
#               <out_dir>/figures/04a_cnv_score_violin.pdf
#               <out_dir>/figures/04b_cnv_score_umap.pdf
#               <out_dir>/logs/seurat_04_sessionInfo.txt
#               <rds_in>                              ← re-saved with cnv_score
#                                                     + cd16_call metadata cols
# Depends on:   Seurat (>=5.0), infercnv (>=1.18), ggplot2, dplyr, optparse
# Notes:
#   * HS-5 stromal / MSC is always excluded — it's a stromal cell line and
#     would dominate the CNV signal as "weird non-haematopoietic".
#   * Reference groups must contain >=50 cells each to give stable baselines;
#     if you have very few B or T cells per sample, pool across samples by
#     using --ref_groups "B cells,T cells" (the default).
#   * The CD16+ Mono/Macrophage classification is heuristic — it compares each
#     cell's CNV burden to the reference (normal) and to known-AML clusters
#     (GMP-like, HSPC/LSC-like). Treat as one line of evidence, not proof.
#   * Heavy run — submit through submit_seurat_04_infer_cnv.sh, not interactive.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 50 * 1024^3)   # 50 GB

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),    type = "character",
                     help = "Seurat .rds from seurat_03 (with cell_type)")
parser <- add_option(parser, c("--gene_order"), type = "character",
                     help = "inferCNV gene_order_file (gene<TAB>chr<TAB>start<TAB>end)")
parser <- add_option(parser, c("--out_dir"),   type = "character",
                     help = "Output directory")
parser <- add_option(parser, c("--ref_groups"),  type = "character",
                     default = "B cells,T cells",
                     help = "Comma-separated normal reference labels [default: %default]")
parser <- add_option(parser, c("--query_groups"), type = "character", default = "",
                     help = "Comma-separated query labels (empty = all malignant clusters)")
parser <- add_option(parser, c("--cutoff"),    type = "double", default = 0.1,
                     help = "Min mean count per gene; 0.1 for 10x [default: %default]")
parser <- add_option(parser, c("--ncores"),    type = "integer", default = 16,
                     help = "Cores for inferCNV [default: %default]")
parser <- add_option(parser, c("--max_cells"), type = "integer", default = 0,
                     help = "Downsample query to this many cells per cluster (0 = off)")
opt <- parse_args(parser)

for (req in c("rds_in", "gene_order", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in))     stop("Input RDS not found: ", opt$rds_in)
if (!file.exists(opt$gene_order)) stop("gene_order file not found: ", opt$gene_order)

run_dir <- file.path(opt$out_dir, "infercnv")
tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(run_dir, tab_dir, fig_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Seurat 04: inferCNV ===")
message("  rds_in     : ", opt$rds_in)
message("  gene_order : ", opt$gene_order)
message("  out_dir    : ", opt$out_dir)
message("  ref_groups : ", opt$ref_groups)
message("  cutoff     : ", opt$cutoff)
message("  ncores     : ", opt$ncores)

# ---- 1. Load and prepare inputs ------------------------------------------

message("--- Loading Seurat object ---")
seu_merged <- readRDS(opt$rds_in)
DefaultAssay(seu_merged) <- "RNA"
if (!"cell_type" %in% colnames(seu_merged@meta.data)) {
  stop("cell_type column required — run seurat_01b / seurat_03 first")
}

ref_groups <- trimws(strsplit(opt$ref_groups, ",")[[1]])
missing_ref <- setdiff(ref_groups, unique(as.character(seu_merged$cell_type)))
if (length(missing_ref) > 0) {
  stop("Reference groups not found in cell_type: ",
       paste(missing_ref, collapse = ", "))
}

# Drop HS-5 stromal / MSC (cell line confounder)
exclude_always <- "HS5 stromal / MSC"
seu_use <- subset(seu_merged, subset = cell_type != exclude_always)

if (nchar(opt$query_groups) > 0) {
  query_groups <- trimws(strsplit(opt$query_groups, ",")[[1]])
  keep_labels  <- union(query_groups, ref_groups)
  seu_use <- subset(seu_use, subset = cell_type %in% keep_labels)
}

# Optional per-cluster downsampling — inferCNV is heavy
if (opt$max_cells > 0) {
  message("--- Downsampling to ", opt$max_cells, " cells per cluster ---")
  Idents(seu_use) <- "cell_type"
  seu_use <- subset(seu_use, downsample = opt$max_cells)
}

message("  cells in inferCNV run: ", ncol(seu_use))
print(table(seu_use$cell_type))

# Annotations file — cell<TAB>label
anno_df <- data.frame(cell = colnames(seu_use),
                      label = as.character(seu_use$cell_type),
                      stringsAsFactors = FALSE)
anno_path <- file.path(run_dir, "annotations.tsv")
write.table(anno_df, anno_path, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Raw counts matrix
counts_mat <- GetAssayData(seu_use, assay = "RNA", layer = "counts")

# ---- 2. Run inferCNV -----------------------------------------------------

message("--- Creating inferCNV object ---")
icnv_obj <- CreateInfercnvObject(
  raw_counts_matrix     = counts_mat,
  annotations_file      = anno_path,
  gene_order_file       = opt$gene_order,
  ref_group_names       = ref_groups,
  delim                 = "\t",
  chr_exclude           = c("chrX", "chrY", "chrM")
)

message("--- Running inferCNV (this is the long step) ---")
icnv_obj <- infercnv::run(
  icnv_obj,
  cutoff                = opt$cutoff,
  out_dir               = run_dir,
  cluster_by_groups     = TRUE,
  denoise               = TRUE,
  HMM                   = TRUE,
  HMM_type              = "i6",
  analysis_mode         = "subclusters",
  num_threads           = opt$ncores,
  output_format         = "pdf",
  save_rds              = TRUE,
  no_prelim_plot        = TRUE
)

# ---- 3. Per-cell CNV score + per-cluster summary -------------------------

message("--- Computing per-cell CNV burden ---")
# The denoised expression matrix sits in icnv_obj@expr.data — deviation from 1
# (the diploid reference baseline). CNV burden = mean squared deviation per cell.
expr_mat <- icnv_obj@expr.data
cnv_score <- colMeans((expr_mat - 1)^2, na.rm = TRUE)

# Map back onto the full Seurat object (NA for excluded cells)
seu_merged$cnv_score <- NA_real_
seu_merged$cnv_score[match(names(cnv_score), colnames(seu_merged))] <- cnv_score

cnv_df <- data.frame(cell = colnames(seu_merged),
                     cell_type = as.character(seu_merged$cell_type),
                     cnv_score = seu_merged$cnv_score)
write.table(cnv_df, file.path(tab_dir, "per_cell_cnv_score.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_summary <- cnv_df |>
  filter(!is.na(cnv_score)) |>
  group_by(cell_type) |>
  summarise(n_cells = n(),
            mean_cnv = mean(cnv_score),
            median_cnv = median(cnv_score),
            sd_cnv = sd(cnv_score),
            .groups = "drop") |>
  arrange(desc(mean_cnv))
write.table(cluster_summary, file.path(tab_dir, "cluster_cnv_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("--- Cluster CNV summary ---")
print(cluster_summary)

# ---- 4. CD16+ Mono/Macrophage classification -----------------------------

message("--- Classifying CD16+ Mono/Macrophage ---")
ref_mean <- cluster_summary |>
  filter(cell_type %in% ref_groups) |>
  pull(mean_cnv) |>
  mean()
aml_known <- c("GMP-like AML blasts", "HSPC / LSC-like blasts", "Mono-like AML",
               "ABCB5+ resistant primitive AML")
aml_mean <- cluster_summary |>
  filter(cell_type %in% aml_known) |>
  pull(mean_cnv) |>
  mean(na.rm = TRUE)
cd16_mean <- cluster_summary |>
  filter(cell_type == "CD16+ Mono / macrophage") |>
  pull(mean_cnv)

cd16_call <- if (length(cd16_mean) == 0 || is.na(cd16_mean)) {
  "not_assessed"
} else {
  midpoint <- (ref_mean + aml_mean) / 2
  if (cd16_mean >= midpoint) "AML-derived" else "Normal-like"
}

cd16_df <- data.frame(
  cluster      = "CD16+ Mono / macrophage",
  ref_mean_cnv = ref_mean,
  aml_mean_cnv = aml_mean,
  cd16_mean_cnv = ifelse(length(cd16_mean) == 0, NA_real_, cd16_mean),
  call         = cd16_call
)
write.table(cd16_df, file.path(tab_dir, "cd16_classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("  CD16+ call: ", cd16_call)

# Annotate Seurat object with the call
seu_merged$cd16_call <- NA_character_
seu_merged$cd16_call[seu_merged$cell_type == "CD16+ Mono / macrophage"] <- cd16_call

# ---- 5. Figures ----------------------------------------------------------

message("--- Plotting ---")
plot_df <- cnv_df |> filter(!is.na(cnv_score))
p_violin <- ggplot(plot_df, aes(x = reorder(cell_type, cnv_score, FUN = median),
                                y = cnv_score, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "CNV burden (mean squared deviation)",
       title = "Per-cell CNV burden by cluster")
pdf(file.path(fig_dir, "04a_cnv_score_violin.pdf"), width = 10, height = 6)
print(p_violin)
dev.off()

if ("umap" %in% Reductions(seu_merged)) {
  pdf(file.path(fig_dir, "04b_cnv_score_umap.pdf"), width = 7, height = 6)
  print(FeaturePlot(seu_merged, features = "cnv_score",
                    order = TRUE, pt.size = 0.3) +
        scale_colour_gradient(low = "grey90", high = "#0B3D5C"))
  dev.off()
}

# ---- 6. Save updated RDS + session info ---------------------------------

message("--- Saving updated RDS ---")
saveRDS(seu_merged, file = opt$rds_in)
message("  updated: ", opt$rds_in)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_04_sessionInfo.txt"))

message("Done!")
