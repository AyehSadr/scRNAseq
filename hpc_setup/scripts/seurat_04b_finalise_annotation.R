#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 04b: Finalise annotation + re-plot
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-24
# Last update:  2026-05-24
# Input:        --rds_in    Seurat .rds after seurat_03 (signatures + vgalen_class)
#                           and seurat_04 (cnv_score). Contains the existing
#                           cell_type column from seurat_01b.
#               --out_dir   Output for the re-rendered figures and tables
# Output:       <rds_in>                                  ← re-saved with:
#                 cell_type_final  (relabelled cluster, ordered factor)
#                 is_aml           (TRUE / FALSE / NA)
#                 is_aml_category  (AML / Normal / Uncertain / CellLine)
#               <out_dir>/figures/04b_final_umap.pdf
#               <out_dir>/figures/04b_final_signature_heatmap.pdf
#               <out_dir>/figures/04b_old_vs_final_confusion.pdf
#               <out_dir>/tables/cell_type_final_mapping.tsv
#               <out_dir>/tables/cell_type_final_counts.tsv
#               <out_dir>/logs/seurat_04b_sessionInfo.txt
# Depends on:   Seurat (>=5.0), pheatmap, ggplot2, dplyr, tibble, optparse
# Notes:
#   * Mapping is fixed in CELL_TYPE_MAPPING below — the result of the Part-3
#     Steps 1+2 validation (transcriptomic + CNV). Edit the mapping table
#     only if a new line of evidence overturns one of the calls.
#   * Re-renders the §1 signature heatmap (matching 03d layout) with the
#     finalised labels on the y-axis, so downstream slides reflect the final
#     annotation.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 20 * 1024^3)

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),  type = "character",
                     help = "Seurat .rds after seurat_03 + seurat_04")
parser <- add_option(parser, c("--out_dir"), type = "character",
                     help = "Output directory")
opt <- parse_args(parser)

for (req in c("rds_in", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(tab_dir, fig_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== Seurat 04b: finalise annotation + re-plot ===")
message("  rds_in  : ", opt$rds_in)
message("  out_dir : ", opt$out_dir)

# ---- 1. Mapping (fixed by Part-3 Steps 1+2 validation) -------------------
# Order matters: factor levels follow this order in plots.
# `category` drives the is_aml_category column.

CELL_TYPE_MAPPING <- tibble::tribble(
  ~cell_type,                          ~cell_type_final,                 ~category,
  "HSPC / LSC-like blasts",            "HSPC / LSC-like AML",            "AML",
  "ABCB5+ resistant primitive AML",    "ABCB5+ resistant LSC AML",       "AML",
  "GMP-like AML blasts",               "GMP-like AML",                   "AML",
  "Promono-like / quiescent myeloid",  "Promono quiescent AML",          "AML",
  "Mono-like AML",                     "Mono-like AML",                  "AML",
  "Stress-response myeloid",           "Stress-response AML",            "AML",
  "Unclassified niche-related",        "Niche-stressed primitive AML",   "AML",
  "Eosinophils / basophils",           "Eosinophils / basophils (uncertain)", "Uncertain",
  "Normal CD16+ monocytes",            "Normal CD16+ monocytes",         "Normal",
  "CD16+ Mono / macrophage",           "Normal CD16+ monocytes",         "Normal",
  "BM endothelial",                    "Normal BM endothelial",          "Normal",
  "B cells",                           "Normal B cells",                 "Normal",
  "T cells",                           "Normal T cells",                 "Normal",
  "Unclassified (low-quality?)",       "Normal-like debris (excl.)",     "Normal",
  "HS5 stromal / MSC",                 "HS-5 stromal (cell line)",       "CellLine"
)

write.table(CELL_TYPE_MAPPING,
            file.path(tab_dir, "cell_type_final_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 2. Load object and remap --------------------------------------------

message("--- Loading Seurat object ---")
seu_merged <- readRDS(opt$rds_in)
DefaultAssay(seu_merged) <- "RNA"

if (!"cell_type" %in% colnames(seu_merged@meta.data)) {
  stop("cell_type column not found — run seurat_01b first")
}

# Build the remap vector keyed by existing cell_type
remap_vec <- setNames(CELL_TYPE_MAPPING$cell_type_final, CELL_TYPE_MAPPING$cell_type)
remap_cat <- setNames(CELL_TYPE_MAPPING$category,        CELL_TYPE_MAPPING$cell_type)

unmapped <- setdiff(unique(as.character(seu_merged$cell_type)), CELL_TYPE_MAPPING$cell_type)
if (length(unmapped) > 0) {
  warning("Unmapped cell_type values (left as NA in cell_type_final): ",
          paste(unmapped, collapse = ", "))
}

new_label    <- remap_vec[as.character(seu_merged$cell_type)]
new_category <- remap_cat[as.character(seu_merged$cell_type)]
names(new_label)    <- colnames(seu_merged)
names(new_category) <- colnames(seu_merged)

# Ordered factor levels — AML first, then Uncertain, then Normal, then CellLine
level_order <- unique(CELL_TYPE_MAPPING$cell_type_final)
seu_merged$cell_type_final  <- factor(new_label, levels = level_order)
seu_merged$is_aml_category  <- factor(new_category,
                                      levels = c("AML", "Uncertain", "Normal", "CellLine"))
seu_merged$is_aml           <- new_category == "AML"

# Count table
counts <- as.data.frame(table(cell_type_final = seu_merged$cell_type_final,
                              is_aml_category = seu_merged$is_aml_category)) |>
  filter(Freq > 0) |>
  arrange(is_aml_category, desc(Freq))
write.table(counts,
            file.path(tab_dir, "cell_type_final_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("--- Final cluster counts ---")
print(counts)

# ---- 3. Re-plot UMAPs (cell_type vs cell_type_final) ---------------------

message("--- Plotting UMAPs ---")
if ("umap" %in% Reductions(seu_merged)) {
  pdf(file.path(fig_dir, "04b_final_umap.pdf"), width = 16, height = 7)
  print(DimPlot(seu_merged,
                group.by = c("cell_type", "cell_type_final"),
                label    = TRUE,
                repel    = TRUE,
                label.size = 3) & NoLegend())
  dev.off()
}

# Category-coloured UMAP — quick visual of AML vs Normal vs Uncertain
if ("umap" %in% Reductions(seu_merged)) {
  cat_cols <- c(AML       = "#D7263D",
                Uncertain = "#C46210",
                Normal    = "#1B98E0",
                CellLine  = "#708090")
  pdf(file.path(fig_dir, "04b_final_umap_category.pdf"), width = 8, height = 7)
  print(DimPlot(seu_merged, group.by = "is_aml_category", cols = cat_cols,
                pt.size = 0.4) +
        ggtitle("AML vs Normal vs Uncertain (post Steps 1+2 validation)"))
  dev.off()
}

# ---- 4. Re-render the signature heatmap with finalised labels -----------
# Re-uses the §1 publication-quality layout from seurat_03 §8b, but on the
# y-axis we now use cell_type_final.

message("--- Re-rendering signature heatmap with finalised labels ---")

ucell_cols <- grep("^UCell_", colnames(seu_merged@meta.data), value = TRUE)
if (length(ucell_cols) == 0) {
  warning("No UCell_* columns found — heatmap skipped")
} else {
  sig_order <- list(
    "vG class (S4A)"      = c("vgalen_hsc_s4a", "vgalen_gmp_s4a", "vgalen_myeloid_s4a"),
    "Stemness / LSC"      = c("lsc17", "eppert_hsc",
                              "van_galen_hsc_like", "van_galen_gmp_like",
                              "van_galen_promono_like"),
    "Quiescence"          = c("cheung_rando_quiescence", "laurenti_hsc_quiescence"),
    "Stress / ISR"        = c("van_den_brink_dissociation", "isr_atf4_targets"),
    "Niche signalling"    = c("baryawno_niche_signal", "tikhonova_msc_response"),
    "Resistance"          = c("arac_metabolism", "mdr_efflux", "aldh_resistance"),
    "Apoptosis"           = c("bcl2_family_pro_survival", "bcl2_family_pro_apoptotic")
  )
  sig_theme <- unlist(lapply(names(sig_order),
                             function(g) setNames(rep(g, length(sig_order[[g]])),
                                                  sig_order[[g]])))
  sig_in_data <- intersect(names(sig_theme),
                           sub("^UCell_", "", ucell_cols))
  sig_theme   <- sig_theme[sig_in_data]

  # Cluster mean UCell scores by cell_type_final
  meta <- seu_merged@meta.data |>
    as_tibble(rownames = "cell") |>
    filter(!is.na(cell_type_final))

  heat_cols <- paste0("UCell_", sig_in_data)
  score_by_cluster <- meta |>
    group_by(cell_type_final, is_aml_category) |>
    summarise(across(all_of(heat_cols), \(x) mean(x, na.rm = TRUE)),
              n_cells = n(), .groups = "drop")

  heat_mat <- as.matrix(score_by_cluster[, heat_cols, drop = FALSE])
  rownames(heat_mat) <- as.character(score_by_cluster$cell_type_final)
  colnames(heat_mat) <- sig_in_data
  heat_z <- scale(heat_mat)

  ann_col <- data.frame(Theme = factor(sig_theme[sig_in_data],
                                       levels = names(sig_order)),
                        row.names = sig_in_data)
  ann_row <- data.frame(Category = score_by_cluster$is_aml_category,
                        row.names = rownames(heat_z))

  theme_colours <- c(
    "vG class (S4A)"     = "#7B2CBF",
    "Stemness / LSC"     = "#D7263D",
    "Quiescence"         = "#2E8B57",
    "Stress / ISR"       = "#E66B00",
    "Niche signalling"   = "#1B98E0",
    "Resistance"         = "#8B4513",
    "Apoptosis"          = "#5A5A5A"
  )
  category_colours <- c(
    AML       = "#D7263D",
    Uncertain = "#C46210",
    Normal    = "#1B98E0",
    CellLine  = "#708090"
  )
  ann_colours <- list(Theme = theme_colours, Category = category_colours)

  zmax <- max(2, ceiling(max(abs(heat_z), na.rm = TRUE)))
  breaks <- seq(-zmax, zmax, length.out = 51)
  heat_palette <- colorRampPalette(c("#1B4F72", "#1B98E0", "white",
                                     "#D7263D", "#6B0F1A"))(50)

  n_per_theme <- table(ann_col$Theme)[unique(ann_col$Theme)]
  col_gaps    <- head(cumsum(n_per_theme), -1)

  # Rows grouped by category — gaps between AML / Uncertain / Normal / CellLine
  row_order <- order(score_by_cluster$is_aml_category,
                     -rowMeans(heat_z[, sig_in_data %in% c("vgalen_hsc_s4a",
                                                           "vgalen_gmp_s4a",
                                                           "vgalen_myeloid_s4a"),
                                       drop = FALSE]))
  heat_z   <- heat_z[row_order, , drop = FALSE]
  ann_row  <- ann_row[row_order, , drop = FALSE]
  cat_seq  <- as.character(score_by_cluster$is_aml_category[row_order])
  row_gaps <- which(diff(as.integer(factor(cat_seq, levels = c("AML","Uncertain","Normal","CellLine")))) != 0)

  pdf(file.path(fig_dir, "04b_final_signature_heatmap.pdf"),
      width  = max(11, 3.0 + ncol(heat_z) * 0.42),
      height = max(7,  2.5 + nrow(heat_z) * 0.42))
  pheatmap(
    heat_z,
    cluster_rows       = FALSE,
    cluster_cols       = FALSE,
    gaps_col           = as.integer(col_gaps),
    gaps_row           = as.integer(row_gaps),
    annotation_col     = ann_col,
    annotation_row     = ann_row,
    annotation_colors  = ann_colours,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    color              = heat_palette,
    breaks             = breaks,
    border_color       = "grey90",
    cellwidth          = 16,
    cellheight         = 16,
    fontsize           = 11,
    fontsize_row       = 11,
    fontsize_col       = 10,
    angle_col          = 45,
    treeheight_row     = 0,
    treeheight_col     = 0,
    main               = "UCell signature scores per cluster (finalised labels, z-scored)"
  )
  dev.off()
}

# ---- 5. Old vs final confusion matrix (sanity check) ---------------------

message("--- Old vs final confusion matrix ---")
conf <- table(old = seu_merged$cell_type, new = seu_merged$cell_type_final)
write.table(conf,
            file.path(tab_dir, "old_vs_final_confusion.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

pdf(file.path(fig_dir, "04b_old_vs_final_confusion.pdf"),
    width  = max(10, ncol(conf) * 0.5 + 4),
    height = max(7,  nrow(conf) * 0.45 + 3))
pheatmap(as.matrix(conf > 0) * 1,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#1B2A4A"))(50),
         border_color = "grey90",
         display_numbers = ifelse(conf > 0, conf, ""),
         number_color = "#FFFFFF",
         fontsize_number = 8,
         main = "Old cell_type (rows) → Finalised cell_type_final (cols)")
dev.off()

# ---- 6. Save updated RDS + session info ----------------------------------

message("--- Saving updated RDS ---")
saveRDS(seu_merged, file = opt$rds_in)
message("  updated: ", opt$rds_in)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_04b_sessionInfo.txt"))

message("Done!")
