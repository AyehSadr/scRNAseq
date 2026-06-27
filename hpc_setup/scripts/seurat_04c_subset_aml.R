#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 04c: Subset to AML-only compartment
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-24
# Last update:  2026-05-24
# Input:        --rds_in    Seurat .rds AFTER seurat_04b — must contain
#                           cell_type_final / is_aml / is_aml_category, plus
#                           the §1 signature scores and §2 cnv_score.
#               --rds_out   Output path for the AML-only Seurat .rds
#               --out_dir   Output directory for QC plots and tables
# Output:       <rds_out>                                  ← AML-only object
#                 reductions:   pca_aml     (AML-specific PCs)
#                               umap_aml    (AML-only UMAP, primary)
#                               umap_full   (subset of full-object UMAP, kept
#                                            for cross-reference only)
#                 meta keeps:   cell_type_final, is_aml_category, condition,
#                               treatment, sample_id, UCell_*, cnv_score
#               <out_dir>/figures/04c_aml_umap.pdf
#               <out_dir>/figures/04c_aml_umap_condition.pdf
#               <out_dir>/figures/04c_aml_umap_treatment.pdf
#               <out_dir>/tables/04c_aml_counts.tsv
#               <out_dir>/logs/seurat_04c_sessionInfo.txt
#               <out_dir>/logs/04c_manifest.tsv
# Depends on:   Seurat (>=5.0), ggplot2, dplyr, tibble, optparse
# Notes:
#   * Subsets to is_aml == TRUE — Normal (CD16+ mono, B, T, BM endo), Uncertain
#     (eos/baso), and CellLine (HS-5) cells are dropped.
#   * Re-runs FindVariableFeatures + ScaleData + RunPCA on the AML subset so
#     the new PCs reflect AML-internal variance rather than lineage variance
#     dominated by HS-5 and B/T cells. Cluster identities (cell_type_final)
#     are NOT recomputed — they stay as the Steps 1+2 validated labels.
#   * The original umap (full-object embedding, subset to AML cells) is kept
#     under the name `umap_full` so any figure that needs the global view
#     for cross-reference still has it.
#   * Downstream consumers:
#       seurat_05_stress_cycling.R   →  AML-only  (this object)
#       seurat_07_pseudobulk_de.R    →  AML-only  (planned)
#       seurat_09_milo.R             →  AML-only  (planned)
#       seurat_10_liana.R            →  FULL      (needs HS-5)
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),  type = "character",
                     help = "Seurat .rds after seurat_04b (full object)")
parser <- add_option(parser, c("--rds_out"), type = "character",
                     help = "Output Seurat .rds for AML-only object")
parser <- add_option(parser, c("--out_dir"), type = "character",
                     help = "Output directory for QC figures / tables / logs")
parser <- add_option(parser, c("--ncores"),  type = "integer", default = 4L,
                     help = "Cores for ScaleData (default 4)")
parser <- add_option(parser, c("--n_pcs"),   type = "integer", default = 30L,
                     help = "Number of PCs for AML-only PCA / UMAP (default 30)")
parser <- add_option(parser, c("--n_hvg"),   type = "integer", default = 2000L,
                     help = "Number of HVGs for AML-only re-embedding (default 2000)")
opt <- parse_args(parser)

for (req in c("rds_in", "rds_out", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(tab_dir, fig_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Seurat 04c: subset to AML-only compartment ===")
message("  rds_in  : ", opt$rds_in)
message("  rds_out : ", opt$rds_out)
message("  out_dir : ", opt$out_dir)
message("  n_pcs   : ", opt$n_pcs)
message("  n_hvg   : ", opt$n_hvg)

# ---- 1. Load full object --------------------------------------------------

message("--- Loading full Seurat object ---")
seu_full <- readRDS(opt$rds_in)
DefaultAssay(seu_full) <- "RNA"

for (req_col in c("cell_type_final", "is_aml", "is_aml_category")) {
  if (!req_col %in% colnames(seu_full@meta.data)) {
    stop("Required metadata column missing: ", req_col,
         " — run seurat_04b first.")
  }
}

# Pre-subset count for the audit log
pre <- as.data.frame(table(is_aml_category = seu_full$is_aml_category)) |>
  arrange(is_aml_category)
message("--- Pre-subset compartment counts ---")
print(pre)

# ---- 2. Subset to AML cells ----------------------------------------------

aml_cells <- colnames(seu_full)[which(seu_full$is_aml)]
if (length(aml_cells) < 100) {
  stop("Suspiciously few AML cells (", length(aml_cells),
       ") — check is_aml column.")
}
message("--- Subsetting to ", length(aml_cells), " AML cells ---")

seu_aml <- subset(seu_full, cells = aml_cells)
DefaultAssay(seu_aml) <- "RNA"

# Drop unused factor levels in cell_type_final so plots are tidy
seu_aml$cell_type_final <- droplevels(seu_aml$cell_type_final)
seu_aml$is_aml_category <- droplevels(seu_aml$is_aml_category)

post <- as.data.frame(table(cell_type_final = seu_aml$cell_type_final)) |>
  filter(Freq > 0) |>
  arrange(desc(Freq))
write.table(post, file.path(tab_dir, "04c_aml_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("--- Post-subset AML cluster counts ---")
print(post)

# ---- 3. Stash full-object UMAP under a new name --------------------------
# The cells.use side of subset() already restricts the existing reductions to
# AML cells; we just rename `umap` → `umap_full` so the AML-specific UMAP we
# build next can live under the canonical `umap` slot.

if ("umap" %in% Reductions(seu_aml)) {
  seu_aml[["umap_full"]] <- seu_aml[["umap"]]
  seu_aml[["umap"]] <- NULL
  message("--- Stashed full-object umap → umap_full ---")
}
if ("pca" %in% Reductions(seu_aml)) {
  seu_aml[["pca_full"]] <- seu_aml[["pca"]]
  seu_aml[["pca"]] <- NULL
  message("--- Stashed full-object pca → pca_full ---")
}

# ---- 4. Re-embed on AML-only variance ------------------------------------
# Re-run HVG / Scale / PCA / Neighbours / UMAP on the AML compartment so the
# embedding reflects AML-internal variance, not contrasts with stroma /
# B / T cells. Clustering is NOT redone — cell_type_final is preserved.

message("--- Re-running variable features + PCA on AML subset ---")
seu_aml <- NormalizeData(seu_aml, verbose = FALSE)
seu_aml <- FindVariableFeatures(seu_aml, selection.method = "vst",
                                nfeatures = opt$n_hvg, verbose = FALSE)
seu_aml <- ScaleData(seu_aml, verbose = FALSE)
seu_aml <- RunPCA(seu_aml, npcs = opt$n_pcs, verbose = FALSE,
                  reduction.name = "pca_aml", reduction.key = "PCaml_")

message("--- Building neighbours + UMAP on AML subset ---")
seu_aml <- FindNeighbors(seu_aml,
                         reduction = "pca_aml",
                         dims = 1:opt$n_pcs,
                         graph.name = c("aml_nn", "aml_snn"),
                         verbose = FALSE)
seu_aml <- RunUMAP(seu_aml,
                   reduction = "pca_aml",
                   dims = 1:opt$n_pcs,
                   reduction.name = "umap",
                   reduction.key  = "UMAPaml_",
                   verbose = FALSE)

# ---- 5. QC plots ----------------------------------------------------------

message("--- Plotting QC UMAPs ---")
cluster_cols <- c(
  "HSPC / LSC-like AML"            = "#7B2CBF",
  "ABCB5+ resistant LSC AML"       = "#5A189A",
  "GMP-like AML"                   = "#D7263D",
  "Promono quiescent AML"          = "#2E8B57",
  "Mono-like AML"                  = "#E66B00",
  "Stress-response AML"            = "#C46210",
  "Niche-stressed primitive AML"   = "#1B98E0"
)
have <- intersect(names(cluster_cols), levels(seu_aml$cell_type_final))
cluster_cols <- cluster_cols[have]

pdf(file.path(fig_dir, "04c_aml_umap.pdf"), width = 9, height = 7)
print(DimPlot(seu_aml,
              reduction = "umap",
              group.by  = "cell_type_final",
              cols      = cluster_cols,
              label     = TRUE, repel = TRUE, pt.size = 0.3) +
        ggtitle("AML compartment — cell_type_final"))
dev.off()

if ("condition" %in% colnames(seu_aml@meta.data)) {
  pdf(file.path(fig_dir, "04c_aml_umap_condition.pdf"), width = 9, height = 7)
  print(DimPlot(seu_aml,
                reduction = "umap",
                group.by  = "condition",
                cols      = c(Conventional = "#1B2A4A", `Co-culture` = "#2E8B57"),
                pt.size   = 0.3) +
          ggtitle("AML compartment — condition"))
  dev.off()
}
if ("treatment" %in% colnames(seu_aml@meta.data)) {
  pdf(file.path(fig_dir, "04c_aml_umap_treatment.pdf"), width = 9, height = 7)
  print(DimPlot(seu_aml,
                reduction = "umap",
                group.by  = "treatment",
                cols      = c(Untreated = "#1B98E0", `Ara-C` = "#D7263D"),
                pt.size   = 0.3) +
          ggtitle("AML compartment — treatment"))
  dev.off()
}
if ("cnv_score" %in% colnames(seu_aml@meta.data)) {
  pdf(file.path(fig_dir, "04c_aml_umap_cnv.pdf"), width = 9, height = 7)
  print(FeaturePlot(seu_aml, features = "cnv_score", reduction = "umap",
                    pt.size = 0.3, order = TRUE) +
          ggtitle("AML compartment — inferCNV burden"))
  dev.off()
}

# Side-by-side: full-object UMAP (subset) vs new AML UMAP
if ("umap_full" %in% Reductions(seu_aml)) {
  pdf(file.path(fig_dir, "04c_umap_full_vs_aml.pdf"), width = 16, height = 7)
  p1 <- DimPlot(seu_aml, reduction = "umap_full", group.by = "cell_type_final",
                cols = cluster_cols, label = TRUE, repel = TRUE, pt.size = 0.3) +
          ggtitle("umap_full (global embedding, AML cells only)") + NoLegend()
  p2 <- DimPlot(seu_aml, reduction = "umap", group.by = "cell_type_final",
                cols = cluster_cols, label = TRUE, repel = TRUE, pt.size = 0.3) +
          ggtitle("umap (AML-only re-embedding)") + NoLegend()
  print(p1 | p2)
  dev.off()
}

# ---- 6. Save AML-only RDS + manifest -------------------------------------

message("--- Saving AML-only RDS ---")
saveRDS(seu_aml, file = opt$rds_out)
message("  written: ", opt$rds_out)

manifest <- tibble::tibble(
  field   = c("rds_in", "rds_out", "n_cells_full", "n_cells_aml",
              "n_clusters_aml", "n_pcs", "n_hvg", "seed", "date"),
  value   = c(opt$rds_in, opt$rds_out,
              ncol(seu_full), ncol(seu_aml),
              length(unique(as.character(seu_aml$cell_type_final))),
              opt$n_pcs, opt$n_hvg, 42L,
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
write.table(manifest, file.path(log_dir, "04c_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_04c_sessionInfo.txt"))

message("Done!")
