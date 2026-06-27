#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 08: Slingshot trajectory inference
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-25
# Last update:  2026-05-25
# Input:        --rds_in    AML-only Seurat .rds from seurat_04c. Must contain
#                           pca_aml + umap reductions, cell_type_final, stroma,
#                           araC, sample_id.
#               --rds_out   Output Seurat .rds with pseudotime columns added.
#               --out_dir   Output directory for figures, tables, logs.
# Output:       <rds_out>                                ← AML-only object with
#                 meta:   slingPseudotime_1 … slingPseudotime_K   (one per lineage)
#                         sling_lineage  (cell's most-assigned lineage)
#               <out_dir>/figures/08a_lineage_umap.pdf           lineage curves
#               <out_dir>/figures/08b_pseudotime_umap.pdf        per-lineage pseudotime
#               <out_dir>/figures/08c_pseudotime_by_cluster.pdf  boxplots
#               <out_dir>/figures/08d_pseudotime_by_condition.pdf conv vs coc per cluster
#               <out_dir>/figures/08e_lineage_topology.pdf       MST + curves
#               <out_dir>/figures/08f_conv_vs_coc_split.pdf      per-condition fits
#               <out_dir>/tables/08_lineage_summary.tsv          terminal/start clusters per lineage
#               <out_dir>/tables/08_pseudotime_by_cluster.tsv    mean/sd per cluster × lineage
#               <out_dir>/tables/08_conv_vs_coc_pseudotime.tsv   Cliff's delta + Wilcoxon
#               <out_dir>/tables/08_tradeseq_top_genes.tsv       (optional, --run_tradeseq)
#               <out_dir>/logs/08_manifest.tsv
#               <out_dir>/logs/seurat_08_sessionInfo.txt
# Depends on:   Seurat (>=5.0), slingshot, SingleCellExperiment, mclust,
#               ggplot2, ggrepel, dplyr, tibble, tidyr, optparse, viridis
#               (optional) tradeSeq for dynamic genes along pseudotime
#
# Notes:
#   * Cohort caveat (carried over from §5a/§5b) ------------------------------
#     One patient (SRAML10), 4 samples — pseudotime is interpretable within
#     this patient but not generalisable across patients. The trajectory
#     topology comparison conv vs coc is descriptive, not inferential at
#     cohort level.
#   * Fitting strategy --------------------------------------------------------
#     - Slingshot is fit on `pca_aml` (30 PCs from seurat_04c). PCA preserves
#       global structure; UMAP distorts it and is used for visualisation only.
#     - clusterLabels = cell_type_final (Steps 1+2 validated labels). Slingshot
#       fits the MST on cluster centroids in PCA space, then smooths principal
#       curves through cells along each branch.
#     - Root cluster: ABCB5+ resistant LSC AML (highest CNV burden + primitive
#       signatures from §1+§2). Falls back to HSPC/LSC-like AML if ABCB5+ has
#       fewer than --min_root_cells cells.
#     - End clusters are NOT specified — letting Slingshot decide which clusters
#       are terminal is informative (e.g. does Promono come out as a terminal
#       branch, or as a midpoint?).
#   * Per-condition fits ------------------------------------------------------
#     A secondary fit is run separately on conv and coc cells to compare
#     topology: same root, same clusters, but each arm gets its own MST.
#     If the topology differs (different number of lineages, different terminal
#     clusters), that is a structural readout of niche effect on differentiation.
#   * Curve embedding on UMAP -------------------------------------------------
#     slingshot::embedCurves projects the PCA-fit curves into UMAP space for
#     visualisation. The pseudotime values themselves are PCA-derived.
#   * tradeSeq (optional) -----------------------------------------------------
#     Behind --run_tradeseq flag (default FALSE). Fits GAMs per gene per
#     lineage; ranks genes by associationTest. Slow — only run when you want
#     the gene-level dynamic atlas.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(optparse)
  library(patchwork)
})

set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),    type = "character",
                     help = "AML-only Seurat .rds (from seurat_04c)")
parser <- add_option(parser, c("--rds_out"),   type = "character",
                     help = "Output Seurat .rds with pseudotime metadata")
parser <- add_option(parser, c("--out_dir"),   type = "character",
                     help = "Output directory")
parser <- add_option(parser, c("--root_cluster"), type = "character",
                     default = "ABCB5+ resistant LSC AML",
                     help = "Starting cluster for Slingshot (default ABCB5+)")
parser <- add_option(parser, c("--root_fallback"), type = "character",
                     default = "HSPC / LSC-like AML",
                     help = "Fallback root if primary has too few cells")
parser <- add_option(parser, c("--min_root_cells"), type = "integer", default = 30L,
                     help = "Minimum cells required in root cluster (default 30)")
parser <- add_option(parser, c("--n_pcs"),     type = "integer", default = 30L,
                     help = "PCs used (must match seurat_04c, default 30)")
parser <- add_option(parser, c("--run_tradeseq"), action = "store_true",
                     default = FALSE,
                     help = "Run tradeSeq dynamic gene analysis (slow)")
parser <- add_option(parser, c("--tradeseq_top_genes"), type = "integer", default = 2000L,
                     help = "HVGs to fit if tradeSeq is on (default 2000)")
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

message("=== Seurat 08: Slingshot trajectory ===")
message("  rds_in     : ", opt$rds_in)
message("  rds_out    : ", opt$rds_out)
message("  out_dir    : ", opt$out_dir)
message("  root       : ", opt$root_cluster)
message("  run_tradeseq: ", opt$run_tradeseq)

# ---- 1. Load AML-only object + validate ----------------------------------

message("--- Loading AML-only Seurat object ---")
seu_aml <- readRDS(opt$rds_in)
DefaultAssay(seu_aml) <- "RNA"

for (req_col in c("cell_type_final", "stroma", "araC", "sample_id")) {
  if (!req_col %in% colnames(seu_aml@meta.data)) {
    stop("Required metadata column missing: ", req_col)
  }
}

# Find the AML-only PCA reduction (case-insensitive — Seurat can save it as
# pca_aml, PCAML, etc., and as.SingleCellExperiment preserves the key)
red_names <- Reductions(seu_aml)
pca_red <- red_names[tolower(red_names) %in% c("pca_aml", "pcaml")]
if (length(pca_red) == 0) {
  message("WARN: pca_aml not found — falling back to 'pca'")
  if (!"pca" %in% red_names) stop("No PCA reduction found in the object")
  pca_red <- "pca"
}
pca_red <- pca_red[1]
umap_red <- if ("umap" %in% red_names) "umap" else red_names[grep("UMAP", red_names, ignore.case = TRUE)[1]]

message("  Using PCA reduction: ", pca_red)
message("  Using UMAP reduction: ", umap_red)

# Drop unused factor levels so MST doesn't include empty clusters
seu_aml$cell_type_final <- droplevels(factor(seu_aml$cell_type_final))
cluster_levels <- levels(seu_aml$cell_type_final)
message("  Clusters in AML object: ", paste(cluster_levels, collapse = ", "))

# Resolve root cluster
root_used <- opt$root_cluster
if (!(root_used %in% cluster_levels) ||
    sum(seu_aml$cell_type_final == root_used) < opt$min_root_cells) {
  message("Root cluster '", root_used, "' missing or has < ",
          opt$min_root_cells, " cells. Using fallback: '",
          opt$root_fallback, "'")
  root_used <- opt$root_fallback
  if (!(root_used %in% cluster_levels)) {
    stop("Fallback root '", root_used, "' not in cluster levels either. Aborting.")
  }
}
message("  Root cluster for Slingshot: ", root_used)

# ---- 2. Helper: convert Seurat → SCE and run Slingshot -------------------

run_slingshot <- function(seu_obj, root_cluster, n_pcs,
                          pca_name = pca_red, umap_name = umap_red,
                          label = "all") {
  message("    Slingshot fit [", label, "] — ",
          ncol(seu_obj), " cells, root = '", root_cluster, "'")
  sce <- as.SingleCellExperiment(seu_obj)
  
  rd_names <- SingleCellExperiment::reducedDimNames(sce)
  message("    Available reductions in SingleCellExperiment: ", 
          paste(rd_names, collapse = ", "))
  
  # Resolve PCA name dynamically (Seurat to SCE conversion can change case/name)
  pca_in_sce <- rd_names[tolower(rd_names) %in% tolower(pca_name)]
  if (length(pca_in_sce) == 0) {
    pca_in_sce <- rd_names[grep("pca", rd_names, ignore.case = TRUE)]
  }
  if (length(pca_in_sce) == 0) {
    stop("Could not find any PCA reduction in SingleCellExperiment object")
  }
  pca_in_sce <- pca_in_sce[1]
  message("    Using PCA in SCE: ", pca_in_sce)
  
  # Resolve UMAP name dynamically
  umap_in_sce <- rd_names[tolower(rd_names) %in% tolower(umap_name)]
  if (length(umap_in_sce) == 0) {
    umap_in_sce <- rd_names[grep("umap", rd_names, ignore.case = TRUE)]
  }
  umap_in_sce <- umap_in_sce[1]
  message("    Using UMAP in SCE: ", umap_in_sce)
  
  # Truncate PCA to n_pcs in case the saved reduction has more
  rd <- SingleCellExperiment::reducedDim(sce, pca_in_sce)
  if (ncol(rd) > n_pcs) {
    SingleCellExperiment::reducedDim(sce, pca_in_sce) <- rd[, seq_len(n_pcs)]
  }
  sce <- tryCatch(
    slingshot::slingshot(sce,
                         clusterLabels = "cell_type_final",
                         reducedDim    = pca_in_sce,
                         start.clus    = root_cluster,
                         approx_points = 150),
    error = function(e) {
      message("    Slingshot error: ", e$message)
      NULL
    }
  )
  if (is.null(sce)) return(NULL)
  # Embed curves on UMAP for plotting
  emb <- tryCatch(
    slingshot::embedCurves(sce, umap_in_sce),
    error = function(e) {
      message("    embedCurves failed: ", e$message)
      NULL
    }
  )
  list(sce = sce, embed = emb, label = label)
}

# ---- 3. Global fit (all AML cells) ---------------------------------------

message("--- Global Slingshot fit ---")
fit_all <- run_slingshot(seu_aml, root_used, opt$n_pcs, label = "all")
if (is.null(fit_all)) stop("Global Slingshot fit failed — see error above.")

pst_all     <- slingshot::slingPseudotime(fit_all$sce)
lineages    <- slingshot::slingLineages(fit_all$sce)
curves_pca  <- slingshot::slingCurves(fit_all$sce, as.df = FALSE)
n_lineages  <- ncol(pst_all)

message("  Lineages discovered: ", n_lineages)
for (i in seq_len(n_lineages)) {
  message("    L", i, ": ", paste(lineages[[i]], collapse = " → "))
}

# Lineage summary table
lineage_tbl <- tibble(
  lineage_id    = paste0("L", seq_along(lineages)),
  n_clusters    = lengths(lineages),
  start_cluster = vapply(lineages, function(x) x[1], character(1)),
  end_cluster   = vapply(lineages, function(x) x[length(x)], character(1)),
  path          = vapply(lineages,
                         function(x) paste(x, collapse = " → "),
                         character(1))
)
write.table(lineage_tbl,
            file.path(tab_dir, "08_lineage_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(lineage_tbl)

# ---- 4. Write pseudotime back to Seurat object metadata ------------------

for (i in seq_len(n_lineages)) {
  seu_aml[[paste0("slingPseudotime_", i)]] <- pst_all[, i]
}

# Most-assigned lineage per cell: pick the lineage with smallest pseudotime
# (i.e. the lineage along which the cell sits closest, by Slingshot's weights)
weights <- slingshot::slingCurveWeights(fit_all$sce)
seu_aml$sling_lineage <- paste0("L", apply(weights, 1, which.max))

# ---- 5. UMAP plots: lineage curves + pseudotime --------------------------

message("--- Plotting Slingshot curves on UMAP ---")

# Pull UMAP coords
umap_df <- as.data.frame(Embeddings(seu_aml, umap_red))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell_type_final <- seu_aml$cell_type_final
umap_df$condition <- ifelse(seu_aml$stroma, "Co-culture", "Conventional")
umap_df$treatment <- ifelse(seu_aml$araC, "Ara-C", "Untreated")

cluster_cols <- c(
  "HSPC / LSC-like AML"           = "#7B2CBF",
  "ABCB5+ resistant LSC AML"      = "#5A189A",
  "GMP-like AML"                  = "#D7263D",
  "Promono quiescent AML"         = "#2E8B57",
  "Mono-like AML"                 = "#E66B00",
  "Stress-response AML"           = "#C46210",
  "Niche-stressed primitive AML"  = "#1B98E0"
)
cluster_cols <- cluster_cols[names(cluster_cols) %in% cluster_levels]

# Curves on UMAP (if embedCurves succeeded)
curve_df <- NULL
if (!is.null(fit_all$embed)) {
  curve_list <- slingshot::slingCurves(fit_all$embed, as.df = TRUE)
  curve_df   <- curve_list  # already a data.frame with Lineage + s1, s2
  # Standardise column names
  colnames(curve_df)[colnames(curve_df) %in% c("Dim.1", "DC1", "UMAP_1")] <- "x"
  if (!"x" %in% colnames(curve_df)) {
    # Fall back: assume first two numeric columns are coords
    num_cols <- which(sapply(curve_df, is.numeric))[1:2]
    colnames(curve_df)[num_cols] <- c("x", "y")
  } else if (!"y" %in% colnames(curve_df)) {
    num_cols <- which(sapply(curve_df, is.numeric))
    y_col    <- setdiff(num_cols, which(colnames(curve_df) == "x"))[1]
    colnames(curve_df)[y_col] <- "y"
  }
}

pdf(file.path(fig_dir, "08a_lineage_umap.pdf"), width = 9, height = 7)
p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = cell_type_final), size = 0.4, alpha = 0.65) +
  scale_colour_manual(values = cluster_cols, name = "Cluster") +
  theme_classic(base_size = 11) +
  ggtitle(paste0("AML lineages — Slingshot (root: ", root_used, ")"))
if (!is.null(curve_df) && all(c("x", "y") %in% colnames(curve_df))) {
  p <- p + geom_path(data = curve_df,
                     aes(x = x, y = y, group = Lineage),
                     colour = "black", linewidth = 0.9)
}
print(p)
dev.off()

# Pseudotime UMAPs — one panel per lineage
pdf(file.path(fig_dir, "08b_pseudotime_umap.pdf"),
    width  = max(7, 4 * min(n_lineages, 3)),
    height = max(5, 4 * ceiling(n_lineages / 3)))
plots <- lapply(seq_len(n_lineages), function(i) {
  d <- umap_df
  d$pst <- pst_all[, i]
  ggplot(d, aes(UMAP_1, UMAP_2, colour = pst)) +
    geom_point(size = 0.4, alpha = 0.65) +
    scale_colour_viridis_c(name = "Pseudotime",
                           option = "magma", na.value = "grey85") +
    theme_classic(base_size = 10) +
    ggtitle(paste0("Lineage ", i, ": ", paste(lineages[[i]], collapse = " → ")))
})
print(patchwork::wrap_plots(plots, ncol = min(3, n_lineages)))
dev.off()

# ---- 6. Pseudotime distribution per cluster ------------------------------

message("--- Pseudotime per cluster ---")
pst_long <- pst_all |>
  as.data.frame() |>
  rownames_to_column("cell") |>
  pivot_longer(cols = -cell,
               names_to = "lineage", values_to = "pst") |>
  mutate(cell_type_final = seu_aml$cell_type_final[match(cell, colnames(seu_aml))],
         condition       = ifelse(seu_aml$stroma[match(cell, colnames(seu_aml))],
                                  "Co-culture", "Conventional")) |>
  filter(!is.na(pst))

# Summary table: mean/median/sd per cluster × lineage
cluster_pst <- pst_long |>
  group_by(cell_type_final, lineage) |>
  summarise(n = n(),
            mean_pst   = mean(pst, na.rm = TRUE),
            median_pst = median(pst, na.rm = TRUE),
            sd_pst     = sd(pst, na.rm = TRUE),
            .groups = "drop") |>
  arrange(lineage, mean_pst)
write.table(cluster_pst,
            file.path(tab_dir, "08_pseudotime_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(fig_dir, "08c_pseudotime_by_cluster.pdf"),
    width  = max(8, 1.5 * length(cluster_levels)),
    height = max(5, 2.5 * n_lineages))
print(
  ggplot(pst_long,
         aes(x = reorder(cell_type_final, pst, FUN = median),
             y = pst, fill = cell_type_final)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    facet_wrap(~ lineage, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = NULL, y = "Slingshot pseudotime",
         title = "Pseudotime per cluster, per lineage")
)
dev.off()

# ---- 7. Per-cluster conv-vs-coc pseudotime comparison --------------------

message("--- Conv vs coc pseudotime, per cluster × lineage ---")

cliff_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 5 || length(y) < 5) return(NA_real_)
  if (requireNamespace("effsize", quietly = TRUE)) {
    return(as.numeric(effsize::cliff.delta(x, y)$estimate))
  }
  # Rank-biserial fallback
  u <- as.numeric(wilcox.test(x, y, exact = FALSE)$statistic)
  1 - (2 * u) / (length(x) * length(y))
}

cond_compare <- pst_long |>
  group_by(cell_type_final, lineage) |>
  summarise(
    n_conv = sum(condition == "Conventional"),
    n_coc  = sum(condition == "Co-culture"),
    delta  = cliff_delta(pst[condition == "Co-culture"],
                         pst[condition == "Conventional"]),
    p_value = ifelse(n_conv >= 5 & n_coc >= 5,
                     wilcox.test(pst[condition == "Co-culture"],
                                 pst[condition == "Conventional"],
                                 exact = FALSE)$p.value,
                     NA_real_),
    .groups = "drop"
  ) |>
  group_by(lineage) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(lineage, desc(abs(delta)))
write.table(cond_compare,
            file.path(tab_dir, "08_conv_vs_coc_pseudotime.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(fig_dir, "08d_pseudotime_by_condition.pdf"),
    width  = max(8, 1.5 * length(cluster_levels)),
    height = max(5, 2.5 * n_lineages))
print(
  ggplot(pst_long,
         aes(x = cell_type_final, y = pst, fill = condition)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.85,
                 position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = c(Conventional = "#1B2A4A",
                                 `Co-culture`  = "#2E8B57")) +
    facet_wrap(~ lineage, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = NULL, y = "Slingshot pseudotime",
         title = "Pseudotime, conventional vs co-culture per cluster")
)
dev.off()

# ---- 8. Per-condition fits (conv only, coc only) -------------------------

message("--- Per-condition Slingshot fits ---")
seu_conv <- seu_aml[, !seu_aml$stroma]
seu_coc  <- seu_aml[,  seu_aml$stroma]

# Need the same root cluster to exist in each subset
have_root <- function(seu, root) {
  root %in% levels(droplevels(factor(seu$cell_type_final))) &&
    sum(seu$cell_type_final == root) >= opt$min_root_cells
}
root_conv <- if (have_root(seu_conv, root_used)) root_used else opt$root_fallback
root_coc  <- if (have_root(seu_coc,  root_used)) root_used else opt$root_fallback

fit_conv <- run_slingshot(seu_conv, root_conv, opt$n_pcs, label = "conv")
fit_coc  <- run_slingshot(seu_coc,  root_coc,  opt$n_pcs, label = "coc")

per_cond_topology <- function(fit, label, root) {
  if (is.null(fit)) {
    return(tibble(arm = label, root = root, status = "failed",
                  n_lineages = NA_integer_, lineages = NA_character_))
  }
  lins <- slingshot::slingLineages(fit$sce)
  tibble(
    arm        = label,
    root       = root,
    status     = "ok",
    n_lineages = length(lins),
    lineages   = paste(vapply(lins, function(x) paste(x, collapse = " → "),
                              character(1)),
                       collapse = " ; ")
  )
}
topology_tbl <- bind_rows(
  per_cond_topology(fit_all,  "all",  root_used),
  per_cond_topology(fit_conv, "conv", root_conv),
  per_cond_topology(fit_coc,  "coc",  root_coc)
)
write.table(topology_tbl,
            file.path(tab_dir, "08_topology_by_arm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(topology_tbl)

# Side-by-side UMAP: conv vs coc curves
plot_arm <- function(fit, arm_label) {
  if (is.null(fit)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(arm_label, " — Slingshot failed")))
  }
  d <- as.data.frame(Embeddings(seu_aml[, colnames(fit$sce)], umap_red))
  colnames(d) <- c("UMAP_1", "UMAP_2")
  d$cell_type_final <- seu_aml$cell_type_final[colnames(fit$sce)]
  p <- ggplot(d, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = cell_type_final), size = 0.4, alpha = 0.65) +
    scale_colour_manual(values = cluster_cols, guide = "none") +
    theme_classic(base_size = 10) +
    ggtitle(paste0(arm_label, " — n_lineages = ",
                   length(slingshot::slingLineages(fit$sce))))
  if (!is.null(fit$embed)) {
    cd <- slingshot::slingCurves(fit$embed, as.df = TRUE)
    if (!"x" %in% colnames(cd)) {
      num_cols <- which(sapply(cd, is.numeric))[1:2]
      colnames(cd)[num_cols] <- c("x", "y")
    } else if (!"y" %in% colnames(cd)) {
      num_cols <- which(sapply(cd, is.numeric))
      y_col    <- setdiff(num_cols, which(colnames(cd) == "x"))[1]
      colnames(cd)[y_col] <- "y"
    }
    p <- p + geom_path(data = cd, aes(x = x, y = y, group = Lineage),
                       colour = "black", linewidth = 0.9)
  }
  p
}
pdf(file.path(fig_dir, "08f_conv_vs_coc_split.pdf"), width = 14, height = 6)
print(plot_arm(fit_conv, "Conventional") | plot_arm(fit_coc, "Co-culture"))
dev.off()

# ---- 9. Optional tradeSeq dynamic gene analysis --------------------------

if (opt$run_tradeseq) {
  message("--- Running tradeSeq dynamic-gene analysis (this is slow) ---")
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    message("    tradeSeq not installed — skipping. Install via BiocManager::install('tradeSeq').")
  } else {
    # Use top HVGs to keep runtime tractable
    seu_aml <- FindVariableFeatures(seu_aml,
                                    selection.method = "vst",
                                    nfeatures = opt$tradeseq_top_genes,
                                    verbose = FALSE)
    hvg <- VariableFeatures(seu_aml)
    counts_mat <- GetAssayData(seu_aml, layer = "counts")[hvg, ]

    sce_ts <- tradeSeq::fitGAM(
      counts        = as.matrix(counts_mat),
      pseudotime    = pst_all,
      cellWeights   = slingshot::slingCurveWeights(fit_all$sce),
      nknots        = 6,
      verbose       = FALSE
    )
    assoc <- tradeSeq::associationTest(sce_ts) |>
      as.data.frame() |>
      rownames_to_column("gene") |>
      arrange(pvalue)
    write.table(assoc,
                file.path(tab_dir, "08_tradeseq_top_genes.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("  tradeSeq complete — top genes written.")
  }
}

# ---- 10. Save updated AML object + manifest ------------------------------

message("--- Saving updated RDS ---")
saveRDS(seu_aml, file = opt$rds_out)

manifest <- tibble::tibble(
  field = c("rds_in", "rds_out",
            "n_cells", "n_clusters", "n_lineages",
            "root_cluster_used", "root_cluster_conv", "root_cluster_coc",
            "pca_reduction", "n_pcs",
            "run_tradeseq", "seed", "date"),
  value = c(opt$rds_in, opt$rds_out,
            ncol(seu_aml), length(cluster_levels), n_lineages,
            root_used, root_conv, root_coc,
            pca_red, opt$n_pcs,
            as.character(opt$run_tradeseq), 42L,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
write.table(manifest, file.path(log_dir, "08_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_08_sessionInfo.txt"))

message("Done!")
