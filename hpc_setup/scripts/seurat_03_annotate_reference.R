#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 03: Reference-based annotation + signature scoring
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-23
# Last update:  2026-05-23
# Input:        --rds_in       Seurat .rds from seurat_01b_annotate_clusters.R
#               --refs_dir     directory populated by 01_download_references.sh
#               --sigs_yml     gene-set YAML (config/signatures.yml)
#               --out_dir      output for tables + figures
# Output:       <rds_in>                              ← re-saved with new metadata
#               <out_dir>/tables/singler_main_calls.tsv
#               <out_dir>/tables/azimuth_calls.tsv
#               <out_dir>/tables/ucell_scores_by_cluster.tsv
#               <out_dir>/tables/cell_type_vs_singler_confusion.tsv
#               <out_dir>/figures/03a_singler_umap.pdf
#               <out_dir>/figures/03b_azimuth_umap.pdf
#               <out_dir>/figures/03c_confusion_existing_vs_singler.pdf
#               <out_dir>/figures/03d_ucell_heatmap.pdf
#               <out_dir>/logs/seurat_03_sessionInfo.txt
# Depends on:   Seurat (>=5.0), SingleR, celldex, SummarizedExperiment,
#               UCell (>=2.6), yaml, pheatmap, ggplot2, dplyr,
#               (optional) Azimuth + SeuratData "bonemarrowref"
# Notes:
#   * Modifies seu_merged in place: adds  singler_main, singler_fine,
#     azimuth_l1, azimuth_l2, azimuth_score, plus one UCell_<sig> column
#     per signature in signatures.yml.
#   * Van Galen 2019 (GSE116256) is built from the GEO RAW tarball on first
#     run, cached as <refs_dir>/van_galen_2019/van_galen_2019_singler_ref.rds.
#   * Azimuth is optional — if the package or "bonemarrowref" is not present
#     the block is skipped with a warning, the rest still runs.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SummarizedExperiment)
  library(UCell)
  library(yaml)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 20 * 1024^3)   # 20 GB for parallel workers

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),  type = "character",
                     help = "Seurat .rds from seurat_01b_annotate_clusters.R")
parser <- add_option(parser, c("--refs_dir"), type = "character", default = NULL,
                     help = "Directory populated by 01_download_references.sh [optional]")
parser <- add_option(parser, c("--sigs_yml"), type = "character",
                     help = "Signatures YAML (config/signatures.yml)")
parser <- add_option(parser, c("--out_dir"),  type = "character",
                     help = "Output directory for tables and figures")
parser <- add_option(parser, c("--ncores"),   type = "integer", default = 8,
                     help = "Cores for SingleR / UCell [default: 8]")
opt <- parse_args(parser)

for (req in c("rds_in", "sigs_yml", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in))   stop("Input RDS not found: ", opt$rds_in)
if (!file.exists(opt$sigs_yml)) stop("Signatures YAML not found: ", opt$sigs_yml)

# Check if refs_dir is provided and valid
has_refs <- !is.null(opt$refs_dir) && dir.exists(opt$refs_dir)
if (has_refs) {
  # Double check that there are files inside
  ref_files <- list.files(opt$refs_dir)
  if (length(ref_files) == 0) {
    message("  refs_dir provided but is empty — reference annotation will be skipped")
    has_refs <- FALSE
  }
} else {
  message("  refs_dir not provided or does not exist — reference annotation will be skipped")
}

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(tab_dir, fig_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== Seurat 03: reference annotation + signature scoring ===")
message("  rds_in   : ", opt$rds_in)
message("  refs_dir : ", ifelse(has_refs, opt$refs_dir, "Not used/skipped"))
message("  sigs_yml : ", opt$sigs_yml)
message("  out_dir  : ", opt$out_dir)

# ---- 1. Load query Seurat object -----------------------------------------

message("--- Loading Seurat object ---")
seu_merged <- readRDS(opt$rds_in)
DefaultAssay(seu_merged) <- "RNA"
if ("cell_type" %in% colnames(seu_merged@meta.data)) {
  message("  existing cell_type column found — will be retained alongside new calls")
} else {
  warning("No cell_type column on input; run seurat_01b first if you want a comparison")
}
message("  ", ncol(seu_merged), " cells × ", nrow(seu_merged), " features")

# ---- 2. Helper: build / load van Galen 2019 SingleR reference ------------

build_van_galen_reference <- function(refs_dir) {
  cache <- file.path(refs_dir, "van_galen_2019", "van_galen_2019_singler_ref.rds")
  if (file.exists(cache)) {
    message("  van Galen reference loaded from cache: ", cache)
    return(readRDS(cache))
  }
  vg_dir <- file.path(refs_dir, "van_galen_2019")
  dem_files <- list.files(vg_dir, pattern = "\\.dem\\.txt\\.gz$", full.names = TRUE)
  ann_files <- list.files(vg_dir, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE)
  if (length(dem_files) == 0 || length(ann_files) == 0) {
    stop("van Galen files not found under ", vg_dir,
         " — run 01_download_references.sh first")
  }
  message("  building van Galen reference from ", length(dem_files), " samples (slow, ~5 min)...")
  exprs_list <- lapply(dem_files, function(f) {
    m <- as.matrix(read.table(f, sep = "\t", header = TRUE, row.names = 1,
                              check.names = FALSE))
    storage.mode(m) <- "double"
    m
  })
  anno_list  <- lapply(ann_files, function(f) {
    read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  # Inner-join genes across samples
  common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
  message("    common genes across samples: ", length(common_genes))
  exprs_mat <- do.call(cbind, lapply(exprs_list, function(m) m[common_genes, , drop = FALSE]))
  anno_df   <- do.call(rbind, anno_list)
  rownames(anno_df) <- anno_df$Cell
  # Drop cells with no PredictionRefined label (the column with cell-type calls)
  label_col <- intersect(c("PredictionRefined", "CellType", "Prediction"), colnames(anno_df))[1]
  if (is.na(label_col)) stop("Could not find a label column in van Galen annotations")
  keep_cells <- intersect(colnames(exprs_mat), rownames(anno_df))
  keep_cells <- keep_cells[!is.na(anno_df[keep_cells, label_col]) &
                             anno_df[keep_cells, label_col] != ""]
  exprs_mat  <- exprs_mat[, keep_cells, drop = FALSE]
  anno_df    <- anno_df[keep_cells, , drop = FALSE]
  vg_ref <- SummarizedExperiment(
    assays  = list(logcounts = log2(exprs_mat + 1)),
    colData = DataFrame(label = anno_df[[label_col]])
  )
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(vg_ref, cache)
  message("  cached: ", cache)
  vg_ref
}

# ---- 3. SingleR — broad reference (HumanPrimaryCellAtlas) ----------------

if (has_refs) {
  message("--- SingleR: HumanPrimaryCellAtlas (broad) ---")
  hpca_ref <- celldex::HumanPrimaryCellAtlasData()
  query_se <- GetAssayData(seu_merged, assay = "RNA", layer = "data")
  sr_hpca <- SingleR(
    test    = query_se,
    ref     = hpca_ref,
    labels  = hpca_ref$label.main,
    fine.tune = TRUE,
    num.threads = opt$ncores
  )
  seu_merged$singler_main <- sr_hpca$labels
} else {
  message("--- SingleR HPCA: skipped (refs_dir not active/empty) ---")
  seu_merged$singler_main <- NA_character_
}

# ---- 4. SingleR — van Galen 2019 (AML-specific) --------------------------

if (has_refs) {
  message("--- SingleR: van Galen 2019 (AML-specific) ---")
  vg_loaded <- tryCatch(build_van_galen_reference(opt$refs_dir),
                        error = function(e) { warning(e$message); NULL })
  if (!is.null(vg_loaded)) {
    query_se <- GetAssayData(seu_merged, assay = "RNA", layer = "data")
    sr_vg <- SingleR(
      test    = query_se,
      ref     = vg_loaded,
      labels  = vg_loaded$label,
      fine.tune = TRUE,
      num.threads = opt$ncores
    )
    seu_merged$singler_fine <- sr_vg$labels
  } else {
    seu_merged$singler_fine <- NA_character_
    message("  skipped — van Galen reference unavailable")
  }
} else {
  message("--- SingleR van Galen: skipped (refs_dir not active/empty) ---")
  seu_merged$singler_fine <- NA_character_
}

# ---- 5. Azimuth — Hao/Stuart human BM reference (optional) ---------------

if (has_refs) {
  message("--- Azimuth: human BM reference (optional) ---")
  azimuth_ok <- requireNamespace("Azimuth", quietly = TRUE) &&
                requireNamespace("SeuratData", quietly = TRUE)
  if (azimuth_ok) {
    bm_installed <- "bonemarrowref" %in% SeuratData::AvailableData()$Dataset
    if (!bm_installed) {
      message("  installing SeuratData::bonemarrowref (first run only)...")
      tryCatch(SeuratData::InstallData("bonemarrowref"),
               error = function(e) { warning(e$message); azimuth_ok <<- FALSE })
    }
  }
  if (azimuth_ok) {
    seu_azi <- Azimuth::RunAzimuth(seu_merged, reference = "bonemarrowref")
    seu_merged$azimuth_l1    <- seu_azi$predicted.celltype.l1
    seu_merged$azimuth_l2    <- seu_azi$predicted.celltype.l2
    seu_merged$azimuth_score <- seu_azi$predicted.celltype.l2.score
  } else {
    message("  Azimuth not available — skipping (install: SeuratData + Azimuth)")
    seu_merged$azimuth_l1    <- NA_character_
    seu_merged$azimuth_l2    <- NA_character_
    seu_merged$azimuth_score <- NA_real_
  }
} else {
  message("--- Azimuth BM: skipped (refs_dir not active/empty) ---")
  seu_merged$azimuth_l1    <- NA_character_
  seu_merged$azimuth_l2    <- NA_character_
  seu_merged$azimuth_score <- NA_real_
}

# ---- 6. UCell — score curated signatures ---------------------------------

message("--- UCell: scoring signatures from ", opt$sigs_yml, " ---")
sigs <- yaml::read_yaml(opt$sigs_yml)
sig_lists <- lapply(sigs, function(x) x$genes)
names(sig_lists) <- names(sigs)
# Drop any signatures whose gene list is empty
sig_lists <- sig_lists[lengths(sig_lists) > 0]
message("  scoring ", length(sig_lists), " signatures: ",
        paste(names(sig_lists), collapse = ", "))

seu_merged <- AddModuleScore_UCell(
  seu_merged,
  features = sig_lists,
  ncores   = opt$ncores,
  name     = ""   # column name = signature name as given
)
# UCell adds columns named "<sig>" — rename to UCell_<sig> for clarity
sig_cols_new <- paste0("UCell_", names(sig_lists))
existing     <- intersect(names(sig_lists), colnames(seu_merged@meta.data))
seu_merged@meta.data[, sig_cols_new] <- seu_merged@meta.data[, existing]
seu_merged@meta.data[, existing]     <- NULL

# ---- 6b. Classify cells based on van Galen S4A signature scores ----------
# Argmax over the three S4A-derived signatures (HSC/Prog, GMP, Myeloid), but
# operating on *z-scored* columns rather than raw UCell scores. UCell baselines
# differ between signatures (the Myeloid list contains more broadly-expressed
# genes, so its floor is higher across all cells); raw argmax therefore biases
# toward Myeloid-like. Per-column z-scoring puts each signature on the same
# scale so argmax compares relative position within the cohort.
#
# Cells whose top z-score < MIN_Z are left "Unassigned" — catches HS-5 stromal,
# endothelial, B/T cells and other populations that should not be force-called
# as one of the three myeloid classes.

message("--- Classifying cells with van Galen S4A signatures (z-scored argmax) ---")
vgalen_sigs    <- c("vgalen_hsc_s4a", "vgalen_gmp_s4a", "vgalen_myeloid_s4a")
class_labels   <- c("HSC-like", "GMP-like", "Myeloid-like")
ucell_sig_cols <- paste0("UCell_", vgalen_sigs)
MIN_Z          <- 0.5   # top z-score must be at least this far above the mean

if (all(ucell_sig_cols %in% colnames(seu_merged@meta.data))) {
  scores_raw <- as.matrix(seu_merged@meta.data[, ucell_sig_cols, drop = FALSE])
  # Per-column (per-signature) z-score across all cells in the cohort
  scores_z   <- scale(scores_raw)
  colnames(scores_z) <- class_labels

  max_idx   <- max.col(scores_z, ties.method = "first")
  row_idx   <- cbind(seq_len(nrow(scores_z)), max_idx)
  max_z     <- scores_z[row_idx]
  max_raw   <- scores_raw[row_idx]

  vgalen_class <- class_labels[max_idx]
  vgalen_class[max_z < MIN_Z] <- "Unassigned"

  seu_merged$vgalen_class       <- vgalen_class
  seu_merged$vgalen_class_z     <- as.numeric(max_z)
  seu_merged$vgalen_class_score <- as.numeric(max_raw)

  message("  Classification composition (MIN_Z = ", MIN_Z, "):")
  print(table(seu_merged$vgalen_class, useNA = "ifany"))
} else {
  warning("Not all vgalen S4A signatures are present in scored metadata — classification skipped.")
}

# ---- 7. Tables -----------------------------------------------------------

message("--- Writing tables ---")
if (has_refs) {
  write.table(
    data.frame(cell = colnames(seu_merged),
               singler_main = seu_merged$singler_main,
               singler_fine = seu_merged$singler_fine),
    file = file.path(tab_dir, "singler_main_calls.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(cell = colnames(seu_merged),
               azimuth_l1 = seu_merged$azimuth_l1,
               azimuth_l2 = seu_merged$azimuth_l2,
               azimuth_score = seu_merged$azimuth_score),
    file = file.path(tab_dir, "azimuth_calls.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# UCell scores by cluster
ucell_cols <- grep("^UCell_", colnames(seu_merged@meta.data), value = TRUE)
score_by_cluster <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  group_by(cell_type) |>
  summarise(across(all_of(ucell_cols), \(x) mean(x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop")
write.table(
  score_by_cluster,
  file = file.path(tab_dir, "ucell_scores_by_cluster.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if ("cell_type" %in% colnames(seu_merged@meta.data) && !all(is.na(seu_merged$singler_main))) {
  conf <- table(existing_cluster = seu_merged$cell_type,
                singler_main     = seu_merged$singler_main)
  write.table(conf, file = file.path(tab_dir, "cell_type_vs_singler_confusion.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
}

if ("vgalen_class" %in% colnames(seu_merged@meta.data)) {
  # Write cell-level calls
  write.table(
    data.frame(cell = colnames(seu_merged),
               vgalen_class = seu_merged$vgalen_class),
    file = file.path(tab_dir, "vgalen_signature_calls.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Confusion table vs manual cell types
  if ("cell_type" %in% colnames(seu_merged@meta.data)) {
    conf_vgalen <- table(existing_cluster = seu_merged$cell_type,
                         vgalen_class     = seu_merged$vgalen_class)
    write.table(conf_vgalen, file = file.path(tab_dir, "cell_type_vs_vgalen_class_confusion.tsv"),
                sep = "\t", quote = FALSE, col.names = NA)
  }
}

# ---- 8. Figures ----------------------------------------------------------

message("--- Plotting ---")
if (!all(is.na(seu_merged$singler_main))) {
  pdf(file.path(fig_dir, "03a_singler_umap.pdf"), width = 14, height = 6)
  print(DimPlot(seu_merged, group.by = c("cell_type", "singler_main"),
                label = TRUE, repel = TRUE) & NoLegend())
  dev.off()
}

if (has_refs && !all(is.na(seu_merged$azimuth_l2))) {
  pdf(file.path(fig_dir, "03b_azimuth_umap.pdf"), width = 14, height = 6)
  print(DimPlot(seu_merged, group.by = c("cell_type", "azimuth_l2"),
                label = TRUE, repel = TRUE) & NoLegend())
  dev.off()
}

if ("cell_type" %in% colnames(seu_merged@meta.data) && !all(is.na(seu_merged$singler_main))) {
  conf_pct <- prop.table(table(seu_merged$cell_type, seu_merged$singler_main), margin = 1)
  pdf(file.path(fig_dir, "03c_confusion_existing_vs_singler.pdf"),
      width = max(8, ncol(conf_pct) * 0.35), height = max(6, nrow(conf_pct) * 0.35))
  pheatmap(as.matrix(conf_pct), cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = FALSE, color = colorRampPalette(c("white", "#0B3D5C"))(50),
           main = "Existing cluster (rows) vs SingleR main label (cols)")
  dev.off()
}

if ("vgalen_class" %in% colnames(seu_merged@meta.data)) {
  # UMAP by vgalen_class
  pdf(file.path(fig_dir, "03e_vgalen_class_umap.pdf"), width = 14, height = 6)
  print(DimPlot(seu_merged, group.by = c("cell_type", "vgalen_class"),
                label = TRUE, repel = TRUE) & NoLegend())
  dev.off()
  
  # Confusion heatmap vs existing cell types
  if ("cell_type" %in% colnames(seu_merged@meta.data)) {
    conf_v_pct <- prop.table(table(seu_merged$cell_type, seu_merged$vgalen_class), margin = 1)
    pdf(file.path(fig_dir, "03f_confusion_existing_vs_vgalen_class.pdf"),
        width = max(8, ncol(conf_v_pct) * 0.5), height = max(6, nrow(conf_v_pct) * 0.35))
    pheatmap(as.matrix(conf_v_pct), cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = FALSE, color = colorRampPalette(c("white", "#0D5F3A"))(50),
             main = "Existing cluster (rows) vs van Galen signature class (cols)")
    dev.off()
  }
}

# ---- 8b. UCell signature heatmap (publication-quality) -------------------
# Signatures are grouped into biological themes (column annotation), cluster
# rows are grouped by cell-type category (row annotation), columns are kept in
# the curated thematic order (no column clustering) so the legend reads left-
# to-right by theme, rows are clustered within groups to reveal which clusters
# share signature profiles.

# Curated signature order — grouped by theme. Only signatures present in the
# scored Seurat object will be plotted; stale UCell columns from previous runs
# (e.g. the deprecated vgalen_*_mmc4 set) are dropped here.
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

# Cluster (row) grouping — broad biological category. Anything unmapped falls
# into "Other" (rare; usually means cluster labels have changed).
cluster_category <- function(x) {
  dplyr::case_when(
    x %in% c("HSPC / LSC-like blasts", "ABCB5+ resistant primitive AML")  ~ "AML primitive",
    x %in% c("GMP-like AML blasts")                                       ~ "AML progenitor",
    x %in% c("Mono-like AML", "Promono-like / quiescent myeloid",
             "Stress-response myeloid", "CD16+ Mono / macrophage")        ~ "AML myeloid",
    x %in% c("Unclassified niche-related", "Unclassified (low-quality?)") ~ "AML ambiguous",
    x %in% c("HS5 stromal / MSC", "BM endothelial")                       ~ "Stromal",
    x %in% c("B cells", "T cells")                                        ~ "Normal lymphoid",
    x %in% c("Eosinophils / basophils")                                   ~ "Normal myeloid",
    TRUE                                                                  ~ "Other"
  )
}

# Build the matrix — cluster means, z-scored per signature column across clusters
heat_cols <- paste0("UCell_", sig_in_data)
heat_mat  <- as.matrix(score_by_cluster[, heat_cols, drop = FALSE])
rownames(heat_mat) <- score_by_cluster$cell_type
colnames(heat_mat) <- sig_in_data           # strip UCell_ prefix for display
heat_z <- scale(heat_mat)

# Annotations — column = signature theme, row = cluster category
ann_col <- data.frame(Theme = factor(sig_theme[sig_in_data],
                                     levels = names(sig_order)),
                      row.names = sig_in_data)
ann_row <- data.frame(Category = factor(cluster_category(rownames(heat_z)),
                                        levels = c("AML primitive", "AML progenitor",
                                                   "AML myeloid",    "AML ambiguous",
                                                   "Stromal", "Normal lymphoid",
                                                   "Normal myeloid", "Other")),
                      row.names = rownames(heat_z))

# Bruna-flavoured palette for annotations
theme_colours <- c(
  "vG class (S4A)"     = "#7B2CBF",   # purple
  "Stemness / LSC"     = "#D7263D",   # ADRN red
  "Quiescence"         = "#2E8B57",   # sea-green
  "Stress / ISR"       = "#E66B00",   # orange
  "Niche signalling"   = "#1B98E0",   # MES blue
  "Resistance"         = "#8B4513",   # saddle brown
  "Apoptosis"          = "#5A5A5A"    # neutral grey
)
category_colours <- c(
  "AML primitive"   = "#7A0F1F",
  "AML progenitor"  = "#D7263D",
  "AML myeloid"     = "#F4A261",
  "AML ambiguous"   = "#C9A86A",
  "Stromal"         = "#708090",
  "Normal lymphoid" = "#1B98E0",
  "Normal myeloid"  = "#4682B4",
  "Other"           = "#BBBBBB"
)
ann_colours <- list(Theme = theme_colours, Category = category_colours)

# Symmetric colour scale so 0 sits at white
zmax <- max(2, ceiling(max(abs(heat_z), na.rm = TRUE)))
breaks <- seq(-zmax, zmax, length.out = 51)
heat_palette <- colorRampPalette(c("#1B4F72", "#1B98E0", "white",
                                   "#D7263D", "#6B0F1A"))(50)

# Gaps between column themes for visual separation
n_per_theme <- table(ann_col$Theme)[unique(ann_col$Theme)]
col_gaps    <- head(cumsum(n_per_theme), -1)

n_rows <- nrow(heat_z); n_cols <- ncol(heat_z)
pdf(file.path(fig_dir, "03d_ucell_heatmap.pdf"),
    width  = max(11, 3.0 + n_cols * 0.42),
    height = max(7,  2.5 + n_rows * 0.42))
pheatmap(
  heat_z,
  cluster_rows       = TRUE,
  cluster_cols       = FALSE,
  gaps_col           = as.integer(col_gaps),
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
  treeheight_row     = 30,
  treeheight_col     = 0,
  main               = "UCell signature scores per cluster (z-scored across clusters)"
)
dev.off()

# ---- 9. Save updated RDS + session info -----------------------------------

message("--- Saving updated RDS ---")
saveRDS(seu_merged, file = opt$rds_in)
message("  updated: ", opt$rds_in)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_03_sessionInfo.txt"))

message("Done!")
