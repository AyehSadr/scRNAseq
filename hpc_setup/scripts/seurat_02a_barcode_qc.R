#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 3a: Barcode QC & Mapping
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),   type="character", help="Path to input general Seurat RDS file")
parser <- add_option(parser, c("--ct_dir"),   type="character", help="Path to Cellecta barcode TSV directory (*.bam.tsv)")
parser <- add_option(parser, c("--out_dir"),  type="character", help="Output directory for RDS and figures")
parser <- add_option(parser, c("--allow_low_recovery"), action="store_true", default=FALSE, help="Bypass the strict 1% barcode recovery threshold (for exploratory analysis)")
opt <- parse_args(parser)

dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- file.path(opt$out_dir, "figures_clonetracker")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# Extract patient prefix from RDS (e.g., AML_SRAML10_seurat_general.RDS -> AML_SRAML10)
rds_base <- basename(opt$rds_in)
patient_prefix <- sub("_seurat_general\\.RDS$", "", rds_base)
if (patient_prefix == rds_base) patient_prefix <- "AML_Patient"

message("=== AML CloneTracker Analysis 3a: Barcode QC ===")
message("Input RDS : ", opt$rds_in)
message("Prefix    : ", patient_prefix)

# 1. Load Seurat Object
if (!file.exists(opt$rds_in)) stop("Input RDS not found")
seu_merged <- readRDS(opt$rds_in)

# 2. Load Cellecta Barcodes
ct_files <- list.files(opt$ct_dir, pattern="\\.bam\\.tsv$", full.names=TRUE)
if (length(ct_files) == 0) stop("No .bam.tsv files found in ", opt$ct_dir)

ct_list <- lapply(ct_files, function(f) {
  sid <- gsub("\\.bam\\.tsv$", "", basename(f))
  df  <- read.delim(f, stringsAsFactors=FALSE)
  if (!all(c("cell_barcode", "clone_barcode", "n_reads") %in% colnames(df))) {
    stop("Missing required columns in ", f)
  }
  df$sample_id <- sid
  df$cell_id   <- paste0(sid, "_", df$cell_barcode)
  df
})
ct_df <- dplyr::bind_rows(ct_list)

# 3. Sanity Check: Barcode length
blens <- table(nchar(ct_df$clone_barcode))
message("\n--- Barcode Lengths (Expected ~44bp for BC14+linker+BC30) ---")
print(blens)

# 4. FIX THE -1 SUFFIX BUG
sample_rn <- head(rownames(seu_merged@meta.data), 3)
sample_id <- head(ct_df$cell_id, 3)
message("\n--- Sanity Check: Barcode Join Suffix ---")
message("Sample Seurat rownames: ", paste(sample_rn, collapse=", "))
message("Sample CT cell_ids:     ", paste(sample_id, collapse=", "))

seurat_has_dash1 <- any(grepl("-1$", rownames(seu_merged@meta.data)))
ct_has_dash1     <- any(grepl("-1$", ct_df$cell_id))

if (seurat_has_dash1 && !ct_has_dash1) {
  message("  [FIX]: Adjusting CT cell_ids to add '-1' suffix to match Seurat rownames")
  ct_df$cell_id <- paste0(ct_df$cell_id, "-1")
}

# 5. Process Barcodes (Purity & Doublet detection)
message("\n--- Processing Barcodes ---")
ct_processed <- ct_df %>%
  arrange(cell_id, desc(n_reads)) %>%
  group_by(cell_id) %>%
  summarise(
    clone_barcode = first(clone_barcode),
    n_reads       = first(n_reads),
    n_reads_2nd   = nth(n_reads, 2, default = 0),
    n_distinct_bc = n_distinct(clone_barcode)
  ) %>%
  mutate(
    barcode_purity = n_reads / (n_reads + n_reads_2nd),
    likely_doublet = n_distinct_bc > 1 & (n_reads_2nd / n_reads > 0.2)
  )

# 6. Inject into Seurat
meta <- seu_merged@meta.data
meta$cell_id <- rownames(meta)
meta <- left_join(meta, ct_processed, by="cell_id")
rownames(meta) <- meta$cell_id
seu_merged <- AddMetaData(seu_merged, metadata=meta[, c("clone_barcode", "n_reads", "n_reads_2nd", "n_distinct_bc", "barcode_purity", "likely_doublet")])

seu_merged$has_clone <- !is.na(seu_merged$clone_barcode) & seu_merged$n_reads >= 1
seu_merged$clone_confident <- !is.na(seu_merged$clone_barcode) & seu_merged$n_reads >= 2

# 7. Recovery Reporting
n_cells    <- ncol(seu_merged)
n_attached <- sum(seu_merged$clone_confident)
pct        <- round(100 * n_attached / n_cells, 2)

message("\n--- Barcode Recovery Summary ---")
message("Total Cells: ", n_cells, " | Confident Clones: ", n_attached, " (", pct, "%)")

per_sample <- seu_merged@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    n_cells = n(),
    n_any_clone = sum(has_clone, na.rm=TRUE),
    n_confident = sum(clone_confident, na.rm=TRUE),
    pct_confident = round(100 * n_confident / n_cells, 2),
    n_unique_clones = length(unique(na.omit(clone_barcode[clone_confident]))),
    n_likely_doublet = sum(likely_doublet, na.rm=TRUE)
  )
print(per_sample)
write.csv(per_sample, file.path(opt$out_dir, paste0(patient_prefix, "_barcode_recovery_per_sample.csv")), row.names=FALSE)

if (pct < 1) {
  msg <- paste0(
          "Barcode recovery is < 1% (", pct, "%).\n",
          "Common causes:\n",
          " 1. Cell Ranger '-1' suffix mismatch — check Seurat rownames vs cell_barcode in TSVs\n",
          " 2. sample_id naming mismatch between TSV filenames and Seurat sample_id metadata\n",
          " 3. Genuinely no enrichment library and barcodes were not captured in GEX\n"
  )
  if (!opt$allow_low_recovery) {
    stop(msg, "\nERROR: Aborting downstream analysis to prevent misinterpretation of underpowered data. Run with --allow_low_recovery to bypass.")
  } else {
    warning(msg, "\nWARNING: --allow_low_recovery flag is active. Proceeding with exploratory analysis on <1% data.")
  }
}

# 8. QC Plots
p_overlay <- DimPlot(seu_merged, cells.highlight=WhichCells(seu_merged, expression=clone_confident==TRUE), cols.highlight="firebrick", cols="lightgrey", reduction="umap") + ggtitle("Confident Clones (n_reads >= 2)")
pdf(file.path(fig_dir, "08a_CloneTracker_Confident_Overlay.pdf"), width=8, height=7)
print(p_overlay)
dev.off()

# Barcode aware cluster table
bc_cluster <- table(seu_merged$seurat_clusters, seu_merged$clone_confident)
colnames(bc_cluster) <- c("No Clone", "Confident Clone")
write.csv(as.data.frame.matrix(bc_cluster), file.path(opt$out_dir, paste0(patient_prefix, "_barcode_by_cluster.csv")))

# UMI Distribution
p_umi <- ggplot(seu_merged@meta.data %>% filter(has_clone), aes(x=sample_id, y=log10(n_reads))) + 
  geom_violin(fill="lightblue") + 
  geom_jitter(height=0, width=0.2, alpha=0.1, size=0.5) + 
  geom_hline(yintercept=log10(2), linetype="dashed", color="red") + # Threshold line
  theme_classic() + 
  labs(title="Clone UMI Distribution per Sample", y="log10(n_reads)") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
pdf(file.path(fig_dir, "08b_CloneTracker_UMI_Distribution.pdf"), width=8, height=6)
print(p_umi)
dev.off()

# 9. Save
out_rds <- file.path(opt$out_dir, paste0(patient_prefix, "_seurat_cellecta_qc.RDS"))
saveRDS(seu_merged, file=out_rds)
message("\nSaved QC RDS: ", out_rds)
