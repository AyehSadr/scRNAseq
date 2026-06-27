#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 3c: Same-Clone DGE
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),   type="character", help="Path to input QC'd Cellecta RDS file")
parser <- add_option(parser, c("--out_dir"),  type="character", help="Output directory")
opt <- parse_args(parser)

dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- file.path(opt$out_dir, "figures_clonetracker")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

rds_base <- basename(opt$rds_in)
patient_prefix <- sub("_seurat_cellecta_qc\\.RDS$", "", rds_base)
if (patient_prefix == rds_base) patient_prefix <- "AML_Patient"

message("=== AML CloneTracker Analysis 3c: Same-Clone DGE ===")
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)
seu_merged <- readRDS(opt$rds_in)

# 1. Filter to confident clones
seu_conf <- subset(seu_merged, subset = clone_confident == TRUE)
if (ncol(seu_conf) < 10) stop("Not enough confident clones for DGE.")

# Ensure we use SCT and proper normalization preparation
DefaultAssay(seu_conf) <- "SCT"

# --- Helper function to run DGE on shared clones ---
run_shared_clone_dge <- function(seu_obj, niche_name, condition_var, cond1, cond2) {
  message("\n--- Analyzing Niche: ", niche_name, " ---")
  
  # Identify clones present in BOTH conditions
  meta <- seu_obj@meta.data
  clones_cond1 <- unique(meta$clone_barcode[meta[[condition_var]] == cond1])
  clones_cond2 <- unique(meta$clone_barcode[meta[[condition_var]] == cond2])
  shared_clones <- intersect(clones_cond1, clones_cond2)
  
  if (length(shared_clones) < 1) {
    message("  No shared clones between ", cond1, " and ", cond2, " in ", niche_name)
    return(NULL)
  }
  
  message("  Found ", length(shared_clones), " shared clones. Subsetting...")
  seu_sub <- subset(seu_obj, subset = clone_barcode %in% shared_clones)
  
  n_cond1 <- sum(seu_sub@meta.data[[condition_var]] == cond1)
  n_cond2 <- sum(seu_sub@meta.data[[condition_var]] == cond2)
  message(sprintf("  Cells per condition: %s=%d, %s=%d", cond1, n_cond1, cond2, n_cond2))
  
  # Must re-prep SCT after subsetting in v5, wrap in tryCatch for low cell counts
  seu_sub <- tryCatch({
    PrepSCTFindMarkers(seu_sub)
  }, error = function(e) {
    warning("  PrepSCTFindMarkers failed (likely too few cells). Proceeding without re-prepping.")
    return(seu_sub)
  })
  
  Idents(seu_sub) <- condition_var
  
  message("  Running FindMarkers (", cond1, " vs ", cond2, ")...")
  res <- tryCatch({
    FindMarkers(seu_sub, ident.1 = cond1, ident.2 = cond2, test.use = "wilcox", logfc.threshold = 0, verbose = FALSE)
  }, error = function(e) {
    warning("  FindMarkers failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res$gene <- rownames(res)
  
  out_csv <- file.path(opt$out_dir, paste0(patient_prefix, "_DGE_SharedClones_", niche_name, "_AraC_", cond1, "_vs_", cond2, ".csv"))
  write.csv(res, out_csv, row.names=FALSE)
  
  if (n_cond1 < 10 || n_cond2 < 10) {
    warning(sprintf("Underpowered: %d vs %d cells. DGE results are not interpretable. ", n_cond1, n_cond2),
            "Skipping volcano output to avoid misleading visualization.")
    return(res)
  }
  
  # Volcano Plot
  res$diffexpressed <- "Not Significant"
  res$diffexpressed[res$avg_log2FC > 0.5 & res$p_val_adj < 0.05] <- "Up-regulated"
  res$diffexpressed[res$avg_log2FC < -0.5 & res$p_val_adj < 0.05] <- "Down-regulated"
  res$diffexpressed <- factor(res$diffexpressed, levels=c("Up-regulated", "Down-regulated", "Not Significant"))
  
  p <- ggplot(res, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=diffexpressed)) +
          geom_point(alpha=0.6, size=1) +
          scale_colour_manual(values=c("Up-regulated"="firebrick", "Down-regulated"="steelblue", "Not Significant"="grey80")) +
          geom_vline(xintercept=c(-0.5, 0.5), linetype="dashed", color="grey40") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
          labs(title=paste0("Same-Clone DGE: ", niche_name, " (AraC ", cond1, " vs ", cond2, ")"), x="log2 Fold Change", y="-log10 adj.p-value", colour="Status") +
          theme_classic()
          
  out_pdf <- file.path(fig_dir, paste0("11_Volcano_SharedClones_", niche_name, ".pdf"))
  pdf(out_pdf, width=8, height=7)
  print(p)
  dev.off()
  
  message("  Done.")
  return(res)
}

# 2. AraC effect on shared clones in NO STROMA
seu_nostroma <- subset(seu_conf, subset = stroma == FALSE)
run_shared_clone_dge(seu_nostroma, "NoStroma", "araC", TRUE, FALSE)

# 3. AraC effect on shared clones WITH STROMA
seu_stroma <- subset(seu_conf, subset = stroma == TRUE)
run_shared_clone_dge(seu_stroma, "WithStroma", "araC", TRUE, FALSE)

message("Same-Clone DGE Complete!")
