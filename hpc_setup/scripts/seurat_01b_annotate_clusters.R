#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 1b: Annotate Clusters
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),   type="character", help="Path to input general Seurat RDS file")
parser <- add_option(parser, c("--out_dir"),  type="character", help="Output directory for RDS and figures")
opt <- parse_args(parser)

dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- file.path(opt$out_dir, "figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

message("=== AML CloneTracker Analysis 1b: Annotate Clusters ===")
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

message("Loading RDS... (this may take a minute)")
seu_merged <- readRDS(opt$rds_in)

# 1. Apply Cluster Annotations
message("--- Annotating cell types ---")
cluster_annotations <- c(
  "0"  = "GMP-like AML blasts",
  "1"  = "HSPC / LSC-like blasts",
  "2"  = "Promono-like / quiescent myeloid",
  "3"  = "Unclassified (low-quality?)",
  "4"  = "Mono-like AML",
  "5"  = "CD16+ Mono / macrophage",
  "6"  = "Stress-response myeloid",
  "7"  = "ABCB5+ resistant primitive AML",
  "8"  = "Unclassified niche-related",
  "9"  = "B cells",
  "10" = "BM endothelial",
  "11" = "T cells",
  "12" = "HS5 stromal / MSC",
  "13" = "Eosinophils / basophils"
)

seu_merged$cell_type <- unname(cluster_annotations[as.character(seu_merged$seurat_clusters)])
seu_merged$cell_type <- factor(seu_merged$cell_type, levels = unique(unname(cluster_annotations)))

# 2. Sanity-check before publishing these labels
message("--- Saving Cell Type composition ---")
cell_counts <- table(seu_merged$cell_type, seu_merged$sample_id)
write.csv(as.data.frame.matrix(cell_counts), file.path(opt$out_dir, "CellType_Counts_per_Sample.csv"))
print(cell_counts)

# 3. Re-plot UMAPs with Cell Types
message("--- Generating new annotated UMAPs ---")

p_umap_type <- DimPlot(seu_merged, group.by = "cell_type", label = TRUE, repel = TRUE) + 
  NoLegend() + 
  ggtitle("UMAP by Cell Type")

pdf(file.path(fig_dir, "04b_UMAP_CellTypes.pdf"), width=10, height=8)
print(p_umap_type)
dev.off()

# Plot split by condition
p_umap_split <- DimPlot(seu_merged, group.by = "cell_type", split.by = "sample_id", ncol=2, label=FALSE) + 
  ggtitle("Cell Types split by Sample") +
  theme(legend.position = "bottom")

pdf(file.path(fig_dir, "05c_UMAP_CellTypes_by_Sample.pdf"), width=12, height=10)
print(p_umap_split)
dev.off()

# 3b. Stacked Bar Plots (Proportions)
message("--- Generating stacked bar plots ---")
meta <- seu_merged@meta.data

# Proportion of cell types per sample
p_prop_sample <- ggplot(meta, aes(x=sample_id, fill=cell_type)) +
  geom_bar(position="fill", color="black", linewidth=0.2) +
  theme_classic() +
  labs(x="Sample", y="Proportion", title="Cell Type Composition per Sample", fill="Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "04c_CellType_Proportions_per_Sample.pdf"), width=8, height=6)
print(p_prop_sample)
dev.off()

# Proportion of samples per cell type
p_prop_celltype <- ggplot(meta, aes(x=cell_type, fill=sample_id)) +
  geom_bar(position="fill", color="black", linewidth=0.2) +
  theme_classic() +
  labs(x="Cell Type", y="Proportion", title="Sample Composition per Cell Type", fill="Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "04d_Sample_Proportions_per_CellType.pdf"), width=10, height=6)
print(p_prop_celltype)
dev.off()

# 4. Save updated RDS
message("--- Saving updated RDS object ---")
saveRDS(seu_merged, file=opt$rds_in)
message("Updated general RDS: ", opt$rds_in)
message("Done!")
