#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Seurat v5 scRNA-seq General Analysis
# Samples: 4 x SRAML10 (Patient 2), GEM-X OCM demultiplexed
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(optparse)
  library(BiocParallel)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# Increase memory limit for parallel workers (needed for large Seurat objects)
options(future.globals.maxSize = 10 * 1024^3) # 10 GB

# --- Options ------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, c("--cr_dir"),   type="character",
                     help="Path to Cell Ranger per_sample_outs directory")
parser <- add_option(parser, c("--out_dir"),  type="character",
                     help="Output directory for RDS and figures")
parser <- add_option(parser, c("--ncores"),   type="integer", default=8,
                     help="Number of cores [default: 8]")
parser <- add_option(parser, c("--mt_max"),   type="double",  default=20,
                     help="Max % mitochondrial reads per cell [default: 20]")
parser <- add_option(parser, c("--min_feat"), type="integer", default=200,
                     help="Min genes per cell [default: 200]")
parser <- add_option(parser, c("--max_feat"), type="integer", default=6000,
                     help="Max genes per cell (doublet proxy) [default: 6000]")
parser <- add_option(parser, c("--sample_sheet"), type="character", default=NULL,
                     help="Path to sample sheet TSV [default: NULL]")
opt <- parse_args(parser)

dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$out_dir, "figures"), recursive=TRUE, showWarnings=FALSE)

message("=== AML CloneTracker Seurat Pipeline (General) ===")
message("Cell Ranger dir : ", opt$cr_dir)
message("Output dir      : ", opt$out_dir)
message("Cores           : ", opt$ncores)

register(MulticoreParam(workers=opt$ncores))

# --- Sample metadata ----------------------------------------------------------
if (!is.null(opt$sample_sheet) && file.exists(opt$sample_sheet)) {
  message("Loading sample metadata dynamically from: ", opt$sample_sheet)
  sample_meta <- read.delim(opt$sample_sheet, stringsAsFactors=FALSE, sep="\t")
  
  # Map columns to script expected names if necessary
  if (!"stroma" %in% colnames(sample_meta) && "stroma_HS5" %in% colnames(sample_meta)) {
    sample_meta$stroma <- sample_meta$stroma_HS5 == "+" | tolower(sample_meta$stroma_HS5) == "true"
  }
  if (!"araC" %in% colnames(sample_meta) && "chemo_AraC" %in% colnames(sample_meta)) {
    sample_meta$araC <- sample_meta$chemo_AraC == "+" | tolower(sample_meta$chemo_AraC) == "true"
  }
  
  # Ensure necessary columns are present
  required_cols <- c("sample_id", "patient", "condition")
  missing_cols <- setdiff(required_cols, colnames(sample_meta))
  if (length(missing_cols) > 0) {
    stop("Sample sheet must contain at least: ", paste(required_cols, collapse=", "))
  }
} else {
  message("No sample sheet provided or file not found. Falling back to default SRAML10 metadata.")
  sample_meta <- data.frame(
    sample_id = c("P2_T11_noStroma_noAraC",
                  "P2_T12_noStroma_AraC",
                  "P2_T15_Stroma_noAraC",
                  "P2_T16_Stroma_AraC"),
    patient   = "SRAML10",
    tube      = c(11, 12, 15, 16),
    stroma    = c(FALSE, FALSE, TRUE,  TRUE),
    araC      = c(FALSE, TRUE,  FALSE, TRUE),
    condition = c("Untreated_NoStromal", "Treated_NoStromal",
                  "Untreated_Stromal",   "Treated_Stromal"),
    stringsAsFactors=FALSE
  )
}

# Define consistent color palette for samples dynamically
unique_sids <- unique(sample_meta$sample_id)
color_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
sample_colors <- setNames(color_palette[1:length(unique_sids)], unique_sids)

# =============================================================================
# STEP 1: Load Cell Ranger filtered matrices
# =============================================================================
message("\n--- Step 1: Loading Cell Ranger matrices ---")

seu_list <- list()
for (i in seq_len(nrow(sample_meta))) {
  sid  <- sample_meta$sample_id[i]
  mdir <- file.path(opt$cr_dir, sid, "sample_filtered_feature_bc_matrix")
  if (!dir.exists(mdir)) {
    stop("Matrix directory not found: ", mdir)
  }
  mat  <- Read10X(data.dir=mdir)
  seu  <- CreateSeuratObject(counts=mat,
                             project=sid,
                             min.cells=3,
                             min.features=opt$min_feat)
  # attach sample metadata
  seu$sample_id <- sid
  seu$patient   <- sample_meta$patient[i]
  seu$tube      <- sample_meta$tube[i]
  seu$stroma    <- sample_meta$stroma[i]
  seu$araC      <- sample_meta$araC[i]
  seu$condition <- sample_meta$condition[i]

  message("  Loaded ", sid, " : ", ncol(seu), " cells / ", nrow(seu), " genes")
  seu_list[[sid]] <- seu
}

# =============================================================================
# STEP 2: Per-sample QC
# =============================================================================
message("\n--- Step 2: QC filtering ---")

# Add MT% to metadata before filtering
seu_list <- lapply(seu_list, function(seu) {
  seu[["pct_mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  seu
})

# Generate combined QC violin plot before filtering
seu_unfiltered <- merge(seu_list[[1]], y=seu_list[-1], add.cell.ids=names(seu_list))
seu_unfiltered$sample_id <- factor(seu_unfiltered$sample_id, levels=names(sample_colors))

p_qc <- VlnPlot(seu_unfiltered,
                features=c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                group.by="sample_id",
                pt.size=0, ncol=3, cols=sample_colors) +
        plot_annotation(title="Pre-filtering QC Metrics")

pdf(file.path(opt$out_dir, "figures", "01_QC_violins.pdf"), width=15, height=6)
print(p_qc)
dev.off()

rm(seu_unfiltered)
gc()

# Now filter the cells
seu_list <- lapply(seu_list, function(seu) {
  subset(seu,
         subset = nFeature_RNA >= opt$min_feat &
                  nFeature_RNA <= opt$max_feat &
                  pct_mt       <  opt$mt_max)
})

# =============================================================================
# STEP 3: Doublet detection with scDblFinder (per sample)
# =============================================================================
message("\n--- Step 3: Doublet detection (scDblFinder) ---")

seu_list <- lapply(names(seu_list), function(sid) {
  seu <- seu_list[[sid]]
  message("  Running scDblFinder on: ", sid, " (", ncol(seu), " cells)")
  result <- tryCatch({
    sce <- as.SingleCellExperiment(seu)
    set.seed(42)
    sce <- scDblFinder(sce, BPPARAM=SerialParam())   # SerialParam avoids multicore issues
    seu$scDblFinder_class <- sce$scDblFinder.class
    seu$scDblFinder_score <- sce$scDblFinder.score
    n_dbl <- sum(seu$scDblFinder_class == "doublet")
    message("  ", sid, ": removed ", n_dbl, " doublets")
    subset(seu, subset=scDblFinder_class == "singlet")
  }, error = function(e) {
    message("  WARNING: scDblFinder failed for ", sid, " — keeping all cells. Error: ", e$message)
    seu$scDblFinder_class <- "singlet"
    seu$scDblFinder_score <- NA
    seu
  })
  result
})
names(seu_list) <- sample_meta$sample_id
message("  Doublet detection complete.")

# =============================================================================
# STEP 4: SCTransform v2 normalization (per sample)
# =============================================================================
message("\n--- Step 4: SCTransform v2 normalization ---")

seu_list <- lapply(seu_list, function(seu) {
  SCTransform(seu,
              vst.flavor  = "v2",
              vars.to.regress = "pct_mt",
              verbose     = FALSE,
              return.only.var.genes = FALSE)
})

# =============================================================================
# STEP 5: Merge + PCA
# =============================================================================
message("\n--- Step 5: Merging and PCA (NO HARMONY) ---")

seu_merged <- merge(seu_list[[1]],
                    y    = seu_list[-1],
                    add.cell.ids = names(seu_list),
                    merge.data   = TRUE)

# Seurat v5: join RNA layers after merge (SCT assay handles this automatically)
DefaultAssay(seu_merged) <- "RNA"
seu_merged <- JoinLayers(seu_merged)
DefaultAssay(seu_merged) <- "SCT"
message("  Layers joined. Cells: ", ncol(seu_merged))

# Select HVGs across samples (correct Seurat v5 function name)
DefaultAssay(seu_merged) <- "SCT"
VariableFeatures(seu_merged) <- SelectIntegrationFeatures(seu_list, nfeatures=3000)
message("  Selected ", length(VariableFeatures(seu_merged)), " variable features")

# Using standard PCA, skipping Harmony based on review of confounding
seu_merged <- RunPCA(seu_merged, npcs=50, verbose=FALSE)
message("  PCA complete")

# =============================================================================
# STEP 6: UMAP + Clustering
# =============================================================================
message("\n--- Step 6: UMAP and clustering ---")

seu_merged <- RunUMAP(seu_merged,
                      reduction = "pca",
                      dims      = 1:30,
                      verbose   = FALSE)

seu_merged <- FindNeighbors(seu_merged,
                            reduction = "pca",
                            dims      = 1:30,
                            verbose   = FALSE)

# Run clustering at multiple resolutions for clustree
for (res in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
  seu_merged <- FindClusters(seu_merged,
                             resolution = res,
                             verbose    = FALSE)
}
# Default clustering for downstream use
seu_merged <- FindClusters(seu_merged, resolution=0.4, verbose=FALSE)
Idents(seu_merged) <- "SCT_snn_res.0.4"

# =============================================================================
# STEP 7: UMAP visualizations
# =============================================================================
message("\n--- Step 7: Generating UMAP plots ---")

p_cluster   <- DimPlot(seu_merged, reduction="umap", label=TRUE)  + ggtitle("Clusters (res 0.4)")
p_sample    <- DimPlot(seu_merged, reduction="umap", group.by="sample_id", cols=sample_colors) + ggtitle("Sample")
p_condition <- DimPlot(seu_merged, reduction="umap", group.by="condition")  + ggtitle("Condition")
p_stroma    <- DimPlot(seu_merged, reduction="umap", group.by="stroma")     + ggtitle("Stroma")
p_araC      <- DimPlot(seu_merged, reduction="umap", group.by="araC")       + ggtitle("AraC")

pdf(file.path(opt$out_dir, "figures", "02a_UMAP_Clusters.pdf"), width=8, height=7)
print(p_cluster)
dev.off()

pdf(file.path(opt$out_dir, "figures", "02b_UMAP_Sample.pdf"), width=8, height=7)
print(p_sample)
dev.off()

pdf(file.path(opt$out_dir, "figures", "02c_UMAP_Condition.pdf"), width=8, height=7)
print(p_condition)
dev.off()

pdf(file.path(opt$out_dir, "figures", "02d_UMAP_Stroma.pdf"), width=8, height=7)
print(p_stroma)
dev.off()

pdf(file.path(opt$out_dir, "figures", "02e_UMAP_AraC.pdf"), width=8, height=7)
print(p_araC)
dev.off()

# Cell count per sample per cluster
cluster_counts <- seu_merged@meta.data %>%
  group_by(sample_id, seurat_clusters) %>%
  tally() %>%
  tidyr::pivot_wider(names_from=seurat_clusters, values_from=n, values_fill=0)
write.csv(cluster_counts,
          file.path(opt$out_dir, "cluster_cell_counts.csv"),
          row.names=FALSE)

# =============================================================================
# STEP 8: Differential gene expression (condition comparisons)
# =============================================================================
message("\n--- Step 8: Differential expression ---")

DefaultAssay(seu_merged) <- "SCT"
seu_merged <- PrepSCTFindMarkers(seu_merged)

# 1. AraC effect (within no-stroma samples) -> T12 vs T11
seu_nostroma <- subset(seu_merged, subset=stroma == FALSE)
seu_nostroma <- PrepSCTFindMarkers(seu_nostroma)
Idents(seu_nostroma) <- "araC"
de_araC_nostroma <- FindMarkers(seu_nostroma, ident.1 = TRUE, ident.2 = FALSE, test.use = "wilcox", logfc.threshold = 0, verbose = FALSE)
write.csv(de_araC_nostroma, file.path(opt$out_dir, "DE_AraC_effect_NoStroma.csv"))

# 2. AraC effect (within stroma samples) -> T16 vs T15
seu_stroma <- subset(seu_merged, subset=stroma == TRUE)
seu_stroma <- PrepSCTFindMarkers(seu_stroma)
Idents(seu_stroma) <- "araC"
de_araC_stroma <- FindMarkers(seu_stroma, ident.1 = TRUE, ident.2 = FALSE, test.use = "wilcox", logfc.threshold = 0, verbose = FALSE)
write.csv(de_araC_stroma, file.path(opt$out_dir, "DE_AraC_effect_WithStroma.csv"))

# 3. Stroma effect (within untreated samples) -> T15 vs T11
seu_notreated <- subset(seu_merged, subset=araC == FALSE)
seu_notreated <- PrepSCTFindMarkers(seu_notreated)
Idents(seu_notreated) <- "stroma"
de_stroma_noarac <- FindMarkers(seu_notreated, ident.1 = TRUE, ident.2 = FALSE, test.use = "wilcox", logfc.threshold = 0, verbose = FALSE)
write.csv(de_stroma_noarac, file.path(opt$out_dir, "DE_Stroma_effect_NoAraC.csv"))

# 4. Stroma effect (within AraC treated samples) -> T16 vs T12
seu_treated <- subset(seu_merged, subset=araC == TRUE)
seu_treated <- PrepSCTFindMarkers(seu_treated)
Idents(seu_treated) <- "stroma"
de_stroma_arac <- FindMarkers(seu_treated, ident.1 = TRUE, ident.2 = FALSE, test.use = "wilcox", logfc.threshold = 0, verbose = FALSE)
write.csv(de_stroma_arac, file.path(opt$out_dir, "DE_Stroma_effect_WithAraC.csv"))

# Cluster markers
all_markers <- FindAllMarkers(seu_merged,
                              only.pos = TRUE,
                              min.pct  = 0.25,
                              logfc.threshold = 0.25,
                              verbose  = FALSE)
write.csv(all_markers, file.path(opt$out_dir, "cluster_markers.csv"), row.names=FALSE)

# --- Top marker dot plot per cluster ---
top5 <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% pull(gene) %>% unique()
pdf(file.path(opt$out_dir, "figures", "04_cluster_marker_dotplot.pdf"), width=16, height=8)
print(DotPlot(seu_merged, features=top5) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        ggtitle("Top 5 markers per cluster"))
dev.off()

# --- Cluster composition bar chart (% per sample) ---
comp_df <- seu_merged@meta.data %>%
  group_by(sample_id, seurat_clusters) %>%
  tally() %>%
  group_by(sample_id) %>%
  mutate(pct_in_sample=n/sum(n)*100) %>%
  ungroup() %>%
  group_by(seurat_clusters) %>%
  mutate(pct_in_cluster=n/sum(n)*100) %>%
  ungroup()

p_comp1 <- ggplot(comp_df, aes(x=sample_id, y=pct_in_sample, fill=seurat_clusters)) +
        geom_bar(stat="identity") +
        labs(x="Sample", y="% of sample", fill="Cluster", title="Cluster composition per sample") +
        theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1))

p_comp2 <- ggplot(comp_df, aes(x=seurat_clusters, y=pct_in_cluster, fill=sample_id)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=sample_colors) +
        labs(x="Cluster", y="% of cluster", fill="Sample", title="Sample composition per cluster") +
        theme_classic()

pdf(file.path(opt$out_dir, "figures", "05a_cluster_composition_per_sample.pdf"), width=8, height=6)
print(p_comp1)
dev.off()

pdf(file.path(opt$out_dir, "figures", "05b_sample_composition_per_cluster.pdf"), width=8, height=6)
print(p_comp2)
dev.off()

# --- Helper function for volcano plots ---
plot_volcano <- function(de_df, title, filename) {
  de_df$gene <- rownames(de_df)
  
  # Categorize genes based on thresholds (logFC > 1 and FDR < 0.05)
  de_df$diffexpressed <- "Not Significant"
  de_df$diffexpressed[de_df$avg_log2FC > 1 & de_df$p_val_adj < 0.05] <- "Up-regulated"
  de_df$diffexpressed[de_df$avg_log2FC < -1 & de_df$p_val_adj < 0.05] <- "Down-regulated"
  
  # Ensure the order of factors for consistent coloring
  de_df$diffexpressed <- factor(de_df$diffexpressed, 
                                levels=c("Up-regulated", "Down-regulated", "Not Significant"))

  p <- ggplot(de_df, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=diffexpressed)) +
          geom_point(alpha=0.6, size=1) +
          scale_colour_manual(values=c("Up-regulated"="firebrick", 
                                       "Down-regulated"="steelblue", 
                                       "Not Significant"="grey80")) +
          geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey40") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
          labs(title=title, x="log2 Fold Change", y="-log10 adj.p-value", colour="Status") +
          theme_classic()
          
  pdf(file.path(opt$out_dir, "figures", filename), width=8, height=7)
  print(p)
  dev.off()
}

plot_volcano(de_araC_nostroma, "AraC effect (No Stroma)", "06a_volcano_AraC_NoStroma.pdf")
plot_volcano(de_araC_stroma,   "AraC effect (With Stroma)", "06b_volcano_AraC_WithStroma.pdf")
plot_volcano(de_stroma_noarac, "Stroma effect (No AraC)", "06c_volcano_Stroma_NoAraC.pdf")
plot_volcano(de_stroma_arac,   "Stroma effect (With AraC)", "06d_volcano_Stroma_WithAraC.pdf")

# =============================================================================
# STEP 9: Pathway / Gene Set Enrichment Analysis (GSEA)
# =============================================================================
message("\n--- Step 9: Pathway Analysis (GSEA) ---")

# Pull MSigDB gene sets (Hallmark, KEGG, GO:BP)
# We use 'collection' as required by msigdbr v10.0.0
m_t2g_h    <- msigdbr(species="Homo sapiens", collection="H") %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g_kegg <- msigdbr(species="Homo sapiens", collection="C2") %>% 
  dplyr::filter(grepl("^KEGG_", gs_name)) %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g_gobp <- msigdbr(species="Homo sapiens", collection="C5") %>% 
  dplyr::filter(grepl("^GOBP_", gs_name)) %>% 
  dplyr::select(gs_name, gene_symbol)

# Helper function to run GSEA and plot
run_gsea <- function(de_results, title_prefix, file_prefix) {
  # Ensure the 'gene' column exists (FindMarkers returns genes as rownames)
  if (!"gene" %in% colnames(de_results)) {
    de_results$gene <- rownames(de_results)
  }
  
  # Clean and sort genes by log2FC for GSEA
  # Remove NA, Inf, and ensure unique names
  de_results <- de_results[!is.na(de_results$avg_log2FC) & is.finite(de_results$avg_log2FC), ]
  gene_list <- de_results$avg_log2FC
  names(gene_list) <- de_results$gene
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing=TRUE)
  
  if (length(gene_list) < 10) {
    message("  Skipping GSEA for ", title_prefix, ": too few genes.")
    return(NULL)
  }
  
  # Hallmark
  set.seed(42)
  gsea_h <- GSEA(gene_list, TERM2GENE=m_t2g_h, pvalueCutoff=0.05, verbose=FALSE)
  if (!is.null(gsea_h) && nrow(gsea_h) > 0) {
    pdf(file.path(opt$out_dir, "figures", paste0(file_prefix, "_GSEA_Hallmark.pdf")), width=12, height=8)
    print(dotplot(gsea_h, showCategory=15, split=".sign") + facet_grid(.~.sign) + ggtitle(paste(title_prefix, "- Hallmark")))
    dev.off()
    write.csv(as.data.frame(gsea_h), file.path(opt$out_dir, paste0(file_prefix, "_GSEA_Hallmark.csv")))
  }
  
  # KEGG
  set.seed(42)
  gsea_kegg <- GSEA(gene_list, TERM2GENE=m_t2g_kegg, pvalueCutoff=0.05, verbose=FALSE)
  if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
    pdf(file.path(opt$out_dir, "figures", paste0(file_prefix, "_GSEA_KEGG.pdf")), width=12, height=8)
    print(dotplot(gsea_kegg, showCategory=15, split=".sign") + facet_grid(.~.sign) + ggtitle(paste(title_prefix, "- KEGG")))
    dev.off()
    write.csv(as.data.frame(gsea_kegg), file.path(opt$out_dir, paste0(file_prefix, "_GSEA_KEGG.csv")))
  }
  
  # GO:BP
  set.seed(42)
  gsea_gobp <- GSEA(gene_list, TERM2GENE=m_t2g_gobp, pvalueCutoff=0.05, verbose=FALSE)
  if (!is.null(gsea_gobp) && nrow(gsea_gobp) > 0) {
    pdf(file.path(opt$out_dir, "figures", paste0(file_prefix, "_GSEA_GOBP.pdf")), width=12, height=8)
    print(dotplot(gsea_gobp, showCategory=15, split=".sign") + facet_grid(.~.sign) + ggtitle(paste(title_prefix, "- GO:BP")))
    dev.off()
    write.csv(as.data.frame(gsea_gobp), file.path(opt$out_dir, paste0(file_prefix, "_GSEA_GOBP.csv")))
  }
}

message("  Running GSEA for AraC effect (No Stroma)...")
run_gsea(de_araC_nostroma, "AraC Effect (No Stroma)", "07a_AraC_NoStroma")

message("  Running GSEA for AraC effect (With Stroma)...")
run_gsea(de_araC_stroma, "AraC Effect (With Stroma)", "07b_AraC_WithStroma")

message("  Running GSEA for Stroma effect (No AraC)...")
run_gsea(de_stroma_noarac, "Stroma Effect (No AraC)", "07c_Stroma_NoAraC")

message("  Running GSEA for Stroma effect (With AraC)...")
run_gsea(de_stroma_arac, "Stroma Effect (With AraC)", "07d_Stroma_WithAraC")

# =============================================================================
# STEP 10: Save final General Seurat object
# =============================================================================
message("\n--- Step 10: Saving General Seurat object ---")

patient_name <- if ("patient" %in% colnames(sample_meta)) unique(sample_meta$patient)[1] else "SRAML10"
if (is.null(patient_name) || is.na(patient_name) || patient_name == "") patient_name <- "Patient"
output_file_name <- paste0("AML_", patient_name, "_seurat_general.RDS")

saveRDS(seu_merged, file=file.path(opt$out_dir, output_file_name))

message("\n=== General Pipeline complete! ===")
message("Output RDS : ", file.path(opt$out_dir, output_file_name))
message("Figures    : ", file.path(opt$out_dir, "figures/"))
