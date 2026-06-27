#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 3b: Clonal Overlap & Sharing
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(ggalluvial)
  library(pheatmap)
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

message("=== AML CloneTracker Analysis 3b: Clonal Overlap ===")
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)
seu_merged <- readRDS(opt$rds_in)

# 1. Filter to confident clones
meta <- seu_merged@meta.data
meta$sample_id <- gsub("^P2_", "", meta$sample_id) # Drop P2_ prefix to save space
meta_conf <- meta %>% filter(clone_confident == TRUE)

if (nrow(meta_conf) == 0) stop("No confident clones found for overlap analysis.")

# 2. Stacked Bar Plot (Proportions)
message("--- Generating Stacked Bar & Alluvial plots ---")
top_n <- 12
clone_freqs <- sort(table(meta_conf$clone_barcode), decreasing=TRUE)
top_clones <- names(clone_freqs)[1:min(length(clone_freqs), top_n)]

meta_conf <- meta_conf %>%
  mutate(clone_group = ifelse(clone_barcode %in% top_clones, clone_barcode, "Other"))
clone_levels <- c(top_clones, "Other")
meta_conf$clone_group <- factor(meta_conf$clone_group, levels=clone_levels)

bar_n <- meta_conf %>% group_by(sample_id) %>% summarise(n = n(), .groups="drop")

p_bar <- ggplot(meta_conf, aes(x=sample_id, fill=clone_group)) +
  geom_bar(position="fill", color="black", linewidth=0.2) +
  geom_text(data=bar_n, aes(x=sample_id, y=1.05, label=paste0("n=", n)), inherit.aes=FALSE, size=3.5) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  theme_classic() +
  labs(x="Sample", y="Proportion of Confident Clones", fill="Clone Barcode", title="Clonal Composition per Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file.path(fig_dir, "09a_CloneTracker_StackedBar.pdf"), width=8, height=7)
print(p_bar)
dev.off()

# 3. Clone x Sample Dot Heatmap
message("--- Generating Dot Heatmap ---")
clone_sample_counts <- meta_conf %>%
  group_by(clone_barcode, sample_id) %>%
  summarise(n_cells = n(), .groups = "drop")

# Shorten clone IDs for readability — use first 8 chars
clone_sample_counts$clone_short <- substr(clone_sample_counts$clone_barcode, 1, 8)

p_dot <- ggplot(clone_sample_counts,
                aes(x = sample_id, y = clone_short, size = n_cells, label = n_cells)) +
  geom_point(colour = "firebrick") +
  geom_text(colour = "white", size = 3) +
  scale_size_area(max_size = 12) +
  theme_classic() +
  labs(title = paste0("Confident clones x samples (n_reads >= 2; ",
                      sum(clone_sample_counts$n_cells), " cells total)"),
       x = "", y = "Clone (first 8 bp of BC14)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "09b_CloneTracker_DotMap.pdf"), width=8, height=7)
print(p_dot)
dev.off()

# 4. Jaccard Overlap Heatmap
message("--- Calculating Clonal Overlap (Jaccard) ---")

# Get unique clones per sample
sample_clones <- split(meta_conf$clone_barcode, meta_conf$sample_id)
sample_clones <- lapply(sample_clones, unique)

samples <- names(sample_clones)
n_samp <- length(samples)
jaccard_mat <- matrix(0, nrow=n_samp, ncol=n_samp, dimnames=list(samples, samples))

for (i in 1:n_samp) {
  for (j in 1:n_samp) {
    if (i == j) {
      jaccard_mat[i,j] <- 1
    } else {
      intersect_len <- length(intersect(sample_clones[[i]], sample_clones[[j]]))
      union_len <- length(union(sample_clones[[i]], sample_clones[[j]]))
      jaccard_mat[i,j] <- intersect_len / union_len
    }
  }
}

# Dynamic column annotations based on Seurat metadata
ann_col <- seu_merged@meta.data %>%
  select(sample_id, stroma, araC) %>%
  mutate(sample_id = gsub("^P2_", "", sample_id)) %>%
  distinct() %>%
  mutate(stroma = as.character(stroma), araC = as.character(araC)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sample_id")

ann_col <- ann_col[samples, , drop=FALSE]

pdf(file.path(fig_dir, "10_CloneTracker_Jaccard_Overlap.pdf"), width=8, height=7)
pheatmap(jaccard_mat, display_numbers = TRUE, main="Jaccard Index of Clonal Overlap", 
         annotation_col = ann_col,
         color = colorRampPalette(c("white", "firebrick"))(50))
dev.off()
write.csv(jaccard_mat, file.path(opt$out_dir, paste0(patient_prefix, "_jaccard_overlap.csv")))

# 4b. Explicit Barcode Survival Metric
message("--- Calculating Barcode Survival (AraC effect) ---")

get_samp <- function(s_val, a_val) {
  # Get sample_id from original metadata (even if 0 clones)
  unique(seu_merged@meta.data$sample_id[seu_merged@meta.data$stroma == s_val & seu_merged@meta.data$araC == a_val])[1]
}

samp_nostroma_noarac <- get_samp(FALSE, FALSE)
samp_nostroma_arac   <- get_samp(FALSE, TRUE)
samp_stroma_noarac   <- get_samp(TRUE, FALSE)
samp_stroma_arac     <- get_samp(TRUE, TRUE)

get_clones <- function(s) {
  if (is.na(s) || is.null(sample_clones[[s]])) return(character(0))
  return(sample_clones[[s]])
}

survival_df <- data.frame(
  pair = c("noStroma_AraC", "Stroma_AraC"),
  pre_clones = c(length(get_clones(samp_nostroma_noarac)),
                 length(get_clones(samp_stroma_noarac))),
  surviving = c(length(intersect(get_clones(samp_nostroma_noarac), get_clones(samp_nostroma_arac))),
                length(intersect(get_clones(samp_stroma_noarac), get_clones(samp_stroma_arac))))
)
survival_df$pct_surviving <- round(100 * survival_df$surviving / ifelse(survival_df$pre_clones == 0, 1, survival_df$pre_clones), 2)
write.csv(survival_df, file.path(opt$out_dir, paste0(patient_prefix, "_AraC_survival.csv")), row.names=FALSE)

# 5. Output 4-way Overlap list
message("--- Exporting 4-way Overlap List ---")
all_clones <- unique(meta_conf$clone_barcode)
overlap_df <- data.frame(clone_barcode = all_clones)
for (s in samples) {
  overlap_df[[s]] <- overlap_df$clone_barcode %in% sample_clones[[s]]
}
overlap_df$n_samples <- rowSums(overlap_df[, -1])
write.csv(overlap_df %>% arrange(desc(n_samples)), file.path(opt$out_dir, paste0(patient_prefix, "_4way_clonal_overlap.csv")), row.names=FALSE)

message("Overlap Analysis Complete!")
