#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 05: Stress signalling + within-cluster
# cell-cycle / quiescence comparisons (Part-3 §4 + §6)
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-24
# Last update:  2026-05-24
# Input:        --rds_in     Seurat .rds after seurat_04b (cell_type_final,
#                            is_aml, vgalen_class, existing UCell scores)
#               --out_dir    Output directory
# Output:       <rds_in>                              ← re-saved with:
#                 S.Score, G2M.Score, Phase
#                 UCell_hallmark_<set>  (UPR, HYPOXIA, ROS, TNFA, INFLAM, P53)
#                 condition, treatment (parsed from sample_id if not present)
#               <out_dir>/figures/05a_stress_violins.pdf
#               <out_dir>/figures/05b_quiescence_promono.pdf
#               <out_dir>/figures/05c_abcb5_arac.pdf
#               <out_dir>/figures/05d_cellcycle_phase_proportions.pdf
#               <out_dir>/figures/05e_effect_size_heatmap.pdf
#               <out_dir>/tables/05_stress_effect_sizes.tsv
#               <out_dir>/tables/05_cell_cycle_summary.tsv
#               <out_dir>/tables/05_promono_quiescence_stats.tsv
#               <out_dir>/tables/05_abcb5_arac_stats.tsv
#               <out_dir>/logs/seurat_05_sessionInfo.txt
# Depends on:   Seurat (>=5.0), UCell (>=2.6), msigdbr, ggplot2, dplyr, tibble,
#               tidyr, ggpubr, pheatmap, optparse
# Notes:
#   * Tests three hypotheses simultaneously:
#       - H1 (conv harsh environment)  → general stress higher in conv
#       - H4 (MSC actively quiesces)   → niche-response higher in co-culture,
#                                        cycling lower in co-culture-Promono
#       - H3b (resistant selection)    → ABCB5+ untreated already low-cycling
#                                        and pro-survival before Ara-C
#   * Effect size: Cliff's delta (preferred) via effsize::cliff.delta, with
#     a fallback to Wilcoxon-based rank-biserial correlation if effsize is
#     not available.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggpubr)
  library(pheatmap)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 20 * 1024^3)

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),  type = "character",
                     help = "Seurat .rds after seurat_04b")
parser <- add_option(parser, c("--out_dir"), type = "character",
                     help = "Output directory")
parser <- add_option(parser, c("--ncores"),  type = "integer", default = 8,
                     help = "Cores for UCell [default: %default]")
opt <- parse_args(parser)

for (req in c("rds_in", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(tab_dir, fig_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== Seurat 05: stress signalling + within-cluster cycling ===")
message("  rds_in  : ", opt$rds_in)
message("  out_dir : ", opt$out_dir)

# ---- 1. Load and parse condition / treatment from sample_id --------------

message("--- Loading Seurat object ---")
seu_merged <- readRDS(opt$rds_in)
DefaultAssay(seu_merged) <- "SCT"

if (!"cell_type_final" %in% colnames(seu_merged@meta.data)) {
  stop("cell_type_final not found — run seurat_04b first")
}
if (!"sample_id" %in% colnames(seu_merged@meta.data)) {
  stop("sample_id column required for condition/treatment parsing")
}

# Parse condition + treatment from sample_id naming (e.g. P2_T11_noStroma_noAraC)
parse_condition <- function(x) {
  ifelse(grepl("noStroma", x, ignore.case = TRUE), "Conventional",
         ifelse(grepl("Stroma",   x, ignore.case = TRUE), "Co-culture",
                NA_character_))
}
parse_treatment <- function(x) {
  ifelse(grepl("noAraC", x, ignore.case = TRUE), "Untreated",
         ifelse(grepl("AraC", x, ignore.case = TRUE), "Ara-C",
                NA_character_))
}
# Overwrite condition and treatment columns unconditionally to ensure consistent labeling
seu_merged$condition <- parse_condition(as.character(seu_merged$sample_id))
seu_merged$treatment <- parse_treatment(as.character(seu_merged$sample_id))
seu_merged$condition <- factor(seu_merged$condition, levels = c("Conventional", "Co-culture"))
seu_merged$treatment <- factor(seu_merged$treatment, levels = c("Untreated", "Ara-C"))

message("--- Condition × treatment matrix ---")
print(table(seu_merged$condition, seu_merged$treatment, useNA = "ifany"))

# ---- 2. Cell-cycle scoring (Tirosh S + G2M) ------------------------------

message("--- Cell-cycle scoring ---")
if (!all(c("S.Score", "G2M.Score", "Phase") %in% colnames(seu_merged@meta.data))) {
  seu_merged <- CellCycleScoring(
    seu_merged,
    s.features  = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    set.ident   = FALSE
  )
} else {
  message("  cell-cycle scores already present — skipping")
}
seu_merged$Phase <- factor(seu_merged$Phase, levels = c("G1", "S", "G2M"))

# ---- 3. Add Hallmark stress signatures via msigdbr + UCell ---------------

message("--- Scoring Hallmark stress signatures with UCell ---")
hallmark_sets <- c(
  hallmark_upr      = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  hallmark_hypoxia  = "HALLMARK_HYPOXIA",
  hallmark_ros      = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  hallmark_tnfa     = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  hallmark_inflam   = "HALLMARK_INFLAMMATORY_RESPONSE",
  hallmark_p53      = "HALLMARK_P53_PATHWAY",
  hallmark_apoptosis = "HALLMARK_APOPTOSIS"
)
msig <- msigdbr(species = "Homo sapiens", category = "H")
sig_lists <- lapply(hallmark_sets, function(gs) {
  unique(msig$gene_symbol[msig$gs_name == gs])
})
names(sig_lists) <- names(hallmark_sets)
sig_lists <- sig_lists[lengths(sig_lists) > 0]

need_score <- !paste0("UCell_", names(sig_lists)) %in% colnames(seu_merged@meta.data)
if (any(need_score)) {
  seu_merged <- AddModuleScore_UCell(
    seu_merged,
    features = sig_lists[need_score],
    ncores   = opt$ncores,
    name     = ""
  )
  new_cols     <- names(sig_lists)[need_score]
  renamed_cols <- paste0("UCell_", new_cols)
  seu_merged@meta.data[, renamed_cols] <- seu_merged@meta.data[, new_cols]
  seu_merged@meta.data[, new_cols]     <- NULL
}

# ---- 4. Effect-size helper ----------------------------------------------
# Cliff's delta if effsize available, else rank-biserial from Wilcoxon U

cliff_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  if (requireNamespace("effsize", quietly = TRUE)) {
    res <- effsize::cliff.delta(x, y)
    return(as.numeric(res$estimate))
  }
  # Rank-biserial via Wilcoxon
  w <- suppressWarnings(wilcox.test(x, y))$statistic
  n1 <- length(x); n2 <- length(y)
  unname(2 * w / (n1 * n2) - 1)
}

wilcox_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, y))$p.value
}

# ---- 5. §4 — Stress comparison: AML clusters, conv vs co-culture (untreated) ----

message("--- §4: stress signalling, conv vs co-culture (untreated only) ---")

aml_clusters <- levels(seu_merged$cell_type_final)[
  levels(seu_merged$cell_type_final) %in%
    as.character(unique(seu_merged$cell_type_final[seu_merged$is_aml]))
]

# Robustness check: filter out clusters with extremely low cell counts (< 5)
# in either condition (Conventional vs Co-culture) in the untreated arm.
# This prevents density estimation errors in geom_violin (which crash PDF generation)
# and avoids statistically meaningless comparisons.
meta_temp <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  filter(treatment == "Untreated",
         !is.na(condition),
         cell_type_final %in% aml_clusters)

cell_counts_by_cond <- table(meta_temp$cell_type_final, meta_temp$condition)
valid_clusters <- rownames(cell_counts_by_cond)[
  rowSums(cell_counts_by_cond >= 5) == 2
]
excluded <- setdiff(aml_clusters, valid_clusters)
if (length(excluded) > 0) {
  message("Excluding low-cell clusters from stress comparisons: ", paste(excluded, collapse = ", "))
}
aml_clusters <- intersect(aml_clusters, valid_clusters)
stress_sigs <- c(
  "UCell_hallmark_upr",
  "UCell_hallmark_hypoxia",
  "UCell_hallmark_ros",
  "UCell_hallmark_tnfa",
  "UCell_hallmark_inflam",
  "UCell_hallmark_p53",
  "UCell_hallmark_apoptosis",
  "UCell_isr_atf4_targets",
  "UCell_van_den_brink_dissociation",
  "UCell_baryawno_niche_signal",
  "UCell_tikhonova_msc_response"
)
stress_sigs <- intersect(stress_sigs, colnames(seu_merged@meta.data))

meta <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  filter(treatment == "Untreated",
         !is.na(condition),
         cell_type_final %in% aml_clusters)

# Effect-size matrix: cluster × signature
effect_tbl <- expand_grid(cell_type_final = aml_clusters, signature = stress_sigs) |>
  rowwise() |>
  mutate(
    n_conv   = sum(meta$cell_type_final == cell_type_final & meta$condition == "Conventional"),
    n_coc    = sum(meta$cell_type_final == cell_type_final & meta$condition == "Co-culture"),
    delta    = cliff_delta(
                 meta[[signature]][meta$cell_type_final == cell_type_final &
                                    meta$condition == "Co-culture"],
                 meta[[signature]][meta$cell_type_final == cell_type_final &
                                    meta$condition == "Conventional"]),
    p_value  = wilcox_p(
                 meta[[signature]][meta$cell_type_final == cell_type_final &
                                    meta$condition == "Co-culture"],
                 meta[[signature]][meta$cell_type_final == cell_type_final &
                                    meta$condition == "Conventional"])
  ) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write.table(effect_tbl, file.path(tab_dir, "05_stress_effect_sizes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Effect-size heatmap (Cliff's delta: positive = higher in co-culture)
es_wide <- effect_tbl |>
  mutate(signature = sub("^UCell_", "", signature)) |>
  select(cell_type_final, signature, delta) |>
  pivot_wider(names_from = signature, values_from = delta) |>
  column_to_rownames("cell_type_final") |>
  as.matrix()

pdf(file.path(fig_dir, "05e_effect_size_heatmap.pdf"),
    width = max(9, ncol(es_wide) * 0.5 + 4),
    height = max(5, nrow(es_wide) * 0.45 + 2))
pheatmap(
  es_wide,
  cluster_rows = TRUE, cluster_cols = TRUE,
  color = colorRampPalette(c("#1B4F72", "#1B98E0", "white", "#D7263D", "#6B0F1A"))(50),
  breaks = seq(-0.6, 0.6, length.out = 51),
  border_color = "grey90",
  cellwidth = 18, cellheight = 18,
  fontsize = 10, fontsize_row = 10, fontsize_col = 9,
  angle_col = 45,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Stress / niche programmes — Cliff's delta (Co-culture vs Conventional, untreated)"
)
dev.off()

# Violin plot per AML cluster
plot_df <- meta |>
  select(cell_type_final, condition, all_of(stress_sigs)) |>
  pivot_longer(cols = all_of(stress_sigs),
               names_to = "signature", values_to = "score") |>
  mutate(signature = sub("^UCell_", "", signature))

p_violins <- ggplot(plot_df, aes(x = condition, y = score, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.size = 0.1, fill = "white") +
  facet_grid(signature ~ cell_type_final, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c("Conventional" = "#3D5A80", "Co-culture" = "#D7263D")) +
  theme_classic(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.text.x = element_text(angle = 30, hjust = 0.5, size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Stress / niche signatures — conventional vs co-culture (untreated)")

pdf(file.path(fig_dir, "05a_stress_violins.pdf"),
    width = max(10, length(aml_clusters) * 1.4),
    height = max(10, length(stress_sigs) * 1.2))
print(p_violins)
dev.off()

# ---- 6. §6a — Quiescence in Promono cluster (conv vs co-culture, untreated) ----

message("--- §6a: quiescence in Promono cluster ---")
quiescence_sigs <- intersect(c("S.Score", "G2M.Score",
                               "UCell_laurenti_hsc_quiescence",
                               "UCell_cheung_rando_quiescence",
                               "UCell_tikhonova_msc_response",
                               "UCell_isr_atf4_targets"),
                             colnames(seu_merged@meta.data))

promono_data <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  filter(treatment == "Untreated",
         cell_type_final == "Promono quiescent AML",
         !is.na(condition))

promono_stats <- tibble(
  signature = quiescence_sigs,
  delta = sapply(quiescence_sigs, function(s) {
    cliff_delta(promono_data[[s]][promono_data$condition == "Co-culture"],
                promono_data[[s]][promono_data$condition == "Conventional"])
  }),
  p_value = sapply(quiescence_sigs, function(s) {
    wilcox_p(promono_data[[s]][promono_data$condition == "Co-culture"],
             promono_data[[s]][promono_data$condition == "Conventional"])
  })
) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write.table(promono_stats, file.path(tab_dir, "05_promono_quiescence_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

promono_long <- promono_data |>
  select(condition, all_of(quiescence_sigs)) |>
  pivot_longer(-condition, names_to = "signature", values_to = "score") |>
  mutate(signature = sub("^UCell_", "", signature))

n_conv <- sum(promono_data$condition == "Conventional", na.rm = TRUE)
n_coc  <- sum(promono_data$condition == "Co-culture", na.rm = TRUE)
if (n_conv >= 3 && n_coc >= 3) {
  pdf(file.path(fig_dir, "05b_quiescence_promono.pdf"),
      width = max(8, length(quiescence_sigs) * 1.4), height = 6)
  print(
    ggplot(promono_long, aes(x = condition, y = score, fill = condition)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white") +
      facet_wrap(~ signature, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Conventional" = "#3D5A80", "Co-culture" = "#D7263D")) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = NULL, y = "Score",
           title = "Promono quiescent AML — quiescence + cycling, conv vs co-culture (untreated)")
  )
  dev.off()
} else {
  message("Skipping 05b Promono plot: insufficient cells (Conventional: ", n_conv, ", Co-culture: ", n_coc, ")")
}

# ---- 7. §6b — ABCB5+ untreated vs Ara-C in conventional arm --------------

message("--- §6b: ABCB5+ untreated vs Ara-C (conventional) ---")
abcb5_data <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  filter(cell_type_final == "ABCB5+ resistant LSC AML",
         condition == "Conventional",
         !is.na(treatment))

abcb5_sigs <- intersect(c("S.Score", "G2M.Score",
                          "UCell_arac_metabolism",
                          "UCell_mdr_efflux",
                          "UCell_aldh_resistance",
                          "UCell_lsc17",
                          "UCell_bcl2_family_pro_survival",
                          "UCell_bcl2_family_pro_apoptotic",
                          "UCell_isr_atf4_targets"),
                        colnames(seu_merged@meta.data))

abcb5_stats <- tibble(
  signature = abcb5_sigs,
  delta = sapply(abcb5_sigs, function(s) {
    cliff_delta(abcb5_data[[s]][abcb5_data$treatment == "Ara-C"],
                abcb5_data[[s]][abcb5_data$treatment == "Untreated"])
  }),
  p_value = sapply(abcb5_sigs, function(s) {
    wilcox_p(abcb5_data[[s]][abcb5_data$treatment == "Ara-C"],
             abcb5_data[[s]][abcb5_data$treatment == "Untreated"])
  })
) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write.table(abcb5_stats, file.path(tab_dir, "05_abcb5_arac_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

abcb5_long <- abcb5_data |>
  select(treatment, all_of(abcb5_sigs)) |>
  pivot_longer(-treatment, names_to = "signature", values_to = "score") |>
  mutate(signature = sub("^UCell_", "", signature))

n_untrt <- sum(abcb5_data$treatment == "Untreated", na.rm = TRUE)
n_arac  <- sum(abcb5_data$treatment == "Ara-C", na.rm = TRUE)
if (n_untrt >= 3 && n_arac >= 3) {
  pdf(file.path(fig_dir, "05c_abcb5_arac.pdf"),
      width = max(10, length(abcb5_sigs) * 1.3), height = 6)
  print(
    ggplot(abcb5_long, aes(x = treatment, y = score, fill = treatment)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white") +
      facet_wrap(~ signature, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Untreated" = "#1B98E0", "Ara-C" = "#D7263D")) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = NULL, y = "Score",
           title = "ABCB5+ resistant LSC AML — untreated vs Ara-C, conventional arm")
  )
  dev.off()
} else {
  message("Skipping 05c ABCB5+ Ara-C plot: insufficient cells (Untreated: ", n_untrt, ", Ara-C: ", n_arac, ")")
}

# ---- 8. Cell-cycle phase proportions: cluster × condition (AML clusters) ----

message("--- Cell-cycle phase proportions ---")
phase_df <- seu_merged@meta.data |>
  as_tibble(rownames = "cell") |>
  filter(treatment == "Untreated",
         cell_type_final %in% aml_clusters,
         !is.na(condition),
         !is.na(Phase)) |>
  group_by(cell_type_final, condition, Phase) |>
  summarise(n = n(), .groups = "drop_last") |>
  mutate(prop = n / sum(n)) |>
  ungroup()

write.table(phase_df, file.path(tab_dir, "05_cell_cycle_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(fig_dir, "05d_cellcycle_phase_proportions.pdf"),
    width = max(10, length(aml_clusters) * 1.4), height = 6)
print(
  ggplot(phase_df, aes(x = condition, y = prop, fill = Phase)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.2) +
    facet_wrap(~ cell_type_final, nrow = 1) +
    scale_fill_manual(values = c(G1 = "#1B98E0", S = "#E66B00", G2M = "#D7263D")) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text = element_text(size = 8)) +
    labs(x = NULL, y = "Cell fraction",
         title = "Cell-cycle phase by cluster — conv vs co-culture (untreated)")
)
dev.off()

# ---- 9. Save updated RDS + session info ---------------------------------

message("--- Saving updated RDS ---")
saveRDS(seu_merged, file = opt$rds_in)
message("  updated: ", opt$rds_in)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_05_sessionInfo.txt"))

message("Done!")
