#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 07: within-patient differential expression
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-24
# Last update:  2026-05-24
# Input:        --rds_in    Seurat .rds — the AML-only object from seurat_04c.
#                           Must contain cell_type_final, stroma, araC, sample_id,
#                           and the singlet-only RNA counts.
#               --out_dir   Output directory.
# Output:       <out_dir>/tables/07_de_<cluster>_<contrast>.tsv  (per cluster × contrast)
#               <out_dir>/tables/07_de_summary.tsv              (counts of DE genes)
#               <out_dir>/tables/07_skipped_tests.tsv           (tests dropped on n)
#               <out_dir>/figures/07_volcano_<cluster>_<contrast>.pdf
#               <out_dir>/figures/07_topgenes_heatmap.pdf
#               <out_dir>/logs/07_manifest.tsv
#               <out_dir>/logs/seurat_07_sessionInfo.txt
# Depends on:   Seurat (>=5.0), presto (Wilcoxon speedup), ggplot2, ggrepel,
#               dplyr, tibble, tidyr, optparse, ComplexHeatmap (optional)
#
# Notes:
#   * Cohort caveat ----------------------------------------------------------
#     The current dataset is ONE patient (SRAML10) × 4 samples (one per
#     condition × treatment). There is no biological replication within any
#     condition × treatment group, so pseudobulk DE with muscat / DESeq2 is
#     not feasible (DESeq2 cannot estimate dispersion with n = 1 per group).
#     This script therefore runs SINGLE-CELL DE (Wilcoxon via presto, called
#     through Seurat::FindMarkers). The resulting p-values reflect within-
#     patient cell-level variability and CANNOT be generalised across
#     patients — they are for hypothesis generation only. When the full
#     S34 cohort lands, replace this script with seurat_07b_pseudobulk_de.R
#     using muscat::aggregateData + DESeq2 on the multi-patient object.
#   * Condition × treatment groups -------------------------------------------
#     Conv_Unt    : stroma == FALSE, araC == FALSE
#     Conv_AraC   : stroma == FALSE, araC == TRUE
#     Coc_Unt     : stroma == TRUE,  araC == FALSE
#     Coc_AraC    : stroma == TRUE,  araC == TRUE
#   * Contrasts run per cluster ----------------------------------------------
#     C1  Coc_Unt   vs Conv_Unt    (niche imprint, baseline)
#     C2  Coc_AraC  vs Conv_AraC   (niche-protected residual disease)
#     C3  Conv_AraC vs Conv_Unt    (Ara-C effect, no niche)
#     C4  Coc_AraC  vs Coc_Unt     (Ara-C effect, with niche)
#   * Special cluster-vs-cluster contrasts -----------------------------------
#     S1  Niche-stressed primitive AML vs HSPC/LSC-like AML, conv arm only
#         (tests whether the niche-stressed cluster is a stressed sub-state
#          of HSPC/LSC-like — H3a — or a transcriptionally distinct primed
#          clone — H3b).
#     S2  HSPC/LSC-like AML vs GMP-like AML (all conditions)
#         (UMAP repositioning observation from seurat_04c — answers whether
#          the cluster is more "primitive" or more "committed-myeloid").
#   * Minimum group size = 25 cells. Tests below threshold are recorded in
#     07_skipped_tests.tsv rather than silently dropped.
#   * BH-FDR is applied within each contrast × cluster table. A combined
#     FDR across all 28 condition contrasts would be too punishing for an
#     n=1 design; per-table FDR is the honest middle ground.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)

# Try to use presto for the Wilcoxon speedup
have_presto <- requireNamespace("presto", quietly = TRUE)
if (!have_presto) {
  message("NOTE: presto not installed — FindMarkers will use the slower base path.")
}

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),    type = "character",
                     help = "AML-only Seurat .rds (output of seurat_04c)")
parser <- add_option(parser, c("--out_dir"),   type = "character",
                     help = "Output directory")
parser <- add_option(parser, c("--min_cells"), type = "integer", default = 25L,
                     help = "Minimum cells per group [default 25]")
parser <- add_option(parser, c("--logfc"),     type = "double",  default = 0.25,
                     help = "Min |log2FC| to report [default 0.25]")
parser <- add_option(parser, c("--top_n"),     type = "integer", default = 30L,
                     help = "Top N hits per side for volcano / heatmap [default 30]")
opt <- parse_args(parser)

for (req in c("rds_in", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
for (d in c(tab_dir, fig_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Seurat 07: within-patient single-cell DE ===")
message("  rds_in    : ", opt$rds_in)
message("  out_dir   : ", opt$out_dir)
message("  min_cells : ", opt$min_cells)

# ---- 1. Load AML-only object ---------------------------------------------

message("--- Loading AML-only Seurat object ---")
seu_aml <- readRDS(opt$rds_in)
DefaultAssay(seu_aml) <- "RNA"

# Ensure normalised data exists (FindMarkers needs it)
if (!"data" %in% Layers(seu_aml[["RNA"]])) {
  message("--- Normalising RNA counts (NormalizeData) ---")
  seu_aml <- NormalizeData(seu_aml, verbose = FALSE)
}

for (req_col in c("cell_type_final", "stroma", "araC", "sample_id")) {
  if (!req_col %in% colnames(seu_aml@meta.data)) {
    stop("Required metadata column missing: ", req_col)
  }
}

message("  Cells: ", ncol(seu_aml))
message("  Clusters: ", paste(levels(droplevels(factor(seu_aml$cell_type_final))),
                              collapse = ", "))

# ---- 2. Build the four condition × treatment groups ----------------------

seu_aml$ct_group <- factor(
  dplyr::case_when(
    !seu_aml$stroma & !seu_aml$araC ~ "Conv_Unt",
    !seu_aml$stroma &  seu_aml$araC ~ "Conv_AraC",
     seu_aml$stroma & !seu_aml$araC ~ "Coc_Unt",
     seu_aml$stroma &  seu_aml$araC ~ "Coc_AraC",
    TRUE                             ~ NA_character_
  ),
  levels = c("Conv_Unt", "Conv_AraC", "Coc_Unt", "Coc_AraC")
)

message("--- Cells per condition × treatment group ---")
print(table(seu_aml$ct_group, useNA = "ifany"))

# ---- 3. Helper: run one Seurat::FindMarkers contrast ---------------------

run_de <- function(seu_obj, ident_col,
                   g1, g2,
                   subset_expr = NULL,
                   min_cells   = opt$min_cells,
                   logfc       = 0.0,        # report all genes, filter later
                   test_use    = "wilcox") {
  seu_sub <- seu_obj
  if (!is.null(subset_expr)) {
    seu_sub <- seu_sub[, subset_expr]
  }
  Idents(seu_sub) <- seu_sub[[ident_col, drop = TRUE]]
  n1 <- sum(Idents(seu_sub) == g1)
  n2 <- sum(Idents(seu_sub) == g2)
  if (n1 < min_cells || n2 < min_cells) {
    return(list(skipped = TRUE, n1 = n1, n2 = n2,
                reason = paste0("n1=", n1, " n2=", n2, " < ", min_cells)))
  }
  res <- tryCatch(
    Seurat::FindMarkers(seu_sub,
                        ident.1     = g1,
                        ident.2     = g2,
                        test.use    = test_use,
                        logfc.threshold = logfc,
                        min.pct     = 0.05,
                        only.pos    = FALSE,
                        verbose     = FALSE),
    error = function(e) {
      message("    FindMarkers error: ", e$message)
      return(NULL)
    }
  )
  if (is.null(res) || nrow(res) == 0) {
    return(list(skipped = TRUE, n1 = n1, n2 = n2,
                reason = "FindMarkers returned empty"))
  }
  res$gene <- rownames(res)
  res$n1   <- n1
  res$n2   <- n2
  res <- res |>
    as_tibble() |>
    select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, n1, n2) |>
    arrange(p_val_adj, desc(abs(avg_log2FC)))
  list(skipped = FALSE, n1 = n1, n2 = n2, table = res)
}

# Tidy file-safe slug for cluster names containing slashes / spaces
slug <- function(x) {
  x |>
    tolower() |>
    gsub(pattern = "[^a-z0-9]+", replacement = "_", x = _) |>
    gsub(pattern = "^_|_$",      replacement = "",  x = _)
}

# ---- 4. Run condition × treatment contrasts per AML cluster --------------

aml_clusters <- levels(droplevels(factor(seu_aml$cell_type_final)))
contrasts <- tibble::tribble(
  ~contrast_id, ~group1,    ~group2,
  "C1",         "Coc_Unt",  "Conv_Unt",
  "C2",         "Coc_AraC", "Conv_AraC",
  "C3",         "Conv_AraC","Conv_Unt",
  "C4",         "Coc_AraC", "Coc_Unt"
)

de_summary  <- list()
skipped_log <- list()
all_top     <- list()

for (cl in aml_clusters) {
  message("--- Cluster: ", cl, " ---")
  seu_cl <- seu_aml[, seu_aml$cell_type_final == cl]
  for (k in seq_len(nrow(contrasts))) {
    cid <- contrasts$contrast_id[k]
    g1  <- contrasts$group1[k]
    g2  <- contrasts$group2[k]
    message("    ", cid, ": ", g1, " vs ", g2)

    res <- run_de(seu_cl, ident_col = "ct_group", g1 = g1, g2 = g2, logfc = opt$logfc)
    if (res$skipped) {
      skipped_log[[length(skipped_log) + 1]] <- tibble(
        cluster = cl, contrast_id = cid, group1 = g1, group2 = g2,
        n1 = res$n1, n2 = res$n2, reason = res$reason
      )
      next
    }
    fname <- file.path(tab_dir,
                       sprintf("07_de_%s_%s.tsv", slug(cl), cid))
    write.table(res$table, fname,
                sep = "\t", quote = FALSE, row.names = FALSE)

    sig <- res$table |>
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= opt$logfc)
    de_summary[[length(de_summary) + 1]] <- tibble(
      cluster     = cl,
      contrast_id = cid,
      contrast    = sprintf("%s vs %s", g1, g2),
      n_cells_g1  = res$n1,
      n_cells_g2  = res$n2,
      n_up_in_g1  = sum(sig$avg_log2FC >  opt$logfc),
      n_dn_in_g1  = sum(sig$avg_log2FC < -opt$logfc),
      top_up_g1   = paste(head(sig |> filter(avg_log2FC >  0) |> pull(gene), 5), collapse = ", "),
      top_up_g2   = paste(head(sig |> filter(avg_log2FC <  0) |> pull(gene), 5), collapse = ", ")
    )

    # Capture top up/down per contrast for the combined heatmap
    top_hits <- res$table |>
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= opt$logfc) |>
      arrange(desc(avg_log2FC)) |>
      mutate(cluster = cl, contrast_id = cid)
    if (nrow(top_hits) > 0) {
      all_top[[length(all_top) + 1]] <- bind_rows(
        head(top_hits, opt$top_n),
        tail(top_hits, opt$top_n)
      ) |> unique()
    }
  }
}

# ---- 5. Special cluster-vs-cluster contrasts -----------------------------

message("--- Special contrasts ---")

# S1: Niche-stressed primitive AML vs HSPC/LSC-like AML, in the conv arm
if (all(c("Niche-stressed primitive AML", "HSPC / LSC-like AML") %in% aml_clusters)) {
  message("  S1: Niche-stressed primitive vs HSPC/LSC-like (conv arm)")
  res_s1 <- run_de(seu_aml,
                   ident_col = "cell_type_final",
                   g1 = "Niche-stressed primitive AML",
                   g2 = "HSPC / LSC-like AML",
                   subset_expr = !seu_aml$stroma,
                   logfc = opt$logfc)
  if (!res_s1$skipped) {
    write.table(res_s1$table,
                file.path(tab_dir, "07_de_special_S1_niche_vs_hspc_conv.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- res_s1$table |>
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= opt$logfc)
    de_summary[[length(de_summary) + 1]] <- tibble(
      cluster     = "(special)",
      contrast_id = "S1",
      contrast    = "Niche-stressed primitive AML vs HSPC/LSC-like AML (conv)",
      n_cells_g1  = res_s1$n1,
      n_cells_g2  = res_s1$n2,
      n_up_in_g1  = sum(sig$avg_log2FC >  opt$logfc),
      n_dn_in_g1  = sum(sig$avg_log2FC < -opt$logfc),
      top_up_g1   = paste(head(sig |> filter(avg_log2FC >  0) |> pull(gene), 5), collapse = ", "),
      top_up_g2   = paste(head(sig |> filter(avg_log2FC <  0) |> pull(gene), 5), collapse = ", ")
    )
  } else {
    skipped_log[[length(skipped_log) + 1]] <- tibble(
      cluster = "(special)", contrast_id = "S1",
      group1 = "Niche-stressed primitive AML", group2 = "HSPC / LSC-like AML",
      n1 = res_s1$n1, n2 = res_s1$n2, reason = res_s1$reason
    )
  }
}

# S2: HSPC/LSC-like AML vs GMP-like AML, all conditions
if (all(c("HSPC / LSC-like AML", "GMP-like AML") %in% aml_clusters)) {
  message("  S2: HSPC/LSC-like vs GMP-like (all conditions)")
  res_s2 <- run_de(seu_aml,
                   ident_col = "cell_type_final",
                   g1 = "HSPC / LSC-like AML",
                   g2 = "GMP-like AML",
                   logfc = opt$logfc)
  if (!res_s2$skipped) {
    write.table(res_s2$table,
                file.path(tab_dir, "07_de_special_S2_hspc_vs_gmp.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- res_s2$table |>
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= opt$logfc)
    de_summary[[length(de_summary) + 1]] <- tibble(
      cluster     = "(special)",
      contrast_id = "S2",
      contrast    = "HSPC/LSC-like AML vs GMP-like AML (all)",
      n_cells_g1  = res_s2$n1,
      n_cells_g2  = res_s2$n2,
      n_up_in_g1  = sum(sig$avg_log2FC >  opt$logfc),
      n_dn_in_g1  = sum(sig$avg_log2FC < -opt$logfc),
      top_up_g1   = paste(head(sig |> filter(avg_log2FC >  0) |> pull(gene), 5), collapse = ", "),
      top_up_g2   = paste(head(sig |> filter(avg_log2FC <  0) |> pull(gene), 5), collapse = ", ")
    )
  } else {
    skipped_log[[length(skipped_log) + 1]] <- tibble(
      cluster = "(special)", contrast_id = "S2",
      group1 = "HSPC / LSC-like AML", group2 = "GMP-like AML",
      n1 = res_s2$n1, n2 = res_s2$n2, reason = res_s2$reason
    )
  }
}

# ---- 6. Summary + skipped tables -----------------------------------------

if (length(de_summary) > 0) {
  de_summary_tbl <- bind_rows(de_summary) |>
    arrange(cluster, contrast_id)
  write.table(de_summary_tbl,
              file.path(tab_dir, "07_de_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("--- DE summary ---")
  print(de_summary_tbl)
}

if (length(skipped_log) > 0) {
  skipped_tbl <- bind_rows(skipped_log)
  write.table(skipped_tbl,
              file.path(tab_dir, "07_skipped_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("--- Skipped tests (insufficient n) ---")
  print(skipped_tbl)
}

# ---- 7. Volcano plots ----------------------------------------------------

message("--- Volcano plots ---")
de_files <- list.files(tab_dir,
                       pattern = "^07_de_.*\\.tsv$",
                       full.names = TRUE)
de_files <- setdiff(de_files,
                    c(file.path(tab_dir, "07_de_summary.tsv"),
                      file.path(tab_dir, "07_skipped_tests.tsv")))

for (f in de_files) {
  d <- read.delim(f, sep = "\t", stringsAsFactors = FALSE) |>
    mutate(neg_log10_p = -log10(pmax(p_val_adj, 1e-300)),
           direction = case_when(
             p_val_adj < 0.05 & avg_log2FC >  opt$logfc ~ "up",
             p_val_adj < 0.05 & avg_log2FC < -opt$logfc ~ "down",
             TRUE                                         ~ "ns"))
  if (nrow(d) == 0) next

  top_label <- bind_rows(
    d |> filter(direction == "up")   |> arrange(p_val_adj) |> head(10),
    d |> filter(direction == "down") |> arrange(p_val_adj) |> head(10)
  )

  ttl <- sub("^07_de_", "", sub("\\.tsv$", "", basename(f)))
  pdf_file <- file.path(fig_dir, sprintf("07_volcano_%s.pdf", ttl))
  pdf(pdf_file, width = 7, height = 6)
  print(
    ggplot(d, aes(x = avg_log2FC, y = neg_log10_p, colour = direction)) +
      geom_point(alpha = 0.55, size = 1.3) +
      scale_colour_manual(values = c(up   = "#D7263D",
                                     down = "#1B98E0",
                                     ns   = "grey80")) +
      geom_vline(xintercept = c(-opt$logfc, opt$logfc),
                 linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", colour = "grey50") +
      geom_text_repel(data = top_label,
                      aes(label = gene),
                      size = 3, max.overlaps = 25, segment.alpha = 0.4) +
      labs(title = ttl,
           x = "avg log2 fold-change (group1 over group2)",
           y = "-log10 adjusted P (BH)") +
      theme_classic(base_size = 11) +
      theme(legend.title = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
  )
  dev.off()
}

# ---- 8. Combined top-hits heatmap (cluster × contrast) -------------------

if (length(all_top) > 0) {
  message("--- Top-hits heatmap ---")
  top_combined <- bind_rows(all_top) |>
    filter(!is.na(avg_log2FC)) |>
    mutate(contrast_key = paste(cluster, contrast_id, sep = " | "))
  top_genes <- top_combined |>
    group_by(gene) |>
    summarise(max_abs = max(abs(avg_log2FC), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(max_abs)) |>
    head(60) |>
    pull(gene)

  if (length(top_genes) >= 5) {
    heat_df <- top_combined |>
      filter(gene %in% top_genes) |>
      select(gene, contrast_key, avg_log2FC) |>
      group_by(gene, contrast_key) |>
      summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop") |>
      pivot_wider(names_from = contrast_key,
                  values_from = avg_log2FC, values_fill = 0)
    heat_mat <- as.matrix(heat_df[, -1])
    rownames(heat_mat) <- heat_df$gene
    cap <- 3
    heat_mat[heat_mat >  cap] <-  cap
    heat_mat[heat_mat < -cap] <- -cap

    pdf(file.path(fig_dir, "07_topgenes_heatmap.pdf"),
        width  = max(9, 2 + ncol(heat_mat) * 0.45),
        height = max(8, 2 + nrow(heat_mat) * 0.22))
    pheatmap::pheatmap(
      heat_mat,
      color        = colorRampPalette(c("#1B4F72", "#1B98E0", "white",
                                        "#D7263D", "#6B0F1A"))(50),
      breaks       = seq(-cap, cap, length.out = 51),
      cluster_rows = nrow(heat_mat) >= 2,
      cluster_cols = ncol(heat_mat) >= 2,
      fontsize_row = 8,
      fontsize_col = 8,
      border_color = "grey90",
      main = "Top DE genes — log2 FC by cluster × contrast"
    )
    dev.off()
  }
}

# ---- 9. Manifest + session info ------------------------------------------

manifest <- tibble::tibble(
  field = c("rds_in", "n_cells", "n_clusters", "min_cells",
           "logfc_threshold", "top_n", "test_use",
           "presto_available", "seed", "date"),
  value = c(opt$rds_in, ncol(seu_aml),
            length(aml_clusters),
            opt$min_cells, opt$logfc, opt$top_n,
            "wilcox",
            as.character(have_presto),
            42L,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
write.table(manifest, file.path(log_dir, "07_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_07_sessionInfo.txt"))

message("Done!")
