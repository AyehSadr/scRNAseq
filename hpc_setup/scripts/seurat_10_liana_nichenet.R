#!/usr/bin/env Rscript
# =============================================================================
# AML CloneTracker XP — Analysis 10: HS-5 ↔ AML cell-cell communication
# -----------------------------------------------------------------------------
# Project:      AML_Cellecta  (conventional vs co-culture, Ara-C response, S34)
# Author:       Ayeh Sadr
# Created:      2026-05-25
# Last update:  2026-05-25
# Input:        --rds_in     FULL Seurat .rds (NOT the AML-only subset) — must
#                            contain HS-5 stromal cells, the validated AML
#                            clusters in cell_type_final, and metadata columns
#                            stroma, araC, sample_id.
#               --out_dir    Output directory.
#               --imprint    Path to YAML / TSV with the niche-imprint gene set
#                            from §5.2 (default: hard-coded list below).
#               --nichenet_data_dir   Directory holding the three NicheNet
#                            prior files (auto-downloaded on first run).
# Output:       <out_dir>/figures/10a_liana_dotplot.pdf                 top LR pairs
#               <out_dir>/figures/10b_nichenet_ligand_activity.pdf      ligand ranking
#               <out_dir>/figures/10c_nichenet_ligand_target.pdf        ligand→target heatmap
#               <out_dir>/figures/10d_nichenet_ligand_receptor.pdf      ligand→receptor heatmap
#               <out_dir>/tables/10_liana_consensus.tsv                 LIANA consensus table
#               <out_dir>/tables/10_nichenet_<cluster>_ligand_activity.tsv  per-target ligand ranking
#               <out_dir>/tables/10_nichenet_top_ligands.tsv            union top-15 ligands
#               <out_dir>/logs/10_manifest.tsv
#               <out_dir>/logs/seurat_10_sessionInfo.txt
# Depends on:   Seurat (>=5.0), liana, nichenetr, OmnipathR, ggplot2, dplyr,
#               tibble, tidyr, optparse, yaml, circlize, ComplexHeatmap
#
# Notes:
#   * Cohort caveat -----------------------------------------------------------
#     Single patient (SRAML10). LIANA / NicheNet do NOT need replication —
#     they work on cell × gene expression patterns and ligand-target priors.
#     The output ligand ranking is patient-specific (because HS-5 expression
#     and AML cluster composition are patient-specific), but the methodology
#     is valid on n=1. This is the single highest-value next analysis given
#     the cohort constraint.
#   * Why the full object, not the AML-only subset --------------------------
#     The whole point is AML ↔ HS-5 signalling. Co-culture wells contain
#     both populations; the script filters to those wells and uses the
#     existing cell_type_final labels (including "HS-5 stromal (cell line)").
#   * Why subset to co-culture wells only -----------------------------------
#     Conventional wells contain no HS-5 cells, so LIANA / NicheNet have
#     no source to draw from. Restricting to coc cells removes the trivial
#     contrast and concentrates on the real signalling environment.
#   * LIANA design ------------------------------------------------------------
#     Runs the consensus pipeline (CellPhoneDB, Connectome, NATMI,
#     SingleCellSignalR, logfc) and aggregates to a robust LR ranking.
#     Output filtered to source = HS-5 → target = each AML cluster.
#   * NicheNet design ---------------------------------------------------------
#     For each AML target cluster, defines ligands as genes expressed by
#     HS-5 above background, defines genes-of-interest as the §5.2 niche-
#     imprint gene set, computes ligand-activity Pearson correlation, and
#     returns the ranked top ligands plus their predicted target genes.
#   * NicheNet prior files (auto-downloaded if missing) ----------------------
#     ligand_target_matrix_nsga2r_final.rds   — ligand × target gene scores
#     lr_network_human_21122021.csv           — curated LR pairs
#     weighted_networks_nsga2r_final.rds      — combined intracellular signalling
#     All from Zenodo doi:10.5281/zenodo.7074291.
#   * Seed = 42.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(optparse)
})

set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)

have_liana    <- requireNamespace("liana",       quietly = TRUE)
have_nichenet <- requireNamespace("nichenetr",   quietly = TRUE)
if (!have_liana) {
  message("NOTE: liana not installed.",
          " Install via remotes::install_github('saezlab/liana').")
}
if (!have_nichenet) {
  message("NOTE: nichenetr not installed.",
          " Install via remotes::install_github('saeyslab/nichenetr').")
}

# ---- 0. CLI ---------------------------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c("--rds_in"),  type = "character",
                     help = "Full Seurat .rds (with HS-5 retained)")
parser <- add_option(parser, c("--out_dir"), type = "character",
                     help = "Output directory")
parser <- add_option(parser, c("--imprint"), type = "character",
                     default = NULL,
                     help = "TSV / YAML with niche-imprint genes (one column 'gene'). Uses §5.2 defaults if missing.")
parser <- add_option(parser, c("--sender"), type = "character",
                     default = "HS-5 stromal (cell line)",
                     help = "Sender cluster name (must match cell_type_final)")
parser <- add_option(parser, c("--nichenet_data_dir"), type = "character",
                     default = NULL,
                     help = "Directory holding NicheNet prior .rds / .csv files")
parser <- add_option(parser, c("--top_n_liana"),    type = "integer", default = 30L,
                     help = "Top N LR pairs to show in LIANA dotplot")
parser <- add_option(parser, c("--top_n_ligands"),  type = "integer", default = 15L,
                     help = "Top N ligands per receiver for NicheNet outputs")
opt <- parse_args(parser)

for (req in c("rds_in", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required --", req)
}
if (!file.exists(opt$rds_in)) stop("Input RDS not found: ", opt$rds_in)

tab_dir <- file.path(opt$out_dir, "tables")
fig_dir <- file.path(opt$out_dir, "figures")
log_dir <- file.path(opt$out_dir, "logs")
nn_dir  <- if (is.null(opt$nichenet_data_dir)) {
  file.path(opt$out_dir, "nichenet_data")
} else {
  opt$nichenet_data_dir
}
for (d in c(tab_dir, fig_dir, log_dir, nn_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Seurat 10: LIANA + NicheNet ===")
message("  rds_in            : ", opt$rds_in)
message("  out_dir           : ", opt$out_dir)
message("  sender            : ", opt$sender)
message("  nichenet_data_dir : ", nn_dir)

# ---- 1. Default niche-imprint gene set (from §5.2) -----------------------

NICHE_IMPRINT_DEFAULT <- c(
  # Iron handling (Ara-C / oxidative response, persist under niche)
  "FTH1", "FTL",
  # NAD salvage / metabolic
  "NAMPT",
  # Redox sensor / fasting-state metabolic regulator
  "TXNIP",
  # MHC / IFN-induced
  "B2M", "IFITM3",
  # Calcium-binding / membrane
  "S100A11",
  # Survival / stem-state regulators (top-up in §5 niche contrasts)
  "EPB41L3", "SLC44A1", "PLSCR1", "AHRR", "TBL1X"
)

if (!is.null(opt$imprint) && file.exists(opt$imprint)) {
  message("--- Reading niche-imprint genes from: ", opt$imprint, " ---")
  if (grepl("\\.ya?ml$", opt$imprint, ignore.case = TRUE) &&
      requireNamespace("yaml", quietly = TRUE)) {
    niche_imprint <- yaml::read_yaml(opt$imprint)$genes
  } else {
    niche_imprint <- read.delim(opt$imprint, stringsAsFactors = FALSE)$gene
  }
  niche_imprint <- unique(niche_imprint[!is.na(niche_imprint) & nzchar(niche_imprint)])
} else {
  niche_imprint <- NICHE_IMPRINT_DEFAULT
}
message("--- Niche-imprint gene set (", length(niche_imprint), " genes) ---")
message("  ", paste(niche_imprint, collapse = ", "))

# ---- 2. Load full object + co-culture subset -----------------------------

message("--- Loading full Seurat object ---")
seu_full <- readRDS(opt$rds_in)
DefaultAssay(seu_full) <- "RNA"

if (!"data" %in% Layers(seu_full[["RNA"]])) {
  message("--- Normalising RNA counts ---")
  seu_full <- NormalizeData(seu_full, verbose = FALSE)
}
for (req_col in c("cell_type_final", "stroma", "sample_id")) {
  if (!req_col %in% colnames(seu_full@meta.data)) {
    stop("Required metadata column missing: ", req_col)
  }
}

if (!(opt$sender %in% levels(droplevels(factor(seu_full$cell_type_final))))) {
  stop("Sender cluster '", opt$sender, "' not found in cell_type_final.\n",
       "  Levels present: ",
       paste(levels(droplevels(factor(seu_full$cell_type_final))), collapse = ", "))
}

message("--- Subsetting to co-culture wells (HS-5 only present there) ---")
seu_coc <- seu_full[, seu_full$stroma]
seu_coc$cell_type_final <- droplevels(factor(seu_coc$cell_type_final))
message("  Cells: ", ncol(seu_coc))
message("  Clusters: ", paste(levels(seu_coc$cell_type_final), collapse = ", "))

aml_clusters <- setdiff(levels(seu_coc$cell_type_final), opt$sender)
aml_clusters <- aml_clusters[grepl("AML", aml_clusters)]
message("  AML receiver clusters: ", paste(aml_clusters, collapse = ", "))
if (length(aml_clusters) == 0) {
  stop("No AML receiver clusters found in the co-culture subset.")
}

Idents(seu_coc) <- seu_coc$cell_type_final

# ---- 3. LIANA — consensus ligand-receptor analysis -----------------------

if (have_liana) {
  message("--- Running LIANA (consensus across CellPhoneDB / NATMI / ...) ---")
  liana_res <- tryCatch(
    liana::liana_wrap(seu_coc,
                      method = c("natmi", "connectome",
                                 "logfc", "sca", "cellphonedb"),
                      resource = "Consensus",
                      verbose = FALSE),
    error = function(e) {
      message("LIANA error: ", e$message); NULL
    }
  )
  if (!is.null(liana_res)) {
    consensus <- liana::liana_aggregate(liana_res, verbose = FALSE) |>
      as_tibble() |>
      filter(source == opt$sender, target %in% aml_clusters) |>
      arrange(aggregate_rank)
    write.table(consensus,
                file.path(tab_dir, "10_liana_consensus.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("  LIANA consensus rows: ", nrow(consensus))

    # Top-N LR pairs per target cluster — dotplot
    top_pairs <- consensus |>
      group_by(target) |>
      slice_min(aggregate_rank, n = opt$top_n_liana) |>
      ungroup() |>
      mutate(lr = paste(ligand.complex, receptor.complex, sep = " → "))

    if (nrow(top_pairs) > 0) {
      pdf(file.path(fig_dir, "10a_liana_dotplot.pdf"),
          width  = max(8, 1.0 * length(aml_clusters) + 5),
          height = max(8, 0.18 * length(unique(top_pairs$lr)) + 3))
      print(
        ggplot(top_pairs,
               aes(x = target,
                   y = reorder(lr, -aggregate_rank),
                   colour = -log10(aggregate_rank),
                   size   = pmin(natmi.edge_specificity, 1, na.rm = TRUE))) +
          geom_point() +
          scale_colour_viridis_c(option = "magma", direction = -1,
                                 name = "-log10 rank") +
          scale_size_continuous(range = c(1, 5), name = "edge specificity") +
          labs(x = NULL, y = NULL,
               title = paste0("LIANA consensus: ", opt$sender, " → AML clusters")) +
          theme_classic(base_size = 9) +
          theme(axis.text.x = element_text(angle = 35, hjust = 1))
      )
      dev.off()
    }
  }
} else {
  message("--- Skipping LIANA (package not available) ---")
}

# ---- 4. NicheNet — ligand activity for niche-imprint genes ---------------

if (have_nichenet) {
  message("--- Loading NicheNet prior networks ---")

  fetch <- function(url, dest) {
    if (!file.exists(dest)) {
      message("  Downloading: ", basename(dest))
      tryCatch({
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)
      }, error = function(e) {
        stop("Failed to download NicheNet resource from Zenodo.\n",
             "  Error details: ", e$message, "\n",
             "  NOTE: If your Falcon compute node does not have outbound internet access,\n",
             "        please run this script once on a login node (which has internet) to cache the files,\n",
             "        or manually download the files and place them in the 'nichenet_data' folder:\n",
             "        - Source:      ", url, "\n",
             "        - Destination: ", dest, "\n")
      })
    }
    dest
  }
  ltm_url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  lr_url  <- "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
  wn_url  <- "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"
  ltm_path <- fetch(ltm_url, file.path(nn_dir, "ligand_target_matrix_nsga2r_final.rds"))
  lr_path  <- fetch(lr_url,  file.path(nn_dir, "lr_network_human_21122021.rds"))
  wn_path  <- fetch(wn_url,  file.path(nn_dir, "weighted_networks_nsga2r_final.rds"))

  ligand_target_matrix <- readRDS(ltm_path)
  lr_network           <- readRDS(lr_path)
  weighted_networks    <- readRDS(wn_path)

  # Expressed genes per cluster (Seurat helper from nichenetr)
  expressed_by <- function(seu_obj, cluster, pct = 0.10) {
    nichenetr::get_expressed_genes(cluster, seu_obj,
                                   pct = pct, assay_oi = "RNA")
  }

  message("--- Sender expressed genes (HS-5) ---")
  sender_expressed <- expressed_by(seu_coc, opt$sender)
  message("  ", length(sender_expressed), " expressed genes in sender")

  all_ligands <- unique(lr_network$from)
  ligands_of_interest <- intersect(sender_expressed, all_ligands)
  message("  Candidate sender ligands: ", length(ligands_of_interest))

  all_top_ligands <- list()
  for (recv in aml_clusters) {
    safe <- gsub("[^A-Za-z0-9]+", "_", recv)
    message("--- NicheNet receiver: ", recv, " ---")
    receiver_expressed <- expressed_by(seu_coc, recv)
    background <- intersect(receiver_expressed, rownames(ligand_target_matrix))
    geneset_oi <- intersect(niche_imprint, background)
    if (length(geneset_oi) < 3) {
      message("  Too few imprint genes expressed in ", recv, " (n=", length(geneset_oi), "). Skipping.")
      next
    }
    expressed_receptors <- intersect(receiver_expressed, unique(lr_network$to))
    potential_ligands <- lr_network |>
      filter(from %in% ligands_of_interest,
             to   %in% expressed_receptors) |>
      pull(from) |>
      unique()
    if (length(potential_ligands) < 3) {
      message("  Too few potential ligands for ", recv, ". Skipping.")
      next
    }

    activity <- nichenetr::predict_ligand_activities(
      geneset             = geneset_oi,
      background_expressed_genes = background,
      ligand_target_matrix       = ligand_target_matrix,
      potential_ligands          = potential_ligands)
    activity <- activity |>
      arrange(desc(aupr_corrected)) |>
      mutate(rank = row_number(),
             receiver = recv)
    write.table(activity,
                file.path(tab_dir,
                          sprintf("10_nichenet_%s_ligand_activity.tsv", safe)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    top_lig <- head(activity$test_ligand, opt$top_n_ligands)
    all_top_ligands[[recv]] <- top_lig

    # Per-receiver ligand activity barplot
    pdf(file.path(fig_dir, sprintf("10b_nichenet_%s_ligand_activity.pdf", safe)),
        width = 7, height = max(4, 0.25 * length(top_lig) + 2))
    print(
      activity |>
        head(opt$top_n_ligands) |>
        ggplot(aes(x = aupr_corrected,
                   y = reorder(test_ligand, aupr_corrected))) +
        geom_col(fill = "#1B2A4A") +
        labs(x = "AUPR (corrected) — ligand activity score",
             y = NULL,
             title = paste0("NicheNet ligand activity → ", recv)) +
        theme_classic(base_size = 10)
    )
    dev.off()
  }

  # Union top ligands across receivers — heatmap
  if (length(all_top_ligands) > 0) {
    union_ligands <- unique(unlist(all_top_ligands))
    lt_sub <- ligand_target_matrix[, intersect(union_ligands, colnames(ligand_target_matrix)),
                                   drop = FALSE]
    target_genes <- intersect(niche_imprint, rownames(lt_sub))
    lt_sub <- lt_sub[target_genes, , drop = FALSE]

    pdf(file.path(fig_dir, "10c_nichenet_ligand_target.pdf"),
        width  = max(8, 0.4 * ncol(lt_sub) + 3),
        height = max(5, 0.3 * nrow(lt_sub) + 2))
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(
        lt_sub,
        color        = colorRampPalette(c("white", "#D7263D", "#6B0F1A"))(50),
        cluster_rows = nrow(lt_sub) >= 2,
        cluster_cols = ncol(lt_sub) >= 2,
        border_color = "grey90",
        fontsize_row = 9, fontsize_col = 9,
        main = "NicheNet ligand × niche-imprint target scores"
      )
    }
    dev.off()

    # Ligand × receptor heatmap (which receptors on AML carry each ligand's signal)
    union_receptors <- lr_network |>
      filter(from %in% union_ligands) |>
      pull(to) |>
      unique()
    lr_long <- lr_network |>
      filter(from %in% union_ligands, to %in% union_receptors) |>
      distinct(from, to) |>
      mutate(value = 1) |>
      pivot_wider(names_from = to, values_from = value, values_fill = 0) |>
      column_to_rownames("from")

    if (nrow(lr_long) > 0 && ncol(lr_long) > 0) {
      pdf(file.path(fig_dir, "10d_nichenet_ligand_receptor.pdf"),
          width  = max(8, 0.3 * ncol(lr_long) + 3),
          height = max(5, 0.3 * nrow(lr_long) + 2))
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(
          as.matrix(lr_long),
          color        = colorRampPalette(c("white", "#1B2A4A"))(50),
          cluster_rows = nrow(lr_long) >= 2,
          cluster_cols = ncol(lr_long) >= 2,
          border_color = "grey90",
          fontsize_row = 9, fontsize_col = 9,
          main = "NicheNet ligand × known receptor (binary)"
        )
      }
      dev.off()
    }

    top_long <- bind_rows(
      lapply(names(all_top_ligands), function(r) {
        tibble(receiver = r, rank = seq_along(all_top_ligands[[r]]),
               ligand = all_top_ligands[[r]])
      })
    )
    write.table(top_long,
                file.path(tab_dir, "10_nichenet_top_ligands.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
} else {
  message("--- Skipping NicheNet (package not available) ---")
}

# ---- 5. Manifest + session info ------------------------------------------

manifest <- tibble::tibble(
  field = c("rds_in", "n_cells_coc", "sender", "n_receivers",
           "n_imprint_genes", "have_liana", "have_nichenet",
           "seed", "date"),
  value = c(opt$rds_in, ncol(seu_coc), opt$sender, length(aml_clusters),
            length(niche_imprint), as.character(have_liana),
            as.character(have_nichenet), 42L,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
write.table(manifest, file.path(log_dir, "10_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "seurat_10_sessionInfo.txt"))

message("Done!")
