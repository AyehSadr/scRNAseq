# scRNA-seq Analysis Brief — ICR Project S34

*Project Manager brief, received 2026-05-04. Source file: `scRNAseq_Bioinformatics_Instructions_S34.docx`. This is a faithful Markdown transcription so the same content is searchable, diffable, and importable to Notion.*

**Project:** Cardiff University · 3 AML patients (SRAML7, SRAML10, SRAML13) · 12 fixed samples · March 2026 · S. Rizzo / A. Tonks / J. Khorashad.

---

## 1. Sample manifest (n=12, all carry CloneTracker XP barcodes)

| Tube | Patient | HS5 stroma | AraC | Cells (×10⁴) | Condition |
|---|---|---|---|---|---|
| 3  | SRAML7  | – | – | 56  | Untreated, no stromal |
| 4  | SRAML7  | – | + | 38  | AraC-treated, no stromal |
| 7  | SRAML7  | + | – | 14  | Untreated, +stromal |
| 8  | SRAML7  | + | + | 14  | AraC-treated, +stromal |
| 11 | SRAML10 | – | – | 28  | Untreated, no stromal |
| 12 | SRAML10 | – | + | 17  | AraC-treated, no stromal |
| 15 | SRAML10 | + | – | 7.5\* | Untreated, +stromal — **low input, use with caution** |
| 16 | SRAML10 | + | + | 6.75\* | AraC-treated, +stromal — **low input, use with caution** |
| 17 | SRAML13 | – | – | 81  | Untreated, no stromal |
| 18 | SRAML13 | – | + | 40  | AraC-treated, no stromal |
| 19 | SRAML13 | + | – | 52.5 | Untreated, +stromal |
| 20 | SRAML13 | + | + | 18.2 | AraC-treated, +stromal |

---

## 2. Analysis 1 — AraC effect (within each niche)

Pairwise comparisons:

- SRAML7: 3 vs 4 (–HS5), 7 vs 8 (+HS5)
- SRAML10: 11 vs 12 (–HS5), 15 vs 16 (+HS5)
- SRAML13: 17 vs 18 (–HS5), 19 vs 20 (+HS5)

Questions to answer per pair:
- Is any cell population / cluster enriched in +AraC vs −AraC, separately for –HS5 and +HS5?
- Is there a transcriptional signature present in +AraC that's absent in −AraC?

Methods: differential cluster abundance (Milo or scCODA), pseudo-bulk DGE (DESeq2 / edgeR; adj. p < 0.05, |log2FC| ≥ 0.5). Deliver volcano plots and top-DEG heatmaps.

---

## 3. Analysis 2 — HS5 stromal niche effect

Pairwise comparisons:

- Untreated leukaemia, ±stroma: 3 vs 7, 11 vs 15, 17 vs 19
- AraC-treated leukaemia, ±stroma (niche-mediated resistance): 4 vs 8, 12 vs 16, 18 vs 20

For each comparison: identify signatures that differ in +HS5 vs −HS5.

Methods: gene set enrichment + pathway analysis using MSigDB Hallmark, KEGG, Reactome, GO BP (fgsea / clusterProfiler / gseapy). Deliver top 20 pathways per comparison as ranked bubble plots (NES vs −log10 adj. p).

---

## 4. Analysis 3 — CloneTracker barcode detection and clonal tracking

**⚠️ Pre-requisite (verbatim from the brief):** "Confirm with ICR that CloneTracker barcode sequences were captured in the sequencing library (feature barcoding / Cell Ranger multi). If not, a separate targeted amplicon-seq run may be required."

This is the open chemistry question we've been flagging. The PM has formally acknowledged it. Email Floriana Manodoro before doing more work on this analysis branch.

Detection metrics (per sample):
- % cells with at least one CloneTracker barcode
- UMI distribution per detected barcode
- Unique barcode count

Clonal overlap:
- Same barcode in −AraC and +AraC of the same patient
  - No-stromal pairs: 3↔4, 11↔12, 17↔18
  - Stromal pairs: 7↔8, 15↔16, 19↔20
- Same barcode across all four conditions of a patient
  - SRAML7: {3, 4, 7, 8} · SRAML10: {11, 12, 15, 16} · SRAML13: {17, 18, 19, 20}
- Report Jaccard similarity, Sankey / alluvial plots, proportion of barcodes surviving AraC.

Clone-specific gene expression:
- For cells sharing the same barcode in −AraC vs +AraC: has expression changed? (DGE on same-barcode cells, per niche.)
- For barcodes spanning all four conditions: trajectory / pseudotime (Monocle3 or scVelo).

---

## 5. Analysis 4 — AML differentiation staging and cell cycle

Cluster all cells by differentiation stage (LSC → MLP → GMP → ProMono → Mono → cDC) and cell cycle phase.

Reference datasets to use:

| Reference | Source | Use |
|---|---|---|
| **van Galen et al. 2019 (Cell)** | AML single-cell hierarchy | **Primary** label transfer (Seurat::TransferData() or scArches/scANVI) |
| **LSC17 signature (Ng et al. 2016, Nature)** | 17-gene LSC signature | Score per cell with AddModuleScore() / AUCell |
| **Tirosh et al. 2016** | G1/S and G2/M gene lists | Cell cycle regression and phase assignment |
| **BeatAML (Tyner et al. 2018, Nature)** | Bulk AML transcriptional subtypes | Cross-reference |

Per-condition / per-patient deliverables: differentiation-state proportions, LSC17 score distributions, cell cycle phase distributions.

Key biological question: does AraC treatment or HS5 co-culture shift the differentiation composition (e.g., quiescent LSC-like cells enriched in surviving / niche-protected populations)?

---

## 6. Pre-processing (uniform across all 12 samples)

| Step | Specification |
|---|---|
| Alignment | Cell Ranger v8+; GRCh38 reference; **include pre-mRNA** for unspliced reads (`--include-introns=true`, default in CR8) |
| QC | genes/cell 200–6,000; UMI ≥ 500; MT% < 20% (AraC samples may have elevated MT% — report distribution before hard cut-off); DoubletFinder or Scrublet per sample |
| Normalisation | Lib-size to 10,000 counts → log1p; top 2,000–3,000 HVGs (VST) |
| Clustering | PCA (30–50 PCs) → SNN k=20 → Leiden (resolution 0.4–1.0) → UMAP |
| Integration | Harmony or scVI on patient ID for cross-patient analyses; **uncorrected data for within-patient pairwise** comparisons (Analyses 1 and 2) |
| Low-input samples | Tubes 15, 16 (SRAML10, +HS5; ~6,750–7,500 cells): relax min. gene threshold to ≥150 if needed; flag in all reports |

---

## 7. Deliverables checklist

| Deliverable | Format |
|---|---|
| QC report | per-sample metrics table; violin plots (genes/cell, UMI, MT%); flag any failing sample |
| UMAP figures | per-patient and integrated; coloured by condition / cluster / patient / diff state / cell cycle |
| Analysis 1 & 2 | DGE tables (CSV); volcano plots; GSEA results (CSV); pathway bubble plots |
| Analysis 3 | barcode detection table; Sankey plots; clonal overlap matrices; clone-specific DGE tables |
| Analysis 4 | differentiation-state UMAP + proportion bar plots; LSC17 score plots; cell cycle plots |
| Data objects | annotated RDS (Seurat) or H5AD (Scanpy) per patient + integrated |

---

## How this changes the pipeline plan

1. **Analysis 3 is now formally pre-requisite-gated** on confirming feature barcoding with ICR. Email Floriana Manodoro this week.
2. **Analyses 1, 2, and 4 are not blocked** on the chemistry question — they only need the GEX library which we already have. The HPC pipeline as scaffolded handles all of them.
3. **New tool requirements** to add to the conda envs (will pick up on next rebuild):
   - R: DoubletFinder, DESeq2, edgeR, fgsea, clusterProfiler, miloR, monocle3, AUCell
   - Python: scrublet, gseapy, scvi-tools, scanpy[leiden]
4. **New reference datasets to download** to scratch (one-off, ~10 GB total): van Galen 2019 processed object (CellxGene or Broad Single Cell Portal), LSC17 gene list, MSigDB collections (msigdbr R package handles this), Tirosh cell cycle (built into Seurat).
5. **Sample sheet** updated to flag tubes 15 and 16 as low-input.
