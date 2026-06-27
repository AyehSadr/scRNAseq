# ICR Project S34 — scRNA-seq Analysis Plan

**Cardiff University · 3 patients (SRAML7, SRAML10, SRAML13) · 12 fixed samples · CloneTracker XP**
*Plan derived from `scRNAseq_Bioinformatics_Instructions_S34.docx` (Rizzo / Tonks / Khorashad, March 2026).*

---

## Phase 0 — Preflight (do before any compute)

**Confirm with ICR:**
- Were CloneTracker XP barcode sequences captured in the GEX library (feature-barcoding read structure in Cell Ranger multi config)? Current understanding: enrichment library was *not* run → expect <10% recovery → raise the option of a targeted amplicon-seq follow-up *now*, not after Analysis 3 fails.
- Cell Ranger multi version, reference build, pre-mRNA inclusion (brief specifies GRCh38, pre-mRNA on for unspliced reads — this is also a prereq for scVelo trajectory work in Analysis 3).

**Lock the sample manifest** (12 samples, 4 conditions × 3 patients) with QC flags for low-input tubes 15 and 16 (SRAML10, +HS5; ~6,750–7,500 cells).

---

## Phase 1 — Pre-processing (uniform, all 12 samples)

| Step | Spec from brief | Notes / deviations from current pipeline |
|---|---|---|
| Alignment | Cell Ranger v8+, GRCh38, pre-mRNA on | Standard |
| QC thresholds | 200–6,000 genes, UMI ≥ 500, MT < 20% | Relax min genes to ≥150 for tubes 15, 16 |
| MT% handling | Report distribution *before* applying hard cut-off (AraC samples may have elevated MT%) | Your current script applies the cut-off directly — modify to plot first, threshold second |
| Doublets | DoubletFinder or Scrublet, per sample | You're using `scDblFinder` — defensible (Germain 2021 benchmark favours it), keep but note in methods |
| Normalisation | Library-size to 10k → log1p, 2,000–3,000 VST HVGs | **Switch from SCTransform** to match brief |
| Clustering | PCA 30–50 → SNN k=20 → Leiden 0.4–1.0 → UMAP | Add `clustree` to justify resolution choice |
| Integration | Harmony/scVI on `patient` for cross-patient only; **uncorrected within-patient** | Drop the `group.by.vars="sample_id"` Harmony call from the per-patient pipeline |

**Output:** per-sample QC report (violins, metrics table, flag failures) + filtered + normalised + clustered RDS per patient.

---

## Phase 2 — Per-patient object construction (×3 patients)

For SRAML7, SRAML10, SRAML13 separately:
- Merge 4 samples, **no batch correction** (sample_id is confounded with condition).
- PCA → SNN → Leiden, UMAP.
- Save annotated RDS: `<patient>_uncorrected.rds`. This is the substrate for Analyses 1 and 2.

---

## Phase 3 — Analysis 1: AraC effect within each niche

**Comparisons** (per patient, 6 total):
- No-stromal: 3 vs 4 (SRAML7), 11 vs 12 (SRAML10), 17 vs 18 (SRAML13)
- +Stromal:    7 vs 8 (SRAML7), 15 vs 16 (SRAML10), 19 vs 20 (SRAML13)

**Per comparison:**
1. **Differential abundance** — Milo (preferred for 2-sample comparisons; scCODA is the alternative if you want compositional null handling). Output: neighbourhood-level beeswarm + DA-coloured UMAP.
2. **Pseudo-bulk DGE** — aggregate counts per-sample-per-cluster, run DESeq2 (`muscat::aggregateData` → `pbDS` is the cleanest path). Thresholds: adj.p < 0.05, |log2FC| ≥ 0.5. Output: per-cluster DE tables + global volcano.
3. **Top-DEG heatmaps** per cluster.
4. **GSEA** with fgsea on MSigDB Hallmark, KEGG, Reactome, GO-BP. Output: ranked bubble plots (NES vs −log10 adj.p), top 20 pathways per comparison.
5. **Clusters unique to treated** — flag any cluster where +AraC frequency >5× −AraC.

---

## Phase 4 — Analysis 2: HS5 stromal niche effect

Same machinery as Phase 3, different comparisons:
- Untreated: 3 vs 7, 11 vs 15, 17 vs 19 (effect of niche on baseline leukaemia)
- AraC-treated: 4 vs 8, 12 vs 16, 18 vs 20 — **niche-mediated AraC resistance**, the headline question of the project.

The treated-pair comparison should get extra interpretive attention: pathway-level convergence across patients (NES sign agreement on Hallmark gene sets) is the kind of result that goes into a figure 1.

---

## Phase 5 — Cross-patient integration

Build a 12-sample integrated object for shared cell-type vocabulary and Analysis 4:
- Merge all 12 samples.
- Harmony or scVI on `patient`.
- Joint Leiden clustering, joint UMAP.
- Save as `S34_integrated.rds` / `.h5ad`.

This object is *not* used for Analyses 1/2 (those stay uncorrected per-patient), but it is required for label transfer in Analysis 4 and for any cross-patient consensus claims.

---

## Phase 6 — Analysis 4: Differentiation staging + cell cycle

1. **Reference assembly:**
   - van Galen et al. 2019 (GSE116256): LSC, MLP, GMP, ProMono, Mono, cDC labels — primary reference.
   - LSC17 signature (Ng et al. 2016) — score per cell.
   - Tirosh et al. 2016 G1/S/G2M lists.
   - Optional: BeatAML (Tyner 2018) for transcriptional subtype cross-reference.

2. **Label transfer** onto integrated object: Seurat `TransferData` or scANVI (scvi-tools). Report mapping confidence per cell; mask cells below threshold.

3. **Scoring:**
   - LSC17 via `AddModuleScore` or AUCell (AUCell is more robust to expression depth).
   - Cell cycle via `CellCycleScoring`; phase assignment.

4. **Reporting:** per patient × per condition → stacked bar of differentiation-state proportions; LSC17 score distribution (violin, with stat test); cell cycle phase composition.

5. **Headline test:** does AraC and/or +HS5 enrich quiescent LSC-like cells (LSC17-high, G1-arrested)? This is the LSC-protection / niche-resistance hypothesis. Test with a pseudo-bulk proportion model (Milo at the cell-state level, or beta-binomial GLM).

---

## Phase 7 — Analysis 3: CloneTracker barcoding

> **Highest-risk component.** Without the enrichment library, recovery is expected to be poor. The plan below assumes the barcodes are detectable in the GEX library at all; if Phase 0 confirms they aren't, this whole phase becomes "request amplicon-seq from ICR" and the rest is parked.

1. **Detection report:** per sample — % cells with ≥1 barcode, UMI distribution per barcode, unique barcode count.
2. **Clonal overlap (within-patient pairs):**
   - No-stromal AraC pairs: 3↔4, 11↔12, 17↔18.
   - +Stromal AraC pairs: 7↔8, 15↔16, 19↔20.
   - Metric: Jaccard similarity. Visual: Sankey / alluvial.
   - Report: proportion of barcodes surviving AraC.
3. **Cross-condition clones (4 conditions per patient):** shared barcode set, table + alluvial across {3,4,7,8} / {11,12,15,16} / {17,18,19,20}.
4. **Clone-specific DGE:**
   - For barcodes shared between −AraC and +AraC: pseudo-bulk DGE on same-barcode cells.
   - Repeat within each niche condition (−HS5 and +HS5).
   - **Caveat:** if recovery is <10%, expect most barcodes to have <5 cells per condition → underpowered. Mitigation: pool clones by trajectory cluster.
5. **Trajectory:** Monocle3 or scVelo on shared multi-condition clones. scVelo requires unspliced reads — only works if Phase 0 pre-mRNA inclusion is confirmed.

---

## Phase 8 — Deliverables (per brief)

| Item | Format |
|---|---|
| QC report | HTML or PDF, per-sample table + violin plots, failure flags |
| UMAP figures | Per-patient (uncorrected) + integrated; coloured by condition / cluster / patient / differentiation state / cycle phase |
| Analysis 1 + 2 | DGE CSVs, volcano plots, GSEA CSVs, pathway bubble plots |
| Analysis 3 | Barcode detection table, Sankeys, clonal overlap matrices, clone-specific DGE CSVs |
| Analysis 4 | Differentiation-state UMAP + proportion bars, LSC17 score plots, cell cycle plots |
| Data objects | `*.RDS` (Seurat) per patient + integrated `*.h5ad` (Scanpy) |

Plus: pin everything to a `renv` lockfile (R) and `environment.yml` (Python), Snakemake or Nextflow wrapper for reproducibility, version-controlled in git.

---

## Suggested execution order (priority)

1. **Phase 0 + 1** — get pre-processing uniform across all 12 samples; this is foundational and unblocks everything.
2. **Phase 2 + 3 + 4** — the headline biology (AraC effect, niche-resistance). These are what the grant cares about.
3. **Phase 5 + 6** — adds biological interpretation and cross-patient consensus.
4. **Phase 7** — CloneTracker. Run the detection report early (after Phase 1) so you know whether to escalate the amplicon-seq request, but full clonal-tracking analysis is downstream.
5. **Phase 8** — deliverables compiled rolling, not at the end.

---

## Gaps in the current pipeline (vs this plan)

Concrete edits to migrate your existing SRAML10 script toward the brief:

- **Drop** `RunHarmony(group.by.vars="sample_id")` from the per-patient pipeline.
- **Replace** `SCTransform(vst.flavor="v2")` with `NormalizeData` (scale.factor=1e4) + `FindVariableFeatures(method="vst", nfeatures=3000)` + `ScaleData(vars.to.regress="pct_mt")`.
- **Replace** Seurat `FindMarkers`/`FindAllMarkers` Wilcoxon with pseudo-bulk DESeq2 via `muscat` or Libra.
- **Add** Milo for differential abundance (Analysis 1/2 part 1).
- **Add** fgsea + MSigDBR for pathway analysis.
- **Add** van Galen reference download + Seurat `TransferData` for Analysis 4.
- **Add** LSC17 scoring (`AUCell`).
- **Modify** the QC step to plot MT% distribution before thresholding (per brief).
- **Restructure** the script as one pipeline-per-patient (×3) + one cross-patient integration script, rather than a single monolithic script.
