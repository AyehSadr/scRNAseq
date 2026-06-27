# AML CloneTracker — Part 3 Analysis Plan

**Project:** AML_Cellecta (ICR S34, 12 AML samples, CloneTracker XP, conventional vs HS-5 co-culture, +/- Ara-C)
**Author:** Ayeh Sadr
**Created:** 2026-05-23
**Source:** `docs/Further_bioinformatic_analysis_part3.docx`

This plan operationalises the 13 questions in *Further_bioinformatic_analysis_part3.docx* on top of what we already have: an integrated Seurat object (`seu_merged`) with marker-based cell-type calls (`cell_type` metadata column) from `seurat_01_general.R` / `seurat_01b_annotate_clusters.R`, and Cellecta barcode QC + clonal overlap from the `seurat_02*` chain.

Each section names the script that will deliver it, the inputs it needs, and the immediate deliverable. The numbering follows the docx; the **execution order** at the bottom is the order to actually run them.

---

## 0. Split AML compartment from Normal / Uncertain / CellLine

**Script:** `seurat_04c_subset_aml.R` (SLURM via `34c_seurat_04c_subset_aml.sbatch`)
**Inputs:** full object after `seurat_04b` (with `is_aml`, `is_aml_category`, `cell_type_final`).
**What it does:** subsets to `is_aml == TRUE`, drops HS-5 / B / T / CD16+ Mono / eosinophils, re-runs HVG + PCA + UMAP on the AML compartment so the embedding reflects AML-internal variance. Cluster identities (`cell_type_final`) are NOT recomputed.
**Why:** every downstream biological question (§4–§7) is about AML biology. Keeping Normal and stromal cells in the object dilutes pseudobulk DE with trivial lineage differences, distorts MiloR neighbourhoods, and corrupts trajectory analysis. The §2 inferCNV step has already run on the full object and §8 LIANA *must* keep the full object because it needs HS-5.
**Output object:** `AML_SRAML10_aml_only.RDS` with reductions `pca_aml`, `umap` (AML-only — primary), and `umap_full` (subset of the global UMAP, kept for cross-reference).
**Downstream consumers:**
- AML-only: `seurat_05_stress_cycling.R`, `seurat_07_pseudobulk_de.R`, `seurat_09_milo.R`, trajectory.
- Full object: `seurat_04_infer_cnv.R` (done), `seurat_10_liana.R`.

## 1. Refine unclassified niche-related and ABCB5+ identity

**Script:** `seurat_03_annotate_reference.R` + `signatures.yml`
**Inputs:** existing Seurat object, van Galen 2019 + Azimuth BM reference (`01_download_references.sh`).
**What it does:** SingleR vs HumanPrimaryCellAtlas (broad) and van Galen 2019 (AML-specific); Azimuth label transfer against Hao/Stuart BM CITE-seq; UCell scoring of LSC17, Eppert HSC, Cheung-Rando + Laurenti quiescence, Baryawno + Tikhonova niche signatures, and the Ara-C / MDR / ALDH resistance lists.
**Deliverables:** `tables/cell_type_vs_singler_confusion.tsv`, `tables/ucell_scores_by_cluster.tsv`, UMAPs and a per-cluster signature heatmap.

## 2. Is CD16+ Mono/Macrophage AML or normal?

**Script:** `seurat_04_infer_cnv.R` (SLURM via `34_seurat_04_infer_cnv.sbatch`)
**Inputs:** Seurat object after seurat_03; inferCNV hg38 gene order file; B + T cells as the diploid reference.
**What it does:** inferCNV in subcluster mode with HMM; per-cell CNV burden score; heuristic CD16+ call by comparing to the reference vs known-AML clusters.
**Deliverables:** chromosome heatmap (`infercnv.pdf`), `tables/per_cell_cnv_score.tsv`, `tables/cd16_classification.tsv`, CNV-score UMAP and violin.
**Companion lines of evidence (separate, optional):** Vartrix if patient mutations are known; PAGA connectivity check from §3.

## 3. Trajectory analysis for conversion hypotheses

**Script:** `seurat_05_trajectory.py` (Scanpy/Scvelo) + `seurat_05b_cellrank.py`
**Inputs:** the integrated Seurat object exported to AnnData (`as.SingleCellExperiment` → `sceasy::convertFormat`).
**What it does:** PAGA on the KNN graph (connectivity sanity check), scVelo dynamical model (direction), CellRank (fate probabilities).
**Targeted edges:** niche↔quiescent, niche↔CD16+, GMP-like↔ABCB5+.
**Deliverables:** PAGA graph, velocity-stream UMAP, terminal-state probabilities table.

## 4. Stress signalling — untreated conventional vs untreated co-culture

**Script:** `seurat_06_stress_scoring.R`
**Inputs:** Seurat object; `signatures.yml` (Hallmark sets, ATF4 targets, van den Brink dissociation).
**What it does:** UCell scoring per cell, Wilcoxon test with Cliff's delta effect size per cluster × programme.
**Deliverables:** effect-size heatmap (cluster × stress programme), per-cell score table.
**Framing:** if co-culture is not generally more stressful but quiescent cells are still enriched, that argues *against* H1 (proliferation arrest from harsh environment) and *toward* H4 (active MSC-induced quiescence).

## 5. Differential expression (within-patient now / pseudobulk later)

**Cohort caveat first.** The current object is one patient (SRAML10) × 4 samples
(one per condition × treatment). There is no biological replication within any
condition × treatment group, so pseudobulk DE with `muscat` / DESeq2 is *not*
feasible — DESeq2 cannot estimate dispersion at n = 1 per group. §5 is therefore
split into a within-patient hypothesis-generation step now and a proper
pseudobulk step gated on cohort expansion.

### 5a. Within-patient single-cell DE — now

**Script:** `seurat_07_de_within_patient.R` (SLURM via `37_*.sbatch`)
**Input:** `AML_SRAML10_aml_only.RDS` (already singlet-only via scDblFinder in §1).
**Test:** Seurat `FindMarkers`, Wilcoxon (presto-accelerated), BH-FDR per table.
**Contrasts per AML cluster:**
1. Coc_Unt vs Conv_Unt — *niche imprint, baseline*
2. Coc_AraC vs Conv_AraC — *niche-protected residual disease (clinical core)*
3. Conv_AraC vs Conv_Unt — *Ara-C effect, no niche*
4. Coc_AraC vs Coc_Unt — *Ara-C effect, with niche*

**Special cluster-vs-cluster contrasts:**
- S1: *Niche-stressed primitive AML vs HSPC/LSC-like AML*, conv arm only —
  is the niche-stressed cluster a stress sub-state of HSPC/LSC (H3a) or a
  transcriptionally distinct primed clone (H3b)?
- S2: *HSPC/LSC-like vs GMP-like AML* (all conditions) — answers the UMAP
  repositioning observation from `seurat_04c`: is HSPC/LSC-like actually a
  primitive cluster or a committed-myeloid one?

**Limit:** p-values reflect cell-level variability within one patient and do
*not* generalise. Listed as hypothesis generation.

### 5b. Pseudobulk DE — when cohort expands

**Script:** `seurat_07b_pseudobulk_de.R` (to write — gated on cohort).
**Backend:** `muscat::aggregateData` + DESeq2.
**Same five comparisons** as §5a plus the `~ condition * treatment` interaction
LRT per cluster.
**Downstream:** fgsea (Hallmark, Reactome), DoRothEA / ChEA3 for TF activity.

## 6. Within-cluster comparisons across conditions

**Script:** `seurat_08_within_cluster.R`
**What it does:** Seurat `CellCycleScoring` + Cheung-Rando + Laurenti + persister signatures, restricted to one cluster at a time; paired-sample Wilcoxon across conditions.
**Two tests:**
- Quiescent myeloid: untreated conv vs untreated co-culture.
- ABCB5+ primitive AML: untreated vs Ara-C in the conventional arm.
**Deliverables:** boxplots, per-cell score tables, statistical summary.

## 7. Differential abundance testing

**Script:** `seurat_09_milo.R`
**What it does:** MiloR on the integrated KNN graph (primary), scCODA Bayesian check (secondary).
**Deliverables:** neighbourhood-level DA plot, FDR-controlled abundance table — replaces the naïve "% per cluster" bar charts.

## 8. Cell–cell communication (HS-5 ↔ AML)

**Script:** `seurat_10_liana.R` + `seurat_10b_nichenet.R`
**Inputs:** co-culture Seurat object *including* HS-5.
**What it does:** LIANA consensus across CellChat/CellPhoneDB/NATMI/SingleCellSignalR for triage; NicheNet to nominate which HS-5 ligands drive the quiescent-cluster transcriptome.
**Deliverables:** LIANA dotplot, NicheNet ligand-activity ranking, top ligand→target heatmap.

## 9. Gene regulatory networks

**Script:** `seurat_11_pyscenic.py` + SLURM wrapper
**What it does:** pySCENIC on the merged object (export to AnnData first).
**Two outputs:** regulon AUC per cluster (validates identities — GATA1/2 in primitive, CEBPA/SPI1 in myeloid); differential regulon activity in quiescent cluster between conv and co-culture (FOXO3, KLF4, BCL6 as a priori candidates).
**Deliverables:** regulon AUC matrix, differential regulon table.

## 10. Absolute cell-number estimation

**Script:** `seurat_12_absolute_counts.R`
**Inputs:** wet-lab viable-cell counts at harvest (needs a small TSV from Jamshid / wet-lab notebook); HS-5 fraction from clustering.
**What it does:** back-calculates absolute AML cells per cluster per condition. Re-draws "% enrichment" plots on absolute counts as a second panel — critical for distinguishing H3a (conversion) from H3b (differential survival).
**Deliverables:** absolute-count tables, paired %-vs-absolute panel.

## 11. Ara-C resistance gene signature

**Lives inside** `seurat_03_annotate_reference.R` already (`arac_metabolism`, `mdr_efflux`, `aldh_resistance` blocks of `signatures.yml`).
**Specific test:** is ABCB5+ coherent across the resistance signature, or is ABCB5 expression isolated?

## 12. Apoptosis vs survival

**Lives inside** `seurat_06_stress_scoring.R` and `seurat_08_within_cluster.R`.
**Specific test:** Hallmark APOPTOSIS + BCL2-family ratio (BCL2 + MCL1 vs BAX + BAK) on the unclassified niche-related cluster, conventional vs co-culture. If conv-niche cells already lean pro-apoptotic, H3b becomes more parsimonious than H3a.

## 13. Clonal / lineage information from natural variation

**Script:** `seurat_13_mt_clones.py` (mgatk → MAESTER pipeline)
**Inputs:** BAM files from CellRanger.
**What it does:** call mitochondrial heteroplasmy, build a clone table per cell, ask whether quiescent-myeloid clones in co-culture overlap with niche-related clones in conventional.
**Deliverables:** clone × cluster × condition contingency table.

---

## Execution order

So that each step feeds the next without rework:

1. **Annotation refinement (§1)** — `seurat_03_annotate_reference.R` *[full object]*
2. **CNV identity check (§2)** — `seurat_04_infer_cnv.R` *[full object — needs B+T reference]*
3. **Finalise annotation (post-§1+§2)** — `seurat_04b_finalise_annotation.R` *[full object]*
4. **Split AML compartment (§0)** — `seurat_04c_subset_aml.R` *[produces AML-only RDS]*
5. **Stress (§4) + within-cluster cycling/quiescence (§6)** — `seurat_05_stress_cycling.R` *[AML-only]*
6. **Differential abundance (§7)** — `seurat_09_milo.R` *[AML-only]*
7. **Within-patient DE (§5a)** — `seurat_07_de_within_patient.R` *[AML-only]* (pseudobulk §5b deferred until cohort expands)
8. **Absolute counts (§10)** — `seurat_12_absolute_counts.R` *[AML-only, needs wet-lab counts]*
9. **Resistance + apoptosis scoring (§11, §12)** — already covered by §1 and §4-6 outputs
10. **Trajectory (§3)** — `seurat_05_trajectory.py` *[AML-only]*
11. **Cell-cell communication (§8)** — `seurat_10_liana.R` + `seurat_10b_nichenet.R` *[full object — needs HS-5]*
12. **Gene regulatory networks (§9)** — `seurat_11_pyscenic.py` *[AML-only]*
13. **MT clonal tracing (§13)** — `seurat_13_mt_clones.py` *[full object — BAM-level]*

---

## Status — what exists vs what still needs writing

| Step | Script(s) | Status |
| --- | --- | --- |
| §1 | `seurat_03_annotate_reference.R` | **written + run** |
| §1 / §11 / §12 | `config/signatures.yml` | **written + audit-corrected** |
| §2 | `seurat_04_infer_cnv.R` + `34_*.sbatch` | **written + run** |
| post-§1+§2 | `seurat_04b_finalise_annotation.R` + `34b_*.sbatch` | **written + run** |
| §0 split | `seurat_04c_subset_aml.R` + `34c_*.sbatch` | **written** (this session) |
| §4 + §6 | `seurat_05_stress_cycling.R` + `35_*.sbatch` | **written + running** |
| §3 | `seurat_05_trajectory.py` | to write |
| §5a | `seurat_07_de_within_patient.R` + `37_*.sbatch` | **written** (this session) |
| §5b | `seurat_07b_pseudobulk_de.R` | deferred (needs multi-patient cohort) |
| §7 | `seurat_09_milo.R` | to write |
| §8 | `seurat_10_liana.R`, `seurat_10b_nichenet.R` | to write |
| §9 | `seurat_11_pyscenic.py` | to write |
| §10 | `seurat_12_absolute_counts.R` | to write (needs viable-cell counts from Jamshid) |
| §13 | `seurat_13_mt_clones.py` | to write |

---

## Open questions / dependencies

- **Multi-patient cohort.** Current object is one patient (SRAML10) × 4 condition × treatment samples. Pseudobulk DE (§5b), MiloR with patient as a covariate (§7), and any cross-patient generalisation are gated on the remaining S34 patients being processed through `seurat_01_general.R`.
- **Patient mutations** (FLT3-ITD, NPM1, DNMT3A, IDH1/2 etc.) for any of the 12 AML samples? If yes, Vartrix gives a clean orthogonal AML/normal call in §2 — request from the clinical metadata sheet.
- **Wet-lab viable-cell counts at harvest** (per well, per condition) for §10 — needed before that step is meaningful.
- **Zeng 2025 reference** — flip `enabled: true` in `configs/references.yml` once the Zenodo DOI is pasted in, for a second-opinion annotation in §1.
