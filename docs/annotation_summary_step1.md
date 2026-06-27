# Annotation Summary — Step 1 (Reference annotation + signature scoring)

**Project:** AML_Cellecta (ICR S34, 12 AML samples, CloneTracker XP)
**Step:** Part-3 §1 — reference-based annotation + UCell signature scoring
**Source script:** `hpc_setup/scripts/seurat_03_annotate_reference.R`
**Inputs read:** existing `seu_merged` from `seurat_01b_annotate_clusters.R`, `config/signatures.yml`, van Galen 2019 reference (if SingleR ran on Falcon)
**Last updated:** 2026-05-23

Each row of the table below collapses what we now know about a cluster from three lines of evidence: existing marker-based `cell_type` annotation, the per-cluster UCell signature heatmap (`03d_ucell_heatmap.pdf`), and the van-Galen S4A z-scored argmax class (`03f_confusion_existing_vs_vgalen_class.pdf`).

"vG class (dominant)" is read from the row of `03f` — the class that most cells in that cluster fall into. "TBD" means the call has to wait for §2 (inferCNV) or a later step.

| Cluster (`cell_type`) | Category | Top elevated signatures | vG class (dominant) | Interpretation | Status / next step |
|---|---|---|---|---|---|
| **GMP-like AML blasts** | AML progenitor | `vgalen_gmp_s4a`, `van_galen_gmp_like`, `arac_metabolism`, `mdr_efflux`, `lsc17` | **GMP-like** (~50%) + some Unassigned | Canonical GMP-like AML blasts. Sit on the HSC→GMP boundary; carry an active drug-metabolism / efflux programme even at baseline. | **Confirmed.** Use as a positive reference for §2 (CNV) and as one arm of the GMP↔ABCB5+ trajectory in §3. |
| **HSPC / LSC-like blasts** | AML primitive | `lsc17`, `eppert_hsc`, `van_galen_hsc_like` (moderate), `arac_metabolism` | Mixed HSC-like + Unassigned | Primitive AML blasts on the HSC↔GMP edge of van Galen's taxonomy. Confidence per cell is lower than for GMP-like, reflecting biological intermediates. | Identity consistent. **§2 will confirm AML** vs residual normal HSPC contamination. |
| **Promono-like / quiescent myeloid** | AML myeloid | `laurenti_hsc_quiescence`, `isr_atf4_targets`, `aldh_resistance`, `tikhonova_msc_response`, `vgalen_myeloid_s4a` | **Myeloid-like** (~80%) | Niche-induced quiescent state — quiescence + ATF4 stress + MSC-response signal co-occur. **Direct evidence for H4** (active MSC-induced quiescence) rather than just "slow-cycling default". | **Confirmed.** Primary cluster for §4 (stress), §6 (within-cluster cycling), §8 (NicheNet). |
| **Mono-like AML** | AML myeloid | `van_galen_promono_like`, `vgalen_myeloid_s4a`, moderate niche | Myeloid-like | Intermediate / mature monocytic AML state. | Confirmed. Use as the "later differentiation" anchor for §3 trajectory. |
| **CD16+ Mono / macrophage** | AML myeloid (TBD) | `vgalen_myeloid_s4a`, `van_galen_promono_like`, `tikhonova_msc_response`, `baryawno_niche_signal` | **Myeloid-like** (~80%) | Mature, niche-engaged monocyte phenotype. Signatures alone cannot distinguish AML-derived monocyte vs contaminating normal monocyte. | **§2 inferCNV is the decisive test.** Look for CD16+ row sharing chr-block CNVs with other AML clusters. |
| **Stress-response myeloid** | AML myeloid | `van_den_brink_dissociation`, `isr_atf4_targets` | Mixed Unassigned / Myeloid-like | Dominantly dissociation-stress signal — likely a partial artifact of sample prep rather than independent biology. | **Reconsider whether to keep as a distinct cluster** or merge into neighbouring myeloid populations after regressing dissociation effects. |
| **ABCB5+ resistant primitive AML** | AML primitive | `lsc17`, `eppert_hsc`, `van_galen_hsc_like`, `arac_metabolism`, `mdr_efflux`, `aldh_resistance` | **Unassigned** (~85%) | Coherent primitive-LSC + drug-resistance programme. ABCB5 is *not* incidental — the whole resistance axis (DCK, SAMHD1, MDR pumps, ALDH) lights up together. Unassigned because the profile straddles canonical vG classes. | **Answers Part-3 §11.** Worth a manuscript footnote that this cluster sits *outside* the standard van Galen taxonomy. |
| **Unclassified niche-related** | AML ambiguous | `lsc17`, `van_galen_hsc_like` (slight), `van_den_brink_dissociation`, `isr_atf4_targets`, `bcl2_family_pro_apoptotic` | Mostly Myeloid-like, some Unassigned | Primitive blasts already leaning pro-apoptotic + stressed under conventional. The pro-apoptotic baseline is the key new finding. | **Strong evidence for H3b** (these cells die in co-culture) over H3a (they convert). Re-test rigorously in §12 once the cluster is isolated. |
| **Unclassified (low-quality?)** | AML ambiguous | Low across nearly all signatures | Unassigned (~70%) | Likely a QC residual — low UMI / low feature / high MT — rather than a distinct biological state. | **Re-examine QC** (revisit `seurat_01_general.R` thresholds). Strongly consider excluding from downstream §5, §6, §7. |
| **HS5 stromal / MSC** | Stromal | `baryawno_niche_signal`, `tikhonova_msc_response`, quiescence | **Unassigned** (~75%) | Stromal cell line, correctly recognised as non-AML by the argmax. Carries niche-signalling and quiescence vocabulary by definition. | Confirmed. **Drop from AML-only analyses** (already excluded in `seurat_04_infer_cnv.R`). Keep in §8 (NicheNet) — it's the ligand source. |
| **BM endothelial** | Stromal | `van_galen_promono_like` (moderate), niche signalling | Mixed HSC-like + Unassigned | Endothelial niche cells — small population. Some bleed-through into HSC-like because of shared niche-cytokine vocabulary in van Galen S4A. | Confirmed identity. Treat as auxiliary in §8 if there are enough cells; otherwise drop. |
| **B cells** | Normal lymphoid | `bcl2_family_pro_survival`; partial HSC-like via `CD52`/`CD74` | HSC-like (~50%) + Unassigned | Normal B cells. The HSC-like bleed-through is a known quirk — `CD52`/`CD74`/`IL2RA` are in the HSC S4A list because they correlate with HSC class in van Galen, but they're broadly expressed in lymphocytes. | Confirmed normal. Used as **diploid reference for §2**. Optional: tighten HSC S4A list by removing CD52/CD74/IL2RA. |
| **T cells** | Normal lymphoid | `bcl2_family_pro_survival`; partial HSC-like via `CD52` | HSC-like (~40%) + Unassigned | Normal T cells. Same `CD52` bleed-through as B cells. | Confirmed normal. **Diploid reference for §2.** |
| **Eosinophils / basophils** | Normal myeloid (rare) | Low across most signatures | Unassigned (~70%) | Rare granulocyte population; transcriptionally distinct from the AML monocytic axis. | Confirmed. Drop from AML-only analyses; small `n` won't support its own statistics anyway. |

## What's answered after Step 1

- **§11 — Ara-C resistance signature in ABCB5+**: yes, coherent. ABCB5 cluster lights up across DCK / SAMHD1 / MDR pumps / ALDH. Not a single-gene incidental.
- **§1 — niche-related and ABCB5+ identity refinement**: ABCB5+ is a primitive-LSC + resistance hybrid that sits outside canonical van Galen classes. Unclassified niche-related is a *stressed pro-apoptotic primitive* state.
- **Early lean on §12 (apoptosis vs survival)**: Unclassified niche-related is the only AML cluster strongly elevated for `bcl2_family_pro_apoptotic` + dissociation stress — **H3b is the simpler hypothesis** for its disappearance under co-culture, ahead of running §12 properly.
- **Early lean on §4 (stress) / §6 (within-cluster quiescence)**: Promono-like / quiescent myeloid carries the canonical niche-induced quiescence programme (`laurenti` + `tikhonova` + `aldh` + `atf4` together) — **H4 evidence**.

## What's still open

- **§2 — CD16+ AML vs normal**: the only remaining identity question. inferCNV is queued via `34_seurat_04_infer_cnv.sbatch`.
- **HSPC / LSC-like blasts AML vs residual normal HSPC**: confirm via CNV; signatures alone are not decisive in clusters whose top class is "Unassigned" or mixed.
- **Stress-response myeloid**: decide whether to retain or merge after regressing dissociation effects (might be partly artifact).
- **Unclassified (low-quality?)**: revisit QC thresholds; likely exclude from §5/§6/§7.

## Caveats / known biases of the Step-1 call

- The van Galen S4A taxonomy is **three classes** (HSC/Prog, GMP, Myeloid). Promonocyte vs Monocyte cannot be resolved by S4A alone; we used the curated `van_galen_promono_like` and `van_galen_mono_mmc4` short lists in the heatmap, but they were not part of the argmax classifier.
- `vgalen_class` uses **per-signature z-scored argmax** with `MIN_Z = 0.5` floor. Cells whose top z-score is below the floor → Unassigned; this is by design but means some primitive AML cells get tagged Unassigned rather than HSC-like (visible for ABCB5+).
- `CD52`, `CD74`, `IL2RA` are in `vgalen_hsc_s4a` because they correlate with HSC class in van Galen 2019 — but they're also expressed in lymphocytes, which is why B/T cells show partial HSC-like calls.

## Files this summary draws from

- `seurat/seurat_03/figures/03d_ucell_heatmap.pdf` — per-cluster signature heatmap
- `seurat/seurat_03/figures/03e_vgalen_class_umap.pdf` — UMAP of cell_type vs vgalen_class
- `seurat/seurat_03/figures/03f_confusion_existing_vs_vgalen_class.pdf` — confusion matrix
- `seurat/seurat_03/tables/ucell_scores_by_cluster.tsv` — numerical scores per cluster
- `seurat/seurat_03/tables/cell_type_vs_singler_confusion.tsv` — if SingleR ran (depends on `has_refs`)
