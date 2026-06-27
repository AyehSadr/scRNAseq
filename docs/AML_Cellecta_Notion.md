# AML × CloneTracker scRNA-seq — Project page

ICR Project S34 · run1930 · 12 AML samples · started March 2026

---

## What this project is

A single-cell transcriptomic study of AML clonal dynamics under stromal co-culture and chemotherapy. Three patient-derived AML samples were lentivirally barcoded with **Cellecta CloneTracker XP**, expanded under four matched conditions, then sequenced on 10x Chromium 3' as fixed cells.

The biological question: how do individual AML clones respond to AraC chemotherapy, with and without protective stromal (HS5) co-culture? Each cell carries a heritable clonal barcode, so we can ask whether specific clones are stroma-protected, intrinsically AraC-resistant, or transcriptionally distinct before treatment.

**Critical caveat:** the libraries were sequenced as **standard 10x 3' GEX only — no separate Cellecta enrichment PCR library was prepared**. Per Cellecta's own analysis (and the slides in this folder), barcode recovery from GEX alone typically reaches <10% of cells. We are running the analysis with that limitation in mind; the pipeline extracts whatever barcode signal is recoverable from R2 of the GEX reads.

---

## Samples (n=12)

| Sample ID | Patient | Tube | HS5 stroma | AraC | Approx cells (×10⁴) |
|---|---|---|---|---|---|
| P1_T03_noStroma_noAraC | SRAML7 | 3 | – | – | 56 |
| P1_T04_noStroma_AraC | SRAML7 | 4 | – | + | 38 |
| P1_T07_Stroma_noAraC | SRAML7 | 7 | + | – | 14 |
| P1_T08_Stroma_AraC | SRAML7 | 8 | + | + | 14 |
| P2_T11_noStroma_noAraC | SRAML10 | 11 | – | – | 28 |
| P2_T12_noStroma_AraC | SRAML10 | 12 | – | + | 17 |
| P2_T15_Stroma_noAraC | SRAML10 | 15 | + | – | 7.5 |
| P2_T16_Stroma_AraC | SRAML10 | 16 | + | + | 6.75 |
| P3_T17_noStroma_noAraC | SRAML13 | 17 | – | – | 81 |
| P3_T18_noStroma_AraC | SRAML13 | 18 | – | + | 40 |
| P3_T19_Stroma_noAraC | SRAML13 | 19 | + | – | 52.5 |
| P3_T20_Stroma_AraC | SRAML13 | 20 | + | + | 18.2 |

3 patients × 4 conditions = 12 samples. Submitted to ICR Genomics on 16/03/2026.

---

## What we are doing — analysis plan

### Phase 1: HPC setup

Stand up the Cardiff Falcon environment so analysis can run reproducibly: directory tree on scratch, conda envs, Cell Ranger install, custom GRCh38 + CloneTracker reference, SLURM submission templates, sample sheet.

**Status (2026-05-04):**
- ✅ Project allocation confirmed (SCWF00196)
- ✅ Directory tree created on Lustre scratch
- ✅ Conda installed on scratch (NOT home — NFS+SQLite incompatibility forces this)
- ✅ Both conda envs built (`scrna`, `cellecta`) with all PM-brief packages
- ⏳ Cell Ranger 10 install (manual download from 10x site)
- ⏳ GRCh38-2024-A reference download (`sbatch slurm/03_download_references.sbatch`)
- ⏳ Real CloneTracker construct FASTA (waiting on Sian Rizzo / kit insert)
- ⏳ Custom CloneTracker reference build (`sbatch slurm/04_build_clonetracker_ref.sbatch`)
- ⏳ Chemistry confirmation from ICR (Floriana Manodoro) — gates Analysis 3

### Phase 2: Per-sample alignment

Run `cellranger count` on each of the 12 samples against a **custom reference**: standard GRCh38-2024-A with the CloneTracker construct appended as an extra contig. Reads that align to the construct retain their cell barcode and UMI, which is the link we need.

### Phase 3: CloneTracker barcode extraction

Two complementary paths run on every sample:

1. **BAM-driven** — parse `possorted_genome_bam.bam`, keep reads with valid `CB:Z` tag that either align to the CloneTracker contig or contain the FBP1 anchor `CCGACCACCGAACGCAACGCACGCA`. Extract the BC14-spacer-BC30 cassette adjacent to the anchor.
2. **Fastq-driven** — parse raw R1+R2 directly as a sanity check; CB+UMI from R1, cassette from R2.

Outputs collapsed with **starcode** (1–2 mismatch tolerance) to fold sequencing errors back into true clones.

### Phase 4: Downstream analysis (per PM brief 2026-05-04)

Four analyses, in priority order:

1. **AraC effect within each niche** — pairwise within-patient (3v4, 7v8, 11v12, 15v16, 17v18, 19v20). Differential cluster abundance (Milo / scCODA) + pseudo-bulk DGE (DESeq2 / edgeR).
2. **HS5 stromal niche effect** — pairwise within-patient (3v7, 4v8, 11v15, 12v16, 17v19, 18v20). GSEA against Hallmark/KEGG/Reactome/GO BP.
3. **CloneTracker clonal tracking** — *gated on confirming feature barcoding with ICR.* Detection metrics, Jaccard overlap, Sankey plots, clone-specific DGE, trajectory analysis.
4. **AML differentiation staging + cell cycle** — label transfer from van Galen 2019; LSC17 score (Ng 2016); Tirosh G1/S/G2M classification. Question: do AraC or HS5 shift the differentiation composition?

QC thresholds: genes/cell 200–6,000 (≥150 for low-input tubes 15/16); UMI ≥ 500; MT% < 20%. Within-patient comparisons use uncorrected data; cross-patient uses Harmony/scVI integration.

Deliverables: per-sample QC, per-patient + integrated UMAPs, DGE CSVs + volcano plots, GSEA CSVs + bubble plots, Sankey clonal-overlap plots, annotated RDS / H5AD per patient + integrated.

### Phase 5: Risk register

- If barcode recovery is <2% (worst case from GEX-only), pivot to clone-aware bulk analysis or request a redo with enrichment library
- Flex chemistry uncertainty — confirm with ICR Genomics before any sequencer time downstream
- Falcon scratch 60-day purge; long-term outputs go to ICR RDS or OneDrive

---

## Where the scripts live

### Local (Mac, source of truth for editing)

```
/Users/ayehsadr/Documents/Jamshid/AML_Cellecta/
├── README.md                                  ← this file's source
├── AML_Cellecta_Notion.md                     ← this Notion-import file
├── Falcon_HPC_Reference.docx                  ← cluster cheat sheet
├── CloneTracker_barcode_workflow.pptx         ← full Cellecta workflow (with enrichment)
├── CloneTracker_scRNAseq_only.pptx            ← failure-modes when no enrichment
├── scRNAseq sample summary Mar2026 ICR Project S34.xlsx
└── hpc_setup/                                 ← rsync this whole folder to Falcon
    ├── README.md
    ├── config/
    │   ├── falcon_env.sh                      ← edit FALCON_PROJECT here
    │   └── samples.tsv
    ├── envs/
    │   ├── scrna.yaml
    │   └── cellecta.yaml
    ├── refs/clonetracker/
    │   ├── clonetracker_construct.fa          ← replace placeholder before step 04
    │   └── clonetracker_construct.gtf
    ├── scripts/
    │   ├── 00_make_dirs.sh
    │   ├── 01_install_miniconda.sh
    │   ├── 02_create_envs.sh
    │   └── extract_clonetracker_barcodes.py   ← BAM + fastq dual extractor
    └── slurm/
        ├── 03_download_references.sbatch
        ├── 04_build_clonetracker_ref.sbatch
        ├── 10_cellranger_count.sbatch         ← array job 1–12
        └── 20_cellecta_extract.sbatch         ← array job 1–12
```

### Falcon HPC (deployment target)

```
~/aml_cellecta_setup/                          ← mirror of hpc_setup/, deployed via rsync
/shared/scratch/SCWF00196/c.medas36/           ← project root on scratch
├── raw/                                       ← Globus drop (S34_*.fastq.gz)
├── refs/
│   ├── refdata-gex-GRCh38-2024-A/
│   └── refdata-GRCh38-CloneTracker/           ← built by step 04
├── envs/
├── software/cellranger-10.0.0/                ← self-installed
├── cellranger/{sample_id}/outs/               ← per-sample CR output
├── cellecta/{sample_id}.bam.tsv               ← extracted barcode tables
├── cellecta/{sample_id}.bam.clones.tsv        ← starcode-collapsed clones
└── logs/                                      ← SLURM .out / .err
```

---

## Cluster reference (Cardiff Falcon)

| Field | Value |
|---|---|
| Login | `ssh c.medas36@falconlogin.cf.ac.uk` (Cardiff VPN required) |
| Open OnDemand | https://falconlogin.cf.ac.uk:8080 |
| Scheduler | SLURM, partition `compute` |
| Project code | SCWF00196 |
| Home | `/shared/home1/c.medas36/` (50 GB, no big data) |
| Scratch | `/shared/scratch/SCWF00196/` (3 TB, no backup, future 60-day purge) |
| ARCCA support | arcca-help@cardiff.ac.uk · wiki.arcca.cf.ac.uk |

---

## Deployment steps

From the Mac, in `~/Documents/Jamshid/AML_Cellecta/hpc_setup/`:

```bash
rsync -av --exclude '.DS_Store' --exclude '__pycache__' \
  ./ c.medas36@falconlogin.cf.ac.uk:/shared/home1/c.medas36/aml_cellecta_setup/
```

Then on Falcon login node:

```bash
cd ~/aml_cellecta_setup
# Edit config/falcon_env.sh — set FALCON_PROJECT to the real code
bash scripts/00_make_dirs.sh
bash scripts/01_install_miniconda.sh
bash scripts/02_create_envs.sh
sbatch slurm/03_download_references.sbatch
# (after step 03 finishes, drop the real CloneTracker FASTA in refs/clonetracker/)
sbatch slurm/04_build_clonetracker_ref.sbatch
```

Once those four steps complete, the cluster is ready for analysis. Per-sample runs:

```bash
sbatch slurm/10_cellranger_count.sbatch
# (after Cell Ranger finishes for at least one sample)
sbatch slurm/20_cellecta_extract.sbatch
```

---

## Open items / blockers

1. ~~**Falcon project code**~~ — confirmed as **SCWF00196** (member via `groups`).
2. **CloneTracker construct sequence** — the cassette FASTA in `refs/clonetracker/` is a placeholder; need the real CloneTracker XP construct from Cellecta kit insert (ask Sian Rizzo if missing).
3. **Chemistry confirmation** — confirm with ICR Genomics (Floriana Manodoro) whether the "fixed cell samples" are standard 10x 3' v3.1 or 10x Flex (probe-based). Flex would invalidate the barcode-recovery half of the pipeline because the CloneTracker transcript would not be captured.
4. **Sample sheet `fastq_prefix` column** — fill in once Globus delivery to `${PROJECT_ROOT}/raw/` completes.

---

## People

| Person | Role | Contact |
|---|---|---|
| Alex Tonks | PI (Cardiff) | TonksA@cardiff.ac.uk |
| Jamshid Khorashad | Project lead (ICR) | jamshid.khorashad@icr.ac.uk |
| Sian Rizzo | Wet lab (Cardiff) | RizzoS3@cardiff.ac.uk |
| Floriana Manodoro | ICR Genomics — chemistry/run details | Floriana.Manodoro@icr.ac.uk |
| Ritika Chauhan | ICR bioinformatician | Ritika.Chauhan@icr.ac.uk |
| ARCCA helpdesk | Cardiff HPC support | arcca-help@cardiff.ac.uk |

---

## How to import this into Notion

1. In Notion, open the workspace where you want the page.
2. Click **+ Add a page** in the sidebar (or open an existing parent page).
3. Click **Import** at the bottom of the empty page (or `⋯` menu → **Import**).
4. Choose **Markdown & CSV**.
5. Select this file: `AML_Cellecta_Notion.md`.

The headings, tables, code blocks, and lists will all map cleanly to Notion blocks. The result is a single page; you can split it into sub-pages later by dragging headings into the sidebar.
