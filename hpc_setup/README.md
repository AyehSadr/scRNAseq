# Cardiff Falcon HPC setup — AML × CloneTracker scRNA-seq

**Project:** ICR Project S34 — 12 AML samples (3 patients × 4 conditions: ±HS5 stroma, ±AraC), submitted 16/03/2026 by S. Rizzo / A. Tonks / J. Khorashad. Sequenced as standard 10x 3' GEX libraries; **no separate Cellecta enrichment library** was prepared, so CloneTracker barcodes must be salvaged from the same R2 reads as the transcriptome.

**Cluster:** Cardiff Falcon (SLURM, Singularity, RHEL-based). User `c.medas36`, project code **SCWF00196** (confirmed via `groups`). Conda lives on Lustre scratch at `/shared/scratch/SCWF00196/c.medas36/software/miniconda3/` (NOT in NFS home — see RUNBOOK Step 4 for why).

**Important caveat (from your own `CloneTracker_scRNAseq_only.pptx`):** with no targeted PCR enrichment, expect <10% of cells to receive a confident CloneTracker call. This setup is built to extract whatever signal is recoverable; it cannot substitute for a proper enrichment library. The pipeline is also designed so that if/when an enrichment library is run later, the same scripts can ingest it with one parameter change.

---

## Looking for the step-by-step?

**Open [`RUNBOOK.md`](RUNBOOK.md)** — that's the actual command-by-command walkthrough with expected outputs and failure-mode advice. The rest of this README is an architecture / reference document.

The very short version:

```bash
# Mac
cd ~/Documents/Jamshid/AML_Cellecta/hpc_setup && bash deploy.sh

# Falcon
ssh c.medas36@falconlogin.cf.ac.uk
cd ~/aml_cellecta_setup
source config/falcon_env.sh
bash scripts/00_make_dirs.sh
bash scripts/01_install_miniconda.sh
bash scripts/02_create_envs.sh
sbatch slurm/03_download_references.sbatch
# (after step 03 finishes, replace placeholder CloneTracker FASTA)
sbatch slurm/04_build_clonetracker_ref.sbatch
```

---

## Pre-flight checklist

1. **VPN connected** (GlobalProtect → ras.cf.ac.uk).
2. **SSH** test: `ssh c.medas36@falconlogin.cf.ac.uk "echo ok && groups"` — must show `SCWF00196`.
3. **Confirm fastq delivery location.** Raw fastqs come via Globus from "ICR RDS — S34 (run1930)" into `${PROJECT_ROOT}/raw/`. Example file: `S34_02_S17_L001_R1_001.fastq.gz` (3.19 GB).
4. **Confirm 10x chemistry.** The samples are listed as "fixed cell samples" — verify with ICR Genomics (Floriana Manodoro) whether they are:
   - **Standard 10x 3' v3.1** — this pipeline works as-is.
   - **10x Flex (Fixed RNA Profiling)** — probe-based; CloneTracker mRNA is **not** captured because no probes target the barcode region. If Flex, the barcode-extraction half is moot. **Resolve before running anything.**

---

## What this setup gives you on Falcon

```
${PROJECT_ROOT}/                          # /shared/scratch/SCWF00196/c.medas36/
├── raw/                       # incoming fastqs from Globus (you populate)
├── refs/
│   ├── refdata-gex-GRCh38-2024-A/        # 10x pre-built human reference
│   └── refdata-GRCh38-CloneTracker/      # custom: GRCh38 + CloneTracker construct contig
├── envs/                       # workspace for env-related scratch (conda itself lives in ~/miniconda3)
├── software/cellranger-X.Y.Z/  # self-installed Cell Ranger
├── cellranger/                 # one subfolder per sample (post-analysis)
├── cellecta/                   # barcode-extraction outputs (post-analysis)
└── logs/                       # SLURM .out / .err
```

---

## Pipeline overview

```
fastq (R1+R2)
    │
    ├──▶  Cell Ranger count                   ──▶  filtered_feature_bc_matrix/
    │     (custom GRCh38 + CloneTracker ref)        possorted_genome_bam.bam
    │                                                  │
    │                                                  ▼
    │                                         reads tagged CB:Z + UB:Z
    │                                         that aligned to the
    │                                         CloneTracker contig
    │                                                  │
    │                                                  ▼
    │                                         parse BC14–sp–BC30 from
    │                                         read sequence using FBP1 anchor
    │
    └──▶  Independent path (sanity check):
          cutadapt on raw R2 → linked adapter
          (FBP1 anchor + N{50}) → starcode/bartab
          → join with R1's CB+UMI
```

Both paths produce a per-cell `(cell_barcode, clone_barcode, n_reads, n_umis)` table. Agreement between them is your QC.

---

## File index

| File | Purpose |
|---|---|
| `RUNBOOK.md` | **Step-by-step instructions to run.** Start here. |
| `deploy.sh` | One-liner Mac→Falcon rsync (no flags needed). |
| `scripts/00_make_dirs.sh` | Create the project directory tree on scratch. |
| `scripts/01_install_miniconda.sh` | Install Miniconda into `~/miniconda3` if missing. |
| `scripts/02_create_envs.sh` | Build the two conda envs from the YAMLs in `envs/`. |
| `scripts/extract_clonetracker_barcodes.py` | BAM- and fastq-mode barcode extractor (used by step 20). |
| `envs/scrna.yaml` | Cell Ranger helpers (samtools, pysam, umi_tools, scanpy, R/Seurat/Harmony/DropletUtils). |
| `envs/cellecta.yaml` | Barcode extraction (cutadapt, starcode, bartab/pycashier via pip, pyfastx). |
| `slurm/03_download_references.sbatch` | Pull 10x's pre-built human ref tarball into `refs/`. |
| `slurm/04_build_clonetracker_ref.sbatch` | Append CloneTracker construct contig + GTF to GRCh38 and run `cellranger mkref`. |
| `slurm/10_cellranger_count.sbatch` | Per-sample `cellranger count` — array job 1–12. |
| `slurm/20_cellecta_extract.sbatch` | Barcode extractor — array job 1–12, runs both BAM and fastq paths + starcode collapse. |
| `refs/clonetracker/clonetracker_construct.fa` | **Placeholder** — replace with real CloneTracker XP cassette before step 04. |
| `refs/clonetracker/clonetracker_construct.gtf` | Matching GTF; `sed` updates length when you replace the FASTA. |
| `config/samples.tsv` | Sample sheet pre-filled from the ICR S34 spreadsheet. Fill `fastq_prefix` once Globus delivery completes. |
| `config/falcon_env.sh` | Sourced at the top of every script: sets `FALCON_PROJECT=SCWF00196`, paths, conda. |

All sbatch scripts have `#SBATCH --account=SCWF00196` baked in, so you don't need to pass `-A` on the command line.

---

## Open items before you can run analysis

1. ~~Falcon project code~~ — confirmed **SCWF00196**.
2. Drop the real CloneTracker XP construct sequence into `refs/clonetracker/clonetracker_construct.fa` (kit insert; ask Sian Rizzo if missing). The GTF length is auto-updated by the runbook step.
3. Fill the `fastq_prefix` column in `config/samples.tsv` once you `ls` the Globus drop in `${PROJECT_ROOT}/raw/`.
4. Confirm 10x chemistry (3' v3.1 vs Flex) with ICR Genomics.
