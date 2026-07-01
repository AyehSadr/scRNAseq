# AML CloneTracker XP — scRNA-seq Analysis Pipeline

This repository contains the complete bioinformatics workflow and pipeline configuration for the single-cell RNA-sequencing (scRNA-seq) analysis of conventional vs. co-culture Ara-C response in AML clone tracking.

The project is designed to run on the **Cardiff Falcon HPC cluster** (using SLURM scheduling and Lustre scratch storage) and sync results locally back to your PC.

**Project Unix Group:** `<your_project_group>` · **Author:** Ayeh Sadr

---

## 📂 Repository Structure

The repository is organized into three major components:

```
├── README.md                      # Project homepage & index
├── hpc_setup/                     # HPC deployment, config, environments & job scripts
│   ├── config/                    # Configuration files (samples.tsv, signatures.yml, paths.yml)
│   ├── envs/                      # Conda environment YAML specifications (scrna.yaml, cellecta.yaml)
│   ├── scripts/                   # Core Python & R pipelines (Seurat, Slingshot, NicheNet, etc.)
│   ├── slurm/                     # SLURM batch submission scripts (.sbatch)
│   ├── deploy.sh                  # Utility script to push local changes to Falcon
│   └── pull_results.sh            # Utility script to pull figures, tables, and logs to local PC
├── docs/                          # Detailed analysis plans, manuals, and validation steps
└── outputs/                       # Clean analysis results (synced back from Falcon scratch)
    └── seurat/                    # Subdivided outputs per analysis step (figures, tables, logs)
```

---

## ⚙️ Pipeline Overview

The pipeline executes the analysis in 10 sequential steps, covering preprocessing to advanced systems biology modeling:

| Step | Script | Description |
|---|---|---|
| **00** | `00_make_dirs.sh` | Sets up the Lustre directory tree on Falcon scratch. |
| **01** | `01_install_miniconda.sh` | Installs Miniconda on scratch to avoid NFS disk locking. |
| **02** | `02_create_envs.sh` | Builds the `scrna` and `cellecta` conda environments. |
| **03** | `seurat_03_annotate_reference.R` | Dual reference-based and publication-grade **Signature-Based (UCell)** cell state annotation. |
| **04** | `seurat_04_infer_cnv.R` | Submits single-cell CNV inference (inferCNV) to distinguish malignant cells from normal. |
| **04b** | `seurat_04b_finalise_annotation.R` | Integrates steps 3 & 4 to finalize labels (e.g. *HSPC-like AML*, *ABCB5+ resistant LSC*, *GMP-like AML*). |
| **05** | `seurat_05_stress_cycling.R` | Stress pathway scoring (UPR, Hypoxia, ROS) and cell-cycle/quiescence profiling. |
| **07** | `seurat_07_de_within_patient.R` | High-performance within-patient single-cell Wilcoxon differential expression. |
| **08** | `seurat_08_trajectory_slingshot.R` | Centroid-based lineage topology and pseudotime fitting using Slingshot. |
| **10** | `seurat_10_liana_nichenet.R` | Ligand-receptor cell-cell communication scoring between HS-5 stromal cells and AML sub-populations. |

---

## 🚀 Getting Started

### 1. Synchronize the Repository to Falcon (from your local Mac)
Run the deploy script to sync code and configuration to Falcon:
```bash
cd hpc_setup
bash deploy.sh
```

### 2. Run the Analysis on Falcon
Log into Falcon, load the environment configurations, and submit your array jobs:
```bash
ssh <your_username>@falconlogin.cf.ac.uk
cd ~/aml_cellecta_setup
source config/falcon_env.sh
sbatch slurm/34b_seurat_04b_finalise_annotation.sbatch
```

### 3. Download Results back to your PC (from your local Mac)
Once jobs complete on the cluster, pull down the publication-ready figures and tables:
```bash
cd hpc_setup
bash pull_results.sh
```
*Note: Large Seurat `.RDS` files are excluded by default to save bandwidth. To include them, run `bash pull_results.sh --include-rds`.*

---

## 🧬 Core Scientific Analysis Details

### Cell State Annotations
Uses single-cell signature scoring via `UCell` using custom marker lists extracted from the van Galen et al. 2019 Cell publication. Cells are classified into key compartments:
*   **HSPC / LSC-like AML**: Primitive, stem-like malignant cells.
*   **ABCB5+ resistant LSC AML**: Stem-like cells showing high CNV burden and drug-resistance signatures.
*   **GMP-like / Promono quiescent / Mono-like AML**: Malignant cells committed to myeloid differentiation states.
*   **Normal Immune**: B-cells, T-cells, and Endothelial cells present in the microenvironment.

### Trajectory Inference
`Slingshot` fits minimum spanning trees (MST) on the `pca_aml` space, smoothing principal curves using `ABCB5+ resistant LSC AML` as the root cluster to reconstruct differentiation kinetics. Fits are split by conventional and co-culture arms to quantify niche-driven lineage bifurcation.

### Cell-Cell Communication (LIANA + NicheNet)
Models the signaling interactome between the HS-5 stromal cells (senders) and leukaemic cells (receivers) in co-culture.
*   **LIANA**: Computes consensus ligand-receptor ranking across five algorithms (CellPhoneDB, NATMI, Connectome, SCA, logFC).
*   **NicheNet**: Ranks stromal ligands by their predicted ability to regulate the niche-imprint gene set (including *FTH1*, *FTL*, *NAMPT*, *TXNIP*, *B2M*) in AML cells.
