#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seurat
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seurat.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seurat.%J.err

set -euo pipefail
SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

CR_DIR="${PROJECT_ROOT}/cellranger/S34_02/S34_02/outs/per_sample_outs"
CT_DIR="${PROJECT_ROOT}/cellecta"
OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "--- Running 01_General Seurat Analysis ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_01_general.R" \
    --cr_dir  "${CR_DIR}"  \
    --out_dir "${OUT_DIR}" \
    --ncores  "${SLURM_CPUS_PER_TASK}" \
    --mt_max  20  \
    --min_feat 200 \
    --max_feat 6000

echo "--- Running 02a_Barcode QC ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02a_barcode_qc.R" \
    --rds_in  "${OUT_DIR}/AML_SRAML10_seurat_general.RDS" \
    --ct_dir  "${CT_DIR}" \
    --out_dir "${OUT_DIR}" \
    --allow_low_recovery

echo "--- Running 02b_Clonal Overlap ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02b_clonal_overlap.R" \
    --rds_in  "${OUT_DIR}/AML_SRAML10_seurat_cellecta_qc.RDS" \
    --out_dir "${OUT_DIR}"

echo "--- Running 02c_Same-Clone DGE ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02c_clone_dge.R" \
    --rds_in  "${OUT_DIR}/AML_SRAML10_seurat_cellecta_qc.RDS" \
    --out_dir "${OUT_DIR}"

echo "Seurat analysis pipeline complete. Outputs: ${OUT_DIR}"

