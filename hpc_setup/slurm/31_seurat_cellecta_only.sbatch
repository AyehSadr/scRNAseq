#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_cellecta
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/seu_cellecta.%J.out
#SBATCH --error=logs/seu_cellecta.%J.err

set -euo pipefail
SCRIPT_DIR="${HOME}/scRNAseq"
if [ ! -d "${SCRIPT_DIR}" ]; then
    SCRIPT_DIR="${HOME}/aml_cellecta_setup"
fi
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

CT_DIR="${PROJECT_ROOT}/cellecta"
OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "=== Running Cellecta Downstream Analysis (Skipping General) ==="

echo "--- Running 02a_Barcode QC ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02a_barcode_qc.R" \
    --rds_in  "${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_general.RDS" \
    --ct_dir  "${CT_DIR}" \
    --out_dir "${OUT_DIR}" \
    --allow_low_recovery

# Force parallel filesystem (Lustre/GPFS) to flush the newly written RDS to disk
sync
sleep 10

echo "--- Running 02b_Clonal Overlap ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02b_clonal_overlap.R" \
    --rds_in  "${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_cellecta_qc.RDS" \
    --out_dir "${OUT_DIR}"

echo "--- Running 02c_Same-Clone DGE ---"
Rscript "${SCRIPT_DIR}/scripts/seurat_02c_clone_dge.R" \
    --rds_in  "${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_cellecta_qc.RDS" \
    --out_dir "${OUT_DIR}"

echo "Cellecta analysis pipeline complete. Outputs: ${OUT_DIR}"
