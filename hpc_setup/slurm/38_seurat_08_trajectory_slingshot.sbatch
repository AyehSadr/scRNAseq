#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_08_sling
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/seu_08_sling.%J.out
#SBATCH --error=logs/seu_08_sling.%J.err

# =============================================================================
# Part-3 §3: Slingshot trajectory inference on the AML-only object.
# Fit on AML-only PCA (pca_aml from seurat_04c), visualise on AML-only UMAP.
# Root: ABCB5+ resistant LSC AML (highest CNV + primitive signatures).
# Global fit + per-condition fits (conv vs coc) for topology comparison.
# Submits: seurat_08_trajectory_slingshot.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="${HOME}/scRNAseq"
if [ ! -d "${SCRIPT_DIR}" ]; then
    SCRIPT_DIR="${HOME}/aml_cellecta_setup"
fi
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_08"

RDS_IN="${OUT_DIR}/AML_${FALCON_PATIENT}_aml_only.RDS"
RDS_OUT="${OUT_DIR}/AML_${FALCON_PATIENT}_aml_only_pst.RDS"

if [[ ! -f "${RDS_IN}" ]]; then
    echo "ERROR: ${RDS_IN} not found — run seurat_04c first."
    exit 1
fi
echo "Input : ${RDS_IN}"
echo "Output: ${RDS_OUT}"

echo "=== Running seurat_08_trajectory_slingshot.R ==="

# To enable tradeSeq dynamic genes (slow, optional), add --run_tradeseq
Rscript "${SCRIPT_DIR}/scripts/seurat_08_trajectory_slingshot.R" \
    --rds_in        "${RDS_IN}" \
    --rds_out       "${RDS_OUT}" \
    --out_dir       "${OUT_DIR}/seurat_08" \
    --root_cluster  "ABCB5+ resistant LSC AML" \
    --root_fallback "HSPC / LSC-like AML" \
    --min_root_cells 30 \
    --n_pcs         30

echo "Done: $(date)"
