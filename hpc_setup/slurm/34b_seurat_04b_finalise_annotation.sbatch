#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_04b_final
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/seu_04b_final.%J.out
#SBATCH --error=logs/seu_04b_final.%J.err

# =============================================================================
# Part-3 finalise-annotation step: relabel clusters using the validated AML /
# Normal / Uncertain / CellLine assignments and re-render the §1 figures.
# Submits: seurat_04b_finalise_annotation.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="${HOME}/scRNAseq"
if [ ! -d "${SCRIPT_DIR}" ]; then
    SCRIPT_DIR="${HOME}/aml_cellecta_setup"
fi
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_04b"

echo "=== Running seurat_04b_finalise_annotation.R ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_04b_finalise_annotation.R" \
    --rds_in  "${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_general.RDS" \
    --out_dir "${OUT_DIR}/seurat_04b"

echo "Done: $(date)"
