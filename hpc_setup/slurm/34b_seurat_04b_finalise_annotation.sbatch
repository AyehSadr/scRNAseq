#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_04b_final
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_04b_final.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_04b_final.%J.err

# =============================================================================
# Part-3 finalise-annotation step: relabel clusters using the validated AML /
# Normal / Uncertain / CellLine assignments and re-render the §1 figures.
# Submits: seurat_04b_finalise_annotation.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_04b"

echo "=== Running seurat_04b_finalise_annotation.R ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_04b_finalise_annotation.R" \
    --rds_in  "${OUT_DIR}/AML_SRAML10_seurat_general.RDS" \
    --out_dir "${OUT_DIR}/seurat_04b"

echo "Done: $(date)"
