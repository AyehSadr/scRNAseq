#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_ref_anno
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_ref_anno.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_ref_anno.%J.err

# =============================================================================
# Step 1 of Part-3 plan: SingleR + Azimuth label transfer + UCell signatures
# Submits: seurat_03_annotate_reference.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
REFS_DIR="${PROJECT_ROOT}/refs/scrna"
SIGS_YML="${SCRIPT_DIR}/config/signatures.yml"

mkdir -p "${OUT_DIR}/seurat_03"

echo "=== Running seurat_03_annotate_reference.R ==="
echo "  PROJECT_ROOT = ${PROJECT_ROOT}"
echo "  REFS_DIR     = ${REFS_DIR}"
echo "  SIGS_YML     = ${SIGS_YML}"

Rscript "${SCRIPT_DIR}/scripts/seurat_03_annotate_reference.R" \
    --rds_in   "${OUT_DIR}/AML_SRAML10_seurat_general.RDS" \
    --refs_dir "${REFS_DIR}" \
    --sigs_yml "${SIGS_YML}" \
    --out_dir  "${OUT_DIR}/seurat_03" \
    --ncores   ${SLURM_CPUS_PER_TASK}

echo "Done: $(date)"
