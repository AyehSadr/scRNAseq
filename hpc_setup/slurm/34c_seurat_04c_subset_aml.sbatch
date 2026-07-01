#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_04c_subset
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/seu_04c_subset.%J.out
#SBATCH --error=logs/seu_04c_subset.%J.err

# =============================================================================
# Part-3 — subset to AML-only compartment.
# Takes the full object after seurat_04b (Steps 1+2 validated labels) and
# produces an AML-only Seurat .rds re-embedded on AML-internal variance.
# Used as input by seurat_05_stress_cycling, seurat_07_pseudobulk_de, and
# seurat_09_milo. LIANA (§8) keeps the full object so HS-5 stays in.
# Submits: seurat_04c_subset_aml.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="${HOME}/scRNAseq"
if [ ! -d "${SCRIPT_DIR}" ]; then
    SCRIPT_DIR="${HOME}/aml_cellecta_setup"
fi
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_04c"

RDS_IN="${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_general.RDS"
RDS_OUT="${OUT_DIR}/AML_${FALCON_PATIENT}_aml_only.RDS"

echo "=== Running seurat_04c_subset_aml.R ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_04c_subset_aml.R" \
    --rds_in  "${RDS_IN}" \
    --rds_out "${RDS_OUT}" \
    --out_dir "${OUT_DIR}/seurat_04c" \
    --ncores  ${SLURM_CPUS_PER_TASK} \
    --n_pcs   30 \
    --n_hvg   2000

echo "Done: $(date)"
