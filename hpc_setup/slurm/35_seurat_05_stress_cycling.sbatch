#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_05_stress
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_05_stress.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_05_stress.%J.err

# =============================================================================
# Part-3 §4 + §6: stress signalling (conv vs co-culture, untreated only) and
# within-cluster cell-cycle / quiescence comparisons (Promono + ABCB5+).
# Submits: seurat_05_stress_cycling.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_05"

echo "=== Running seurat_05_stress_cycling.R ==="

# NOTE: consumes the AML-only object produced by seurat_04c_subset_aml.R.
# If that file is missing, fall back to the full object — the script already
# filters to aml_clusters internally, so the result is biologically the same
# (just slower and noisier on global QC).
RDS_IN="${OUT_DIR}/AML_SRAML10_aml_only.RDS"
if [[ ! -f "${RDS_IN}" ]]; then
    echo "WARN: ${RDS_IN} not found — falling back to full object"
    RDS_IN="${OUT_DIR}/AML_SRAML10_seurat_general.RDS"
fi
echo "Input: ${RDS_IN}"

Rscript "${SCRIPT_DIR}/scripts/seurat_05_stress_cycling.R" \
    --rds_in  "${RDS_IN}" \
    --out_dir "${OUT_DIR}/seurat_05" \
    --ncores  ${SLURM_CPUS_PER_TASK}

echo "Done: $(date)"
