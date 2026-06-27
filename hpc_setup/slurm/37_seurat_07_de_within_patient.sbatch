#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_07_de
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_07_de.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_07_de.%J.err

# =============================================================================
# Part-3 §5: within-patient single-cell DE on the AML-only object.
# COHORT CAVEAT: only one patient (SRAML10) × 4 samples — pseudobulk DE
# (muscat / DESeq2) is not feasible. This step runs single-cell Wilcoxon DE
# via Seurat::FindMarkers for hypothesis generation. When the full S34
# cohort lands, swap in seurat_07b_pseudobulk_de.R.
# Submits: seurat_07_de_within_patient.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_07"

# Consume the AML-only object from seurat_04c; fall back to the full object
# if it isn't there (FindMarkers will still work — the script filters to
# cell_type_final levels that exist).
RDS_IN="${OUT_DIR}/AML_SRAML10_aml_only.RDS"
if [[ ! -f "${RDS_IN}" ]]; then
    echo "WARN: ${RDS_IN} not found — falling back to full object"
    RDS_IN="${OUT_DIR}/AML_SRAML10_seurat_general.RDS"
fi
echo "Input: ${RDS_IN}"

echo "=== Running seurat_07_de_within_patient.R ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_07_de_within_patient.R" \
    --rds_in    "${RDS_IN}" \
    --out_dir   "${OUT_DIR}/seurat_07" \
    --min_cells 25 \
    --logfc     0.25 \
    --top_n     30

echo "Done: $(date)"
