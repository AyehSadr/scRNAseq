#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_annot
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/seu_annot.%J.out
#SBATCH --error=logs/seu_annot.%J.err

set -euo pipefail
SCRIPT_DIR="${HOME}/scRNAseq"
if [ ! -d "${SCRIPT_DIR}" ]; then
    SCRIPT_DIR="${HOME}/aml_cellecta_setup"
fi
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"

echo "=== Running Cluster Annotation ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_01b_annotate_clusters.R" \
    --rds_in  "${OUT_DIR}/AML_${FALCON_PATIENT}_seurat_general.RDS" \
    --out_dir "${OUT_DIR}"

echo "Annotation complete. Updated RDS and generated new UMAPs in ${OUT_DIR}/figures"
