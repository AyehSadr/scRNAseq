#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_annot
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_annot.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_annot.%J.err

set -euo pipefail
SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"

echo "=== Running Cluster Annotation ==="

Rscript "${SCRIPT_DIR}/scripts/seurat_01b_annotate_clusters.R" \
    --rds_in  "${OUT_DIR}/AML_SRAML10_seurat_general.RDS" \
    --out_dir "${OUT_DIR}"

echo "Annotation complete. Updated RDS and generated new UMAPs in ${OUT_DIR}/figures"
