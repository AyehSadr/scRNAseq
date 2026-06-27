#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_infercnv
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_infercnv.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_infercnv.%J.err

# =============================================================================
# Step 2 of Part-3 plan: inferCNV — is the CD16+ Mono/Macrophage cluster
# AML-derived or normal? Uses B + T cells as the diploid reference.
# Submits: seurat_04_infer_cnv.R
#
# Notes:
#   * inferCNV is heavy — 24h walltime and 128G memory are conservative for
#     ~30k cells. Halve --max_cells (in the R call below) if the run OOMs.
#   * Requires hpc_setup/refs/infercnv/gene_order_file.hg38.txt — download
#     once on a login node:
#       wget -O ${PROJECT_ROOT}/refs/infercnv/gene_order_file.hg38.txt \\
#         https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
GENE_ORDER="${PROJECT_ROOT}/refs/infercnv/gene_order_file.hg38.txt"

mkdir -p "${OUT_DIR}/seurat_04"

if [[ ! -f "${GENE_ORDER}" ]]; then
  echo "ERROR: gene order file missing — see header of this script for the wget command"
  exit 1
fi

echo "=== Running seurat_04_infer_cnv.R ==="
echo "  PROJECT_ROOT = ${PROJECT_ROOT}"
echo "  GENE_ORDER   = ${GENE_ORDER}"

Rscript "${SCRIPT_DIR}/scripts/seurat_04_infer_cnv.R" \
    --rds_in     "${OUT_DIR}/AML_SRAML10_seurat_general.RDS" \
    --gene_order "${GENE_ORDER}" \
    --out_dir    "${OUT_DIR}/seurat_04" \
    --ref_groups "B cells,T cells" \
    --cutoff     0.1 \
    --ncores     ${SLURM_CPUS_PER_TASK} \
    --max_cells  3000

echo "Done: $(date)"
