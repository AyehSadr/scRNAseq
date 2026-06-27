#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=seu_10_comm
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/shared/scratch/SCWF00196/c.medas36/logs/seu_10_comm.%J.out
#SBATCH --error=/shared/scratch/SCWF00196/c.medas36/logs/seu_10_comm.%J.err

# =============================================================================
# Part-3 §8: HS-5 ↔ AML cell-cell communication via LIANA + NicheNet.
# CONSUMES THE FULL OBJECT (with HS-5) — not the AML-only subset.
# Restricts to co-culture wells where HS-5 is physically present.
# LIANA gives consensus ligand-receptor pairs across five methods;
# NicheNet ranks ligands by their predicted ability to drive the §5.2
# niche-imprint gene set in each AML receiver cluster.
# Submits: seurat_10_liana_nichenet.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate scrna

OUT_DIR="${PROJECT_ROOT}/seurat"
mkdir -p "${OUT_DIR}/seurat_10"

# Consume the FULL object (HS-5 retained), NOT the AML-only subset
RDS_IN="${OUT_DIR}/AML_SRAML10_seurat_general.RDS"
if [[ ! -f "${RDS_IN}" ]]; then
    echo "ERROR: ${RDS_IN} not found"
    exit 1
fi

# Cache the NicheNet prior matrices in the project space so subsequent runs
# don't re-download. ~500MB total.
NN_DATA_DIR="${PROJECT_ROOT}/refs/nichenet"
mkdir -p "${NN_DATA_DIR}"

echo "=== Running seurat_10_liana_nichenet.R ==="
echo "Input              : ${RDS_IN}"
echo "NicheNet data dir  : ${NN_DATA_DIR}"

Rscript "${SCRIPT_DIR}/scripts/seurat_10_liana_nichenet.R" \
    --rds_in            "${RDS_IN}" \
    --out_dir           "${OUT_DIR}/seurat_10" \
    --sender            "HS-5 stromal (cell line)" \
    --nichenet_data_dir "${NN_DATA_DIR}" \
    --top_n_liana       30 \
    --top_n_ligands     15

echo "Done: $(date)"
