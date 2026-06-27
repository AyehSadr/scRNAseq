#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=ref_download
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ref_download.%J.out
#SBATCH --error=logs/ref_download.%J.err

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config/falcon_env.sh"
mkdir -p "${REF_DIR}" "${LOG_DIR}"
cd "${REF_DIR}"

# 10x's pre-built human GEX reference (latest 2024-A as of writing).
# If you need a different build, swap the URL — list at:
# https://www.10xgenomics.com/support/software/cell-ranger/downloads
TARBALL="refdata-gex-GRCh38-2024-A.tar.gz"
URL="https://cf.10xgenomics.com/supp/cell-exp/${TARBALL}"

if [ -d "refdata-gex-GRCh38-2024-A" ]; then
    echo "Reference already extracted; skipping."
else
    if [ ! -f "${TARBALL}" ]; then
        echo "Downloading ${URL} ..."
        wget --tries=3 --continue -q "${URL}"
    fi
    echo "Verifying tarball is non-empty ..."
    [ -s "${TARBALL}" ] || { echo "Empty tarball" >&2; exit 1; }
    echo "Extracting ${TARBALL} ..."
    tar -xzf "${TARBALL}"
    rm -f "${TARBALL}"
fi

echo "Reference contents:"
ls -la "${REF_GEX}"
echo "Genome size:"
du -sh "${REF_GEX}"
