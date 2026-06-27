#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=ct_extract
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/ct_extract.%A_%a.out
#SBATCH --error=logs/ct_extract.%A_%a.err

set -euo pipefail
SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"
conda activate cellecta

# Skip header
ROW=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(awk -v r="${ROW}" 'NR==r' "${SAMPLE_SHEET}")

POOL_ID=$(echo "${LINE}" | cut -f1)
SAMPLE_ID=$(echo "${LINE}" | cut -f2)

if [ -z "${POOL_ID}" ] || [ -z "${SAMPLE_ID}" ]; then
    echo "ERROR: Invalid line in samples.tsv at row ${ROW}" >&2
    exit 1
fi

BAM="${CR_OUT}/${POOL_ID}/${POOL_ID}/outs/per_sample_outs/${SAMPLE_ID}/sample_alignments.bam"
OUT_BAM="${CC_OUT}/${SAMPLE_ID}.bam.tsv"
mkdir -p "${CC_OUT}"

# --- BAM-driven extraction ---
if [ -f "${BAM}" ]; then
    python "${SCRIPT_DIR}/scripts/extract_clonetracker_barcodes.py" bam \
        --bam "${BAM}" \
        --out "${OUT_BAM}" \
        --construct-contig CloneTracker_construct
else
    echo "BAM not found for ${SAMPLE_ID} (expected ${BAM}); skipping extraction."
    exit 1
fi

# --- starcode error-correction on clone barcodes ---
# Collapse 1-2 mismatch sequencing errors so duplicate reads of the same
# clone don't fragment into many "near-clones".
if [ -f "${OUT_BAM}" ]; then
    awk -F'\t' 'NR>1 {print $3}' "${OUT_BAM}" | sort -u \
        | starcode -d 2 --print-clusters > "${CC_OUT}/${SAMPLE_ID}.bam.clones.tsv" || true
fi

echo "Done extracting CloneTracker barcodes for ${SAMPLE_ID}"
ls -la "${CC_OUT}" | grep "${SAMPLE_ID}" || true
