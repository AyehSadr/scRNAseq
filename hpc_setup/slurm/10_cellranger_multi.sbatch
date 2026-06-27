#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=cr_multi
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-12:00:00
#SBATCH --array=1-1              # Change to 1-3 if 3 pools exist
#SBATCH --output=logs/cr_multi.%A_%a.out
#SBATCH --error=logs/cr_multi.%A_%a.err

set -euo pipefail
SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"

# Extract unique pool_ids from samples.tsv, skipping the header and empty lines
mapfile -t POOLS < <(tail -n +2 "${SAMPLE_SHEET}" | cut -f1 | sort -u | grep -v '^$')

# Check if array task ID is valid
if [ "${SLURM_ARRAY_TASK_ID}" -gt "${#POOLS[@]}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds number of unique pools (${#POOLS[@]})." >&2
    exit 1
fi

# Array is 1-based, bash array is 0-based
POOL_ID="${POOLS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "Processing Pool: ${POOL_ID}"

OUT="${CR_OUT}/${POOL_ID}"
mkdir -p "${OUT}"
cd "${OUT}"

CSV_FILE="${OUT}/multi_${POOL_ID}.csv"

# --- GENERATE MULTI.CSV DYNAMICALLY ---
cat << EOF > "${CSV_FILE}"
[gene-expression]
reference,${REF_GEX}
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
${POOL_ID},${RAW_DIR},Gene Expression

[samples]
sample_id,ocm_barcode_ids,description
EOF

# Extract samples belonging to this pool (col2=sample_id, col9=cmo_id/ocm_id, col8=condition)
tail -n +2 "${SAMPLE_SHEET}" | awk -F'\t' -v pool="${POOL_ID}" '$1 == pool {print $2 "," $9 "," $8}' >> "${CSV_FILE}"

echo "Generated multi.csv:"
cat "${CSV_FILE}"
echo "----------------------------------------"

# Run cellranger multi
cellranger multi \
    --id="${POOL_ID}" \
    --csv="${CSV_FILE}" \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=120

echo "Cell Ranger Multi done for ${POOL_ID}"
ls -la "${POOL_ID}/outs/per_sample_outs" || true
