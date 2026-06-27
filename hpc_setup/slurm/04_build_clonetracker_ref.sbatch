#!/bin/bash --login
# Account is set dynamically via SBATCH_ACCOUNT env var (sourced from config/falcon_env.sh)
#SBATCH --job-name=ct_mkref
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ct_mkref.%J.out
#SBATCH --error=logs/ct_mkref.%J.err

set -euo pipefail
SCRIPT_DIR="/shared/home1/c.medas36/aml_cellecta_setup"
source "${SCRIPT_DIR}/config/falcon_env.sh"

if [ ! -d "${REF_GEX}" ]; then
    echo "ERROR: 10x reference not found at ${REF_GEX}. Run 03_download_references.sbatch first." >&2
    exit 1
fi

CT_FA="${SCRIPT_DIR}/refs/clonetracker/clonetracker_construct.fa"
CT_GTF="${SCRIPT_DIR}/refs/clonetracker/clonetracker_construct.gtf"

# Sanity: refuse to build a reference against the placeholder sequence
if grep -q "PLACEHOLDER" "${CT_FA}"; then
    echo "ERROR: clonetracker_construct.fa is still the placeholder." >&2
    echo "Replace it with the real CloneTracker XP cassette sequence." >&2
    exit 1
fi

# Check construct length matches GTF
CT_LEN=$(awk 'BEGIN{l=0} /^>/{next} {l+=length($0)} END{print l}' "${CT_FA}")
GTF_END=$(awk -F'\t' '!/^##/{print $5; exit}' "${CT_GTF}")
if [ "${CT_LEN}" != "${GTF_END}" ]; then
    echo "ERROR: FASTA length (${CT_LEN}) does not match GTF end coordinate (${GTF_END})." >&2
    echo "Edit refs/clonetracker/clonetracker_construct.gtf so the exon spans 1..${CT_LEN}." >&2
    exit 1
fi

# Combine genome + construct
WORK="${REF_DIR}/_ct_build"
mkdir -p "${WORK}"
cd "${WORK}"

cp "${REF_GEX}/fasta/genome.fa" combined.fa
cat "${CT_FA}" >> combined.fa

cp "${REF_GEX}/genes/genes.gtf" combined.gtf
cat "${CT_GTF}" >> combined.gtf

# Run cellranger mkref
which cellranger || { echo "cellranger not on PATH; check CELLRANGER_BIN in falcon_env.sh" >&2; exit 1; }

cellranger mkref \
    --genome="$(basename "${REF_CT}")" \
    --fasta=combined.fa \
    --genes=combined.gtf \
    --nthreads="${SLURM_CPUS_PER_TASK}" \
    --memgb=60

# Move into place
mv "$(basename "${REF_CT}")" "${REF_DIR}/"
cd "${REF_DIR}"
rm -rf "${WORK}"

echo "Custom reference built at ${REF_CT}"
ls -la "${REF_CT}"
