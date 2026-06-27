#!/usr/bin/env bash
# =============================================================================
# generate_multi_csv.sh
#
# Build a cellranger-multi config for one ICR-S34 OCM pool from a small
# manifest. Designed so once Floriana confirms the per-pool OB mapping and
# fastq prefix for SRAML7 and SRAML13, you can produce both configs in seconds
# (and re-emit SRAML10 for sanity-check).
#
# Usage:
#   ./generate_multi_csv.sh manifest.tsv > multi_<pool>.csv
#
# Manifest format — TAB-separated, header line required:
#   pool_id    ref    chemistry    fastq_prefix    fastq_dir
#   patient    sample_id_1    ob_1    description_1
#   patient    sample_id_2    ob_2    description_2
#   patient    sample_id_3    ob_3    description_3
#   patient    sample_id_4    ob_4    description_4
#
# Header keys ('ref', 'chemistry', etc.) are taken from a meta block at the top
# (rows whose first field is one of: pool_id, ref, chemistry, fastq_prefix,
# fastq_dir, create_bam). Sample rows come after, no key name in column 1.
#
# Example invocation:
#   ./generate_multi_csv.sh manifests/SRAML7.tsv > multi_SRAML7.csv
#   ./generate_multi_csv.sh manifests/SRAML10.tsv > multi_SRAML10.csv
#   ./generate_multi_csv.sh manifests/SRAML13.tsv > multi_SRAML13.csv
# =============================================================================

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <manifest.tsv>" >&2
    exit 1
fi
manifest="$1"
[[ -f "$manifest" ]] || { echo "manifest not found: $manifest" >&2; exit 1; }

# defaults
pool_id=""
ref=""
chemistry="SC3Pv4-OCM"
fastq_prefix=""
fastq_dir=""
create_bam="true"

declare -a sample_lines

while IFS=$'\t' read -r f1 f2 f3 f4 _; do
    # skip blank / comment lines
    [[ -z "${f1// }" ]] && continue
    [[ "$f1" =~ ^# ]] && continue

    case "$f1" in
        pool_id)       pool_id="$f2" ;;
        ref)           ref="$f2" ;;
        chemistry)     chemistry="$f2" ;;
        fastq_prefix)  fastq_prefix="$f2" ;;
        fastq_dir)     fastq_dir="$f2" ;;
        create_bam)    create_bam="$f2" ;;
        *)
            # sample row: sample_id<TAB>ob_id<TAB>description
            # (column 1 is sample_id, not a meta key)
            sample_lines+=("${f1},${f2},${f3:-}")
            ;;
    esac
done < "$manifest"

# sanity
for v in pool_id ref fastq_prefix fastq_dir; do
    if [[ -z "${!v}" ]]; then
        echo "manifest missing required key: $v" >&2
        exit 1
    fi
done
if (( ${#sample_lines[@]} == 0 )); then
    echo "no sample rows in $manifest" >&2
    exit 1
fi

# emit the config
cat <<HEADER
# =============================================================================
# Cell Ranger multi config — pool: ${pool_id}
# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ) by generate_multi_csv.sh
# Source manifest: ${manifest}
# =============================================================================

[gene-expression]
ref,${ref}
chemistry,${chemistry}
create-bam,${create_bam}

[libraries]
fastq_id,fastqs,lanes,feature_types
${fastq_prefix},${fastq_dir},any,Gene Expression

[samples]
sample_id,ocm_barcode_ids,description
HEADER

for line in "${sample_lines[@]}"; do
    echo "${line}"
done
