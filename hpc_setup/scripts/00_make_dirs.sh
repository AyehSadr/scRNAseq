#!/bin/bash
# Run on the Falcon login node.
# Creates the project directory tree on scratch.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config/falcon_env.sh"

if [ "${FALCON_PROJECT}" = "SCWF000XX" ] || [ -z "${FALCON_PROJECT}" ]; then
    echo "ERROR: edit hpc_setup/config/falcon_env.sh and set FALCON_PROJECT first." >&2
    exit 1
fi

if [ ! -d "/shared/scratch/${FALCON_PROJECT}" ]; then
    echo "ERROR: /shared/scratch/${FALCON_PROJECT} does not exist." >&2
    echo "Confirm in Coldfront that you are added to the project, then re-run." >&2
    exit 1
fi

mkdir -p \
    "${RAW_DIR}" \
    "${REF_DIR}" \
    "${REF_DIR}/clonetracker" \
    "${ENV_DIR}" \
    "${CR_OUT}" \
    "${CC_OUT}" \
    "${LOG_DIR}" \
    "${PROJECT_ROOT}/software"

echo "Directory tree under ${PROJECT_ROOT}:"
ls -la "${PROJECT_ROOT}"
echo
echo "Disk usage of project root:"
du -sh "${PROJECT_ROOT}"
echo
echo "Quota check (scratch):"
df -h "${PROJECT_ROOT}" || true
