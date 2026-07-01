#!/bin/bash
# Push hpc_setup/ from this Mac to ~/scRNAseq on Falcon.
# Run from your Mac, in the hpc_setup folder:
#     bash deploy.sh
#
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Configuration ---
# Set this to your Falcon HPC username (or export FALCON_USER in your environment)
FALCON_USER="${FALCON_USER:-your_username}"
FALCON_HOST="falconlogin.cf.ac.uk"

if [ "${FALCON_USER}" = "your_username" ]; then
    echo "ERROR: Please edit hpc_setup/deploy.sh and set your FALCON_USER username first." >&2
    exit 1
fi

REMOTE_HOST="${FALCON_USER}@${FALCON_HOST}"
REMOTE_PATH="/shared/home1/${FALCON_USER}/scRNAseq/"

echo "Pushing ${HERE}/  →  ${REMOTE_HOST}:${REMOTE_PATH}"
rsync -av --delete \
    --exclude '.DS_Store' \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    "${HERE}/" "${REMOTE_HOST}:${REMOTE_PATH}"

echo
echo "Done. Next, on Falcon:"
echo "  ssh ${REMOTE_HOST}"
echo "  cd ~/aml_cellecta_setup"
echo "  source config/falcon_env.sh"
echo "  cat RUNBOOK.md           # step-by-step from here"
