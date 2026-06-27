#!/bin/bash
# Push hpc_setup/ from this Mac to ~/aml_cellecta_setup on Falcon.
# Run from your Mac, in the hpc_setup folder:
#     bash deploy.sh
#
# Falcon login (c.medas36@falconlogin.cf.ac.uk) requires the Cardiff VPN.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REMOTE_HOST="c.medas36@falconlogin.cf.ac.uk"
REMOTE_PATH="/shared/home1/c.medas36/aml_cellecta_setup/"

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
