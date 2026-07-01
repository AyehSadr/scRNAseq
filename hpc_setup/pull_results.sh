#!/bin/bash
# =============================================================================
# AML CloneTracker XP — Pull Results Script
# -----------------------------------------------------------------------------
# Pulls Seurat analysis results (figures, tables, logs) from Cardiff Falcon HPC
# to this local Mac.
#
# Run from your Mac, in the hpc_setup folder:
#     bash pull_results.sh
#
# Options:
#     bash pull_results.sh --include-rds    (Includes the massive Seurat .RDS objects)
#
# Author:       Ayeh Sadr
# Created:      2026-05-25
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(dirname "${HERE}")"

# --- Configuration ---
# Set this to your Falcon HPC username and project group (or export FALCON_USER / FALCON_PROJECT)
FALCON_USER="${FALCON_USER:-your_username}"
FALCON_PROJECT="${FALCON_PROJECT:-your_project_group}"
FALCON_HOST="falconlogin.cf.ac.uk"

if [ "${FALCON_USER}" = "your_username" ]; then
    echo "ERROR: Please edit hpc_setup/pull_results.sh and set your FALCON_USER username first." >&2
    exit 1
fi

REMOTE_HOST="${FALCON_USER}@${FALCON_HOST}"
REMOTE_SCRATCH="/shared/scratch/${FALCON_PROJECT}/${FALCON_USER}/seurat/"
LOCAL_OUT_DIR="${WORKSPACE_ROOT}/outputs/seurat"

INCLUDE_RDS=false
if [[ "${1:-}" == "--include-rds" ]]; then
    INCLUDE_RDS=true
fi

mkdir -p "${LOCAL_OUT_DIR}"

echo "====================================================================="
echo "Pulling Seurat analysis results from Falcon HPC..."
echo "  Remote: ${REMOTE_HOST}:${REMOTE_SCRATCH}"
echo "  Local:  ${LOCAL_OUT_DIR}"
if [ "$INCLUDE_RDS" = true ]; then
    echo "  Mode:   Downloading EVERYTHING (including massive Seurat .RDS datasets)"
else
    echo "  Mode:   Downloading FIGURES, TABLES, LOGS (excluding massive .RDS files)"
fi
echo "====================================================================="
echo

# Build rsync arguments
RSYNC_ARGS=("-avh" "--progress")

# Exclusions
RSYNC_ARGS+=("--exclude" ".DS_Store")
RSYNC_ARGS+=("--exclude" "*/.DS_Store")

if [ "$INCLUDE_RDS" = false ]; then
    # Exclude all RDS files (both uppercase and lowercase extensions)
    RSYNC_ARGS+=("--exclude" "*.RDS")
    RSYNC_ARGS+=("--exclude" "*.rds")
fi

# Execute rsync
rsync "${RSYNC_ARGS[@]}" "${REMOTE_HOST}:${REMOTE_SCRATCH}/" "${LOCAL_OUT_DIR}/"

echo
echo "====================================================================="
echo "Success! All results have been downloaded successfully to your PC:"
echo "  ${LOCAL_OUT_DIR}/"
echo "====================================================================="
