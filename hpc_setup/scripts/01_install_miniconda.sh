#!/bin/bash
# Install Miniconda to ${PROJECT_ROOT}/software/miniconda3 (on scratch).
# Falcon home is 50 GB with a tight inode cap; conda envs do not fit there.
# Safe to re-run; idempotent for the install AND for the channel/solver config.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config/falcon_env.sh"

# --- Detect and clean up an old home-directory install --------------------
OLD_HOME_INSTALL="${HOME}/miniconda3"
if [ -d "${OLD_HOME_INSTALL}" ] && [ "${OLD_HOME_INSTALL}" != "${CONDA_BASE}" ]; then
    echo "============================================================"
    echo "Found a previous Miniconda install at ${OLD_HOME_INSTALL}."
    echo "We're moving conda to scratch (${CONDA_BASE}) to avoid home"
    echo "quota / inode issues. Removing the old install in 5s — Ctrl-C"
    echo "to abort if you want to keep it."
    echo "============================================================"
    sleep 5
    rm -rf "${OLD_HOME_INSTALL}" "${HOME}/.conda" "${HOME}/.condarc"
    # Strip out the old conda init block from .bashrc
    if grep -q '>>> conda initialize >>>' "${HOME}/.bashrc"; then
        sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "${HOME}/.bashrc"
    fi
    echo "Old conda removed."
fi

CONDA_BIN="${CONDA_BASE}/bin/conda"

if [ -d "${CONDA_BASE}" ] && [ -x "${CONDA_BIN}" ]; then
    echo "Miniconda already installed at ${CONDA_BASE} — skipping bootstrap."
    "${CONDA_BIN}" --version
else
    mkdir -p "$(dirname "${CONDA_BASE}")"
    cd "$(dirname "${CONDA_BASE}")"
    INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
    if [ ! -f "${INSTALLER}" ]; then
        echo "Downloading Miniconda installer..."
        wget -q "https://repo.anaconda.com/miniconda/${INSTALLER}"
    fi
    echo "Installing Miniconda to ${CONDA_BASE} ..."
    bash "${INSTALLER}" -b -p "${CONDA_BASE}"
    "${CONDA_BIN}" init bash
    rm -f "${INSTALLER}"
fi

# --- Channel config (idempotent) ------------------------------------------
# bioconda + conda-forge only. NOT defaults — conda 25+ requires ToS for that
# channel, and bioconda is designed to work without it.
# See https://bioconda.github.io/#usage
echo "Configuring conda channels..."

# Wipe defaults from BOTH possible condarc locations (user-level + base-level)
for rc in "${HOME}/.condarc" "${CONDA_BASE}/.condarc"; do
    [ -f "${rc}" ] && sed -i '/^[[:space:]]*-[[:space:]]*defaults[[:space:]]*$/d' "${rc}"
done

"${CONDA_BIN}" config --remove-key channels 2>/dev/null || true
"${CONDA_BIN}" config --add channels bioconda
"${CONDA_BIN}" config --add channels conda-forge
"${CONDA_BIN}" config --set channel_priority strict

# --- Solver / pkgs_dirs settings -------------------------------------------
# Use the classic solver (not libmamba) — conda 26's libmamba shard cache
# uses SQLite which is fragile on Lustre. Classic is slower but reliable.
# Pkgs cache stays in CONDA_BASE/pkgs (on scratch), which has plenty of room.
echo "Applying solver settings..."
"${CONDA_BIN}" config --set solver classic

# --- Final state ----------------------------------------------------------
echo
echo "=== Final conda config ==="
"${CONDA_BIN}" config --show channels
"${CONDA_BIN}" config --show solver
"${CONDA_BIN}" config --show envs_dirs 2>/dev/null || true
echo
echo "Done. Open a new shell (or 'source ~/.bashrc') to pick up conda."
echo "Conda lives at: ${CONDA_BASE}"
