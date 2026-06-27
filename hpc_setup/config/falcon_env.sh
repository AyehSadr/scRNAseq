# Sourced at the top of every SLURM script and login-node script.
# Edit the two values below once your Falcon project is approved.

# --- Falcon project ---------------------------------------------------------
export FALCON_PROJECT="SCWF00196"            # unix group name, confirmed via `groups` on 2026-05-04
export FALCON_USER="c.medas36"

# SLURM accounting account name — usually matches FALCON_PROJECT but can differ.
# Verify with: sacctmgr show user "$USER" -s format=user,defaultaccount,account%30 -p
# If you get "Invalid account or account/partition combination" on sbatch, set this
# to whatever account `sacctmgr` shows next to your username.
export FALCON_SLURM_ACCOUNT="scwf00196_a_tonks_192"

# Tell SLURM to use this account by default — SBATCH_ACCOUNT is honored by sbatch
# without any --account directive in the script. Single source of truth.
export SBATCH_ACCOUNT="${FALCON_SLURM_ACCOUNT}"

# --- Paths derived from the above ------------------------------------------
export PROJECT_ROOT="/shared/scratch/${FALCON_PROJECT}/${FALCON_USER}"
export RAW_DIR="${PROJECT_ROOT}/raw"
export REF_DIR="${PROJECT_ROOT}/refs"
export ENV_DIR="${PROJECT_ROOT}/envs"
export CR_OUT="${PROJECT_ROOT}/cellranger"
export CC_OUT="${PROJECT_ROOT}/cellecta"
export LOG_DIR="${PROJECT_ROOT}/logs"
export SCRIPT_DIR="${HOME}/aml_cellecta_setup"   # where this repo lives on Falcon
export SAMPLE_SHEET="${SCRIPT_DIR}/config/samples.tsv"

# --- Conda -----------------------------------------------------------------
# Miniconda lives on SCRATCH (not home). Falcon home is 50 GB with a tight
# inode cap; conda envs have ~10K+ small files each and quickly trip "disk
# I/O error" sqlite failures when installed to home. Scratch has 3 TB / 1M
# inodes per project, plenty of room.
export CONDA_BASE="${PROJECT_ROOT}/software/miniconda3"
# shellcheck disable=SC1091
if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi

# --- 10x Cell Ranger -------------------------------------------------------
# Falcon currently does not ship Cell Ranger as a module. Two install options:
#   (a) module load (if/when ARCCA adds it):
#       module load cellranger/<version>   # if/when ARCCA adds it
#   (b) self-installed under ${PROJECT_ROOT}/software/cellranger-${CELLRANGER_VERSION}/:
export CELLRANGER_VERSION="10.0.0"
export CELLRANGER_BIN="${PROJECT_ROOT}/software/cellranger-${CELLRANGER_VERSION}"
export PATH="${CELLRANGER_BIN}:${PATH}"

# --- Reference layout ------------------------------------------------------
export REF_GEX="${REF_DIR}/refdata-gex-GRCh38-2024-A"
export REF_CT="${REF_DIR}/refdata-GRCh38-CloneTracker"

# --- Sanity check ----------------------------------------------------------
if [ "${FALCON_PROJECT}" = "SCWF000XX" ] || [ -z "${FALCON_PROJECT}" ]; then
    echo "WARNING: FALCON_PROJECT is unset or still the placeholder." >&2
    echo "Edit hpc_setup/config/falcon_env.sh before submitting jobs." >&2
fi
