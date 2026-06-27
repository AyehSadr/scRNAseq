#!/bin/bash
# Build the two conda envs declared in envs/*.yaml.
# Idempotent — uses `mamba env update` if the env already exists.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config/falcon_env.sh"

if ! command -v mamba >/dev/null 2>&1; then
    echo "mamba not found; falling back to conda (slower)."
    SOLVER=conda
else
    SOLVER=mamba
fi

create_or_update () {
    local name="$1" yaml="$2"
    if conda env list | awk '{print $1}' | grep -qx "${name}"; then
        echo ">>> Updating env ${name} from ${yaml}"
        ${SOLVER} env update -n "${name}" -f "${yaml}"
    else
        echo ">>> Creating env ${name} from ${yaml}"
        ${SOLVER} env create -n "${name}" -f "${yaml}"
    fi
}

create_or_update scrna     "${HERE}/../envs/scrna.yaml"
create_or_update cellecta  "${HERE}/../envs/cellecta.yaml"

echo
echo "Envs ready:"
conda env list | grep -E '^(scrna|cellecta)\b' || true
