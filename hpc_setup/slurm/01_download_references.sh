#!/bin/bash
# ============================================================
# Title:        Download scRNA-seq reference atlases for AML annotation
# Project:      AML_Cellecta — conventional vs co-culture, Ara-C response
# Author:       Ayeh Sadr
# Created:      2026-05-23
# Last update:  2026-05-23
# Input:        config/references.yml  — accessions, URLs, target subdirs
#               config/paths.yml       — REFS_DIR (per-machine, gitignored)
# Output:       ${REFS_DIR}/van_galen_2019/GSE116256_RAW.tar  + extracted
#               ${REFS_DIR}/baryawno_2019/GSE128423_RAW.tar   + extracted
#               ${REFS_DIR}/tikhonova_2019/GSE108892_RAW.tar  + extracted
#               ${REFS_DIR}/msigdb/h.all.v2023.2.Hs.symbols.gmt
#               ${REFS_DIR}/manifest_references.tsv           (size + md5 per file)
#               logs/01_download_references_<timestamp>.log
# Depends on:   bash 4+, wget, tar, gzip, md5sum, yq (v4+, https://github.com/mikefarah/yq)
# Notes:        Idempotent — re-running skips any file already present with non-zero
#               size. md5 is recorded in the manifest; if expected_md5 is set in
#               references.yml it is verified after download.
#
#               Triana 2021, Zeng 2025, and the Azimuth BM reference are intentionally
#               left out of this shell script — they are easier and more reliable to
#               install from inside R (SeuratData / Azimuth / direct figshare links
#               that change). They will be handled in 01b_install_r_references.R.
#
#               Zeng 2025 specifically: confirm the Zenodo DOI from the paper's
#               data-availability statement and add to config/references.yml.
#
#               Run on a Falcon login/transfer node — most clusters block outbound
#               HTTP from compute partitions, so sbatch is unnecessary and may fail.
# ============================================================

set -euo pipefail

# ---- 0. Config -----------------------------------------------------------

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
TIMESTAMP="$(date +%Y-%m-%d_%H%M%S)"

PATHS_YML="${PROJECT_DIR}/config/paths.yml"
REFS_YML="${PROJECT_DIR}/config/references.yml"

LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"

# Mirror everything to the log file
exec > >(tee -a "${LOG_FILE}") 2>&1

# Hard-fail early if tooling is missing rather than mid-download
for cmd in wget tar md5sum yq; do
  command -v "${cmd}" >/dev/null || { echo "ERROR: '${cmd}' not on PATH"; exit 1; }
done

[[ -f "${PATHS_YML}" ]] || { echo "ERROR: missing ${PATHS_YML} — copy paths.yml.example"; exit 1; }
[[ -f "${REFS_YML}"  ]] || { echo "ERROR: missing ${REFS_YML}"; exit 1; }

REFS_DIR="$(yq -r '.refs_dir' "${PATHS_YML}")"
[[ -z "${REFS_DIR}" || "${REFS_DIR}" == "null" ]] && {
  echo "ERROR: refs_dir not set in ${PATHS_YML}"; exit 1; }

mkdir -p "${REFS_DIR}"
MANIFEST="${REFS_DIR}/manifest_references.tsv"
[[ -f "${MANIFEST}" ]] || \
  printf "reference\tfile\tbytes\tmd5\tdownloaded_at\n" > "${MANIFEST}"

echo "[$(date)] === ${SCRIPT_NAME} start ==="
echo "[$(date)] PROJECT_DIR = ${PROJECT_DIR}"
echo "[$(date)] REFS_DIR    = ${REFS_DIR}"
echo "[$(date)] LOG_FILE    = ${LOG_FILE}"

# ---- 1. Helpers ----------------------------------------------------------

download_if_missing () {
  # $1 = URL
  # $2 = absolute target path
  # $3 = reference key (for manifest)
  # $4 = expected md5 (optional, "" to skip verification)
  local url="$1" target="$2" ref_key="$3" expected_md5="${4:-}"
  local dir; dir="$(dirname "${target}")"
  mkdir -p "${dir}"

  if [[ -s "${target}" ]]; then
    echo "  [skip] ${target} ($(du -h "${target}" | cut -f1))"
    return 0
  fi

  echo "  [get ] ${url}"
  # --continue resumes a half-downloaded *.part rather than restarting
  wget --continue --tries=3 --timeout=60 --quiet --show-progress \
       -O "${target}.part" "${url}"
  mv "${target}.part" "${target}"

  local actual_md5 bytes
  actual_md5="$(md5sum "${target}" | awk '{print $1}')"
  bytes="$(stat -c%s "${target}" 2>/dev/null || stat -f%z "${target}")"

  if [[ -n "${expected_md5}" && "${actual_md5}" != "${expected_md5}" ]]; then
    echo "ERROR: md5 mismatch for ${target}"
    echo "  expected ${expected_md5}"
    echo "  got      ${actual_md5}"
    exit 1
  fi

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${ref_key}" "${target}" "${bytes}" "${actual_md5}" "$(date -Iseconds)" \
    >> "${MANIFEST}"
}

extract_tar_once () {
  # $1 = archive
  # $2 = output dir
  # $3 = marker filename inside outdir — if it exists, skip extraction
  local archive="$1" outdir="$2" marker="$3"
  mkdir -p "${outdir}"
  if [[ -e "${outdir}/${marker}" ]]; then
    echo "  [skip] ${archive} already extracted (${marker} present)"
    return 0
  fi
  echo "  [tar ] ${archive} -> ${outdir}"
  tar -xf "${archive}" -C "${outdir}"
}

fetch_block () {
  # $1 = top-level key in references.yml (e.g. 'van_galen_2019')
  local key="$1"
  local enabled url subdir filename expected_md5 extract_marker
  enabled="$(yq -r ".${key}.enabled // true" "${REFS_YML}")"
  if [[ "${enabled}" != "true" ]]; then
    echo "[$(date)] --- ${key}: disabled in references.yml — skipped ---"
    return 0
  fi
  url="$(yq -r ".${key}.url"      "${REFS_YML}")"
  subdir="$(yq -r ".${key}.subdir" "${REFS_YML}")"
  filename="$(yq -r ".${key}.filename // \"\"" "${REFS_YML}")"
  expected_md5="$(yq -r ".${key}.md5 // \"\"" "${REFS_YML}")"
  extract_marker="$(yq -r ".${key}.extract_marker // \"\"" "${REFS_YML}")"

  [[ -z "${filename}" || "${filename}" == "null" ]] && filename="$(basename "${url}")"

  local outdir="${REFS_DIR}/${subdir}"
  local target="${outdir}/${filename}"

  echo "[$(date)] --- ${key} ---"
  download_if_missing "${url}" "${target}" "${key}" "${expected_md5}"

  if [[ -n "${extract_marker}" && "${extract_marker}" != "null" ]]; then
    extract_tar_once "${target}" "${outdir}" "${extract_marker}"
  fi
}

# ---- 2. Fetch each reference --------------------------------------------

for ref_key in van_galen_2019 baryawno_2019 tikhonova_2019 msigdb_hallmark; do
  fetch_block "${ref_key}"
done

# ---- 3. Summary ---------------------------------------------------------

echo ""
echo "[$(date)] === Download summary ==="
du -sh "${REFS_DIR}"/*/ 2>/dev/null | sort -h || true
echo ""
echo "[$(date)] Manifest: ${MANIFEST}"
echo "[$(date)] Log:      ${LOG_FILE}"
echo "[$(date)] === ${SCRIPT_NAME} done ==="
