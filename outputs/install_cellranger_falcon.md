# Installing Cell Ranger on Falcon for the SRAML10/SRAML7/SRAML13 OCM pools

## Why version matters

GEM-X Universal On-Chip Multiplexing (OCM) — the chemistry the ICR Genomics
email confirmed — is **only supported from Cell Ranger v9.0 onwards**. v8.0 added
GEM-X 3' v4, but `cellranger multi` cannot do OB-bead demultiplexing on those
4-plex pools until v9.0. So the install target is **Cell Ranger v9.0 or later**
(latest at time of writing — check the [release notes](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-release-notes)
before downloading).

The chemistry string for the `[gene-expression]` block is `SC3Pv4-OCM`. The
sample sheet's `ocm_barcode_ids` accepts the IDs `OB1`, `OB2`, `OB3`, `OB4` —
Cell Ranger 9 ships the underlying barcode whitelists internally so you don't
have to download them.

## Sizing

Cell Ranger 9 install footprint: ~10 GB (binary + bundled binaries).
GRCh38-2024-A reference: ~20 GB unpacked.
Working directory per `cellranger multi` run: budget ~150–250 GB scratch
per OCM pool (565 M reads in your SRAML10 case → BAMs + intermediate per-OB
matrices; expect SRAML7/13 to be similar).

So three pools concurrently could need ~750 GB scratch headroom. If
`/shared/scratch/SCWF00196/c.medas36/` is tight, run pools serially.

## Step 1 — get the download link

Cell Ranger downloads are gated by a time-limited signed URL — you can't
`wget` a fixed URL. From a browser **on your laptop**:

1. Go to <https://www.10xgenomics.com/support/software/cell-ranger/downloads>.
2. Sign in with your 10x Cloud account (the one you've used for previous
   downloads — Jamshid or Sian's group account if you don't have your own).
3. Pick the Cell Ranger 9.x Linux tarball (~600 MB, name is
   `cellranger-9.X.X.tar.gz`).
4. Right-click the green "Download" button → "Copy Link Address". This is the
   short-lived signed URL.
5. Do the same for the GRCh38 reference under "References" → "Human (GRCh38)
   2024-A". The reference tarball link is on `cf.10xgenomics.com` and is
   non-gated (you can copy it once and reuse).

## Step 2 — install on Falcon

```bash
ssh c.medas36@hawklogin.cf.ac.uk    # or your usual Falcon login

# project scratch
SCRATCH=/shared/scratch/SCWF00196/c.medas36
mkdir -p ${SCRATCH}/software ${SCRATCH}/refs ${SCRATCH}/fastqs ${SCRATCH}/cellranger_runs
cd ${SCRATCH}/software

# 1) Cell Ranger — paste the signed URL from Step 1 in single quotes
#    (signed URLs contain &, ?, =; quoting prevents shell expansion)
wget -O cellranger-9.tar.gz 'PASTE_SIGNED_URL_HERE'
tar -xzvf cellranger-9.tar.gz   # creates cellranger-9.X.X/
rm cellranger-9.tar.gz          # save 600 MB

# 2) GRCh38 reference (URL is stable, no signing)
cd ${SCRATCH}/refs
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
rm refdata-gex-GRCh38-2024-A.tar.gz
```

## Step 3 — make Cell Ranger discoverable

Don't put it in `~/.bashrc` — that runs on every login including SLURM job
nodes which may not need it and slows logins. Create a sourceable env file:

```bash
cat > ${SCRATCH}/software/env_cellranger.sh <<'EOF'
# source this at the top of any SLURM script that calls cellranger
SCRATCH=/shared/scratch/SCWF00196/c.medas36
export PATH="${SCRATCH}/software/cellranger-9.X.X/bin:${PATH}"
export CELLRANGER_REF="${SCRATCH}/refs/refdata-gex-GRCh38-2024-A"
EOF
# update X.X to whatever version you actually downloaded
```

## Step 4 — verify

```bash
source ${SCRATCH}/software/env_cellranger.sh
cellranger --version          # should print 9.X.X
cellranger sitecheck > sitecheck.log    # checks gcc, java, perl, ulimits
cellranger testrun --id=tiny_test       # ~30 min sanity run, ~10 GB scratch
```

If `sitecheck` complains about `ulimit -u` (process limit) being low — common
on shared HPC — add `ulimit -u 8192` to your env file. If it complains about
locale, add `export LC_ALL=C.UTF-8`.

## Step 5 — wire into the multi config

In each `multi_SRAMLx.csv`, the `ref` line should now point at:

```
ref,/shared/scratch/SCWF00196/c.medas36/refs/refdata-gex-GRCh38-2024-A
```

Which is what's already in the manifests I generated. So no further edits
needed if you used `${SCRATCH}/refs/`.

## Step 6 — submit a SLURM job

`cellranger multi` is single-node parallel. A starter SLURM script:

```bash
#!/bin/bash
#SBATCH --job-name=cr_multi_SRAML10
#SBATCH --partition=compute            # or whatever Falcon's GEX partition is
#SBATCH --account=scwf00196
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cr_multi_SRAML10_%j.out

set -euo pipefail
SCRATCH=/shared/scratch/SCWF00196/c.medas36
source ${SCRATCH}/software/env_cellranger.sh

cd ${SCRATCH}/cellranger_runs
cellranger multi \
    --id=SRAML10_multi \
    --csv=${SCRATCH}/configs/multi_SRAML10.csv \
    --localcores=16 \
    --localmem=120
```

For 565 M reads / 4 OB samples you can expect ~12–20 hours wall on a 16-core
node. Halve `--localmem` from `--mem` to leave headroom for Cell Ranger's
forking subprocesses; running OOM late in the pipeline is the most common
Cell Ranger failure mode.

## Common gotchas

- **Don't** run `cellranger` from `/home` on Falcon — the SQLite databases it
  writes for STAR indices choke on the home directory's NFS. Always run from
  `/shared/scratch/`. (You already discovered this with miniconda; same root
  cause.)
- **Don't** auto-detect chemistry. For OCM pools you must set
  `chemistry,SC3Pv4-OCM` explicitly — auto-detect can mistakenly call them
  plain `SC3Pv4` and skip OB demultiplexing.
- The `[libraries]` section's `fastqs` path **must contain only the fastqs
  for that one pool**. If the run delivers all 12 samples' fastqs in a single
  flat directory, sub-symlink each pool's fastqs into a per-pool subdir, e.g.:

  ```bash
  mkdir -p ${SCRATCH}/fastqs/run1930/SRAML10_pool
  ln -s ${SCRATCH}/fastqs/run1930_raw/S34_02_S14_*.fastq.gz \
       ${SCRATCH}/fastqs/run1930/SRAML10_pool/
  ```

## References

- [Cell Ranger v9.0 release notes — GEM-X 3'v4 4-plex / OCM support added](https://www.10xgenomics.com/support/software/cell-ranger/9.0/release-notes/cr-release-notes)
- [Cell Ranger downloads (signed-URL gated)](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
- [Cell Ranger multi config CSV options reference](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts)
- [GEM-X Universal Multiplex (OCM) Barcode lists — KB article (the link Floriana sent)](https://kb.10xgenomics.com/s/article/35032129227405-GEM-X-Universal-Multiplex-or-OCM-Barcode-lists-and-understanding-its-usage)
