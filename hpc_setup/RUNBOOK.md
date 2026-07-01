# Runbook — Cardiff Falcon HPC setup for AML × CloneTracker

Step-by-step. Run the commands in order; do not skip ahead. Each step ends with **Expected output** (what success looks like) and **If it fails** (what to do).

**Project:** <your_project_group> (Tonks lab) · **User:** <your_username> · **Cluster:** Falcon (SLURM, Lustre)

---

## Step 0 — Prerequisites (one-time)

On your **Mac**:

- Cardiff VPN connected (GlobalProtect → ras.cf.ac.uk).
- SSH key set up to Falcon (or you'll be prompted for your COGs password each time).

Test:

```bash
ssh <your_username>@falconlogin.cf.ac.uk "echo ok && groups"
```

**Expected:** `ok` and a groups list including `<your_project_group>`.

---

## Step 1 — Push the setup folder to Falcon (Mac side)

From your Mac, in the project folder:

```bash
cd ~/Documents/Jamshid/AML_Cellecta/hpc_setup
bash deploy.sh
```

**Expected:** rsync output ending with "Done. Next, on Falcon: ...".

**If it fails:** check VPN; check the SSH connectivity test above; re-run `deploy.sh`.

---

## Step 2 — SSH into Falcon and source the env

```bash
ssh <your_username>@falconlogin.cf.ac.uk
cd ~/aml_cellecta_setup
source config/falcon_env.sh
echo "FALCON_PROJECT=${FALCON_PROJECT}"
echo "PROJECT_ROOT=${PROJECT_ROOT}"
```

**Expected:**
```
FALCON_PROJECT=<your_project_group>
PROJECT_ROOT=/shared/scratch/<your_project_group>/<your_username>
```

No "WARNING: FALCON_PROJECT is unset..." message.

**If you see the warning:** open `config/falcon_env.sh` in nano/vim and confirm the first line has `FALCON_PROJECT="<your_project_group>"`, not `SCWF000XX`.

---

## Step 3 — Make the project directory tree

```bash
bash scripts/00_make_dirs.sh
```

**Expected:** prints the directory tree under `${PROJECT_ROOT}` (raw, refs, envs, cellranger, cellecta, logs, software). Disk usage line shows ~28K (just inodes).

**If it fails with "Permission denied":** confirm `groups` still includes `<your_project_group>`. If yes, contact ARCCA — it should never deny.

---

## Step 4 — Install Miniconda on scratch

> ⚠️ **Important — conda must live on scratch, not home.** Falcon home is NFS, and conda 26's default libmamba solver writes a SQLite shard cache. NFS doesn't support SQLite's fcntl locking, so any conda operation against NFS throws `sqlite3.OperationalError: disk I/O error`. The install script handles this automatically: it installs to `${PROJECT_ROOT}/software/miniconda3` (Lustre scratch), drops the `defaults` channel (avoids Anaconda ToS prompt), and sets `solver=classic` (avoids the libmamba shard cache entirely).

```bash
bash scripts/01_install_miniconda.sh
source ~/.bashrc          # or open a new shell
which conda                # /shared/scratch/<your_project_group>/<your_username>/software/miniconda3/bin/conda
conda --version
conda config --show channels   # only conda-forge + bioconda
conda config --show solver     # classic
```

If you previously installed conda to `~/miniconda3` (home), the script will detect it and remove it after a 5-second warning before installing fresh on scratch.

**If `which conda` shows a home-dir path:** `~/.bashrc` still has the old init block. Run `bash scripts/01_install_miniconda.sh` again — it will fix the bashrc — then `source ~/.bashrc` and re-check.

---

## Step 5 — Build the two conda environments

```bash
bash scripts/02_create_envs.sh
conda env list
```

**Expected (after ~20–40 min on classic solver):**
```
cellecta   /shared/scratch/<your_project_group>/<your_username>/software/miniconda3/envs/cellecta
scrna      /shared/scratch/<your_project_group>/<your_username>/software/miniconda3/envs/scrna
```

The `scrna` env is large (R 4.4 + Seurat 5.5 + Harmony + DropletUtils + DESeq2 + edgeR + fgsea + clusterProfiler + AUCell + SingleR + celldex + Scanpy + Scrublet + gseapy + scvi-tools). Allow plenty of time. The `cellecta` env is small (cutadapt + starcode + samtools + pysam + biopython + bartab + pycashier) — a few minutes.

Smoke test:

```bash
conda activate cellecta
cutadapt --version
starcode --version
samtools --version | head -1
conda deactivate

conda activate scrna
python -c "import scanpy, scrublet, gseapy, anndata; print('scrna python OK')"
R -e 'suppressPackageStartupMessages({library(Seurat); library(DESeq2); library(fgsea); library(SingleR); library(AUCell); library(harmony); library(DropletUtils)}); cat("scrna R OK\n")'
conda deactivate
```

**If you re-run `02_create_envs.sh`** with envs already present, it falls through to `conda env update -n <name> -f <yaml>` — safe and idempotent. To force a clean rebuild: `conda env remove -n scrna && bash scripts/02_create_envs.sh`. When `conda env remove` asks twice — first to confirm the package list, second to confirm directory deletion — answer `y` to both for a true clean rebuild, or `y` then `n` if you only want to drop the packages but keep the directory.

---

## Step 6 — Install Cell Ranger (manual, ~5 GB)

10x's Cell Ranger isn't a Falcon module. Register and grab the link from https://www.10xgenomics.com/support/software/cell-ranger/downloads (free academic account). Copy the v10.0.0 (or whichever version is in `CELLRANGER_VERSION` in `config/falcon_env.sh`) tarball URL — it's a long signed URL valid for ~1 hour.

On Falcon:

```bash
cd ${PROJECT_ROOT}/software
wget -O cellranger-${CELLRANGER_VERSION}.tar.gz "<paste signed URL here>"
tar -xzf cellranger-${CELLRANGER_VERSION}.tar.gz
rm cellranger-${CELLRANGER_VERSION}.tar.gz
ls cellranger-${CELLRANGER_VERSION}/

# Confirm it's on PATH (falcon_env.sh adds it)
source ~/aml_cellecta_setup/config/falcon_env.sh
which cellranger
cellranger --version
```

**Expected:** the version that matches `${CELLRANGER_VERSION}` in `config/falcon_env.sh` (currently `10.0.0`).

**If `which cellranger` shows nothing:** check `${CELLRANGER_VERSION}` in `config/falcon_env.sh` matches the actual extracted folder name. If you grabbed a different version, edit one line in `falcon_env.sh` and re-source it.

---

## Step 7 — Submit the reference download

```bash
cd ~/aml_cellecta_setup
sbatch slurm/03_download_references.sbatch
squeue -u <your_username>
```

**Expected:** a JobID is printed. `squeue` shows the job in PD (pending) then R (running). Runs ~30–60 min depending on Falcon load. Check progress with:

```bash
tail -f ${PROJECT_ROOT}/logs/ref_download.<JobID>.out
```

When `squeue -u <your_username>` shows the job is gone:

```bash
ls -la ${REF_GEX}
du -sh ${REF_GEX}      # should be ~15 GB
```

**If the job fails:** check `logs/ref_download.<JobID>.err`. Most common cause is wget timeout — re-submit; the script resumes the partial download.

---

## Step 8 — Replace the placeholder CloneTracker construct

While step 7 runs in the background, you have ~30–60 min to get the real CloneTracker XP cassette sequence and update two files. Two ways to obtain it:

1. **Cellecta kit insert / certificate of analysis** — Sian Rizzo or Jamshid will have the PDF that came with the Tonks-lab CloneTracker XP order. The full plasmid map and cassette FASTA are on it.
2. **Cellecta website / customer portal** — the Barcode-3' construct sequence is on the product page (`CloneTracker-XP_pBA438CMV-3UTR-RFP-2A-PuroR.fa` or similar; exact name depends on which reporter+selection you ordered).

Once you have the FASTA, on Falcon:

```bash
cd ~/aml_cellecta_setup/refs/clonetracker

# Backup the placeholder so you can compare format
cp clonetracker_construct.fa clonetracker_construct.fa.placeholder

# Replace it (paste in nano, or scp from your Mac)
nano clonetracker_construct.fa

# Get the new length
NEW_LEN=$(awk 'BEGIN{l=0} /^>/{next} {l+=length($0)} END{print l}' clonetracker_construct.fa)
echo "Construct length: ${NEW_LEN}"

# Update the GTF end coordinate to match
sed -i "s/\texon\t1\t720\t/\texon\t1\t${NEW_LEN}\t/" clonetracker_construct.gtf
grep "^CloneTracker_construct" clonetracker_construct.gtf
```

**Expected:** the GTF line now ends with the real construct length, not 720. The fasta header should be `>CloneTracker_construct` (don't include any spaces or special chars).

**If you don't have the sequence yet:** stop here, email Sian/Jamshid, and resume when you do. Step 9 will refuse to run against the placeholder.

---

## Step 9 — Build the custom Cell Ranger reference

After step 7 has finished and the construct FASTA is real:

```bash
cd ~/aml_cellecta_setup
sbatch slurm/04_build_clonetracker_ref.sbatch
squeue -u <your_username>
```

Runs ~45 min. Watch:

```bash
tail -f ${PROJECT_ROOT}/logs/ct_mkref.<JobID>.out
```

When done:

```bash
ls -la ${REF_CT}
du -sh ${REF_CT}      # ~15 GB
${CELLRANGER_BIN}/cellranger count --help | head -3   # sanity that cellranger is happy
```

**Expected:** `${REF_CT}` exists with `fasta/`, `genes/`, `star/` subfolders. STAR index built without error.

**If it fails with "FASTA length mismatch":** the GTF wasn't updated to match the new construct length — re-run Step 8's `sed` line.

---

## Step 10 — End-of-setup checkpoint

You should now have:

```bash
df -h ${PROJECT_ROOT}                  # plenty of room
du -sh ${PROJECT_ROOT}/refs/*           # two refs, each ~15 GB
ls -la ${PROJECT_ROOT}/software/cellranger-*  # cellranger present
conda env list                          # scrna + cellecta
ls ${PROJECT_ROOT}/raw/                 # empty for now (Globus drop pending)
```

If all four are good, **HPC setup is complete**. Next phase (per-sample analysis) starts only when:

- ICR Genomics (Floriana Manodoro) confirms 10x chemistry (3' v3.1 vs Flex). Email her.
- Globus drop completes — fastqs land in `${PROJECT_ROOT}/raw/`. Then fill the `fastq_prefix` column in `config/samples.tsv`.

---

## Phase 2 (later — once data arrives)

**Don't run these yet.** Listed for completeness so you know what's next.

```bash
# Per-sample Cell Ranger (12 samples, parallelism limited to 4 at once)
sbatch slurm/10_cellranger_count.sbatch

# After at least one Cell Ranger sample finishes, run the barcode extractor on it
sbatch slurm/20_cellecta_extract.sbatch
```

---

## Troubleshooting cheat sheet

| Symptom | Cause | Fix |
|---|---|---|
| `mkdir: Permission denied` under `/shared/scratch/<your_project_group>/` | Not in the <your_project_group> group | `groups` to verify; ask ARCCA if missing |
| `sbatch: error: invalid account specified` | Account directive mismatch | Confirm `#SBATCH --account=<your_project_group>` matches your `groups` |
| `cellranger: command not found` | Path issue or version mismatch | Re-source `config/falcon_env.sh`; confirm `${CELLRANGER_BIN}` matches actual install dir |
| `ImportError` from python tools | Wrong env activated | `conda activate scrna` (or cellecta) then re-run |
| Job sits in PD with `(QOSMaxJobsPerUser)` | Too many submitted | Wait or reduce array concurrency from `%4` to `%2` |
| Job killed with `OOM` | Hit memory limit | Bump `#SBATCH --mem=` in the relevant sbatch file |
| Reference download timed out | wget 4-hour limit | Re-submit; the script resumes partial files |
| `sqlite3.OperationalError: disk I/O error` from conda | Conda installed on NFS home; SQLite + NFS is broken | Re-run `01_install_miniconda.sh` — it auto-detects an old `~/miniconda3` and reinstalls on scratch |
| `sbatch: Invalid account or account/partition combination specified` | SLURM accounting account name differs from the unix group `<your_project_group>` | Run `sacctmgr show user "$USER" -s format=user,defaultaccount,account%30 -p` to find the real account name; set `FALCON_SLURM_ACCOUNT` in `config/falcon_env.sh` to that value, source the file, resubmit |
| `CondaToSNonInteractiveError` for `pkgs/main` or `pkgs/r` | `defaults` channel sneaked into `~/.condarc` or an env yaml | `conda config --remove channels defaults`; verify yamls have only `conda-forge` and `bioconda`; `01_install_miniconda.sh` strips this on every run |

---

## What to send back if you get stuck

When asking for help on any step, paste:

1. The exact command you ran.
2. The full error message (last ~30 lines of stderr).
3. Output of `squeue -u <your_username>` (if SLURM-related).
4. Output of `source config/falcon_env.sh && env | grep -E 'FALCON|PROJECT|REF|CELLRANGER'`.

That makes it ~10x faster to diagnose.
