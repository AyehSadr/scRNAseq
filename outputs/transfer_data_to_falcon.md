# Transferring run1930 fastqs from ICR RDS to Falcon

The run1930 data sits on the ICR RDS Globus endpoint (the same place Jamshid /
Sian have shared previous sequencing batches from). We need it on Falcon at
`/shared/scratch/SCWF00196/c.medas36/fastqs/run1930/`. Globus is the right tool
for this: ~565 M reads × 3 pools is on the order of 200–400 GB; rsync over SSH
will time out, scp will be slow, and the transfer needs to survive the link
flapping overnight.

## Endpoints to identify

| Side | Likely endpoint | Where to confirm |
|------|-----------------|------------------|
| **Source** | ICR RDS / EDS (External Data Sharing) | Email from Jamshid or ICR Genomics — they share the endpoint name and the dataset path when delivery is ready. ICR's [Scientific Computing Service](https://www.icr.ac.uk/research-and-discoveries/other-facilities-and-technology/scientific-computing-service) operates EDS over Globus. |
| **Destination** | Cardiff ARCCA Falcon (or the older Hawk RDS access VM if Falcon's endpoint isn't yet provisioned) | [Supercomputing Wales portal](https://portal.supercomputing.wales/) and `arcca-help@cardiff.ac.uk` — ask for "the Globus collection / endpoint for Falcon scratch". The previous Hawk endpoint exposed `/shared/scratch/` via the [RDS access VM](https://portal.supercomputing.wales/index.php/hawk-cardiff-research-datastore-rds-access-vm/); the Falcon equivalent should be similarly named. |

If Falcon's Globus collection isn't yet exposed, ARCCA may ask you to use
**Globus Connect Personal** on a Cardiff submit node — that turns Falcon's
filesystem into an ad-hoc personal endpoint. Email `arcca-help@cardiff.ac.uk`
first; they'll tell you which.

## Option A — Globus web UI (recommended for a one-off transfer)

1. Browser → <https://app.globus.org/>, log in with **Cardiff University** SSO
   (CRSiD or the same credentials you use for Falcon).
2. Two-pane File Manager:
   - Left pane: search the ICR endpoint name (e.g. "ICR EDS" / "ICR RDS").
     You'll need ICR to have **shared** the run1930 path with you specifically
     — Jamshid / ICR Genomics typically do this when they email "data is
     available". If you can't see it, that share hasn't been set up yet.
   - Right pane: search the Cardiff Falcon endpoint and navigate to
     `/shared/scratch/SCWF00196/c.medas36/fastqs/`. Create `run1930/` if it
     doesn't exist (the UI has a "New folder" button).
3. Select the run1930 directory in the source pane and click "Transfer or Sync
   to..." — under the gear menu enable:
   - "verify file integrity after transfer" (md5)
   - "preserve source file modification times"
4. Click "Start". You'll get an email per transfer task when it finishes.
5. Expect 200–400 GB. At a typical 100 MB/s sustained between research
   networks (JISC + Janet + Globus tuning) that's 30–60 minutes; in practice
   plan for 2–4 hours and let it run overnight if needed.

## Option B — Globus CLI on Falcon (better if you'll do this repeatedly)

Useful if you want to script it or kick off the transfer from a SLURM batch
job. The `globus` CLI is a small Python package.

```bash
# install once into a small dedicated conda env on Falcon
SCRATCH=/shared/scratch/SCWF00196/c.medas36
source ${SCRATCH}/software/miniconda3/etc/profile.d/conda.sh
conda create -n globus python=3.12 -y
conda activate globus
pip install globus-cli

# log in (opens a URL — copy to laptop browser, paste the auth code back)
globus login

# find endpoints
globus endpoint search 'ICR'         # find the source
globus endpoint search 'Cardiff'     # or 'ARCCA' / 'Falcon' / 'Hawk RDS'

# capture UUIDs once you see them
SRC_EP=<icr-endpoint-uuid>
DST_EP=<falcon-endpoint-uuid>

# kick off transfer (recursive, sync mode = checksum, skip already-transferred)
globus transfer "${SRC_EP}:/path/to/run1930/" \
                "${DST_EP}:/shared/scratch/SCWF00196/c.medas36/fastqs/run1930/" \
                --recursive --sync-level=checksum \
                --label "S34 run1930 to Falcon" \
                --notify on,off,inactive

# the command prints a Task ID; check status with:
globus task show <task-id>
globus task wait <task-id>            # blocks until done — useful in scripts
```

## Verifying the transfer

Globus' built-in checksum verify (Option A step 3 / `--sync-level=checksum` in
Option B) catches in-flight corruption. Belt and braces — once it lands on
Falcon, also confirm the file count and total size match what ICR sent:

```bash
SCRATCH=/shared/scratch/SCWF00196/c.medas36
cd ${SCRATCH}/fastqs/run1930

# expected: 12 samples × N lanes × 2 reads (R1, R2) per sample
ls *.fastq.gz | wc -l
du -sh .

# if ICR sent an md5sums.txt alongside the fastqs, check it:
md5sum -c md5sums.txt
```

For S34 with 12 fixed-cell samples (3 OCM pools × 4 OBs each) you'd typically
see:

| Pool | Illumina sample | Approx fastq count (R1+R2 across lanes) |
|------|----------------|-----------------------------------------|
| SRAML7 (Patient 1) | `S34_01_S?` | 8 files (4 lanes × R1/R2) typical |
| SRAML10 (Patient 2) | `S34_02_S14` ✓ | 8 files |
| SRAML13 (Patient 3) | `S34_03_S?` | 8 files |

Total ≈ 24 fastq files if the run was 4-lane. Confirm against ICR's manifest.

## Layout once landed

After the transfer, separate per-pool fastqs into subdirectories so the
`[libraries]` section of each multi config can point to a clean directory.
The manifest TSVs already assume this layout:

```
/shared/scratch/SCWF00196/c.medas36/fastqs/run1930/
├── run1930_raw/                    # original Globus drop
│   ├── S34_01_S??_L001_R1_001.fastq.gz
│   ├── S34_01_S??_L001_R2_001.fastq.gz
│   ├── S34_02_S14_L001_R1_001.fastq.gz
│   ├── ...
├── SRAML7_pool/                    # symlinks to S34_01_*
├── SRAML10_pool/                   # symlinks to S34_02_*
└── SRAML13_pool/                   # symlinks to S34_03_*
```

A throwaway script to do the symlinking once you know the prefixes:

```bash
cd ${SCRATCH}/fastqs/run1930
mkdir -p SRAML7_pool SRAML10_pool SRAML13_pool
ln -sf "$(pwd)"/run1930_raw/S34_01_*.fastq.gz SRAML7_pool/
ln -sf "$(pwd)"/run1930_raw/S34_02_*.fastq.gz SRAML10_pool/
ln -sf "$(pwd)"/run1930_raw/S34_03_*.fastq.gz SRAML13_pool/
```

Symlinks (rather than copies) keep your scratch usage flat and let
`cellranger multi` find each pool's fastqs in isolation.

## What to ask ARCCA / ICR if anything is missing

Email `arcca-help@cardiff.ac.uk`:

> "Hi — I need to transfer ~400 GB of sequencing fastqs from the ICR RDS
> Globus endpoint to my project scratch on Falcon
> (`/shared/scratch/SCWF00196/c.medas36/`). Could you confirm (a) the Globus
> collection / endpoint name for Falcon scratch, and (b) whether I should use
> the built-in collection or run Globus Connect Personal on a submit node?
> My Falcon username is c.medas36, project SCWF00196."

If the ICR side hasn't shared the run1930 directory yet, ping Jamshid or ICR
Genomics:

> "Has run1930 been deposited to ICR RDS for sharing yet? If so, could the
> share be granted to my Globus identity (the one tied to my Cardiff
> University account)? I'd like to pull it directly to Falcon."
