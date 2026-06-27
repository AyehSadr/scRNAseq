#!/usr/bin/env python3
"""
Extract CloneTracker XP barcodes from 10x scRNA-seq R2 reads (no enrichment library).

Two complementary modes:

  1) bam   — parse a Cell Ranger BAM, keep reads with a valid CB:Z and either:
             (a) aligned to the CloneTracker_construct contig, or
             (b) unmapped/multimapped but whose read sequence contains the
                 FBP1 anchor.
  2) fastq — parse raw R1+R2 directly, pair on read name, extract CB+UMI from
             R1 and the BC14-spacer-BC30 cassette from R2.

Outputs a TSV: cell_barcode, umi, clone_barcode, n_reads (after collapse).

Usage examples
--------------
# Mode 1: BAM-driven
python extract_clonetracker_barcodes.py bam \
    --bam   ${PROJECT_ROOT}/cellranger/P1_T03_noStroma_noAraC/outs/possorted_genome_bam.bam \
    --out   ${PROJECT_ROOT}/cellecta/P1_T03_noStroma_noAraC.tsv \
    --construct-contig CloneTracker_construct

# Mode 2: fastq-driven sanity check
python extract_clonetracker_barcodes.py fastq \
    --r1   ${PROJECT_ROOT}/raw/S34_02_S17_L001_R1_001.fastq.gz \
    --r2   ${PROJECT_ROOT}/raw/S34_02_S17_L001_R2_001.fastq.gz \
    --out  ${PROJECT_ROOT}/cellecta/P1_T03_noStroma_noAraC.fastq.tsv

The FBP1 anchor (constant region adjacent to the cassette) defaults to the
sequence shown in the CloneTracker workflow slides:

    CCGACCACCGAACGCAACGCACGCA

The cassette structure assumed: BC14 - 6nt spacer - BC30. Override with --layout.
"""
from __future__ import annotations
import argparse
import gzip
import sys
import re
from collections import defaultdict
from pathlib import Path

# ---- Defaults (matching CloneTracker XP) ---------------------------------
DEFAULT_ANCHOR = "CCGACCACCGAACGCAACGCACGCA"     # FBP1 site, immediately 3' of BC30
DEFAULT_BC14_LEN = 14
DEFAULT_SPACER_LEN = 6
DEFAULT_BC30_LEN = 30
DEFAULT_CB_LEN = 16
DEFAULT_UMI_LEN = 12

# Reverse complement helper
_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def find_cassette(seq: str, anchor: str, bc14: int, sp: int, bc30: int):
    """Return (clone_barcode, strand) if cassette found, else (None, None).
    The cassette in the cDNA orientation is:  BC14 - spacer - BC30 - <anchor>
    so we look for the anchor and take the bc14+sp+bc30 nt immediately 5' of it.
    Also tries the reverse complement.
    """
    span = bc14 + sp + bc30
    # forward
    idx = seq.find(anchor)
    if idx >= span:
        cassette = seq[idx - span: idx]
        return f"{cassette[:bc14]}_{cassette[bc14+sp:]}", "+"
    # reverse complement
    rc_anchor = revcomp(anchor)
    idx = seq.find(rc_anchor)
    if idx >= 0 and idx + len(rc_anchor) + span <= len(seq):
        cassette_rc = seq[idx + len(rc_anchor): idx + len(rc_anchor) + span]
        cassette = revcomp(cassette_rc)
        return f"{cassette[:bc14]}_{cassette[bc14+sp:]}", "-"
    return None, None


# ---- BAM mode ------------------------------------------------------------
def run_bam(args):
    import pysam  # noqa: WPS433 — lazy import so fastq mode works without pysam
    bam = pysam.AlignmentFile(args.bam, "rb")
    counts = defaultdict(int)
    seen = 0; kept = 0

    target_tid = None
    if args.construct_contig:
        try:
            target_tid = bam.get_tid(args.construct_contig)
        except KeyError:
            print(f"Warning: contig {args.construct_contig} not in BAM header — skipping aligned-mode filter", file=sys.stderr)

    for rec in bam.fetch(until_eof=True):
        seen += 1
        if not rec.has_tag("CB"):
            continue
        cb = rec.get_tag("CB")
        if cb == "-":
            continue
        ub = rec.get_tag("UB") if rec.has_tag("UB") else "-"
        seq = rec.query_sequence or ""
        # Two ways in: aligned to construct OR sequence contains anchor
        on_construct = (target_tid is not None and rec.reference_id == target_tid)
        clone, _strand = find_cassette(seq, args.anchor, args.bc14, args.spacer, args.bc30)
        if not (on_construct or clone):
            continue
        if clone is None and on_construct:
            # try anchor on the reverse strand of the BAM record
            clone, _ = find_cassette(revcomp(seq), args.anchor, args.bc14, args.spacer, args.bc30)
        if clone is None:
            continue
        counts[(cb, ub, clone)] += 1
        kept += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("cell_barcode\tumi\tclone_barcode\tn_reads\n")
        for (cb, ub, clone), n in sorted(counts.items(), key=lambda kv: -kv[1]):
            fh.write(f"{cb}\t{ub}\t{clone}\t{n}\n")
    print(f"BAM records seen: {seen}; barcode-bearing: {kept}; unique (CB,UMI,clone): {len(counts)}", file=sys.stderr)


# ---- fastq mode ----------------------------------------------------------
def run_fastq(args):
    cb_len, umi_len = args.cb_len, args.umi_len
    counts = defaultdict(int)
    seen = 0; kept = 0
    with open_maybe_gz(args.r1) as f1, open_maybe_gz(args.r2) as f2:
        while True:
            h1 = f1.readline();  s1 = f1.readline();  _ = f1.readline();  _ = f1.readline()
            h2 = f2.readline();  s2 = f2.readline();  _ = f2.readline();  _ = f2.readline()
            if not h1 or not h2:
                break
            seen += 1
            s1 = s1.strip(); s2 = s2.strip()
            if len(s1) < cb_len + umi_len:
                continue
            cb  = s1[:cb_len]
            umi = s1[cb_len: cb_len + umi_len]
            clone, _ = find_cassette(s2, args.anchor, args.bc14, args.spacer, args.bc30)
            if clone is None:
                continue
            counts[(cb, umi, clone)] += 1
            kept += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("cell_barcode\tumi\tclone_barcode\tn_reads\n")
        for (cb, ub, clone), n in sorted(counts.items(), key=lambda kv: -kv[1]):
            fh.write(f"{cb}\t{ub}\t{clone}\t{n}\n")
    print(f"Read pairs seen: {seen}; barcode-bearing: {kept}; unique (CB,UMI,clone): {len(counts)}", file=sys.stderr)


# ---- CLI -----------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--anchor",  default=DEFAULT_ANCHOR, help=f"Constant region 3' of cassette (default: {DEFAULT_ANCHOR})")
    common.add_argument("--bc14",    type=int, default=DEFAULT_BC14_LEN)
    common.add_argument("--spacer",  type=int, default=DEFAULT_SPACER_LEN)
    common.add_argument("--bc30",    type=int, default=DEFAULT_BC30_LEN)

    pb = sub.add_parser("bam", parents=[common])
    pb.add_argument("--bam", required=True)
    pb.add_argument("--out", required=True)
    pb.add_argument("--construct-contig", default="CloneTracker_construct")
    pb.set_defaults(func=run_bam)

    pf = sub.add_parser("fastq", parents=[common])
    pf.add_argument("--r1", required=True)
    pf.add_argument("--r2", required=True)
    pf.add_argument("--out", required=True)
    pf.add_argument("--cb-len",  type=int, default=DEFAULT_CB_LEN)
    pf.add_argument("--umi-len", type=int, default=DEFAULT_UMI_LEN)
    pf.set_defaults(func=run_fastq)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)
