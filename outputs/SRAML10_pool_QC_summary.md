# SRAML10 pool — pre-alignment QC summary

Source: `multiqc_report.html` (FastQC v0.12.1 + FastQ Screen v0.15.3, Illumina sample `S34_02_S14`).
Pool composition: 4 samples (tubes 11/12/15/16) demultiplexed by bead barcode OB1–OB4.
Note: this MultiQC covers only this one Illumina sample, so all numbers below are pool-level — per-tube metrics will only exist after `cellranger multi`.

## Read structure

| Read | Length | M reads | % unique | % duplicates | FastQC dup level |
|------|--------|---------|----------|--------------|-------------------|
| R1   | 28 bp  | 564.8   | ~51%     | ~49%         | warn              |
| R2   | 90 bp  | 564.8   | ~12%     | ~88%         | fail              |

R1 = 28 bp confirms **GEM-X 3' v4** (16 bp 10x cell barcode + 12 bp UMI). R2 dup level "fail" is expected for scRNA-seq (high pre-collapse duplication is mostly UMI-coverage of expressed transcripts, not library failure) — this is not a real concern.

## FastQC module statuses

| Module                       | R1   | R2   | Action |
|------------------------------|------|------|--------|
| Basic Statistics             | pass | pass | — |
| Per Base Sequence Quality    | pass | pass | — |
| Per Tile Sequence Quality    | pass | warn | benign |
| Per Sequence Quality Scores  | pass | pass | — |
| Per Base Sequence Content    | warn | fail | expected — 10x BC structure on R1, random-hexamer bias on R2 |
| Per Sequence GC Content      | pass | warn | check `human` reference cleanliness post-align |
| Per Base N Content           | pass | pass | — |
| Sequence Length Distribution | pass | pass | — |
| Sequence Duplication Levels  | warn | fail | expected for scRNA-seq |
| Overrepresented Sequences    | warn | warn | likely poly-A / mtRNA — re-check post-align |
| Adapter Content              | pass | pass | — |

Nothing here suggests the pool is unusable. Everything flagged is the standard pattern for 10x 3' GEX: structured R1 + transcript-content R2 always trips a few FastQC modules.

## FastQ Screen — contamination check (R2 only — R1 is 28 bp BC+UMI and unmappable by design)

R2 hits across 8 reference databases (% of reads):

| Category                      | Adapters | E.coli | Human  | Mouse  | No hits | PhiX | Primers | rRNA |
|-------------------------------|---------:|-------:|-------:|-------:|--------:|-----:|--------:|-----:|
| Multiple Hits, Multiple Genomes | 0.12   | 0.00   | 51.56  | 43.92  | n/a     | 0.00 | 0.12    | 2.84 |
| One Hit, Multiple Genomes       | 0.00   | 0.00   | 1.98   | 7.07   | n/a     | 0.00 | 0.00    | 0.00 |
| Multiple Hits, One Genome       | 0.00   | 0.00   | 18.67  | 0.02   | n/a     | 0.00 | 0.01    | 0.00 |
| One Hit, One Genome             | 0.00   | 0.00   | **26.75** | 0.01 | 0.89    | 0.00 | 0.00    | 0.00 |

Read this as: **27% of R2 reads map cleanly to human only**, 19% are repetitive but human-specific, ~50% map to multi-position regions shared between human and mouse (conserved sequence — rRNA, repeat elements, mitochondrial), and **only 0.01% map cleanly to mouse** — i.e. there is no mouse contamination. The high "Multiple Hits, Multiple Genomes" stack against both human and mouse looks alarming but is a known artefact of FastQ Screen + conserved repetitive content; Cell Ranger will collapse this against the GRCh38 reference and the conserved-repeat reads will be excluded by multimapping filters.

The **0.89% no-hits** rate is excellent — suggests low adapter / primer-dimer load.

## Implications for `cellranger multi`

- Set `chemistry,SC3Pv4-OCM` (don't auto-detect — auto-detect can confuse OCM with SC3Pv4).
- 565 M reads across 4 samples ≈ 141 M reads/tube — comfortably above the 20K reads/cell × ~10K cells target for each.
- No reason to filter or trim before alignment; Cell Ranger will handle adapters internally.

## Open questions to send back to ICR Genomics

1. Are the SRAML7 (tubes 3/4/7/8) and SRAML13 (tubes 17/18/19/20) pools each a separate Illumina sample (e.g. `S34_01_S?`, `S34_03_S?`), and does the OB1–OB4 mapping follow the same numerical order as in SRAML10?
2. Is there a separate Cellecta/CloneTracker XP feature-barcoding library for any of these pools, or is this GEX-only?
3. Can we get MultiQC reports for the other two pools as well (or the full run-level report)?

## Notes on per-tube metrics

True per-tube counts (cells, mean reads/cell, median genes/cell, sequencing saturation) will only be available after `cellranger multi` runs and produces a `per_sample_outs/` directory for each of OB1–OB4. The MultiQC numbers here represent the entire pool combined and cannot tell you whether any one of the 4 samples is dominant or under-represented — that's an OB-distribution question that comes out of Cell Ranger.
