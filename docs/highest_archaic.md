# Highest-archaic AADR module

`archaic.highest_archaic` performs a credibility-aware all-sample scan on the
validated Phase 2-4 output, then optionally rereads the AADR EIGENSTRAT panel for
genome-wide sensitivity tests on selected extremes.

## Reused pipeline components

- `archaic.anno`: AADR metadata parsing (extended here for master ID,
  publication/DOI, and library type).
- `archaic.panel.Panel`: memory-safe packed EIGENSTRAT reads and reference
  frequencies.
- `archaic.stats`: the validated Neanderthal f4-ratio, D-statistics, and
  block-jackknife implementation.
- `archaic.refs`: Altai, Vindija, Denisova, Chimp, Mbuti, and Yoruba definitions.
- Phase 2 global exclusions and Phase 4 all-sample estimates/QC covariates.
- The existing 200,000 informative-SNP high-confidence threshold, coverage,
  contamination, and damage thresholds.

## New components

- Broad, high-confidence, and distribution-derived elite analysis sets.
- Raw, lower-bound, Denisovan-affinity, geographic, temporal, and residual ranks.
- Explicit artifact-risk components and five-level credibility classifications.
- Direct candidate recomputation with transversions, Yoruba instead of Mbuti,
  swapped Altai/Vindija roles, leave-one-chromosome/block analyses, block
  bootstrap, and SNP-count-matched repeated subsampling.
- Twelve figures, top-10 candidate reports, a main report, manifests,
  checkpoints, subset mode, and dry-run mode.
- Separate transversion-only, damage-proxy, and geotemporal-residual rankings in
  addition to the twelve required TSV tables.

## Important estimator boundary

The validated Neanderthal f4-ratio is a proportion. `D_Den` is a relative
Denisovan-affinity D-statistic and is not a percentage. Consequently this module
does not add Neanderthal percentage to `D_Den`, and leaves Denisovan and combined
percentages null. A combined percentage requires a separately validated
Denisovan calibration model or second independent Denisovan scale reference.

Transversion-only estimates are a postmortem-damage robustness proxy. Packed
AADR genotypes cannot be re-filtered by terminal-base position, base quality, or
mapping quality; those controls require BAM/CRAM reads.

## Commands

Full run with candidate sensitivity:

```bash
python -m archaic.highest_archaic \
  --aadr-data C:/Users/USER/aadr_v66 \
  --metadata results/phase4_1240k_global_analysis.csv \
  --excluded results/phase2_1240k_global_excluded.csv \
  --config configs/highest_archaic.yaml \
  --output results/highest_archaic \
  --threads AUTO \
  --resume
```

Fast summary-only scan:

```bash
python -m archaic.highest_archaic --skip-sensitivity --resume
```

Dry run and subset examples:

```bash
python -m archaic.highest_archaic --dry-run
python -m archaic.highest_archaic --subset Oase1_d.AG.BY.AA,F6-620.AG.BY.AA
```

Equivalent installed entry points are `archaic-highest` and
`archaic-pipeline highest-archaic`.

Optional segment-table follow-up (summarization only; no unsupported calls):

```bash
python -m archaic.highest_archaic_segments \
  --candidates candidate_ids.txt \
  --segments validated_segments.csv \
  --output results/highest_archaic/segment_followup.tsv
```

The random seed is fixed in the YAML configuration. Candidate sensitivity is
checkpointed in `.sensitivity.json`; `--resume` reuses it only when the config
and relevant input file metadata match.
