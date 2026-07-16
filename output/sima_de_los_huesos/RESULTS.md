# Sima de los Huesos: laptop-safe AADR-compatible scan

## Verdict

The published nuclear BAMs are public and small enough for laptop work. The
five BAMs from ENA study `PRJEB10597` total **98,726,265 bytes** and all five
archive MD5 checks passed.

This run is **not eligible for the repository's full Phase 2-9 ancient-Eurasian
pipeline**. Even the combined pseudo-haploid call set has only 25,765 callable
1240K sites, below the pipeline's 30,000-site entry floor, and the Sima
individuals are ~430,000 years old rather than members of the pipeline's
modern-human target cohort. Treat this as a targeted compatibility scan, not a
whole-genome ancestry estimate.

## Compatible Phase-3 statistics

Calls were made at local AADR v66.1 1240K sites from the published,
duplicate-removed L35/MQ30 BAMs, retaining bases with base quality >=30. One
base per site was deterministically selected (pseudo-haploid) and passed to the
same f4-ratio / D-statistic functions used in `phase3_estimate.py`.

| Pooled Sima call set | All SNPs | Transversions only |
| --- | ---: | ---: |
| Callable 1240K sites | 25,765 | 5,204 |
| f4-ratio alpha (Altai/Vindija scale) | 6.73% +/- 2.56% | 9.85% +/- 5.25% |
| f4-ratio informative sites | 11,427 | 3,123 |
| Neanderthal D Z-score | 2.79 | 1.62 |
| Denisovan D Z-score | 2.55 | 0.47 |

The all-SNP values are transition-sensitive; the transversion-only repeat has
far fewer sites and neither D-statistic reaches the conventional |Z| >=3
threshold. Therefore these numbers do **not** establish a percentage of
Neanderthal or Denisovan admixture. They are compatible with the published
conclusion that the nuclear data place the Sima hominins on the Neanderthal
lineage, while avoiding an overclaim from a sparse, very ancient sample.

## Reproduction

```powershell
$env:PYTHONPATH='tmp/sima_deps'
& 'C:\Users\benne\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  tools\sima_de_los_huesos_scan.py `
  --bam-dir data\sima_de_los_huesos_prjeb10597 `
  --output output\sima_de_los_huesos
```

Full numerical output is in `sima_aadr_f4_results.csv`; source checksums and
filters are recorded in `run_manifest.json`.
