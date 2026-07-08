# Changelog

## [0.5.0] — Positional & Denisovan extensions

Three AADR-native additions that move beyond the single genome-wide Neanderthal
number, plus honest documentation of where the ascertained pseudo-haploid AADR
data cannot support an analysis (an outcome, not a bug).

### Added
- **`archaic/cohort.py`** (`pooled_freq`): RAM-safe *streaming* pooled allele
  frequency, so a cohort of thousands can be pooled over 1.2M SNPs without the
  ~18 GB one-shot read `panel.frequencies` would need. Verified bit-identical to
  `panel.frequencies` for any chunk size (`tests/test_cohort.py`).
- **`archaic/windows.py`**: pure, tested genomic sliding-window aggregator
  (`window_scan`, `robust_z`, `empirical_p`) for position-resolved scans
  (`tests/test_windows.py`).
- **`local_archaic_scan.py`**: genome-wide windowed archaic-affinity landscape
  (desert/peak map) built on the archaic-informative-allele machinery. Chosen
  *instead of* a haplotype HMM caller (Sprime/hmmix/IBDmix), which cannot run on
  ascertained pseudo-haploid AADR data. Honest outcome: the **peak** direction
  recovers known adaptive-introgression loci (BNC2, OCA2-HERC2, FADS, HLA, KRT in
  the top percentiles) as a positive control, but published Neanderthal
  **deserts** are *not* recovered (Mann-Whitney p≈0.47) — windowed archaic-allele
  frequency is ILS-dominated (cf. FADS_REPORT.md). Reported accordingly
  (`LOCAL_ARCHAIC_REPORT.md`, `results/figures/fig_local_archaic_scan.png`).
- **`xchrom_depletion.py`**: X-vs-autosome Neanderthal f4-ratio depletion. The
  autosomal arm reproduces the expected ~2–3%, but the analysis is **inconclusive
  on AADR by construction**: 1240K has *zero* outgroup X genotypes (Chimp.REF,
  Gorilla.REF, Ancestor.REF all blank on X), and HO — the only panel whose Chimp
  covers the X — yields just ~1.7k usable X SNPs, so α(X)'s SE is ~17x the
  autosomal one and floored by X SNP count (pooling cannot help). Refuses to run
  on a panel whose outgroup lacks X. Documented as a data limitation
  (`XCHROM_REPORT.md`).
- **`denisovan_survey.py`** + **`archaic/loci.py::denisovan_informative`**:
  Denisovan-ancestry survey with a *validated positive control*. The pooled
  D_Den statistic recovers the known Denisovan ancestry of Oceanians (Papuan
  Z≈+6), grading to ~0 in West Eurasians/Africans; against that calibrated scale,
  ancient Eurasians are a **controlled null** (no gradient, no Bonferroni
  outliers) — the Denisovan counterpart to `FINDINGS.md`, strengthened by an
  explicit power anchor. Includes an EPAS1 vignette. Outputs
  `results/denisovan_1240k_{survey,outliers}.csv`,
  `results/figures/fig_denisovan_survey.png`, `DENISOVAN_REPORT.md`.
- **Tests** (`tests/test_cohort.py`, `tests/test_windows.py`,
  `tests/test_loci_denisovan.py`): 9 new unit tests, all passing.

## [0.4.0] — West-Eurasian ancestry admixture (Steppe/Yamnaya, WHG, EHG, CHG, …)

### Added
- **`archaic/qpadm.py`**: a constrained ("supervised admixture") qpAdm solver
  (`qpadm_constrained`, simplex-bounded via SciPy SLSQP) alongside the existing
  unconstrained "rotating outgroup" f4-ratio qpAdm (Haak et al. 2015), plus
  `compete_models`/`archaic.ancestry.decompose_best` to rank several candidate
  source models per target. **Performance rewrite**: every f4 term needed by
  the linear system, and by each of the 50 leave-one-block-out jackknife
  replicates, is now computed via one block-sum table per term
  (`_build_system`) instead of being recomputed from scratch inside the
  jackknife loop — ~50x fewer full-genome scans, verified bit-identical to
  the naive per-block recomputation (`tests/test_ancestry.py`). This is what
  makes an 18-target × 5-model competition (90 fits, ~476k SNPs each) run in
  minutes instead of hours.
- **`archaic/ancestry.py`**: a verified West-Eurasian source-population
  library (WHG, EHG, CHG, Anatolia_N, Iran_N, Levant_N/Natufian,
  Steppe_Yamnaya, ANE — each an AADR `group_id` predicate checked against the
  local 1240K release for correctness and zero source/target overlap),
  cohort resolution (kinship-pruned mean-genome profiles), five candidate
  admixture models (`west3` through `deep5`), and `decompose`/`decompose_best`
  high-level entry points.
- **`ancestry_decomposition.py`**: applies the engine to a chronological
  transect of 9 ancient European cohorts (Balkans Neolithic → Imperial Roman)
  plus 9 modern populations, competing all 5 models per target and reporting
  a fixed reference model (`west3`) for direct cross-target comparison,
  alongside each cohort's group-level Neanderthal ancestry. Reproduces the
  Steppe-migration signal with no manual tuning (Corded Ware 73% Steppe,
  Bronze-Age steppe cultures 86%, Sardinia lowest at 11%, Finland/Russia
  highest at ~62-66%). Outputs `results/ancestry/`, `reports/ancestry/`
  (4 figures, `PAPER_ancestry.md`, `Ancestry_admixture_survey.pdf`).
- **`tests/test_ancestry.py`** (7 tests): synthetic-data correctness for the
  new qpAdm/ancestry modules, including a regression guard proving the
  vectorised block-sum jackknife exactly matches a naive brute-force
  recomputation, and a model-competition test confirming the true generating
  model is ranked first.
- `archaic-pipeline ancestry` CLI subcommand (`archaic/cli.py`).

### Changed
- README: new sub-study section for the ancestry-admixture survey.

### Removed
- The personal consumer-DNA analysis capability (archaic and ancestry
  estimates for an individual's own direct-to-consumer genotype file) and all
  associated results/reports have been removed from this repository and its
  history for privacy.

## [0.3.0] — Reproducibility, CLI, and housekeeping

### Added
- **Synthetic-data smoke test** (`archaic/synthetic.py`, `archaic/smoke.py`,
  `tests/test_smoke_synthetic.py`): builds a small, fully synthetic
  AADR-shaped panel (.ind/.snp/.geno) with a known Neanderthal-admixture
  fraction and runs it through the real `Panel`/`PackedGeno`/`stats` code
  path, checking that the estimator recovers a sane, correctly-signed
  archaic-ancestry signal. This is secondary to, and does not replace, the
  AADR-based pipeline or the coalescent-simulation validation
  (`validate_simulation.py`) — its only job is to catch reader/estimator
  plumbing regressions in seconds, in CI or for anyone without a copy of the
  (non-redistributable) AADR data. Run via `pytest` or
  `archaic-pipeline smoke-test`.
- **Unit tests for `archaic/lib_eigenstrat.py`** (`tests/test_lib_eigenstrat.py`,
  11 tests): the packed TGENO/GENO binary reader had zero coverage despite
  being the most fragile, bit-arithmetic-heavy code in the pipeline. Covers
  exact round-trip for both packed layouts, arbitrary/unsorted SNP+individual
  selection, chunked vs. unchunked reads, missing-genotype coding, and
  corrupt-file-size / unknown-magic error paths.
- **CLI entry point** (`archaic/cli.py`, `archaic-pipeline` console script):
  `archaic-pipeline validate|prepare|estimate|...|all --panel 1240k` or
  `archaic-pipeline smoke-test`, instead of `python phase3_estimate.py`.
  Pure dispatch to the existing phase scripts — no behaviour change, just a
  shorter invocation once the package is installed (`pip install -e .`).
- **Logging** (`archaic/log_utils.py`): timestamped, level-controlled
  (`ARCHAIC_LOG_LEVEL` env var) progress logging for the orchestrator
  (`run_pipeline.py`) and the long-running Phase 3 estimation loop
  (`phase3_estimate.py`, ETA per chunk). Phase/report scripts keep using
  plain `print()` for their tabular results and markdown bodies by design —
  that's their actual output, not a diagnostic log.

### Fixed
- `phase1_validate.py`: validation gate G6 (Denisovan check) was computing a
  proper margin test (`den_ok`) but then checking a weaker inline condition
  instead of using it. Now uses the intended check. Re-run confirms the gate
  still passes 7/7.
- `archaic/__init__.py`: `__version__` was stuck at `0.1.0` while
  `pyproject.toml`/`CITATION.cff` had already moved to `0.2.0`. Synced.

### Changed
- Removed dead code flagged by `pyflakes`: unused imports (`os`, `json`,
  `pandas`, `sys`, `csv`, `numpy`, `scipy.stats`, `matplotlib.lines.Line2D`)
  and unused local variables across `etruscan_study.py`, `export_plink.py`,
  `fads_report.py`, `generate_report.py`, `phase2_prepare.py`,
  `phase6_outliers.py`, `phase7_reports.py`, `validate_published.py`,
  `high_archaic_survey.py`, `tools/compare_admixtools.py`.
- Dropped stray `f` prefixes on string literals with no placeholders
  (`fads_report.py`, `phase6_outliers.py`, `validate_published.py`,
  `oase1_bam_pipeline/summarize_segments.py`).
- CI (`.github/workflows/ci.yml`): added Python 3.10 to the test matrix (to
  match the `requires-python = ">=3.10"` floor in `pyproject.toml`), enabled
  pip caching, and added a `pyflakes` lint step so this class of issue is
  caught automatically going forward.
- Split dev-only tooling (`pytest`, `pyflakes`) out of `requirements.txt` into
  a new `requirements-dev.txt`.
- README: added CI/license/Python-version badges.

## [0.2.0] — Etruscan case study
- ADMIXTOOLS 2 concordance validation; PLINK export; qpAdm; kinship
  robustness; Etruscan methods paper (`reports/Etruscan_paper.pdf`).

## [0.1.0] — Initial pipeline
- 9-phase genome-wide archaic-ancestry estimation pipeline; validated against
  published estimates and coalescent simulation.
