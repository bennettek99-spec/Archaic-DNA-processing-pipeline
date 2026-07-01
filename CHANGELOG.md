# Changelog

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
