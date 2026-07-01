# Changelog

## [Unreleased]

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
