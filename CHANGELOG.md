# Changelog

## Unreleased

### Changed

- Repository layout. The 32 flat entry-point scripts moved from the repository
  root into `scripts/`, and the 12 report and policy markdown files moved into
  `docs/` (authored) and `docs/studies/` (generated). Each moved script anchored
  its paths with `os.path.dirname(os.path.abspath(__file__))` meaning "the
  repository root", so all 51 anchors were re-based one level; the report
  generators now write into `docs/studies/`. `run_pipeline.py` deliberately
  keeps a `scripts/`-level anchor because it locates its sibling phase scripts.
  Commands change from `python phase3_estimate.py` to
  `python scripts/phase3_estimate.py`; the five console entry points are
  unaffected. The repository root now holds nine files instead of forty-four.
- Upstream phase intermediates are no longer tracked: phase-2 metadata and
  sample manifests, phase-3 estimate tables, the phase-5 PCA table, phase-6
  residuals, and the highest-archaic outlier dump (27 MB). The phase-4 analysis
  tables stay tracked because every published study reads them. The boundary is
  documented in `docs/DATA.md`.
- README leads with the two papers — the global AADR survey and the Iron Age
  Etruscan study — and carries a quickstart that ends at a passing smoke test.

### Fixed

- Simulated introgression truth was over-attributed. Reading true tracts from
  msprime migration records takes each record's interval to be inherited by
  every descendant sample, but that interval is the span of the ancestral
  lineage when it moved and is far wider than the segment any one modern sample
  receives. The measured archaic fraction therefore depended on how much
  sequence was simulated — 0.158 at 10 Mb against 0.093 at 30 Mb — and a
  1400-generation pulse came back as 921 and 1180 generations. Truth now comes
  from a census placed above every archaic pulse plus `link_ancestors`, which
  recovers the simulated proportion at 0.0400 ± 0.0089 against a target of
  0.0400 with no dependence on sequence length. The same defect is present in
  `msprime_backend._extract`, which is left in place because the tract-level
  M1–M10 outputs were validated against it, but is now documented as unusable
  as truth.
- The same extraction was also quadratic — it rebuilt a tree per migration
  record, and a 50 Mb replicate never finished. The census path is a single
  pass.
- `posterior_decode` accepted a `threshold` argument and ignored it, so callers
  could believe they were thresholding a decode when they were not.
  Thresholding belongs to `extract_runs`, which lets one decode be
  re-thresholded cheaply — the mechanism the specificity check relies on.
- `caller_calibration.invert` only rejected an exactly-zero slope, so a
  near-degenerate calibration curve produced a confident nonsense estimate
  instead of an error. It now tests the span the curve covers against the scale
  of the decay values.
- `fit_hmm` raised `NameError` rather than a useful error for `max_iter=0`, and
  accepted inputs too short to fit.
- `run_replicate` and `run_scenario` defaulted to different variant densities
  (1.0 and 0.40) and different sequence lengths, so a run could silently sit at
  a different operating point from the curve it inverts on. Both now share
  named constants.
- `pyflakes` lint failure on `main`: 26 fragments of placeholder-free f-strings
  in the Neanderthal source-contrast report text.

### Added

- New `archaic_admixture_dating/skov_hmm.py`: a Skov-type two-state Poisson HMM
  over outgroup-private variant counts, with Baum-Welch fitting, posterior
  decoding, run extraction, and truncated-exponential decay estimation. The
  forward-backward recursion is numba-JIT'd, which takes a whole-genome
  individual from 43 hours to about 5 seconds and is what makes calibration
  laptop-feasible at all. Reports the HMM's own admixture-time parameter and
  the decoded run decay separately, because their ratio is the decoder
  inflation the calibration exists to measure.
- New `archaic_admixture_dating/genotype_simulation.py`: genotype-level
  simulation of the Skov observation process under the published Jacobs et al.
  (2019) demography from stdpopsim, with the Denisovan pulse time as a free
  parameter. The published model encodes a young pulse as a finding, so running
  it unmodified would assume what is under test; everything else stays at its
  published value. Pulses older than the Papuan/Ghost merge are retargeted onto
  the Papuan ancestor.
- New `archaic_admixture_dating/caller_calibration.py` and the
  `papuan_denisovan_v1/run_calibration.py` entry point: sweep the true pulse
  time, run the whole chain to a decoded decay, fit the calibration curve, and
  invert it on the real measurement with a replicate bootstrap. This replaces
  the tract-level observation model, whose error was length-preserving and
  therefore structurally unable to reproduce posterior-decoding run inflation.
  20 unit tests.
- New `archaic/source_contrast.py`: a reusable layer for asking *which*
  Neanderthal population a cohort descends from, rather than how much
  Neanderthal ancestry it has. Provides `D_VA = D(X, Yoruba; Vindija, Altai)`, a
  Vindija/Altai-symmetric normaliser `D_NEA = D(X, Yoruba; NeaAvg, Chimp)` that
  measures Neanderthal quantity without responding to its source, per-cohort
  block tables, a **paired** block jackknife that cancels the archaic genomes'
  shared sampling noise (empirically 1.9x tighter than quadrature), an
  origin-constrained single-source fit with jackknifed residuals, a detection-
  limit calculator, and technical-covariate diagnostics. 15 unit tests, no AADR
  data needed.
- New `archaic.cohort.pooled_freq_multi`: pooled allele frequencies for many
  overlapping cohorts in a single streaming pass over the genotypes, so a study
  spanning most of the AADR traverses the `.geno` once instead of once per
  cohort.
- New `neanderthal_source.py` study (`archaic-pipeline neanderthal-source`):
  Altai-versus-Vindija affinity across 10,954 unique Eurasian ancient genomes in
  41 dated cohorts plus 18 present-day anchors, with a pooled-frequency cache so
  reruns take seconds. Recovers both anchors — every Neanderthal-carrying cohort
  is displaced towards Vindija (Prufer et al. 2017), and the Denisovan genome is
  pulled towards Altai (`-0.147`, Z = -6.7) as its ~1% Denisovan-related
  ancestry predicts. Upper Palaeolithic Europeans and Palaeolithic north-east
  Asians are indistinguishable (`-0.0009 +/- 0.0060`), as are present-day French
  and Han (`-0.0008 +/- 0.0056`); a single proportional relation
  `D_VA = 2.24 x D_NEA` fits 49 of 53 cohorts within 2 SE and none departs after
  Bonferroni. **Stated detection limit: 0.0098 in D_VA units, 13% of the total
  Vindija-over-Altai signal (2% for the best-powered pairs).** An age-correlated
  deficit in the oldest cohorts is reported explicitly as an unresolved
  candidate, not a finding, with four checks against it. Oase1 is admitted below
  the SNP floor and reported as unplaceable on this axis at 25,838 usable sites.


- New `archaic/transect.py`: a reusable pooled time-transect layer (date
  binning, pooled cohort archaic statistics, archaic-free f4-ratio ancestry
  fractions, mixture prediction with full error propagation, coverage-matched
  SNP sets, endpoint contrasts and inverse-variance trends), with 14 unit tests
  that need no AADR data.
- New `oceania_transect.py` study: a 3,000-year Remote Oceania archaic-ancestry
  time transect (`archaic-pipeline oceania-transect`). Pooled Denisovan
  affinity in Vanuatu rises from indistinguishable-from-zero in the founding
  Lapita horizon to 75% of the present-day Papuan level, tracking an
  independently measured Papuan-related ancestry influx; the Guam/Marianas
  control stays flat, giving a difference-in-differences of +0.030 +/- 0.012
  (Z = 2.58). Includes coverage-matched recomputation, a transversions-only
  damage sensitivity, and an explicit finding that the parallel Neanderthal
  rise appears in the control too and is therefore *not* attributed to the
  influx.
- Published a dedicated EV16A.SG Etruscan example with a portable HTML report,
  machine-readable sensitivity results, provenance, and explicit comparison to
  the segment-supported Oase1 evidence standard.

### Fixed

- Made `highest-archaic --subset` robust when a focused subset contains no
  high/elite-confidence sample, and prevented subset reports from overwriting
  the canonical all-sample report.

## [0.6.0] - 2026-07-17

### Added

- Added `archaic.denisovan_genome` / `archaic-denisovan`, a reference-aware
  Denisova 3 workflow with genotype QC, same-specimen replicate concordance,
  Denisovan/Neanderthal lineage tests, Denisova 11 F1 controls, transversion and
  chromosome sensitivity, reference-defined marker fingerprints, machine-readable
  provenance, and a self-contained HTML report. The module deliberately does not
  turn reference-dependent archaic statistics into a Denisovan percentage.
- Added `python -m archaic.highest_archaic`, a resumable all-AADR extreme scan
  with distribution-aware confidence sets, lower-bound ranking, explicit
  artifact scoring, top-candidate genotype sensitivity controls, twelve figures,
  per-candidate reports, deterministic seeds, dry-run/subset modes, and tests.
- Extended AADR annotation parsing with master ID, publication/DOI, and library
  type fields. Denisovan and combined percentages remain deliberately null
  because the validated `D_Den` statistic is a relative affinity, not a calibrated
  ancestry proportion.
- Replaced the monolithic project README with a concise onboarding page and a
  structured documentation set covering setup, methods, interpretation, studies,
  reports, and evidence labels.
- Added a data/artifact policy documenting the intentionally retained PRJEB10597
  BAM files, accessions, publisher MD5 values, repository SHA-256 values, and
  exploratory evidence boundary.
- Added contribution, security, conduct, issue, pull-request, and dependency
  maintenance files for the public GitHub repository.
- Added repository documentation/version checks and Windows smoke coverage to CI.

### Robust outlier model and portability hardening

#### Added
- **`archaic/neighborhood.py`**: shared local-neighbour residual model for Phase 6
  and Phase 9. The expected-value uncertainty now uses the precision-weighted
  effective neighbour count instead of the nominal K, so a few high-coverage
  neighbours cannot make the expected mean look more certain than it is.
- **Phase 9 robustness perturbations**: duplicate-library-root pruning, feature
  weight perturbations (ancestry-only, equal geo/time, weaker PCs), optional
  alternate PCA files (`phase5_pca.py --offset-frac ... --out ...`), and focused
  local READ-style kinship pruning around the nominal top candidates.
- **High-confidence technical gates** in `config.yaml` / `phase4_normalize.py`:
  SNP count remains the main floor, with available coverage, contamination lower
  bound, and damage metadata now used to exclude technically risky genomes from
  individual-level outlier calls.
- **Continent-audit metadata** in Phase 2: retained samples now record
  `continent_uncertain`, and unrecognized-country Eurasia fallbacks are surfaced
  in the sample flags.

#### Changed
- Removed developer-machine AADR paths from built-in defaults. A fresh checkout now
  requires `aadr_dir` in config, `ARCHAIC_CONFIG`, or `ARCHAIC_AADR_DIR`, and fails
  with a setup message if none is configured. `config.local.yaml` is auto-detected
  and ignored for machine-local paths.
- `phase5_pca.py` can write alternate PCA files from shifted even-SNP grids for
  PCA-subset sensitivity checks.

#### Fixed
- Removed lingering unused bindings flagged by `pyflakes` in `denisovan_survey.py`,
  `etruscan_study.py`, and `phase8_figures.py`.

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
