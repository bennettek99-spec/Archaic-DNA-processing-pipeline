# Papuan Denisovan V1 implementation plan

Status: implemented and smoke-validated; authorized real-data input pending
Created: 2026-07-25
Primary profile: laptop (Windows-compatible, restartable blocks under 55 minutes)

## Scientific objective and interpretation boundary

V1 will test whether an observed Papuan Denisovan tract-length distribution is
more compatible with one older pulse, one relatively recent pulse, two pulses,
prolonged flow, or an older pulse plus demographic and measurement processes
that can mimic recent gene flow. Dates will always be conditional on tract
calling, recombination map, generation time, truncation model, and tested
demographies. V1 will not translate a relative Denisovan-affinity statistic into
an ancestry percentage and will not claim Denisovan survival to a particular
date unless late direct admixture is distinguishable from the alternatives.

## Repository components to reuse

- `archaic.config`: YAML configuration conventions and explicit local-path
  handling.
- `archaic.log_utils`: timestamped progress logging for long stages.
- `archaic.stats`: array-oriented estimates and linked-data uncertainty
  conventions.
- `archaic.panel`, `archaic.lib_eigenstrat`, and `archaic.refs`: optional
  AADR/EIGENSTRAT population and reference integration without duplicating
  genotype readers.
- `archaic.simulate`: existing `msprime` dependency and deterministic
  simulation conventions.
- `archaic.synthetic`: compact, redistribution-safe synthetic fixtures.
- Existing repository rules: no large genomic inputs in Git, explicit
  interpretation status, focused changes, reproducible commands, and
  scientific fail-closed behavior.

The new package will not change existing estimators or generated result files.

## Files to create

```text
archaic_admixture_dating/
├── README.md
├── __init__.py
├── __main__.py
├── cli.py
├── config.py
├── logging_utils.py
├── checkpointing.py
├── downloads.py
├── manifests.py
├── tract_schema.py
├── tract_import.py
├── tract_filtering.py
├── tract_summary.py
├── dating_single_pulse.py
├── dating_two_pulse.py
├── dating_continuous.py
├── simulations.py
├── model_comparison.py
├── bootstrap.py
├── diagnostics.py
├── plotting.py
├── reporting.py
├── configs/papuan_denisovan_v1.yaml
├── projects/papuan_denisovan_v1/
│   ├── IMPLEMENTATION_PLAN.md
│   ├── README.md
│   ├── DATA_SOURCES.md
│   ├── DATA_GOVERNANCE.md
│   ├── METHODS.md
│   ├── LIMITATIONS.md
│   └── run_v1.py
├── tests/
│   ├── fixtures/
│   ├── test_checkpointing.py
│   ├── test_dating.py
│   ├── test_downloads.py
│   ├── test_simulations.py
│   ├── test_tract_schema.py
│   └── test_workflow_smoke.py
└── outputs/.gitkeep
```

Root packaging metadata will be changed only as needed to include the new
package and console entry point.

## Dependencies

Required dependencies already present: Python 3.10+, NumPy, pandas, SciPy,
Matplotlib, and PyYAML.

Optional dependencies:

- `msprime>=1.3` for coalescent pilot/production simulations (already exposed
  by the repository's `sim` extra).
- `pyarrow>=14` for Parquet input/output. TSV remains fully supported when
  `pyarrow` is unavailable.

Snakemake is deferred until the Python checkpoint orchestrator is stable; the
V1 CLI itself provides dependency checks, atomic checkpoints, and safe
stop/resume behavior without making smoke mode depend on a separate workflow
runtime.

## Likely data sources and governance

Candidate modern Papuan sources include SGDP public/authorized releases,
HGDP/1000 Genomes comparison panels where suitable, and published tract tables
whose reuse terms permit this analysis. African references will preferentially
use Yoruba or Mbuti and East Asian comparisons Han or Dai. The high-coverage
Altai Denisovan or an approved derivative panel is the archaic reference.
Raw Papuan genomes, controlled-access data, and redistribution-restricted
tracts will never be bundled or downloaded without explicit authorization.

The primary V1 route is `import-tracts` from an established caller or published
validated tract table. IBDmix is the leading primary caller adapter because it
can work without phasing and has established archaic-reference workflows.
admixfrog is the secondary candidate for sensitivity where input and compute
allow. V1 will not implement a novel tract caller.

## Storage estimates

- Source metadata and manifests: under 10 MB.
- Imported compressed tract tables: typically 10 MB to 2 GB.
- Recombination maps and masks: roughly 0.5-2 GB combined.
- Public phased WGS panels: tens to hundreds of GB and therefore manual,
  component-wise, or controlled-access only in laptop mode.
- Pilot simulations: under 2 GB, with summary-only retention by default.
- Full simulations: potentially tens of GB; disabled by the laptop defaults.
- Reports, tables, checkpoints, and plots: typically under 500 MB per run.

The storage estimator will require expected sizes before large transfers,
reserve working headroom, refuse insufficient disk space, and never silently
start a multi-gigabyte download.

## Computational bottlenecks

- Whole-genome tract inference and phasing.
- Genetic-map interpolation for millions of interval endpoints.
- Bootstrap refits of mixture models.
- Coalescent simulation across many chromosomes and competing demographies.
- Applying realistic caller error and ascertainment to simulated tracts.

Mitigation: import validated tracts, stream/chunk tables, process by chromosome,
use deterministic per-unit checkpoints, start with chromosomes 21 and 22,
retain summary statistics instead of tree sequences, and stop before the
configured deadline.

## Dating and model-comparison strategy

- Single pulse: left-truncated exponential tract-length likelihood in Morgans,
  with generations since pulse as the decay parameter.
- Two pulses: ordered two-component mixture of left-truncated exponentials,
  fitted by bounded likelihood optimization with collapse/separation warnings.
- Prolonged flow: an interpretable uniform-time interval approximation whose
  tract density is a mixture of exponentials integrated numerically.
- Uncertainty: chromosome-block and sample bootstrap, deterministic seeds,
  generation-time sensitivity, and likelihood-profile/fit diagnostics.
- Model comparison: log likelihood, AIC, BIC, held-out score where feasible,
  bootstrap stability, parameter recovery, posterior predictive summaries, and
  warnings. Complexity alone will never determine the conclusion.

## Simulation strategy

All ten prompt models (M1-M10) will be configuration-defined. Smoke mode uses a
fast tract-level generative approximation so the whole workflow runs in
minutes. With `msprime` installed, pilot and full profiles add coalescent
replicates on configured chromosome subsets. Seeds are derived from the master
seed, model ID, and replicate ID. Each replicate is an atomic checkpoint unit.
Observed and simulated summaries use the same length threshold, masks, and
configurable detection-error approximation.

Demographic alternatives explicitly include a severe Papuan bottleneck, recent
modern-human mixing between groups with different Denisovan ancestry, selected
long tracts, phasing/detection error, recombination-map error, and two divergent
Denisovan sources. These models are sensitivity scenarios, not literal claims
about unobserved populations.

## Exact laptop-safe checkpoint strategy

1. Every run receives a normalized configuration hash and immutable config
   snapshot.
2. State is stored in `outputs/<run_id>/checkpoints/<task>.json`.
3. A task records software version, config hash, input fingerprints, master
   seed, completed units, failed units, timestamps, and output paths.
4. Units are chromosomes, download components, tract batches, simulation
   replicates, bootstrap replicates, model fits, or plot groups.
5. Checkpoints and result tables use temporary sibling files plus atomic
   replacement.
6. A checkpoint is updated only after output validation.
7. `--resume` skips a unit only when its checkpoint metadata and output
   fingerprints remain valid.
8. Each long loop checks a deadline before starting another unit. At the
   default 55-minute block with a five-minute stop buffer, it exits with a
   distinct resumable status after the current atomic unit.
9. Partial HTTP downloads remain as `.part` files with byte counts and ETag or
   Last-Modified metadata where available; validation precedes final rename.
10. Random seeds are deterministic per unit, so interrupted work is exactly
    reproducible.

## Unresolved scientific risks

- Available public data may not include enough legally reusable, phased Papuan
  genomes for de novo caller validation.
- Short-tract detection thresholds can erase the signal from older pulses and
  bias dates younger.
- Two-pulse and prolonged-flow models may be non-identifiable at realistic
  sample sizes.
- Bottlenecks, recent modern-human population mixing, selection, caller error,
  and map error can preserve or create apparent long-tract excess.
- A single Altai Denisovan reference does not uniquely identify the true donor
  population or geography.
- Tracts from divergent Denisovan-related sources may be merged by a caller.
- Indigenous data governance can prohibit an analysis even when files are
  technically reachable.
- Published tract tables may encode caller-specific filters that cannot be
  reconstructed fully.

Fail-closed behavior: when sample size, callability, recovery, identifiability,
or model separation is inadequate, the report will say `not distinguishable`
or `not estimable`.

## Phased implementation checklist

### Phase 1 - scaffolding

- [x] Inspect repository and reusable components.
- [x] Write this implementation plan.
- [x] Create package/module layout and conservative configuration.
- [x] Add CLI skeleton and status/dry-run behavior.
- [x] Add atomic logging/checkpoint helpers.
- [x] Add compact synthetic tract fixtures.

### Phase 2 - safe data handling

- [x] Implement manifests and storage estimates.
- [x] Implement resumable range downloads and checksum verification.
- [x] Add disk-space and access/governance guards.
- [x] Test interruption, partial preservation, and atomic manifests.

### Phase 3 - tract ingestion and QC

- [x] Implement standard TSV/Parquet tract schema.
- [x] Add generic/IBDmix import adapters and provenance preservation.
- [x] Add filtering, masks, overlap diagnostics, and summaries.
- [x] Add QC tables and plots.

### Phase 4 - dating

- [x] Implement truncated single-pulse fit.
- [x] Implement ordered two-pulse mixture and warnings.
- [x] Implement prolonged-flow approximation.
- [x] Add chromosome/sample bootstrap and sensitivity surfaces.

### Phase 5 - simulations

- [x] Define M1-M10 in configuration.
- [x] Implement deterministic tract-level smoke simulations.
- [x] Implement `msprime` migration-tract pilot/full interface.
- [x] Add recovery and confusion diagnostics.

### Phase 6 - comparison and reporting

- [x] Produce model-comparison and sensitivity tables.
- [x] Generate cautious portable HTML and Markdown reports.
- [x] Record exact commands, environment, Git state, manifests, and hashes.
- [x] Add Snakemake wrapper with per-rule resources and resumable CLI units.
- [x] Complete smoke workflow and test suite.

## Validation status

- Synthetic smoke workflow: complete end to end.
- Focused plus repository regression tests: 76 passed.
- Static checks: `pyflakes archaic_admixture_dating` clean.
- Wheel build: successful in isolated PEP 517 mode.
- Real Papuan inference: not run because no authorized tract/genome input was
  supplied. Synthetic outputs are workflow/calibration evidence only.

## Deferred beyond V1

- GPU or large-cloud processing.
- A new tract-calling algorithm.
- Definitive source-population geography.
- ARG-scale whole-genome analysis.
- Automatic controlled-access acquisition.
- Claims about the terminal survival date of Denisovans.
