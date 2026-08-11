# Papuan Denisovan admixture dating V1

## Question

Can Papuan Denisovan tract lengths distinguish a relatively recent pulse near
25-35 kya from older 40-55 kya admixture, two pulses, prolonged flow, or older
admixture plus demographic and measurement effects?

## V1 route

V1 imports established tract calls, validates a standard schema, filters with
explicit exclusion reasons, summarizes QC, fits three tract-length models, and
calibrates interpretation with M1-M10 simulations. Source classification
defaults to `unresolved`; reference similarity is not treated as unique
geographic identification.

## Required input

Use only data whose access, consent, and reuse terms permit the analysis.
Genetic lengths must come from an appropriate recombination map. Do not pass
raw base-pair lengths as the dating unit.

## Smoke workflow

```bash
python -m archaic_admixture_dating.cli run-all \
  --profile smoke \
  --run-id papuan_v1_smoke \
  --resume
```

## Laptop workflow

```bash
python -m archaic_admixture_dating.cli run-all \
  --profile laptop \
  --run-id papuan_v1_authorized \
  --input /approved/path/tracts.tsv \
  --caller ibdmix \
  --max-block-minutes 55 \
  --resume
```

The HTML report is written to
`archaic_admixture_dating/outputs/<run_id>/report/report.html`.

## Laptop failure-anatomy workflow

The follow-up workflow uses the same 89-person S4/S5 source files and the
cached chromosome-specific GRCh37 maps. It measures threshold, long-tail,
chromosome, individual, source-affinity, published-parameter, and candidate-
selection-locus sensitivity before running a bounded M1-M10 observation-
process calibration:

```bash
python -m archaic_admixture_dating.failure_anatomy \
  --s5 /approved/path/skov_2018_s5_segments.tsv \
  --s4 /approved/path/skov_2018_s4_parameters.xlsx \
  --genetic-map-dir /path/to/grch37_hapmap_phase2 \
  --selected-loci-bed archaic_admixture_dating/projects/papuan_denisovan_v1/references/gower_2021_melanesian_ai_candidates_grch37.bed \
  --output archaic_admixture_dating/outputs/papuan_s5_failure_anatomy

python -m archaic_admixture_dating.observation_calibration \
  --observed-overall archaic_admixture_dating/outputs/papuan_s5_failure_anatomy/overall.tsv \
  --config archaic_admixture_dating/configs/papuan_denisovan_v1.yaml \
  --replicates 20 \
  --tracts-per-replicate 5000 \
  --output archaic_admixture_dating/outputs/papuan_observation_calibration
```

The result remains `not estimable`; see `REAL_DATA_RESULTS.md` and the compact
committed summary tables.
