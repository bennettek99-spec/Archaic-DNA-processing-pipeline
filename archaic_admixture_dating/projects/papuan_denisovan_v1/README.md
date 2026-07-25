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
