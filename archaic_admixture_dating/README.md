# Archaic admixture dating

This repository module provides restartable, uncertainty-aware tract-length
dating for the Papuan Denisovan V1 project. It compares a single pulse, two
ordered pulses, and an interpretable prolonged-flow approximation, then
calibrates those fits against configured M1-M10 demographic and measurement
scenarios.

The module does not call novel Denisovan tracts. Its primary V1 route imports
tracts from an established caller or a legally reusable validated tract table.
Raw genomes, controlled-access inputs, and large simulations remain outside
Git.

## Quick smoke test

From the repository root:

```bash
python -m archaic_admixture_dating.cli run-all \
  --profile smoke \
  --run-id smoke_validation \
  --resume
```

The command creates synthetic tract data, applies the same schema and QC used
for imported data, fits all three dating models, runs M1-M10 smoke simulations,
bootstraps chromosomes, and generates a portable HTML report.

## Import an established tract table

```bash
python -m archaic_admixture_dating.cli import-tracts \
  --profile laptop \
  --run-id papuan_authorized_v1 \
  --input /approved/path/tracts.tsv \
  --caller ibdmix \
  --population Papuan
```

Generic tables can map caller-specific columns repeatedly with
`--map-column source=standard`. Imported tables must include sample,
population, chromosome, base-pair start/end, and genetic length in cM (or
genetic start/end).

## Safe one-hour laptop workflow

```bash
python -m archaic_admixture_dating.cli run-all \
  --profile laptop \
  --run-id papuan_authorized_v1 \
  --input /approved/path/tracts.tsv \
  --caller ibdmix \
  --max-block-minutes 55 \
  --resume
```

Each long loop works in atomic units, checks the deadline before starting a new
unit, writes an atomic checkpoint only after output validation, and exits with
code 75 when it pauses safely. Run the identical command again to continue.
Valid partial downloads are retained as `.part` files and resumed with HTTP
range requests.

Inspect progress:

```bash
python -m archaic_admixture_dating.cli status \
  --profile laptop \
  --run-id papuan_authorized_v1
```

## Runtime profiles

- `smoke`: synthetic fixtures, one thread, minute-scale end-to-end validation.
- `laptop`: four threads or fewer, chromosome/replicate checkpoints, limited
  bootstrap/simulation counts, 55-minute blocks.
- `full`: larger bootstrap and simulation counts for a workstation or cloud;
  still restartable.

## Snakemake wrapper

Install the optional workflow extra and set `workflow.tract_input`,
`workflow.profile`, and `workflow.run_id` in a machine-local copy of the YAML:

```bash
python -m pip install -e ".[workflow]"
snakemake \
  --snakefile archaic_admixture_dating/workflows/Snakefile \
  --cores 4 \
  --rerun-incomplete
```

Snakemake tracks the input, configuration, log, report, provenance, and output
manifest. The invoked Python workflow owns finer-grained chromosome,
simulation-replicate, bootstrap, model-fit, and download checkpoints and
voluntarily pauses before the configured block deadline.

## Core commands

Run `python -m archaic_admixture_dating.cli --help` and each subcommand's
`--help` for the full interface. The CLI includes project initialization, data
inspection, storage estimation, guarded downloads, tract import, QC, all three
fits, simulations, model comparison, sensitivity, report generation, the
end-to-end workflow, and status.

## Interpretation boundary

A fitted recent date is not evidence by itself that Denisovans survived until
that date. Bottlenecks, recent mixing between modern-human populations with
different Denisovan ancestry, selection, map error, tract-detection error, and
multiple Denisovan-related sources can generate a long-tract excess. Reports
must retain those alternatives unless calibrated analyses reject them.

## Cleanup

Run outputs live below `archaic_admixture_dating/outputs/<run_id>` and are
ignored by Git. A completed run can be removed by deleting that one explicitly
named run directory after preserving its report, configuration snapshot,
provenance, and any non-reproducible authorized inputs. Never delete a shared
download cache or an unresolved `.part` file during cleanup.
