# Contributing

Contributions are welcome when they preserve the repository's scientific and
reproducibility boundaries.

## Before opening a change

1. Open an issue for a new estimator, dataset, or interpretation change.
2. Keep AADR and unapproved BAM/CRAM files outside Git.
3. Reuse the existing packed-panel reader, statistics, configuration, and
   provenance machinery where possible.
4. State whether a result is validated, supported, exploratory, or
   inconclusive/data-limited.

## Development setup

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[test,sim]"
```

Run the checks used in CI:

```bash
pyflakes archaic *.py tests tools oase1_bam_pipeline
python -m pytest -q
archaic-pipeline smoke-test
python tools/check_repo_docs.py
```

## Scientific changes

A new or changed statistic should include:

- a precise definition and reference populations;
- synthetic or simulation tests with known expected direction;
- a null or negative control and a positive control when available;
- uncertainty appropriate to linked genomic data;
- missingness, coverage, contamination, damage, and ascertainment limits;
- explicit behavior under transversions, alternate references/outgroups, or
  chromosome/block sensitivity when relevant.

Do not convert `D_Den` to a Denisovan percentage or combined archaic percentage
without an independently justified calibration and corresponding validation.

## Generated results

Stage files explicitly. Do not use `git add -A` in a working tree containing
local analyses. Retain a generated result only when its command, input release,
provenance, and interpretation status are documented in
[docs/DATA.md](docs/DATA.md) or the study report.

## Publishing a focused change

GitHub authentication and the `origin` remote are already configured for this
checkout. Existing tracked branches can be published with `git push`. For a new
branch, the repository's `push.autoSetupRemote` setting makes its first
`git push` create and track the matching `origin` branch automatically.

For a guarded commit-and-push workflow, use the PowerShell helper with an
explicit file list:

```powershell
.\tools\publish.ps1 -Message "Describe the focused change" -Path README.md, CONTRIBUTING.md
```

It refuses a pre-staged index, checks whether the branch is behind its remote,
rejects ignored and out-of-repository paths, shows the exact staged files,
requires `PUBLISH` confirmation, then verifies that the remote SHA matches and
reports `0 0` local/remote parity. It deliberately does not support an implicit
"all files" mode.

## Pull requests

Keep pull requests focused. Include:

- what changed and why;
- scientific and user impact;
- commands used for validation;
- any outputs deliberately excluded;
- whether documentation, citation metadata, or the changelog changed.
