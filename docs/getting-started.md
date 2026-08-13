# Getting started

## Requirements

- Python 3.10, 3.11, or 3.12
- Windows PowerShell, Linux, macOS, or WSL
- An AADR v66.p1 panel for real-data workflows
- Enough memory and disk for the selected panel and retained outputs

The synthetic smoke test does not need AADR data.

## Install

```bash
git clone https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline.git
cd Archaic-DNA-processing-pipeline
python -m venv .venv
```

Activate the environment:

```powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
```

```bash
# Linux, macOS, or WSL
source .venv/bin/activate
```

Then install:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
```

The command-line entry points are designed for an editable checkout because
the phase scripts live at the repository root.

## Run the synthetic smoke test

```bash
archaic-pipeline smoke-test
```

This creates a temporary AADR-shaped packed panel, exercises the real
`PackedGeno`, `Panel`, f4-ratio, and D-statistic code, then removes the temporary
data. It checks plumbing and estimator direction; it is not a replacement for
the real-data or coalescent validation.

## Configure AADR

Obtain AADR from the Reich Lab / Harvard Dataverse under its own terms. Do not
copy the panel into this repository.

Copy `config.yaml` to `config.local.yaml` and edit:

```yaml
aadr_dir: "/path/to/aadr"
```

The configured directory should contain files such as:

```text
v66.p1_1240K.geno
v66.p1_1240K.snp
v66.p1_1240K.ind
v66.p1_1240K.anno
```

`config.local.yaml` is ignored by Git. Environment alternatives are:

```bash
export ARCHAIC_AADR_DIR=/path/to/aadr
# or
export ARCHAIC_CONFIG=/path/to/config.yaml
```

In PowerShell:

```powershell
$env:ARCHAIC_AADR_DIR = "C:\path\to\aadr"
```

## Validate before running the full workflow

```bash
archaic-pipeline validate --panel 1240k
```

The validation step checks known African, Eurasian, Oceanian, and archaic
controls. Review the output rather than treating command completion alone as a
scientific pass.

## Run Phases 2-9

```bash
archaic-pipeline all --panel 1240k
```

Equivalent direct invocation:

```bash
python scripts/run_pipeline.py --panel 1240k
```

Progress messages include timestamps and can be made more verbose with
`ARCHAIC_LOG_LEVEL=DEBUG`.

The main HTML report is written to `reports/archaic_report_1240k.html`.
Intermediate phase outputs are written under `results/`.

## Development checks

```bash
python -m pip install -e ".[test,sim]"
python -m pytest -q
pyflakes archaic *.py tests tools oase1_bam_pipeline
python tools/check_repo_docs.py
```

## Common problems

- **No AADR path configured:** create `config.local.yaml` or set
  `ARCHAIC_AADR_DIR`.
- **AADR files not found:** confirm the panel prefix and that `.geno`, `.snp`,
  `.ind`, and `.anno` share the same prefix.
- **A low-coverage candidate looks extreme:** check SNP count, standard error,
  contamination/damage flags, transversion sensitivity, chromosome deletion,
  and block stability before interpretation.
- **A result is absent from a report:** reports are assembled from existing
  result artifacts; confirm the producing phase completed successfully.
