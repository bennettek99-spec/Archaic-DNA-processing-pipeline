# Environment setup

## Tested environment

- Windows 11 x86-64
- Python 3.12.10
- Exact tested packages in `requirements-papuan-v1.lock`

Create an isolated environment from the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements-papuan-v1.lock
.\.venv\Scripts\python.exe -m pip install -e .
```

The normal project extras remain available:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test,sim]"
```

`pyarrow` is optional and required only for Parquet. TSV remains the portable
default.

## Validation

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pyflakes archaic_admixture_dating
.\.venv\Scripts\python.exe -m build --wheel
.\.venv\Scripts\python.exe -m archaic_admixture_dating.cli run-all `
  --profile smoke `
  --run-id papuan_v1_smoke `
  --resume
```

Every run stores the resolved configuration, exact command, Git identity,
Python/platform information, deterministic seeds, input manifest, output
manifest, checkpoints, and hashes.
