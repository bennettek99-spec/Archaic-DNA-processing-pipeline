"""Portable Markdown/HTML reporting with explicit scientific guardrails."""

from __future__ import annotations

import base64
import html
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .checkpointing import atomic_write_json, atomic_write_text, sha256_file, utc_now
from .diagnostics import interpretation_status


def _dating_estimability(
    model_table: pd.DataFrame,
    fits: dict[str, dict[str, Any]],
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    single_warnings = set(fits.get("single_pulse", {}).get("warning_flags", []))
    if "estimate_at_parameter_bound" in single_warnings:
        reasons.append("single-pulse estimate hit its configured parameter bound")
    if "poor_single_exponential_fit" in single_warnings:
        reasons.append("single-pulse exponential goodness-of-fit failed")
    for model_id in ("two_pulse", "continuous_flow"):
        fit = fits.get(model_id, {})
        if float(fit.get("younger_generations", float("inf"))) <= 50.01:
            reasons.append(f"{model_id} younger component hit its configured parameter bound")
    if not model_table.empty:
        best = model_table.iloc[0]
        if str(best.get("parameter_recovery_quality", "")).lower() == "poor":
            reasons.append("best-BIC model had poor simulation parameter recovery")
        accuracy = pd.to_numeric(
            pd.Series([best.get("simulation_classification_accuracy")]),
            errors="coerce",
        ).iloc[0]
        if pd.notna(accuracy) and float(accuracy) < 0.70:
            reasons.append("best-BIC model failed the 70% simulation classification floor")
    return not reasons, reasons


def _markdown_table(frame: pd.DataFrame, columns: list[str] | None = None, limit: int = 20) -> str:
    view = frame.loc[:, columns] if columns else frame
    view = view.head(limit)
    if view.empty:
        return "_No rows._"
    headers = [str(value) for value in view.columns]
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in view.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |")
    return "\n".join(lines)


def git_info(repository: str | Path) -> dict[str, str | None]:
    root = Path(repository)

    def run(*args: str) -> str | None:
        result = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "remote": run("remote", "get-url", "origin"),
    }


def write_report(
    output_dir: str | Path,
    *,
    config: dict[str, Any],
    tracts: pd.DataFrame,
    summaries: dict[str, pd.DataFrame],
    model_table: pd.DataFrame,
    fits: dict[str, dict[str, Any]],
    warnings: list[str],
    command: str,
    repository: str | Path,
    figure: str | Path | None = None,
    simulation_summary: pd.DataFrame | None = None,
    sensitivity_table: pd.DataFrame | None = None,
    calibration_table: pd.DataFrame | None = None,
    confusion_table: pd.DataFrame | None = None,
    uncertainty_table: pd.DataFrame | None = None,
    published_benchmark: dict[str, Any] | None = None,
) -> tuple[Path, Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    estimable, estimability_reasons = _dating_estimability(model_table, fits)
    status = interpretation_status(warnings) if estimable else "not-estimable"
    best = model_table.iloc[0]
    plausible = model_table.loc[model_table["delta_bic"] < 6, "model_name"].tolist()
    figure_markdown = f"![QC and model comparison]({Path(figure).name})" if figure else "_No figure generated._"
    warning_lines = "\n".join(f"- {warning}" for warning in warnings)
    estimability_lines = "\n".join(f"- {reason}" for reason in estimability_reasons)
    if estimable:
        executive_summary = f"""Under the tested tract-length models, **{best['model_name']}** has the lowest
BIC. Models remaining within ΔBIC < 6 are: {", ".join(plausible)}. This result is
conditional on the imported tract calls, masks, detection threshold,
recombination map, generation-time assumptions, and simulated alternatives."""
        conclusion = f"""Under the tested models, the observed tract distribution is more compatible
with {best['model_name']} than models with substantially larger BIC, but all
models within ΔBIC < 6 and the listed demographic/caller alternatives remain
plausible."""
    else:
        executive_summary = f"""**No event date is estimable from this S5 interval-length analysis.**
Although {best['model_name']} has the lowest BIC among the fitted models, the
fits reached lower parameter bounds and/or failed simulation-recovery checks.
The BIC ranking and fitted dates are diagnostic outputs only and must not be
interpreted as a late admixture event."""
        conclusion = """The S5 interval export is unsuitable for direct event dating with the current
independent-tract exponential models. The fitted dates are rejected rather
than reported as biological estimates. A caller-aware transition model or
independently calibrated tract calls with a genomic recombination map is
required before event-time inference."""
    benchmark_section = "_No matched published parameter benchmark was supplied._"
    if published_benchmark:
        summary = published_benchmark["summary"]
        bootstrap = published_benchmark["bootstrap_median"]
        benchmark_section = f"""The matched S4 workbook contains published HMM transition-parameter fits
for {summary['n_individuals']} Papuan individuals, the same sample set found in
S5. The per-individual admixture-time parameter has median
**{summary['median_generations']:.1f} generations** (range
{summary['minimum_generations']:.1f}-{summary['maximum_generations']:.1f}).
At 29 years per generation this is **{summary['median_kya_29']:.2f} kya**.
A deterministic 10,000-replicate bootstrap of the cross-individual median is
{bootstrap['ci_low_kya_29']:.2f}-{bootstrap['ci_high_kya_29']:.2f} kya; this
describes uncertainty in the median of published per-person estimates, not
the full uncertainty of an admixture event date."""
    report = f"""# Papuan Denisovan admixture dating V1 report

Generated: {utc_now()}
Interpretation status: **{status}**

## Executive summary

{executive_summary}

It does **not** establish that Denisovans survived until any fitted or
published date.

## Research question

Can the observed Denisovan tract distribution distinguish a 25-35 kya event
from older admixture around 40-55 kya, two pulses, prolonged flow, or an older
pulse plus demographic and measurement effects that mimic a younger signal?

## Data and permissions

This run contains {len(tracts):,} retained tracts from
{tracts['sample_id'].nunique():,} individuals. The workflow does not infer
ethical permission from technical availability; raw genomes and
redistribution-restricted data are excluded from Git.

## Methods

Lengths are analyzed in centimorgans/Morgans. The single-pulse model is a
left-truncated exponential. The two-pulse model is an ordered two-component
mixture with collapse and separation warnings. Prolonged flow is approximated
as a uniform interval mixture. Dates are conditional estimates, not direct
fossil or survival dates.

The fast single-pulse uncertainty uses
{config['dating']['bootstrap_replicates']} full-data bootstrap replicates. The
exploratory two-pulse and continuous-flow uncertainty uses
{config['dating']['complex_bootstrap_replicates']} replicates, with at most
{config['dating'].get('complex_bootstrap_max_tracts') or 'all'} tracts per
replicate; point estimates use all retained tracts.

## Matched published S4 benchmark

{benchmark_section}

## Quality-control results

{_markdown_table(summaries['overall'])}

{figure_markdown}

## Observed tract distributions

### By population

{_markdown_table(summaries['population'])}

### By chromosome

{_markdown_table(summaries['chromosome'], limit=30)}

### Sample missingness/callability

{_markdown_table(summaries.get('sample_missingness', pd.DataFrame()), limit=40)}

### Caller/batch diagnostics

{_markdown_table(summaries.get('caller_batch', pd.DataFrame()), limit=20)}

### Source-class summaries

{_markdown_table(summaries.get('source_class', pd.DataFrame()), limit=20)}

## Single-pulse fit

```json
{json.dumps(fits['single_pulse'], indent=2)}
```

## Two-pulse fit

```json
{json.dumps(fits['two_pulse'], indent=2)}
```

## Prolonged-flow fit

```json
{json.dumps(fits['continuous_flow'], indent=2)}
```

## Simulation calibration

{_markdown_table(simulation_summary if simulation_summary is not None else pd.DataFrame())}

Tract-level smoke simulations test plumbing and approximate identifiability.
They do not substitute for calibrated caller-aware coalescent simulations on
appropriately governed genomic data.

### Parameter recovery

{_markdown_table(calibration_table if calibration_table is not None else pd.DataFrame(), limit=40)}

### Model confusion

{_markdown_table(confusion_table if confusion_table is not None else pd.DataFrame(), limit=30)}

## Model comparison

{_markdown_table(model_table.drop(columns=['fit_json'], errors='ignore'))}

The best-fitting model is not automatically uniquely supported. Models marked
competitive or plausible remain live explanations under the tested
assumptions.

## Sensitivity analyses

{_markdown_table(sensitivity_table if sensitivity_table is not None else pd.DataFrame(), limit=80)}

The table evaluates generation time, minimum length, confidence,
recombination-map scale, leave-one-chromosome-out fits, longest-tract
exclusion, and selected-locus exclusion where annotations exist. Chromosome
and sample bootstrap are stored as machine-readable tables. M7 and configured
caller-error simulations cover alternate bottleneck and detection assumptions.
Unavailable rows remain explicit rather than being silently omitted.

### Bootstrap uncertainty

{_markdown_table(uncertainty_table if uncertainty_table is not None else pd.DataFrame(), limit=30)}

## Alternative explanations

M7-M10 cover an older pulse plus a severe bottleneck, modern-human population
mixing, two divergent Denisovan-related sources, and selection on a tract
subset. Phasing error, tract detection, and map error remain caller- and
dataset-dependent.

## Automatic warnings

{warning_lines}

### Estimability failures

{estimability_lines if estimability_lines else "_None._"}

## Limitations

- Tract dates depend on an appropriate genetic map and tract-call calibration.
- Short-tract censoring can bias pulse dates younger.
- Two pulses and prolonged flow may be non-identifiable.
- Altai similarity does not uniquely locate or name the donor population.
- Relative Denisovan affinity is not an absolute Denisovan percentage.
- Technical availability does not override Indigenous data governance.

## Conclusions

{conclusion}

Direct late Denisovan admixture is not demonstrated by a young tract-length
estimate alone.

## Reproducibility information

Exact command:

```text
{command}
```

Software: archaic-admixture-dating {__version__}; Python
{sys.version.split()[0]}; {platform.platform()}.

The run directory contains a configuration snapshot, model tables, tract
tables, checkpoints, input/output manifests, hashes, and Git metadata.
"""
    markdown_path = output / "report.md"
    atomic_write_text(markdown_path, report)

    image_html = ""
    if figure:
        figure_path = Path(figure)
        encoded = base64.b64encode(figure_path.read_bytes()).decode("ascii")
        image_html = f'<img alt="QC and model comparison" src="data:image/png;base64,{encoded}">'
    html_body = html.escape(report)
    html_body = f"<pre>{html_body}</pre>"
    html_report = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Papuan Denisovan V1 report</title>
<style>
body{{font:16px/1.55 system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;color:#202124}}
pre{{white-space:pre-wrap;font:inherit}} img{{max-width:100%;height:auto;border:1px solid #ddd}}
</style></head><body>{image_html}{html_body}</body></html>"""
    html_path = output / "report.html"
    atomic_write_text(html_path, html_report)

    provenance = {
        "generated_at": utc_now(),
        "command": command,
        "software_version": __version__,
        "python_version": sys.version,
        "git": git_info(repository),
        "outputs": {
            "report_md": {"path": str(markdown_path), "sha256": sha256_file(markdown_path)},
            "report_html": {"path": str(html_path), "sha256": sha256_file(html_path)},
        },
        "interpretation_status": status,
    }
    atomic_write_json(output / "provenance.json", provenance)
    return markdown_path, html_path
