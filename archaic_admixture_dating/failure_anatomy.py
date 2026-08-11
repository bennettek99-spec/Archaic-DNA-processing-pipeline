"""Laptop-safe diagnostics for Papuan S4/S5 tract-dating failure modes."""

from __future__ import annotations

import argparse
import base64
import html
import json
from io import BytesIO
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .checkpointing import atomic_write_json, atomic_write_text, sha256_file, utc_now
from .dating_single_pulse import fit_single_pulse
from .tract_filtering import filter_tracts, load_bed, overlap_mask
from .tract_import import import_tracts


AUTOSOMES = [str(value) for value in range(1, 23)]
GENERATION_TIME_YEARS = 29.0
MINIMUM_LENGTH_CM = 0.02


def _portable_fingerprint(path: str | Path) -> dict[str, object]:
    source = Path(path)
    return {
        "filename": source.name,
        "bytes": source.stat().st_size,
        "sha256": sha256_file(source),
    }


def _effective_decay(lengths: Iterable[float], threshold_cm: float) -> dict[str, object]:
    values = np.asarray(list(lengths), dtype=float)
    values = values[np.isfinite(values) & (values >= threshold_cm)]
    excess = (values - threshold_cm) / 100.0
    positive_sum = float(excess.sum())
    if len(values) < 10 or positive_sum <= 0:
        return {
            "n_tracts": len(values),
            "effective_generations_unbounded": np.nan,
            "ks_statistic": np.nan,
            "ks_pvalue": np.nan,
        }
    generations = len(values) / positive_sum
    ks = stats.kstest(excess, "expon", args=(0.0, 1.0 / generations))
    return {
        "n_tracts": len(values),
        "effective_generations_unbounded": float(generations),
        "ks_statistic": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
    }


def _length_summary(frame: pd.DataFrame) -> dict[str, object]:
    values = frame["length_cm"].to_numpy(float)
    if len(values) == 0:
        return {"tract_count": 0}
    result = {
        "tract_count": int(len(values)),
        "individuals": int(frame["sample_id"].nunique()),
        "total_length_cm": float(values.sum()),
        "mean_length_cm": float(values.mean()),
        "median_length_cm": float(np.median(values)),
        "q75_length_cm": float(np.quantile(values, 0.75)),
        "q90_length_cm": float(np.quantile(values, 0.90)),
        "q95_length_cm": float(np.quantile(values, 0.95)),
        "q99_length_cm": float(np.quantile(values, 0.99)),
        "fraction_ge_0_5cm": float(np.mean(values >= 0.5)),
        "fraction_ge_1cm": float(np.mean(values >= 1.0)),
        "fraction_ge_2cm": float(np.mean(values >= 2.0)),
    }
    decay = _effective_decay(values, MINIMUM_LENGTH_CM)
    result.update(decay)
    fit = fit_single_pulse(
        values,
        minimum_length_cm=MINIMUM_LENGTH_CM,
        generation_time_years=GENERATION_TIME_YEARS,
    )
    result.update(
        {
            "bounded_single_generations": fit["generations"],
            "bounded_single_kya": fit["kya"],
            "fit_warning_flags": ";".join(fit["warning_flags"]),
        }
    )
    return result


def _analysis_sets(tracts: pd.DataFrame) -> dict[str, pd.DataFrame]:
    source = tracts["source_class"].astype(str)
    return {
        "all_archaic": tracts,
        "denisovan_broad": tracts.loc[
            source.isin(["denisovan_affinity", "denisovan_affinity_strict"])
        ].copy(),
        "denisovan_strict": tracts.loc[source.eq("denisovan_affinity_strict")].copy(),
        "neanderthal_affinity": tracts.loc[source.eq("neanderthal_affinity")].copy(),
        "unresolved_affinity": tracts.loc[source.eq("unresolved")].copy(),
    }


def _group_summary(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    grouper = columns[0] if len(columns) == 1 else columns
    for key, group in frame.groupby(grouper, dropna=False, sort=True):
        values = key if isinstance(key, tuple) else (key,)
        rows.append({**dict(zip(columns, values)), **_length_summary(group)})
    return pd.DataFrame(rows)


def _survival_tables(sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for label, frame in sets.items():
        values = frame["length_cm"].to_numpy(float)
        if len(values) == 0:
            continue
        upper = max(MINIMUM_LENGTH_CM * 1.01, float(np.quantile(values, 0.999)))
        grid = np.geomspace(MINIMUM_LENGTH_CM, upper, 80)
        for threshold in grid:
            rows.append(
                {
                    "analysis_set": label,
                    "length_cm": threshold,
                    "survival_probability": float(np.mean(values >= threshold)),
                }
            )
    return pd.DataFrame(rows)


def _threshold_sensitivity(sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    thresholds = [0.02, 0.03, 0.05, 0.10, 0.20, 0.50, 1.0]
    for label, frame in sets.items():
        for threshold in thresholds:
            rows.append(
                {
                    "analysis_set": label,
                    "minimum_length_cm": threshold,
                    **_effective_decay(frame["length_cm"], threshold),
                }
            )
    return pd.DataFrame(rows)


def _trim_sensitivity(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    ordered = np.sort(frame["length_cm"].to_numpy(float))
    for fraction in [0.0, 0.001, 0.005, 0.01, 0.05]:
        keep = len(ordered) if fraction == 0 else max(10, int(len(ordered) * (1 - fraction)))
        values = ordered[:keep]
        cutoff = float(values[-1])
        rows.append(
            {
                "longest_fraction_removed": fraction,
                "upper_length_cutoff_cm": cutoff,
                **_effective_decay(values, MINIMUM_LENGTH_CM),
            }
        )
    return pd.DataFrame(rows)


def _leave_one_group_out(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for value in sorted(frame[column].dropna().astype(str).unique()):
        retained = frame.loc[frame[column].astype(str).ne(value)]
        rows.append(
            {
                "group_column": column,
                "excluded_group": value,
                **_effective_decay(retained["length_cm"], MINIMUM_LENGTH_CM),
                "median_length_cm": float(retained["length_cm"].median()),
            }
        )
    return pd.DataFrame(rows)


def _tail_concentration(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    cutoff = float(frame["length_cm"].quantile(0.99))
    tail = frame.loc[frame["length_cm"].ge(cutoff)]
    result = (
        tail.groupby(column, dropna=False)
        .agg(tail_tracts=("length_cm", "size"), tail_length_cm=("length_cm", "sum"))
        .reset_index()
    )
    result["tail_tract_fraction"] = result["tail_tracts"] / len(tail)
    result["tail_length_fraction"] = result["tail_length_cm"] / tail["length_cm"].sum()
    result["tail_cutoff_cm"] = cutoff
    return result.sort_values("tail_length_fraction", ascending=False).reset_index(drop=True)


def _class_comparisons(sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    labels = ["denisovan_broad", "denisovan_strict", "neanderthal_affinity"]
    rows: list[dict[str, object]] = []
    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            a = sets[left]["length_cm"].to_numpy(float)
            b = sets[right]["length_cm"].to_numpy(float)
            test = stats.ks_2samp(a, b, alternative="two-sided", method="asymp")
            rows.append(
                {
                    "left": left,
                    "right": right,
                    "left_n": len(a),
                    "right_n": len(b),
                    "left_median_cm": float(np.median(a)),
                    "right_median_cm": float(np.median(b)),
                    "median_difference_cm": float(np.median(a) - np.median(b)),
                    "ks_statistic": float(test.statistic),
                    "ks_pvalue": float(test.pvalue),
                }
            )
    return pd.DataFrame(rows)


def _sample_features(sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    for label in ["denisovan_broad", "denisovan_strict", "neanderthal_affinity"]:
        summary = _group_summary(sets[label], ["sample_id"])
        keep = [
            "sample_id",
            "tract_count",
            "total_length_cm",
            "mean_length_cm",
            "median_length_cm",
            "q90_length_cm",
            "q95_length_cm",
            "fraction_ge_1cm",
            "effective_generations_unbounded",
        ]
        summary = summary[keep].rename(
            columns={column: f"{label}_{column}" for column in keep if column != "sample_id"}
        )
        pieces.append(summary)
    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece, on="sample_id", how="outer", validate="one_to_one")
    return merged


def _bootstrap_spearman(x: np.ndarray, y: np.ndarray, *, seed: int, replicates: int = 1000):
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(replicates):
        index = rng.integers(0, len(x), len(x))
        if np.unique(x[index]).size < 2 or np.unique(y[index]).size < 2:
            continue
        values.append(float(stats.spearmanr(x[index], y[index]).statistic))
    if not values:
        return np.nan, np.nan
    return tuple(np.quantile(values, [0.025, 0.975]))


def build_s4_bridge(
    s4_path: str | Path,
    sample_features: pd.DataFrame,
    *,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    s4 = pd.read_excel(s4_path, sheet_name="Human population parameters")
    published = s4.loc[s4["Outgroup"].astype(str).eq("Whole~world")].copy()
    published = published.rename(
        columns={
            "name": "sample_id",
            "Dataset": "dataset",
            "Admixture_time": "published_admixture_generations",
            "Admixture_proportion": "published_admixture_proportion",
            "Coalescence_with_Human": "published_coalescence_human_generations",
            "Coalescence_with_Archaic": "published_coalescence_archaic_generations",
        }
    )
    columns = [
        "sample_id",
        "dataset",
        "published_admixture_generations",
        "published_admixture_proportion",
        "published_coalescence_human_generations",
        "published_coalescence_archaic_generations",
    ]
    bridge = published[columns].merge(
        sample_features,
        on="sample_id",
        how="inner",
        validate="one_to_one",
    )
    if len(bridge) != 89:
        raise ValueError(f"Expected 89 matched S4/S5 individuals, found {len(bridge)}")

    features = [
        column
        for column in bridge.columns
        if column.startswith(("denisovan_broad_", "denisovan_strict_", "neanderthal_affinity_"))
    ]
    rows: list[dict[str, object]] = []
    target = bridge["published_admixture_generations"].to_numpy(float)
    for index, feature in enumerate(features):
        values = bridge[feature].to_numpy(float)
        finite = np.isfinite(values) & np.isfinite(target)
        if finite.sum() < 10 or np.unique(values[finite]).size < 2:
            continue
        pearson = stats.pearsonr(target[finite], values[finite])
        spearman = stats.spearmanr(target[finite], values[finite])
        low, high = _bootstrap_spearman(
            target[finite], values[finite], seed=seed + index, replicates=1000
        )
        rows.append(
            {
                "feature": feature,
                "n": int(finite.sum()),
                "pearson_r": float(pearson.statistic),
                "pearson_pvalue": float(pearson.pvalue),
                "spearman_rho": float(spearman.statistic),
                "spearman_pvalue": float(spearman.pvalue),
                "spearman_bootstrap_low": float(low),
                "spearman_bootstrap_high": float(high),
            }
        )
    correlations = pd.DataFrame(rows).sort_values(
        "spearman_rho", key=lambda value: value.abs(), ascending=False
    )

    cohort_rows: list[dict[str, object]] = []
    feature = "denisovan_broad_median_length_cm"
    for dataset, group in bridge.groupby("dataset"):
        result = stats.spearmanr(
            group["published_admixture_generations"], group[feature]
        )
        cohort_rows.append(
            {
                "dataset": dataset,
                "n": len(group),
                "feature": feature,
                "spearman_rho": float(result.statistic),
                "spearman_pvalue": float(result.pvalue),
            }
        )
    return bridge, correlations, pd.DataFrame(cohort_rows)


def _map_constant_summary(tracts: pd.DataFrame) -> pd.DataFrame:
    work = tracts.copy()
    work["constant_length_cm"] = work["length_bp"] * 1.2e-8 * 100.0
    work["map_constant_ratio"] = work["length_cm"] / work["constant_length_cm"]
    rows: list[dict[str, object]] = []
    for columns in [["chromosome"], ["dataset"], ["source_class"]]:
        grouper = columns[0]
        for value, group in work.groupby(grouper, dropna=False):
            ratio = group["map_constant_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
            rows.append(
                {
                    "group_type": grouper,
                    "group": value,
                    "tract_count": len(group),
                    "median_map_constant_ratio": float(ratio.median()),
                    "q10_map_constant_ratio": float(ratio.quantile(0.10)),
                    "q90_map_constant_ratio": float(ratio.quantile(0.90)),
                    "median_map_length_cm": float(group["length_cm"].median()),
                    "median_constant_length_cm": float(group["constant_length_cm"].median()),
                }
            )
    return pd.DataFrame(rows)


def _selected_locus_sensitivity(frame: pd.DataFrame, bed_path: str | Path | None) -> pd.DataFrame:
    if bed_path is None:
        return pd.DataFrame(
            [
                {
                    "status": "not_run_no_validated_grch37_selection_bed",
                    "tracts_overlapping_selected_loci": np.nan,
                    "tracts_retained": len(frame),
                }
            ]
        )
    intervals = load_bed(bed_path)
    indexed = frame.reset_index(drop=True)
    selected = overlap_mask(indexed, intervals)
    retained = indexed.loc[~selected]
    return pd.DataFrame(
        [
            {
                "status": "complete",
                "tracts_overlapping_selected_loci": int(selected.sum()),
                "tracts_retained": len(retained),
                **_effective_decay(retained["length_cm"], MINIMUM_LENGTH_CM),
            }
        ]
    )


def _figure(
    survival: pd.DataFrame,
    map_constant: pd.DataFrame,
    bridge: pd.DataFrame,
    threshold: pd.DataFrame,
) -> str:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for label, group in survival.groupby("analysis_set"):
        if label == "all_archaic":
            continue
        axes[0, 0].plot(group["length_cm"], group["survival_probability"], label=label)
    axes[0, 0].set(xscale="log", yscale="log", xlabel="Length (cM)", ylabel="Survival")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_title("Empirical tract-length survival")

    chrom = map_constant.loc[map_constant["group_type"].eq("chromosome")].copy()
    chrom["order"] = pd.to_numeric(chrom["group"], errors="coerce")
    chrom = chrom.sort_values("order")
    axes[0, 1].plot(chrom["order"], chrom["median_map_constant_ratio"], marker="o")
    axes[0, 1].axhline(1.0, color="black", linewidth=1, linestyle="--")
    axes[0, 1].set(xlabel="Chromosome", ylabel="Median map/constant length ratio")
    axes[0, 1].set_title("Effect of local recombination map")

    for dataset, group in bridge.groupby("dataset"):
        axes[1, 0].scatter(
            group["published_admixture_generations"],
            group["denisovan_broad_median_length_cm"],
            label=dataset,
            alpha=0.8,
        )
    axes[1, 0].set(
        xlabel="Published S4 HMM admixture parameter (generations)",
        ylabel="S5 broad median tract length (cM)",
    )
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_title("Per-person S4-S5 bridge")

    broad = threshold.loc[threshold["analysis_set"].eq("denisovan_broad")]
    axes[1, 1].plot(
        broad["minimum_length_cm"],
        broad["effective_generations_unbounded"],
        marker="o",
    )
    axes[1, 1].set(xscale="log", xlabel="Minimum tract length (cM)", ylabel="Effective decay rate")
    axes[1, 1].set_title("Single-exponential threshold instability")
    fig.tight_layout()
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def _table_html(frame: pd.DataFrame, rows: int = 15) -> str:
    return frame.head(rows).to_html(index=False, border=0, classes="data", float_format=lambda x: f"{x:.5g}")


def _markdown_table(frame: pd.DataFrame, rows: int = 20) -> str:
    shown = frame.head(rows).copy()
    columns = [str(column) for column in shown.columns]

    def render(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in shown.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(render(value) for value in row) + " |")
    return "\n".join(lines)


def run_failure_anatomy(
    *,
    s5_path: str | Path,
    s4_path: str | Path,
    genetic_map_directory: str | Path,
    output_directory: str | Path,
    selected_loci_bed: str | Path | None = None,
    seed: int = 20260807,
) -> dict[str, object]:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)

    imported, malformed = import_tracts(
        s5_path,
        caller="skov_hmm",
        chromosomes=AUTOSOMES,
        genetic_map_directory=genetic_map_directory,
    )
    retained, excluded = filter_tracts(
        imported,
        minimum_length_cm=MINIMUM_LENGTH_CM,
        minimum_confidence=0.8,
        minimum_callable_fraction=0.8,
    )

    s4 = pd.read_excel(s4_path, sheet_name="Human population parameters")
    cohort = (
        s4.loc[s4["Outgroup"].astype(str).eq("Whole~world"), ["name", "Dataset"]]
        .drop_duplicates("name")
        .rename(columns={"name": "sample_id", "Dataset": "dataset"})
    )
    retained = retained.merge(cohort, on="sample_id", how="left", validate="many_to_one")
    if retained["dataset"].isna().any():
        missing = retained.loc[retained["dataset"].isna(), "sample_id"].unique()
        raise ValueError(f"Missing S4 cohort labels for {len(missing)} S5 individuals")

    sets = _analysis_sets(retained)
    overall = pd.DataFrame(
        [{"analysis_set": label, **_length_summary(frame)} for label, frame in sets.items()]
    )
    by_chromosome = pd.concat(
        [
            _group_summary(frame, ["chromosome"]).assign(analysis_set=label)
            for label, frame in sets.items()
        ],
        ignore_index=True,
    )
    by_dataset = pd.concat(
        [
            _group_summary(frame, ["dataset"]).assign(analysis_set=label)
            for label, frame in sets.items()
        ],
        ignore_index=True,
    )
    survival = _survival_tables(sets)
    threshold = _threshold_sensitivity(sets)
    trim = _trim_sensitivity(sets["denisovan_broad"])
    loo_chromosome = _leave_one_group_out(sets["denisovan_broad"], "chromosome")
    loo_individual = _leave_one_group_out(sets["denisovan_broad"], "sample_id")
    tail_chromosome = _tail_concentration(sets["denisovan_broad"], "chromosome")
    tail_individual = _tail_concentration(sets["denisovan_broad"], "sample_id")
    comparisons = _class_comparisons(sets)
    sample_features = _sample_features(sets)
    bridge, correlations, cohort_correlations = build_s4_bridge(
        s4_path, sample_features, seed=seed
    )
    map_constant = _map_constant_summary(retained)
    selected = _selected_locus_sensitivity(sets["denisovan_broad"], selected_loci_bed)

    tables = {
        "overall": overall,
        "by_chromosome": by_chromosome,
        "by_dataset": by_dataset,
        "survival": survival,
        "threshold_sensitivity": threshold,
        "longest_tract_trim": trim,
        "leave_one_chromosome_out": loo_chromosome,
        "leave_one_individual_out": loo_individual,
        "tail_by_chromosome": tail_chromosome,
        "tail_by_individual": tail_individual,
        "source_class_comparisons": comparisons,
        "sample_features": sample_features,
        "s4_s5_bridge": bridge,
        "s4_s5_correlations": correlations,
        "s4_s5_cohort_correlations": cohort_correlations,
        "map_constant_comparison": map_constant,
        "selected_locus_sensitivity": selected,
    }
    for name, frame in tables.items():
        atomic_write_text(output / f"{name}.tsv", frame.to_csv(sep="\t", index=False))

    malformed_counts = malformed.get("_exclusion_reason", pd.Series(dtype=str)).value_counts().to_dict()
    excluded_counts = excluded.get("_exclusion_reason", pd.Series(dtype=str)).value_counts().to_dict()
    image = _figure(survival, map_constant, bridge, threshold)
    strongest = correlations.iloc[0].to_dict() if not correlations.empty else {}
    broad = overall.loc[overall["analysis_set"].eq("denisovan_broad")].iloc[0].to_dict()
    conclusion = (
        "The map-aware S5 distribution remains incompatible with a stable single "
        "exponential: the effective decay rate changes with the detection threshold, "
        "and the bounded fit is not an event-date estimate."
    )
    summary = {
        "status": "complete_not_estimable",
        "created_at": utc_now(),
        "papuan_individuals": int(retained["sample_id"].nunique()),
        "valid_map_intervals": int(len(imported)),
        "malformed_intervals": int(len(malformed)),
        "qc_retained_all_archaic": int(len(retained)),
        "qc_excluded": int(len(excluded)),
        "malformed_reasons": malformed_counts,
        "qc_exclusion_reasons": excluded_counts,
        "broad_denisovan": broad,
        "strongest_s4_s5_correlation": strongest,
        "selected_locus_sensitivity_status": str(selected.iloc[0]["status"]),
        "conclusion": conclusion,
        "interpretation_guardrail": (
            "Effective decay rates and optimizer outputs are diagnostics, not biological "
            "admixture dates, unless goodness-of-fit and simulation recovery pass."
        ),
    }
    atomic_write_json(output / "summary.json", summary)
    manifest = {
        "created_at": utc_now(),
        "inputs": {
            "s4": _portable_fingerprint(s4_path),
            "s5": _portable_fingerprint(s5_path),
        },
        "genetic_map": {
            "build": "GRCh37",
            "source": "HapMap Phase II map lifted by Adam Auton",
            "files": [
                _portable_fingerprint(path)
                for path in sorted(Path(genetic_map_directory).glob("genetic_map_GRCh37_chr*.txt.gz"))
            ],
        },
        "parameters": {
            "generation_time_years": GENERATION_TIME_YEARS,
            "minimum_length_cm": MINIMUM_LENGTH_CM,
            "minimum_archaic_posterior": 0.8,
            "autosomes": AUTOSOMES,
            "random_seed": seed,
        },
    }
    atomic_write_json(output / "provenance.json", manifest)

    report_markdown = f"""# Papuan Denisovan S5 failure anatomy

Status: **not estimable**

{conclusion}

## Core counts

{_markdown_table(overall)}

## S4-S5 bridge

The bridge matched all **{len(bridge)}** published `Whole~world` S4 individuals to
their map-aware S5 tract features. Correlations quantify concordance only; the S4
HMM parameter is not treated as an independent tract-length estimate.

{_markdown_table(correlations, 12)}

## Detection-threshold stability

{_markdown_table(threshold.loc[threshold['analysis_set'].eq('denisovan_broad')])}

## Long-tail sensitivity

{_markdown_table(trim)}

## Interpretation

- `MeanProb` is a generic archaic posterior, not a Denisovan probability.
- Denisovan and Neanderthal labels are affinity rules based on relative sharing.
- Boundary fits and threshold-dependent decay rates are rejected as event dates.
- Published S4 HMM parameters are compared for concordance, not re-estimated.
- Selected-locus result: `{selected.iloc[0]['status']}`.
"""
    atomic_write_text(output / "report.md", report_markdown)
    report_html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Papuan Denisovan failure anatomy</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;line-height:1.5}}h1,h2{{color:#24324a}}.verdict{{padding:1rem;background:#fff2cc;border-left:5px solid #c58b00}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}th,td{{border:1px solid #ddd;padding:.35rem;text-align:right}}th:first-child,td:first-child{{text-align:left}}img{{max-width:100%}}code{{background:#f3f3f3;padding:.1rem .3rem}}</style></head>
<body><h1>Papuan Denisovan S5 failure anatomy</h1>
<div class="verdict"><strong>Verdict: not estimable.</strong> {html.escape(conclusion)}</div>
<p>All 89 S4 individuals were matched to S5. This report diagnoses why the exported
tract distribution cannot support an independent biological event date.</p>
<img alt="Failure anatomy diagnostics" src="data:image/png;base64,{image}">
<h2>Analysis sets</h2>{_table_html(overall, 20)}
<h2>S4-S5 correlations</h2>{_table_html(correlations, 15)}
<h2>Cohort summaries</h2>{_table_html(by_dataset, 20)}
<h2>Threshold sensitivity</h2>{_table_html(threshold.loc[threshold['analysis_set'].eq('denisovan_broad')], 20)}
<h2>Longest-tract trimming</h2>{_table_html(trim, 10)}
<h2>Tail contributors</h2>{_table_html(tail_individual, 15)}
<h2>Scientific boundary</h2><p><code>MeanProb</code> remains a generic archaic
posterior. Relative affinity is not an ancestry percentage. Boundary estimates,
failed goodness-of-fit, and poor recovery remain <strong>not estimable</strong>.</p>
</body></html>"""
    atomic_write_text(output / "report.html", report_html)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s5", required=True)
    parser.add_argument("--s4", required=True)
    parser.add_argument("--genetic-map-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--selected-loci-bed")
    parser.add_argument("--seed", type=int, default=20260807)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = run_failure_anatomy(
        s5_path=args.s5,
        s4_path=args.s4,
        genetic_map_directory=args.genetic_map_dir,
        output_directory=args.output,
        selected_loci_bed=args.selected_loci_bed,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
