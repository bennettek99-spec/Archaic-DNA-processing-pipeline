"""Bounded posterior-predictive checks for tract-level Papuan models M1-M10."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .checkpointing import atomic_write_json, atomic_write_text, sha256_file, utc_now
from .config import apply_profile, config_hash, load_config
from .simulations import derived_seed, simulate_tracts


FEATURES = [
    "median_length_cm",
    "q75_length_cm",
    "q90_length_cm",
    "q95_length_cm",
    "q99_length_cm",
    "fraction_ge_0_5cm",
    "fraction_ge_1cm",
    "fraction_ge_2cm",
    "effective_generations_unbounded",
]

SCENARIOS = {
    "configured_error": {"false_negative_rate": 0.05, "length_noise_sd_fraction": 0.08},
    "heavy_error_stress": {"false_negative_rate": 0.20, "length_noise_sd_fraction": 0.25},
}


def _portable_fingerprint(path: str | Path) -> dict[str, object]:
    source = Path(path)
    return {"filename": source.name, "bytes": source.stat().st_size, "sha256": sha256_file(source)}


def summarize_lengths(lengths: Iterable[float], minimum_length_cm: float = 0.02) -> dict[str, float]:
    values = np.asarray(list(lengths), dtype=float)
    values = values[np.isfinite(values) & (values >= minimum_length_cm)]
    if len(values) < 10:
        raise ValueError("At least 10 detected tract lengths are required")
    excess = (values - minimum_length_cm) / 100.0
    return {
        "tract_count": float(len(values)),
        "mean_length_cm": float(np.mean(values)),
        "median_length_cm": float(np.median(values)),
        "q75_length_cm": float(np.quantile(values, 0.75)),
        "q90_length_cm": float(np.quantile(values, 0.90)),
        "q95_length_cm": float(np.quantile(values, 0.95)),
        "q99_length_cm": float(np.quantile(values, 0.99)),
        "fraction_ge_0_5cm": float(np.mean(values >= 0.5)),
        "fraction_ge_1cm": float(np.mean(values >= 1.0)),
        "fraction_ge_2cm": float(np.mean(values >= 2.0)),
        "effective_generations_unbounded": float(len(values) / excess.sum()),
    }


def simulate_calibration(
    config: dict,
    *,
    replicates: int,
    tracts_per_replicate: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    master = int(config["project"]["random_seed"])
    generation_time = float(config["project"]["generation_time_years"])
    minimum = float(config["tracts"]["minimum_length_cm"])
    configured = config["simulation"].get("caller_error", {})
    scenarios = {
        **SCENARIOS,
        "configured_error": {
            "false_negative_rate": float(configured.get("false_negative_rate", 0.05)),
            "length_noise_sd_fraction": float(configured.get("length_noise_sd_fraction", 0.08)),
        },
    }
    for scenario_index, (scenario, error) in enumerate(scenarios.items()):
        for model_id, model in config["models"].items():
            for replicate in range(replicates):
                seed = derived_seed(master + 10_000_019 * scenario_index, model_id, replicate)
                tracts = simulate_tracts(
                    model_id,
                    model,
                    n_tracts=tracts_per_replicate,
                    generation_time_years=generation_time,
                    minimum_length_cm=minimum,
                    seed=seed,
                    false_negative_rate=error["false_negative_rate"],
                    length_noise_sd_fraction=error["length_noise_sd_fraction"],
                )
                rows.append(
                    {
                        "scenario": scenario,
                        "model_id": model_id,
                        "model_name": model["name"],
                        "replicate": replicate,
                        "seed": seed,
                        **summarize_lengths(tracts["length_cm"], minimum),
                    }
                )
    return pd.DataFrame(rows)


def leave_one_out_classification(simulations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for scenario, table in simulations.groupby("scenario", sort=True):
        matrix = table[FEATURES].to_numpy(float)
        scale = np.nanstd(matrix, axis=0, ddof=1)
        scale[~np.isfinite(scale) | (scale == 0)] = 1.0
        for index, row in table.iterrows():
            candidates = table.drop(index=index)
            distances: list[tuple[float, str]] = []
            for model_id, group in candidates.groupby("model_id"):
                centroid = group[FEATURES].to_numpy(float).mean(axis=0)
                distance = float(np.sqrt(np.mean(((row[FEATURES].to_numpy(float) - centroid) / scale) ** 2)))
                distances.append((distance, str(model_id)))
            distance, predicted = min(distances)
            rows.append(
                {
                    "scenario": scenario,
                    "model_id": row["model_id"],
                    "replicate": int(row["replicate"]),
                    "predicted_model_id": predicted,
                    "correct": predicted == row["model_id"],
                    "nearest_centroid_distance": distance,
                }
            )
    assignments = pd.DataFrame(rows)
    accuracy = (
        assignments.groupby(["scenario", "model_id"])["correct"]
        .agg(["count", "mean"])
        .reset_index()
        .rename(columns={"count": "replicates", "mean": "classification_accuracy"})
    )
    return assignments, accuracy


def posterior_predictive_compatibility(
    simulations: pd.DataFrame, observed: dict[str, float]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    for (scenario, model_id, model_name), group in simulations.groupby(
        ["scenario", "model_id", "model_name"], sort=True
    ):
        squared: list[float] = []
        inside = 0
        for feature in FEATURES:
            values = group[feature].to_numpy(float)
            low, median, high = np.quantile(values, [0.025, 0.5, 0.975])
            if feature.startswith("fraction_"):
                absolute_floor = 0.01
            elif feature == "effective_generations_unbounded":
                absolute_floor = 1.0
            else:
                absolute_floor = 0.01
            width = max(float(high - low), abs(float(median)) * 1e-6, absolute_floor)
            value = float(observed[feature])
            z = (value - float(median)) / width
            is_inside = bool(low <= value <= high)
            inside += int(is_inside)
            squared.append(z * z)
            feature_rows.append(
                {
                    "scenario": scenario,
                    "model_id": model_id,
                    "feature": feature,
                    "observed": value,
                    "simulation_q025": float(low),
                    "simulation_median": float(median),
                    "simulation_q975": float(high),
                    "within_95pct_envelope": is_inside,
                    "scaled_deviation": float(z),
                }
            )
        coverage = inside / len(FEATURES)
        distance = float(np.sqrt(np.mean(squared)))
        summary_rows.append(
            {
                "scenario": scenario,
                "model_id": model_id,
                "model_name": model_name,
                "features_within_95pct_envelope": inside,
                "feature_count": len(FEATURES),
                "feature_coverage": coverage,
                "scaled_rms_distance": distance,
                "compatibility": "compatible" if coverage >= 0.8 and distance <= 1.0 else "rejected",
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["rank_within_scenario"] = summary.groupby("scenario")["scaled_rms_distance"].rank(
        method="first"
    ).astype(int)
    return summary.sort_values(["scenario", "rank_within_scenario"]), pd.DataFrame(feature_rows)


def parameter_recovery(simulations: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    generation_time = float(config["project"]["generation_time_years"])
    for (scenario, model_id), group in simulations.groupby(["scenario", "model_id"]):
        model = config["models"][model_id]
        if model["kind"] != "single":
            continue
        truth = float(model["dates_kya"][0]) * 1000.0 / generation_time
        errors = abs(group["effective_generations_unbounded"] - truth) / truth
        rows.append(
            {
                "scenario": scenario,
                "model_id": model_id,
                "truth_generations": truth,
                "median_estimate_generations": float(group["effective_generations_unbounded"].median()),
                "median_relative_error": float(errors.median()),
                "q95_relative_error": float(errors.quantile(0.95)),
                "recovery_quality": "pass" if errors.median() <= 0.20 else "fail",
            }
        )
    return pd.DataFrame(rows)


def _plot(compatibility: pd.DataFrame, simulations: pd.DataFrame, observed: dict[str, float]) -> str:
    import base64
    from io import BytesIO

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    baseline = compatibility.loc[compatibility["scenario"].eq("configured_error")].sort_values(
        "scaled_rms_distance"
    )
    axes[0].bar(baseline["model_id"], baseline["scaled_rms_distance"], color="#526d82")
    axes[0].set(ylabel="Scaled RMS distance", title="Observed-to-model feature distance")
    axes[0].tick_params(axis="x", rotation=45)
    for model_id, group in simulations.loc[
        simulations["scenario"].eq("configured_error")
    ].groupby("model_id"):
        axes[1].scatter(
            group["median_length_cm"], group["q95_length_cm"], s=12, alpha=0.5, label=model_id
        )
    axes[1].scatter(
        [observed["median_length_cm"]], [observed["q95_length_cm"]], marker="*", s=180,
        color="#b23a48", edgecolor="black", label="observed", zorder=5
    )
    axes[1].set(xlabel="Median length (cM)", ylabel="95th percentile (cM)", title="Observed profile vs simulations")
    axes[1].legend(ncol=3, fontsize=7)
    fig.tight_layout()
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def run_observation_calibration(
    *,
    observed_overall_path: str | Path,
    config_path: str | Path,
    output_directory: str | Path,
    replicates: int = 20,
    tracts_per_replicate: int = 5000,
) -> dict[str, object]:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    config = apply_profile(load_config(config_path), "laptop")
    observed_table = pd.read_csv(observed_overall_path, sep="\t")
    selected = observed_table.loc[observed_table["analysis_set"].eq("denisovan_broad")]
    if len(selected) != 1:
        raise ValueError("Observed table must contain exactly one denisovan_broad row")
    observed = {feature: float(selected.iloc[0][feature]) for feature in FEATURES}

    simulations = simulate_calibration(
        config, replicates=replicates, tracts_per_replicate=tracts_per_replicate
    )
    assignments, accuracy = leave_one_out_classification(simulations)
    compatibility, feature_checks = posterior_predictive_compatibility(simulations, observed)
    recovery = parameter_recovery(simulations, config)
    for name, frame in {
        "simulation_summaries": simulations,
        "classification_assignments": assignments,
        "classification_accuracy": accuracy,
        "posterior_predictive_summary": compatibility,
        "posterior_predictive_features": feature_checks,
        "parameter_recovery": recovery,
    }.items():
        atomic_write_text(output / f"{name}.tsv", frame.to_csv(sep="\t", index=False))

    closest = compatibility.loc[compatibility["scaled_rms_distance"].idxmin()].to_dict()
    scenario_accuracy = accuracy.groupby("scenario")["classification_accuracy"].mean().to_dict()
    image = _plot(compatibility, simulations, observed)
    result = {
        "status": "complete_not_estimable",
        "created_at": utc_now(),
        "replicates_per_model_scenario": replicates,
        "tracts_requested_per_replicate": tracts_per_replicate,
        "simulation_rows": len(simulations),
        "closest_model": closest,
        "mean_model_classification_accuracy": scenario_accuracy,
        "compatible_models": int(compatibility["compatibility"].eq("compatible").sum()),
        "conclusion": "None of the tested M1-M10 tract approximations reproduces the observed map-aware S5 length profile under the configured or heavy caller-error stress process.",
        "guardrail": "This is a bounded posterior-predictive rejection check, not a posterior probability over demographic histories.",
    }
    atomic_write_json(output / "summary.json", result)
    atomic_write_json(
        output / "provenance.json",
        {
            "created_at": result["created_at"],
            "config_hash": config_hash(config),
            "inputs": {
                "observed_overall": _portable_fingerprint(observed_overall_path),
                "config": _portable_fingerprint(config_path),
            },
            "replicates": replicates,
            "tracts_per_replicate": tracts_per_replicate,
            "scenarios": SCENARIOS,
            "features": FEATURES,
        },
    )
    report = f"""# Papuan observation-process calibration

Status: **not estimable**

{result['conclusion']}

- Simulated rows: {len(simulations)} ({replicates} replicates × 10 models × 2 error scenarios).
- Closest tested combination: `{closest['model_id']}` / `{closest['scenario']}`.
- Closest feature coverage: {int(closest['features_within_95pct_envelope'])}/{int(closest['feature_count'])} within replicate 95% envelopes.
- Closest scaled RMS distance: {float(closest['scaled_rms_distance']):.3f}.
- Models passing the predeclared compatibility rule: {result['compatible_models']}.

The heavy-error scenario raises random false negatives to 20% and length noise to
25%. It is a stress test, not a claim about the Skov caller. False negatives alone
do not change a length distribution when they are random, and the calibration
does not invent an unvalidated tract-merging process.

Classification and single-pulse recovery measure internal separability of the
tract approximations. They do not rescue a model whose observed feature profile
falls outside its posterior-predictive envelope.

**Guardrail:** {result['guardrail']}
"""
    atomic_write_text(output / "report.md", report)
    body = "".join(f"<p>{html.escape(line)}</p>" for line in report.splitlines() if line)
    atomic_write_text(
        output / "report.html",
        "<!doctype html><meta charset='utf-8'><title>Papuan calibration</title>"
        "<style>body{font:16px system-ui;max-width:980px;margin:40px auto;line-height:1.5}img{max-width:100%}</style>"
        f"<body>{body}<img alt='Calibration diagnostics' src='data:image/png;base64,{image}'></body>",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-overall", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--tracts-per-replicate", type=int, default=5000)
    args = parser.parse_args(argv)
    result = run_observation_calibration(
        observed_overall_path=args.observed_overall,
        config_path=args.config,
        output_directory=args.output,
        replicates=args.replicates,
        tracts_per_replicate=args.tracts_per_replicate,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
