"""Simulation classification, confusion, and parameter-recovery summaries."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .model_comparison import compare_models
from .tract_schema import read_tracts


def expected_family(model: dict) -> str:
    if model["kind"] in {"single", "bottleneck", "selection"}:
        return "single_pulse"
    if model["kind"] in {"two", "modern_mixing", "divergent_sources"}:
        return "two_pulse"
    if model["kind"] == "continuous":
        return "continuous_flow"
    return "unknown"


def calibrate_simulations(output_dir: str | Path, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    simulation_dir = Path(output_dir) / "simulations"
    rows: list[dict[str, object]] = []
    generation_time = float(config["project"]["generation_time_years"])
    minimum_length = float(config["tracts"]["minimum_length_cm"])
    for path in sorted(simulation_dir.glob("M*_replicate_*.tsv")):
        stem = path.stem
        model_id = stem.split("_", 1)[0]
        replicate = int(stem.rsplit("_", 1)[1])
        tracts = read_tracts(path)
        if len(tracts) < 10:
            rows.append(
                {
                    "model_id": model_id,
                    "replicate": replicate,
                    "status": "insufficient_detected_tracts",
                }
            )
            continue
        table, fits = compare_models(
            tracts["length_cm"],
            minimum_length_cm=minimum_length,
            generation_time_years=generation_time,
            single_bounds=tuple(config["dating"]["single_pulse"]["bounds_generations"]),
            two_minimum_separation_generations=config["dating"]["two_pulse"]["minimum_separation_generations"],
        )
        predicted = str(table.iloc[0]["model_id"])
        model = config["models"][model_id]
        expected = expected_family(model)
        recovery_error = np.nan
        if model["kind"] == "single":
            truth = float(model["dates_kya"][0])
            recovery_error = abs(float(fits["single_pulse"]["kya"]) - truth) / truth
        elif model["kind"] in {"two", "divergent_sources"}:
            truth_old, truth_young = max(model["dates_kya"]), min(model["dates_kya"])
            fit = fits["two_pulse"]
            recovery_error = (
                abs(float(fit["older_kya"]) - truth_old) / truth_old
                + abs(float(fit["younger_kya"]) - truth_young) / truth_young
            ) / 2.0
        elif model["kind"] == "continuous":
            truth_old, truth_young = max(model["dates_kya"]), min(model["dates_kya"])
            fit = fits["continuous_flow"]
            recovery_error = (
                abs(float(fit["older_kya"]) - truth_old) / truth_old
                + abs(float(fit["younger_kya"]) - truth_young) / truth_young
            ) / 2.0
        rows.append(
            {
                "model_id": model_id,
                "replicate": replicate,
                "status": "ok",
                "expected_family": expected,
                "predicted_family": predicted,
                "correct": predicted == expected,
                "relative_parameter_error": recovery_error,
                "parameter_recovery_quality": (
                    "good"
                    if np.isfinite(recovery_error) and recovery_error <= 0.20
                    else "poor"
                    if np.isfinite(recovery_error)
                    else "not_applicable"
                ),
            }
        )
    recovery = pd.DataFrame(rows)
    valid = recovery.loc[recovery.get("status", pd.Series(dtype=str)).eq("ok")]
    if valid.empty:
        confusion = pd.DataFrame(
            columns=["expected_family", "predicted_family", "count", "fraction_within_expected"]
        )
    else:
        confusion = (
            valid.groupby(["expected_family", "predicted_family"])
            .size()
            .rename("count")
            .reset_index()
        )
        totals = confusion.groupby("expected_family")["count"].transform("sum")
        confusion["fraction_within_expected"] = confusion["count"] / totals
    return recovery, confusion
