"""Model fitting and complexity-aware comparison for observed tract lengths."""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np
import pandas as pd

from .dating_continuous import continuous_log_likelihood, fit_continuous_flow
from .dating_single_pulse import (
    fit_single_pulse,
    prepare_excess_morgans,
    single_log_likelihood,
)
from .dating_two_pulse import fit_two_pulse, mixture_log_likelihood


def _row(
    model_id: str,
    name: str,
    fit: dict[str, Any],
    parameters: int,
    heldout_score: float = np.nan,
) -> dict[str, Any]:
    n = int(fit["n_tracts"])
    ll = float(fit["log_likelihood"])
    return {
        "model_id": model_id,
        "model_name": name,
        "parameters": parameters,
        "log_likelihood": ll,
        "aic": 2 * parameters - 2 * ll,
        "bic": math.log(max(n, 1)) * parameters - 2 * ll,
        "cross_validated_score": heldout_score,
        "bootstrap_support": np.nan,
        "parameter_recovery_quality": "not_yet_calibrated",
        "warning_flags": ";".join(fit.get("warning_flags", [])),
        "fit_json": json.dumps(fit, sort_keys=True),
    }


def compare_models(
    lengths_cm,
    *,
    minimum_length_cm: float,
    generation_time_years: float,
    single_bounds: tuple[float, float] = (100, 4000),
    two_minimum_separation_generations: float = 100,
) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    values = np.asarray(lengths_cm, dtype=float)
    values = values[np.isfinite(values)]
    fits = {
        "single_pulse": fit_single_pulse(
            values,
            minimum_length_cm=minimum_length_cm,
            generation_time_years=generation_time_years,
            bounds_generations=single_bounds,
        ),
        "two_pulse": fit_two_pulse(
            values,
            minimum_length_cm=minimum_length_cm,
            generation_time_years=generation_time_years,
            minimum_separation_generations=two_minimum_separation_generations,
        ),
        "continuous_flow": fit_continuous_flow(
            values,
            minimum_length_cm=minimum_length_cm,
            generation_time_years=generation_time_years,
        ),
    }
    heldout = {"single_pulse": np.nan, "two_pulse": np.nan, "continuous_flow": np.nan}
    if len(values) >= 20:
        test_mask = np.arange(len(values)) % 5 == 0
        train, test = values[~test_mask], values[test_mask]
        test_excess = prepare_excess_morgans(test, minimum_length_cm)
        train_fits = {
            "single_pulse": fit_single_pulse(
                train,
                minimum_length_cm=minimum_length_cm,
                generation_time_years=generation_time_years,
                bounds_generations=single_bounds,
            ),
            "two_pulse": fit_two_pulse(
                train,
                minimum_length_cm=minimum_length_cm,
                generation_time_years=generation_time_years,
                minimum_separation_generations=two_minimum_separation_generations,
            ),
            "continuous_flow": fit_continuous_flow(
                train,
                minimum_length_cm=minimum_length_cm,
                generation_time_years=generation_time_years,
            ),
        }
        heldout["single_pulse"] = single_log_likelihood(
            test_excess, train_fits["single_pulse"]["generations"]
        ) / len(test_excess)
        heldout["two_pulse"] = mixture_log_likelihood(
            test_excess,
            train_fits["two_pulse"]["older_generations"],
            train_fits["two_pulse"]["younger_generations"],
            train_fits["two_pulse"]["weight_older"],
        ) / len(test_excess)
        heldout["continuous_flow"] = continuous_log_likelihood(
            test_excess,
            train_fits["continuous_flow"]["older_generations"],
            train_fits["continuous_flow"]["younger_generations"],
        ) / len(test_excess)
    rows = [
        _row("single_pulse", "Single pulse", fits["single_pulse"], 1, heldout["single_pulse"]),
        _row("two_pulse", "Two discrete pulses", fits["two_pulse"], 3, heldout["two_pulse"]),
        _row(
            "continuous_flow",
            "Uniform prolonged flow",
            fits["continuous_flow"],
            2,
            heldout["continuous_flow"],
        ),
    ]
    table = pd.DataFrame(rows).sort_values(["bic", "aic"]).reset_index(drop=True)
    table["delta_bic"] = table["bic"] - table["bic"].min()
    table["bic_plausibility"] = np.where(
        table["delta_bic"] < 2,
        "competitive",
        np.where(table["delta_bic"] < 6, "plausible", "disfavored_under_tested_assumptions"),
    )
    if len(table) > 1 and table.loc[1, "delta_bic"] < 2:
        table.loc[:1, "warning_flags"] = table.loc[:1, "warning_flags"].map(
            lambda value: ";".join(filter(None, [value, "models_not_distinguishable"]))
        )
    return table, fits
