"""Automatic quality and interpretation warnings."""

from __future__ import annotations

import pandas as pd


def collect_warnings(
    tracts: pd.DataFrame,
    model_table: pd.DataFrame | None = None,
    bootstrap_table: pd.DataFrame | None = None,
    calibration_table: pd.DataFrame | None = None,
) -> list[str]:
    warnings: list[str] = []
    if len(tracts) < 50:
        warnings.append("inadequate sample size: fewer than 50 retained tracts")
    if tracts["sample_id"].nunique() < 5:
        warnings.append("inadequate individual coverage: fewer than five individuals")
    callable_values = tracts["callable_fraction"].dropna()
    if len(callable_values) and callable_values.mean() < 0.8:
        warnings.append("low mean callable fraction")
    by_chromosome = tracts.groupby("chromosome")["length_cm"].mean()
    if len(by_chromosome) > 2 and by_chromosome.std() > by_chromosome.mean():
        warnings.append("strong disagreement across chromosomes")
    if len(tracts) and tracts.nlargest(max(1, int(len(tracts) * 0.01)), "length_cm")["length_cm"].sum() > 0.25 * tracts["length_cm"].sum():
        warnings.append("result may depend on a few long tracts")
    if model_table is not None and len(model_table) > 1 and model_table.iloc[1]["delta_bic"] < 2:
        warnings.append("tested models are not distinguishable by BIC")
    if bootstrap_table is not None and len(bootstrap_table):
        failed = (bootstrap_table["status"] != "ok").mean()
        if failed > 0.1:
            warnings.append("unstable pulse estimates: over 10% of bootstrap fits failed")
    if calibration_table is not None and len(calibration_table):
        valid = calibration_table.loc[calibration_table["status"] == "ok"]
        correct = valid["correct"]
        if correct.dtype != bool:
            correct = correct.astype(str).str.lower().map({"true": True, "false": False})
        if len(valid) and correct.dropna().mean() < 0.7:
            warnings.append("poor simulation classification accuracy under configured scenarios")
        measured = valid["relative_parameter_error"].dropna()
        if len(measured) and (measured > 0.2).mean() > 0.3:
            warnings.append("poor parameter recovery in over 30% of calibrated simulation fits")
    warnings.append(
        "direct late Denisovan admixture is not established unless modern-human mixing, "
        "bottlenecks, selection, map error, and tract-caller error are rejected"
    )
    return warnings


def interpretation_status(warnings: list[str]) -> str:
    severe = (
        "inadequate sample size",
        "not distinguishable",
        "unstable pulse",
        "low mean callable",
        "poor simulation",
        "poor parameter recovery",
    )
    return "inconclusive/data-limited" if any(any(term in warning for term in severe) for warning in warnings) else "exploratory"
