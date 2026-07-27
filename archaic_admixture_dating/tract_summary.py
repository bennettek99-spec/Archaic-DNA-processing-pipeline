"""Observed tract summaries used by QC, simulation, and model comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _summarize(group: pd.DataFrame) -> pd.Series:
    lengths = group["length_cm"].to_numpy(dtype=float)
    return pd.Series(
        {
            "tract_count": len(group),
            "total_length_cm": float(np.sum(lengths)),
            "mean_length_cm": float(np.mean(lengths)) if len(lengths) else np.nan,
            "median_length_cm": float(np.median(lengths)) if len(lengths) else np.nan,
            "q75_length_cm": float(np.quantile(lengths, 0.75)) if len(lengths) else np.nan,
            "q90_length_cm": float(np.quantile(lengths, 0.90)) if len(lengths) else np.nan,
            "q95_length_cm": float(np.quantile(lengths, 0.95)) if len(lengths) else np.nan,
            "tracts_ge_0_1cm": int(np.sum(lengths >= 0.1)),
            "tracts_ge_0_2cm": int(np.sum(lengths >= 0.2)),
            "mean_callable_fraction": float(group["callable_fraction"].mean()),
            "mean_posterior_denisovan": float(group["posterior_denisovan"].mean()),
        }
    )


def summarize_tracts(tracts: pd.DataFrame) -> dict[str, pd.DataFrame]:
    if tracts.empty:
        empty = pd.DataFrame()
        return {
            "overall": empty,
            "individual": empty,
            "population": empty,
            "chromosome": empty,
            "sample_missingness": empty,
            "caller_batch": empty,
            "source_class": empty,
        }
    overall = _summarize(tracts).to_frame().T
    individual = (
        tracts.groupby(["sample_id", "population"], dropna=False)
        .apply(_summarize)
        .reset_index()
    )
    population = (
        tracts.groupby(["population"], dropna=False)
        .apply(_summarize)
        .reset_index()
    )
    chromosome = (
        tracts.groupby(["chromosome"], dropna=False)
        .apply(_summarize)
        .reset_index()
    )
    sample_missingness = (
        tracts.groupby(["sample_id", "population"], dropna=False)
        .agg(
            mean_callable_fraction=("callable_fraction", "mean"),
            minimum_callable_fraction=("callable_fraction", "min"),
            tract_count=("length_cm", "size"),
        )
        .reset_index()
    )
    sample_missingness["estimated_missing_fraction"] = (
        1.0 - sample_missingness["mean_callable_fraction"]
    )
    caller_batch = (
        tracts.groupby(["caller"], dropna=False)
        .apply(_summarize)
        .reset_index()
    )
    source_class = (
        tracts.groupby(["source_class"], dropna=False)
        .apply(_summarize)
        .reset_index()
    )
    return {
        "overall": overall,
        "individual": individual,
        "population": population,
        "chromosome": chromosome,
        "sample_missingness": sample_missingness,
        "caller_batch": caller_batch,
        "source_class": source_class,
    }


def survival_curve(lengths_cm: np.ndarray, points: int = 100) -> pd.DataFrame:
    values = np.sort(np.asarray(lengths_cm, dtype=float))
    if len(values) == 0:
        return pd.DataFrame(columns=["length_cm", "survival_probability"])
    grid = np.linspace(values.min(), values.max(), points)
    survival = np.asarray([(values >= threshold).mean() for threshold in grid])
    return pd.DataFrame({"length_cm": grid, "survival_probability": survival})
