"""Left-truncated single-pulse tract-length dating."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats


def prepare_excess_morgans(lengths_cm, minimum_length_cm: float) -> np.ndarray:
    values = np.asarray(lengths_cm, dtype=float)
    values = values[np.isfinite(values)]
    values = values[values >= minimum_length_cm]
    if len(values) == 0:
        raise ValueError("No finite tracts remain above the minimum length")
    excess = (values - minimum_length_cm) / 100.0
    if np.sum(excess) <= 0:
        raise ValueError("Tracts contain no length information above the truncation threshold")
    return excess


def single_log_likelihood(excess_morgans: np.ndarray, generations: float) -> float:
    if generations <= 0:
        return -math.inf
    return float(len(excess_morgans) * math.log(generations) - generations * np.sum(excess_morgans))


def fit_single_pulse(
    lengths_cm,
    *,
    minimum_length_cm: float = 0.02,
    generation_time_years: float = 29.0,
    bounds_generations: tuple[float, float] = (100.0, 4000.0),
) -> dict[str, Any]:
    excess = prepare_excess_morgans(lengths_cm, minimum_length_cm)
    raw = len(excess) / float(np.sum(excess))
    generations = float(np.clip(raw, *bounds_generations))
    log_likelihood = single_log_likelihood(excess, generations)
    se_generations = generations / math.sqrt(len(excess))
    ci_low = max(bounds_generations[0], generations - 1.96 * se_generations)
    ci_high = min(bounds_generations[1], generations + 1.96 * se_generations)
    ks = stats.kstest(excess, "expon", args=(0.0, 1.0 / generations))
    warnings: list[str] = []
    if generations in bounds_generations:
        warnings.append("estimate_at_parameter_bound")
    if len(excess) < 50:
        warnings.append("inadequate_sample_size")
    if ks.pvalue < 0.01:
        warnings.append("poor_single_exponential_fit")
    return {
        "model_id": "single_pulse",
        "n_tracts": len(excess),
        "minimum_length_cm": minimum_length_cm,
        "generations": generations,
        "years": generations * generation_time_years,
        "kya": generations * generation_time_years / 1000.0,
        "ci_low_generations": ci_low,
        "ci_high_generations": ci_high,
        "ci_low_kya": ci_low * generation_time_years / 1000.0,
        "ci_high_kya": ci_high * generation_time_years / 1000.0,
        "log_likelihood": log_likelihood,
        "ks_statistic": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
        "warning_flags": warnings,
    }
