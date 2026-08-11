"""Uniform-interval prolonged-flow approximation."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from .dating_single_pulse import prepare_excess_morgans


def continuous_log_likelihood(
    excess_morgans: np.ndarray,
    older_generations: float,
    younger_generations: float,
    grid_points: int = 64,
) -> float:
    if older_generations <= younger_generations or younger_generations <= 0:
        return -math.inf
    rates = np.linspace(younger_generations, older_generations, grid_points)
    log_density = np.log(rates)[:, None] - rates[:, None] * excess_morgans[None, :]
    return float(np.sum(logsumexp(log_density, axis=0) - math.log(grid_points)))


def fit_continuous_flow(
    lengths_cm,
    *,
    minimum_length_cm: float = 0.02,
    generation_time_years: float = 29.0,
) -> dict[str, Any]:
    excess = prepare_excess_morgans(lengths_cm, minimum_length_cm)

    def unpack(parameters: np.ndarray) -> tuple[float, float]:
        younger = math.exp(float(parameters[0]))
        older = younger + math.exp(float(parameters[1]))
        return older, younger

    def objective(parameters: np.ndarray) -> float:
        older, younger = unpack(parameters)
        if younger < 50 or older > 6000:
            return 1e12
        value = -continuous_log_likelihood(excess, older, younger)
        return value if np.isfinite(value) else 1e12

    result = minimize(
        objective,
        np.log([900.0, 700.0]),
        method="L-BFGS-B",
        bounds=[
            (math.log(50), math.log(4000)),
            (math.log(1), math.log(5950)),
        ],
    )
    older, younger = unpack(result.x)
    warnings: list[str] = []
    if not result.success:
        warnings.append("optimizer_not_converged")
    if older - younger < 100:
        warnings.append("flow_interval_collapsed")
    if len(excess) < 100:
        warnings.append("inadequate_sample_size_for_continuous_flow")
    return {
        "model_id": "continuous_flow",
        "n_tracts": len(excess),
        "minimum_length_cm": minimum_length_cm,
        "older_generations": older,
        "younger_generations": younger,
        "older_kya": older * generation_time_years / 1000.0,
        "younger_kya": younger * generation_time_years / 1000.0,
        "duration_generations": older - younger,
        "log_likelihood": -float(result.fun),
        "converged": bool(result.success),
        "warning_flags": warnings,
    }
