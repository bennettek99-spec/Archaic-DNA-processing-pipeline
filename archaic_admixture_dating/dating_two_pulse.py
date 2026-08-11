"""Ordered two-component mixture of truncated tract-length exponentials."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

from .dating_single_pulse import prepare_excess_morgans


def _unpack(parameters: np.ndarray) -> tuple[float, float, float]:
    younger = math.exp(float(parameters[0]))
    older = younger + math.exp(float(parameters[1]))
    weight_older = float(expit(parameters[2]))
    return older, younger, weight_older


def mixture_log_likelihood(
    excess_morgans: np.ndarray,
    older_generations: float,
    younger_generations: float,
    weight_older: float,
) -> float:
    if older_generations <= younger_generations or not 0 < weight_older < 1:
        return -math.inf
    components = np.vstack(
        [
            math.log(weight_older) + math.log(older_generations) - older_generations * excess_morgans,
            math.log1p(-weight_older) + math.log(younger_generations) - younger_generations * excess_morgans,
        ]
    )
    return float(np.sum(logsumexp(components, axis=0)))


def fit_two_pulse(
    lengths_cm,
    *,
    minimum_length_cm: float = 0.02,
    generation_time_years: float = 29.0,
    minimum_separation_generations: float = 100.0,
) -> dict[str, Any]:
    excess = prepare_excess_morgans(lengths_cm, minimum_length_cm)

    def objective(parameters: np.ndarray) -> float:
        older, younger, weight = _unpack(parameters)
        if younger < 50 or older > 6000:
            return 1e12
        value = -mixture_log_likelihood(excess, older, younger, weight)
        return value if np.isfinite(value) else 1e12

    starts = [
        (1700, 1050, 0.65),
        (1400, 850, 0.5),
        (2200, 1100, 0.8),
        (1100, 600, 0.35),
    ]
    fits = []
    for older, younger, weight in starts:
        encoded = np.array(
            [math.log(younger), math.log(older - younger), math.log(weight / (1 - weight))]
        )
        fits.append(
            minimize(
                objective,
                encoded,
                method="L-BFGS-B",
                bounds=[
                    (math.log(50), math.log(4000)),
                    (math.log(1), math.log(5950)),
                    (-8, 8),
                ],
            )
        )
    best = min(fits, key=lambda result: result.fun)
    older, younger, weight = _unpack(best.x)
    separation = older - younger
    warnings: list[str] = []
    if not best.success:
        warnings.append("optimizer_not_converged")
    if weight < 0.05 or weight > 0.95:
        warnings.append("component_weight_near_zero")
    if separation < minimum_separation_generations:
        warnings.append("pulse_dates_not_separable")
    if len(excess) < 100:
        warnings.append("inadequate_sample_size_for_two_pulse_model")
    return {
        "model_id": "two_pulse",
        "n_tracts": len(excess),
        "minimum_length_cm": minimum_length_cm,
        "older_generations": older,
        "younger_generations": younger,
        "older_kya": older * generation_time_years / 1000.0,
        "younger_kya": younger * generation_time_years / 1000.0,
        "weight_older": weight,
        "separation_generations": separation,
        "log_likelihood": -float(best.fun),
        "converged": bool(best.success),
        "warning_flags": warnings,
    }
