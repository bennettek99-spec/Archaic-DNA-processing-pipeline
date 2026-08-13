from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from archaic_admixture_dating.caller_calibration import (
    REAL_DECODED_DECAY,
    fit_curve,
    invert,
    required_inflation,
)


def _synthetic_curve(slope=0.6, intercept=400.0, noise=20.0, seed=1):
    """A calibration table with a known linear truth."""
    rng = np.random.default_rng(seed)
    rows = []
    for t in (600.0, 900.0, 1200.0, 1500.0, 1800.0):
        for r in range(4):
            decoded = slope * t + intercept + rng.normal(0, noise)
            rows.append(
                {
                    "pulse_generations": t,
                    "decoded_decay": decoded,
                    "fitted_generations": decoded / 0.78,
                    "decoded_over_fitted": 0.78,
                }
            )
    return pd.DataFrame(rows)


def test_fit_curve_recovers_known_slope_and_intercept():
    table = _synthetic_curve(slope=0.6, intercept=400.0, noise=5.0)
    slope, intercept = fit_curve(table)
    assert abs(slope - 0.6) < 0.02
    assert abs(intercept - 400.0) < 25.0


def test_inversion_recovers_a_planted_date():
    table = _synthetic_curve(slope=0.6, intercept=400.0, noise=5.0)
    planted = 0.6 * 1300.0 + 400.0
    result = invert(table, planted, n_boot=200)
    assert abs(result["point_estimate_generations"] - 1300.0) < 60.0


def test_inversion_interval_brackets_the_point_estimate():
    table = _synthetic_curve(noise=40.0)
    result = invert(table, 1200.0, n_boot=400)
    low, high = result["ci95_generations"]
    assert low <= result["point_estimate_generations"] <= high


def test_inversion_rejects_a_degenerate_curve():
    flat = pd.DataFrame(
        {
            "pulse_generations": [600.0, 900.0, 1200.0],
            "decoded_decay": [800.0, 800.0, 800.0],
            "decoded_over_fitted": [0.78, 0.78, 0.78],
            "fitted_generations": [1025.0, 1025.0, 1025.0],
        }
    )
    with pytest.raises(ValueError):
        invert(flat, 800.0)


def test_required_inflation_reports_the_shortfall():
    """An observation below the curve needs extra compression, not less."""
    table = _synthetic_curve(slope=0.6, intercept=400.0, noise=5.0)
    result = required_inflation(table, observed_decay=655.3, candidate_generations=1550.0)
    assert result["predicted_decoded_decay"] > result["observed_decoded_decay"]
    assert result["extra_compression_needed"] < 1.0
    assert result["required_total_inflation"] > result["simulated_inflation"]


def test_required_inflation_is_neutral_when_the_curve_already_fits():
    table = _synthetic_curve(slope=0.6, intercept=400.0, noise=1.0)
    on_curve = 0.6 * 1550.0 + 400.0
    result = required_inflation(table, on_curve, candidate_generations=1550.0)
    assert abs(result["extra_compression_needed"] - 1.0) < 0.02


def test_real_anchor_constant_is_the_documented_value():
    assert REAL_DECODED_DECAY == pytest.approx(655.3)
