from __future__ import annotations

import numpy as np
import pandas as pd

from archaic_admixture_dating.bootstrap import bootstrap_fit
from archaic_admixture_dating.dating_single_pulse import fit_single_pulse
from archaic_admixture_dating.dating_two_pulse import fit_two_pulse
from archaic_admixture_dating.model_comparison import compare_models


def _lengths(rate, n, seed, threshold=0.02):
    rng = np.random.default_rng(seed)
    return threshold + rng.exponential(1 / rate, size=n) * 100


def test_single_pulse_recovers_approximate_true_date():
    true_generations = 1500
    fit = fit_single_pulse(_lengths(true_generations, 6000, 1), generation_time_years=29)
    assert abs(fit["generations"] - true_generations) / true_generations < 0.05
    assert abs(fit["years"] - fit["generations"] * 29) < 1e-9


def test_truncation_correction_is_memoryless_and_sensible():
    values = _lengths(1200, 10000, 2, threshold=0.02)
    high_threshold = 0.08
    retained = values[values >= high_threshold]
    fit = fit_single_pulse(retained, minimum_length_cm=high_threshold)
    assert abs(fit["generations"] - 1200) / 1200 < 0.08


def test_two_pulse_data_do_not_always_collapse():
    rng = np.random.default_rng(3)
    older = _lengths(1900, 7000, 4)
    younger = _lengths(500, 3000, 5)
    values = np.concatenate([older, younger])
    rng.shuffle(values)
    fit = fit_two_pulse(values, minimum_separation_generations=100)
    assert fit["separation_generations"] > 100
    assert "pulse_dates_not_separable" not in fit["warning_flags"]
    assert 0.05 < fit["weight_older"] < 0.95


def test_one_pulse_does_not_prefer_false_second_pulse_by_bic():
    values = _lengths(1300, 5000, 6)
    table, _ = compare_models(
        values,
        minimum_length_cm=0.02,
        generation_time_years=29,
    )
    single = table.set_index("model_id").loc["single_pulse"]
    two = table.set_index("model_id").loc["two_pulse"]
    assert single["bic"] < two["bic"]


def test_bootstrap_is_reproducible_with_fixed_seed():
    values = _lengths(1200, 1000, 7)
    frame = pd.DataFrame(
        {
            "length_cm": values,
            "chromosome": np.tile(np.arange(1, 21), 50).astype(str),
        }
    )
    kwargs = {"minimum_length_cm": 0.02, "generation_time_years": 29}
    first = bootstrap_fit(frame, fit_single_pulse, replicates=20, seed=99, fitter_kwargs=kwargs)
    second = bootstrap_fit(frame, fit_single_pulse, replicates=20, seed=99, fitter_kwargs=kwargs)
    pd.testing.assert_frame_equal(first, second)
