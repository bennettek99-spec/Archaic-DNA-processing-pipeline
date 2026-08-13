from __future__ import annotations

import numpy as np

from archaic_admixture_dating.skov_hmm import (
    call_individual,
    decay_generations,
    extract_runs,
    fit_hmm,
    poisson_emissions,
    run_lengths_morgans,
)

WINDOW_BP = 1000
RECOMBINATION_RATE = 1.2e-8
# The two Poisson rates measured on the real 89 Papuan individuals.
MODERN_RATE = 0.0256
ARCHAIC_RATE = 0.2245


def _simulate_chain(generations, archaic_fraction, n_windows, seed):
    """Draw a two-state Markov chain with a known admixture parameter."""
    rng = np.random.default_rng(seed)
    p_leave = generations * RECOMBINATION_RATE * WINDOW_BP
    p_enter = p_leave * archaic_fraction / (1.0 - archaic_fraction)
    state = np.empty(n_windows, dtype=np.int8)
    current = 0
    draws = rng.random(n_windows)
    for i in range(n_windows):
        if current == 0:
            if draws[i] < p_enter:
                current = 1
        elif draws[i] < p_leave:
            current = 0
        state[i] = current
    rates = np.array([MODERN_RATE, ARCHAIC_RATE])
    return state, rng.poisson(rates[state])


def test_poisson_emissions_match_scipy():
    from scipy.stats import poisson

    obs = np.array([0, 1, 2, 5, 9])
    rates = np.array([0.03, 0.25])
    emissions = poisson_emissions(obs, rates)
    for j, rate in enumerate(rates):
        assert np.allclose(emissions[:, j], poisson.pmf(obs, rate))


def test_fit_recovers_known_emission_rates():
    _, obs = _simulate_chain(1500.0, 0.07, 400_000, seed=1)
    fit = fit_hmm(obs, window_bp=WINDOW_BP, recombination_rate=RECOMBINATION_RATE)
    assert abs(fit.rates[0] - MODERN_RATE) / MODERN_RATE < 0.15
    assert abs(fit.rates[1] - ARCHAIC_RATE) / ARCHAIC_RATE < 0.15
    assert fit.rates[1] > fit.rates[0]


def test_fit_recovers_known_admixture_parameter():
    _, obs = _simulate_chain(1500.0, 0.07, 600_000, seed=2)
    fit = fit_hmm(obs, window_bp=WINDOW_BP, recombination_rate=RECOMBINATION_RATE)
    assert abs(fit.fitted_generations - 1500.0) / 1500.0 < 0.20


def test_archaic_state_is_always_index_one():
    """Label switching must not leak into downstream code."""
    for seed in (3, 4, 5):
        _, obs = _simulate_chain(1200.0, 0.05, 120_000, seed=seed)
        fit = fit_hmm(obs, window_bp=WINDOW_BP, recombination_rate=RECOMBINATION_RATE)
        assert fit.rates[1] > fit.rates[0]


def test_decoding_inflates_run_lengths():
    """The distortion the calibration exists to measure must be present."""
    state, obs = _simulate_chain(1500.0, 0.07, 600_000, seed=6)
    result = call_individual(
        obs, window_bp=WINDOW_BP, recombination_rate=RECOMBINATION_RATE
    )
    padded = np.concatenate(([0], state, [0]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    true_median = np.median(edges[1::2] - edges[0::2])
    decoded_median = np.median(result["ends"] - result["starts"])
    assert decoded_median > true_median
    assert result["decoded_over_fitted"] < 1.0


def test_extract_runs_handles_edges_and_gaps():
    posterior = np.array([0.9, 0.9, 0.1, 0.8, 0.1, 0.99])
    starts, ends = extract_runs(posterior)
    assert starts.tolist() == [0, 3, 5]
    assert ends.tolist() == [2, 4, 6]


def test_extract_runs_empty_when_nothing_called():
    starts, ends = extract_runs(np.zeros(50))
    assert starts.size == 0 and ends.size == 0


def test_decay_generations_matches_truncated_exponential():
    rng = np.random.default_rng(9)
    threshold = 5e-4
    lengths = threshold + rng.exponential(1 / 800.0, size=40_000)
    assert abs(decay_generations(lengths, threshold) - 800.0) / 800.0 < 0.05


def test_decay_generations_returns_nan_when_nothing_survives():
    assert np.isnan(decay_generations(np.array([1e-5, 2e-5]), 1e-3))


def test_run_lengths_convert_windows_to_morgans():
    lengths = run_lengths_morgans(
        np.array([0, 10]), np.array([5, 20]), WINDOW_BP, RECOMBINATION_RATE
    )
    assert np.allclose(lengths, [5 * 1000 * 1.2e-8, 10 * 1000 * 1.2e-8])
