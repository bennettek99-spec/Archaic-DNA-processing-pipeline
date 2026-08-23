import numpy as np
import pandas as pd

from archaic.admixture_dating import (
    ExponentialFit,
    calendar_interval,
    combine_aggregates,
    covariance_curve,
    derived_dosage,
    fit_exponential_curve,
    pair_covariance_aggregates,
)


def test_derived_dosage_uses_chimp_polarization_and_preserves_missing():
    genotype = np.array([0, 1, 2, -1, 0, 1, 2, -1], dtype=np.int8)
    chimp = np.array([0, 0, 0, 0, 2, 2, 2, 2], dtype=np.int8)
    observed = derived_dosage(genotype, chimp)
    assert observed.tolist() == [0, 1, 2, -1, 2, 1, 0, -1]


def test_pair_covariance_matches_direct_bin_calculation():
    chrom = np.array(["1", "1", "1", "1", "2", "2", "2"])
    # Positions in Morgans: chr1 has three pairs in [0.02, 1.02) cM.
    gpos = np.array([0.0, 0.003, 0.007, 0.015, 0.0, 0.004, 0.02])
    genotype = np.array([0, 2, 2, 0, 0, 2, 2], dtype=float)
    centers, per_chrom = pair_covariance_aggregates(
        chrom, gpos, genotype, min_cm=0.02, max_cm=1.02, bin_cm=1.0)
    curve = covariance_curve(centers, combine_aggregates(per_chrom))

    # Included pairs are four on chr1 and one on chr2.
    x = np.array([0.0, 0.0, 2.0, 2.0, 0.0])
    y = np.array([2.0, 2.0, 2.0, 0.0, 2.0])
    expected = np.sum((x - x.mean()) * (y - y.mean())) / (len(x) - 1)
    assert curve.loc[0, "n_pairs"] == 5
    assert np.isclose(curve.loc[0, "covariance"], expected)


def test_exponential_fit_recovers_known_generations():
    distance_cm = np.arange(0.0225, 1.0, 0.005)
    true_generations = 120.0
    covariance = 0.2 * np.exp(-true_generations * distance_cm / 100.0) + 0.01
    curve = pd.DataFrame({
        "distance_cm": distance_cm,
        "covariance": covariance,
        "n_pairs": np.full(len(distance_cm), 1000),
    })
    fit = fit_exponential_curve(curve, min_pairs=50)
    assert fit.converged
    assert abs(fit.generations - true_generations) < 0.1
    assert fit.r_squared > 0.999999


def test_calendar_interval_is_reproducible_and_contains_point():
    first = calendar_interval(80, 8, 45_000, 500, seed=17, draws=20_000)
    second = calendar_interval(80, 8, 45_000, 500, seed=17, draws=20_000)
    assert first == second
    point, low, high = first
    assert point == 47_320
    assert low < point < high


def test_pair_covariance_rejects_bad_distance_bounds():
    chrom = np.array(["1", "1"])
    gpos = np.array([0.0, 0.005])
    g = np.array([0, 2], dtype=float)
    for min_cm, max_cm in ((0.5, 0.1), (-0.1, 1.0)):
        try:
            pair_covariance_aggregates(chrom, gpos, g,
                                       min_cm=min_cm, max_cm=max_cm, bin_cm=1.0)
            raised = False
        except ValueError:
            raised = True
        assert raised


def test_pair_covariance_rejects_nonpositive_bin():
    chrom = np.array(["1", "1"])
    gpos = np.array([0.0, 0.005])
    g = np.array([0, 2], dtype=float)
    try:
        pair_covariance_aggregates(chrom, gpos, g,
                                   min_cm=0.0, max_cm=1.0, bin_cm=0.0)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_combine_aggregates_excluding_every_chromosome_raises():
    agg = np.zeros((4, 2))
    per_chrom = {"1": agg.copy(), "2": agg.copy()}
    try:
        combine_aggregates(per_chrom, exclude="1")
        # excluding one still leaves "2", so this must NOT raise
    except ValueError:
        assert False, "excluding one of two chromosomes should be valid"
    try:
        combine_aggregates({"1": agg.copy()}, exclude="1")
        assert False, "excluding the only chromosome should raise"
    except ValueError:
        pass


def test_covariance_curve_marks_single_pair_as_nan():
    centers = np.array([0.5])
    agg = np.array([[1.0], [1.0], [1.0], [1.0]])   # sum_x, sum_y, sum_xy, n=1
    curve = covariance_curve(centers, agg)
    assert curve.loc[0, "n_pairs"] == 1
    assert np.isnan(curve.loc[0, "covariance"])     # needs n > 1


def test_fit_exponential_curve_rejects_too_few_bins():
    dist = np.arange(0.0225, 0.06, 0.005)
    curve = pd.DataFrame({"distance_cm": dist,
                          "covariance": np.zeros(len(dist)),
                          "n_pairs": np.full(len(dist), 100)})
    try:
        fit_exponential_curve(curve, min_pairs=50)
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_fit_exponential_curve_weighted_recovers_same_rate():
    distance_cm = np.arange(0.0225, 1.0, 0.005)
    true_generations = 80.0
    covariance = 0.3 * np.exp(-true_generations * distance_cm / 100.0) + 0.005
    curve = pd.DataFrame({"distance_cm": distance_cm,
                          "covariance": covariance,
                          "n_pairs": np.full(len(distance_cm), 1000)})
    fit = fit_exponential_curve(curve, min_pairs=50, weighted=True)
    assert fit.converged
    assert abs(fit.generations - true_generations) < 1.0
