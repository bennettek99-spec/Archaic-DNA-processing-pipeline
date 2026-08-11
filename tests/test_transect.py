"""
Unit tests for archaic.transect — the pooled time-transect layer.

These run without AADR data. The load-bearing checks are the algebraic ones:
an f4-ratio ancestry fraction must recover an exactly-constructed mixture
proportion, and a mixture prediction must be exact for a linear statistic.
Run with: pytest -q
"""
import os, shutil, sys, tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import stats as st
from archaic import transect as tr
from archaic.panel import Panel
from archaic.synthetic import write_synthetic_panel


def test_assign_bins_half_open_and_first_match():
    bins = [tr.TimeBin("late", 0, 1000), tr.TimeBin("early", 1000, 3000)]
    lab = tr.assign_bins([500, 1000, 1001, 3000, 3001, np.nan], bins)
    # hi_bp >= date > lo_bp, so 1000 belongs to "late" and 1001 to "early"
    assert list(lab[:5]) == ["late", "late", "early", "early", ""]
    assert lab[5] == ""                                  # non-finite never matches


def test_assign_bins_overlapping_bins_take_the_first():
    bins = [tr.TimeBin("a", 0, 2000), tr.TimeBin("b", 1000, 3000)]
    lab = tr.assign_bins([1500], bins)
    assert lab[0] == "a"


def test_timebin_mid():
    assert tr.TimeBin("x", 100, 300).mid_bp == 200


def test_common_rows_and_floor():
    c1 = np.array([0, 5, 9, 2])
    c2 = np.array([3, 0, 7, 4])
    assert list(tr.common_rows([c1, c2])) == [False, False, True, True]
    assert list(tr.common_rows([c1, c2], min_calls=4)) == [False, False, True, False]
    assert list(tr.snp_floor_rows(c1, 5)) == [False, True, True, False]
    assert len(tr.common_rows([])) == 0


def test_ancestry_fraction_recovers_exact_mixture():
    """POOL built as an exact 0.35/0.65 frequency mixture must read back 0.35.

    f4 is linear in the POOL frequency, so this is an algebraic identity, not a
    statistical approximation — any deviation is a bug.
    """
    rng = np.random.default_rng(11)
    n = 4000
    ref = {k: rng.uniform(0.05, 0.95, n) for k in
           ("Australian", "Mbuti", "Papuan", "Ami")}
    frac = 0.35
    pool = frac * ref["Papuan"] + (1 - frac) * ref["Ami"]
    block = st.assign_blocks(n, 20)
    r = tr.ancestry_fraction(pool, ref, block, sahul="Australian", out="Mbuti",
                             source="Papuan", base="Ami", n_blocks=20)
    assert abs(r["theta"] - frac) < 1e-9
    assert r["n_used"] == n


def test_ancestry_fraction_zero_and_one_endpoints():
    rng = np.random.default_rng(12)
    n = 2000
    ref = {k: rng.uniform(0.05, 0.95, n) for k in
           ("Australian", "Mbuti", "Papuan", "Ami")}
    block = st.assign_blocks(n, 20)
    for pool, want in [(ref["Ami"], 0.0), (ref["Papuan"], 1.0)]:
        r = tr.ancestry_fraction(pool, ref, block, "Australian", "Mbuti",
                                 "Papuan", "Ami", n_blocks=20)
        assert abs(r["theta"] - want) < 1e-9


def test_ancestry_fraction_mask_restricts_snp_set():
    rng = np.random.default_rng(13)
    n = 1000
    ref = {k: rng.uniform(0.05, 0.95, n) for k in
           ("Australian", "Mbuti", "Papuan", "Ami")}
    pool = 0.5 * ref["Papuan"] + 0.5 * ref["Ami"]
    block = st.assign_blocks(n, 10)
    mask = np.zeros(n, dtype=bool)
    mask[:400] = True
    r = tr.ancestry_fraction(pool, ref, block, "Australian", "Mbuti", "Papuan",
                             "Ami", n_blocks=10, mask=mask)
    assert r["n_used"] == 400
    assert abs(r["theta"] - 0.5) < 1e-9        # still exact on the subset


def test_predict_mixture_and_residual_z():
    assert abs(tr.predict_mixture(0.25, 0.032, -0.004) - 0.005) < 1e-12
    # observed sits 2 combined-SEs above prediction
    z = tr.mixture_residual_z(0.010, 0.003, 0.004, 0.004)
    assert abs(z - 0.006 / np.hypot(0.003, 0.004)) < 1e-12
    assert np.isnan(tr.mixture_residual_z(0.01, 0.0, 0.01, 0.0))


def test_predict_mixture_se_terms():
    # pure fraction uncertainty: endpoints known exactly
    se = tr.predict_mixture_se(0.5, 0.02, 0.03, 0.0, -0.01, 0.0)
    assert abs(se - 0.04 * 0.02) < 1e-12
    # pure anchor uncertainty at the endpoints of the mixture
    assert abs(tr.predict_mixture_se(1.0, 0.0, 0.03, 0.005, -0.01, 0.004)
               - 0.005) < 1e-12
    assert abs(tr.predict_mixture_se(0.0, 0.0, 0.03, 0.005, -0.01, 0.004)
               - 0.004) < 1e-12
    # terms combine in quadrature and never shrink below any single term
    both = tr.predict_mixture_se(0.5, 0.02, 0.03, 0.005, -0.01, 0.004)
    assert both > max(0.04 * 0.02, 0.5 * 0.005, 0.5 * 0.004)


def test_cohort_contrast():
    c = tr.cohort_contrast(-0.008, 0.0075, 0.0258, 0.0054)
    assert abs(c["diff"] - 0.0338) < 1e-9
    assert abs(c["se"] - np.hypot(0.0075, 0.0054)) < 1e-12
    assert c["z"] > 3
    # a flat control reads ~0 with a finite SE
    flat = tr.cohort_contrast(-0.0087, 0.006, -0.0052, 0.004)
    assert abs(flat["z"]) < 1
    assert np.isnan(tr.cohort_contrast(0.1, 0.0, 0.2, 0.0)["z"])


def test_weighted_trend_exact_on_noiseless_line():
    x = np.array([100.0, 500.0, 1200.0, 2500.0])
    y = 0.003 + 2e-6 * x
    se = np.full(4, 0.001)
    slope, slope_se, z = tr.weighted_trend(x, y, se)
    assert abs(slope - 2e-6) < 1e-15
    assert slope_se > 0 and abs(z) > 1


def test_weighted_trend_downweights_noisy_point():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 2.0, 30.0])          # last point wildly off
    tight = np.array([0.01, 0.01, 0.01, 10.0])   # ...and wildly uncertain
    slope, _, _ = tr.weighted_trend(x, y, tight)
    assert abs(slope - 1.0) < 0.05               # noisy point barely moves the fit
    loose = np.full(4, 0.01)
    slope2, _, _ = tr.weighted_trend(x, y, loose)
    assert slope2 > 5                            # equal weights: it dominates


def test_weighted_trend_degenerate_inputs():
    assert np.isnan(tr.weighted_trend([1, 2], [1, 2], [1, 1])[0])          # <3 pts
    assert np.isnan(tr.weighted_trend([1, 1, 1], [1, 2, 3], [1, 1, 1])[0])  # no x spread


def test_pooled_archaic_stats_on_synthetic_panel():
    """End-to-end through the real packed reader: a cohort with known archaic
    ancestry reads positive, and a mask genuinely restricts the SNP set."""
    tmp = tempfile.mkdtemp(prefix="archaic_transect_")
    try:
        prefix, _ = write_synthetic_panel(tmp, n_snp=6000, n_test=30,
                                          alpha_true=0.05, missing_rate=0.1, seed=7)
        panel = Panel(prefix, autosomes_only=True)
        ref, _ = panel.frequencies({k: dict(pops=[k]) for k in
                                    ("Altai", "Vindija", "Denisova", "Chimp", "Mbuti")})
        block = st.assign_blocks(panel.n_snp, 20)
        cols = panel.cols_for(pops=["Test"])
        out = tr.pooled_archaic_stats(panel, cols, ref, block, n_blocks=20, chunk=8)
        assert out["n_ind"] == len(cols)
        assert out["alpha"] > 0.01                      # recovers the seeded signal
        assert out["alpha_nsnp"] > 1000
        assert len(out["_count"]) == panel.n_snp

        mask = np.zeros(panel.n_snp, dtype=bool)
        mask[:1500] = True
        sub = tr.pooled_archaic_stats(panel, cols, ref, block, n_blocks=20,
                                      chunk=8, mask=mask)
        assert sub["alpha_nsnp"] < out["alpha_nsnp"]
        assert sub["alpha_nsnp"] <= 1500
    finally:
        panel = None
        shutil.rmtree(tmp, ignore_errors=True)
