"""
Unit tests for archaic.source_contrast — the Vindija-vs-Altai source contrast.

These check the properties the study's conclusions actually lean on: that the
block-table route reproduces the ordinary D-statistic, that the paired jackknife
really does cancel shared noise (and is not merely smaller by accident), that the
symmetric normaliser is symmetric, and that a known source difference planted in
synthetic data is recovered with the right sign and size. Run with: pytest -q
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import source_contrast as sc
from archaic import stats as st
from archaic.cohort import pooled_freq, pooled_freq_multi

NB = 20


def _refs(n, rng):
    """A plausible reference set: pseudo-haploid archaics, an African baseline."""
    return dict(
        V=rng.integers(0, 2, n).astype(float),
        A=rng.integers(0, 2, n).astype(float),
        C=rng.integers(0, 2, n).astype(float),
        B=rng.random(n),
    )


def test_block_table_matches_dstat():
    """D_VA from a block table equals st.dstat on the same inputs."""
    rng = np.random.default_rng(0)
    n = 5000
    r = _refs(n, rng)
    X = rng.random(n)
    blk = st.assign_blocks(n, NB)
    t = sc.build_block_table("x", X, r["B"], r["V"], r["A"], r["C"], blk, NB)
    got = sc.cohort_stat(t, "D_VA")
    want = st.dstat({"X": X, "B": r["B"], "V": r["V"], "A": r["A"]},
                    "X", "B", "V", "A", blk, NB)
    # dstat uses every SNP finite for its own four pops; the block table also
    # requires the D_NEA terms to be finite. With no NaNs the sets coincide.
    assert np.isclose(got["theta"], want["theta"], rtol=1e-10)
    assert np.isclose(got["se"], want["se"], rtol=1e-8)


def test_neanderthal_average_is_symmetric_and_nan_safe():
    v = np.array([0.0, 1.0, np.nan, 0.5])
    a = np.array([1.0, 1.0, 0.0, np.nan])
    m = sc.neanderthal_average(v, a)
    assert np.allclose(m[:2], [0.5, 1.0])
    assert np.isnan(m[2]) and np.isnan(m[3])
    assert np.allclose(sc.neanderthal_average(v, a)[:2],
                       sc.neanderthal_average(a, v)[:2])


def test_D_NEA_ignores_which_neanderthal():
    """Swapping the two archaic genomes must not move D_NEA (it moves D_VA)."""
    rng = np.random.default_rng(1)
    n = 6000
    r = _refs(n, rng)
    X = rng.random(n)
    blk = st.assign_blocks(n, NB)
    t1 = sc.build_block_table("a", X, r["B"], r["V"], r["A"], r["C"], blk, NB)
    t2 = sc.build_block_table("b", X, r["B"], r["A"], r["V"], r["C"], blk, NB)
    assert np.isclose(t1.d_ne_theta, t2.d_ne_theta, rtol=1e-12)
    assert np.isclose(t1.d_va_theta, -t2.d_va_theta, rtol=1e-12)


def test_paired_difference_of_a_cohort_with_itself_is_zero():
    rng = np.random.default_rng(2)
    n = 4000
    r = _refs(n, rng)
    X = rng.random(n)
    blk = st.assign_blocks(n, NB)
    t = sc.build_block_table("x", X, r["B"], r["V"], r["A"], r["C"], blk, NB)
    d = sc.paired_difference(t, t, "D_VA")
    assert np.isclose(d["diff"], 0.0, atol=1e-12)
    assert np.isclose(d["se"], 0.0, atol=1e-12)
    # and the naive combination would have claimed a large error
    assert d["se_independent"] > 0


def test_paired_jackknife_cancels_shared_archaic_noise():
    """Two cohorts measured against the SAME archaics: pairing must be tighter.

    This is the property the detection limit depends on. The shared term is the
    archaic genomes' own sampling noise, which is identical for both cohorts.
    """
    rng = np.random.default_rng(3)
    n = 20000
    r = _refs(n, rng)
    blk = st.assign_blocks(n, NB)
    base = rng.random(n)
    X1 = np.clip(base + 0.02 * rng.standard_normal(n), 0, 1)
    X2 = np.clip(base + 0.02 * rng.standard_normal(n), 0, 1)
    t1 = sc.build_block_table("1", X1, r["B"], r["V"], r["A"], r["C"], blk, NB)
    t2 = sc.build_block_table("2", X2, r["B"], r["V"], r["A"], r["C"], blk, NB)
    d = sc.paired_difference(t1, t2, "D_VA")
    assert d["se"] < d["se_independent"]


def _neanderthal_sources(r, lean=0.30):
    """A Vindija-leaning and an Altai-leaning source with EQUAL Neanderthal-ness.

    Both are built as the symmetric Neanderthal average plus or minus a multiple
    of the difference between the two genomes. Because (V-A) is non-zero only
    where the average is 0.5, the two sources have the same D_NEA and opposite
    D_VA — which is exactly the situation the study has to be able to tell apart,
    and exactly what a source built only from (V-A) would fail to reproduce (it
    would have no Neanderthal affinity at all, leaving the fit degenerate).
    """
    nea = 0.5 * (r["V"] + r["A"])
    return nea + lean * (r["V"] - r["A"]), nea - lean * (r["V"] - r["A"])


def test_proportional_fit_recovers_planted_slope():
    """Cohorts built as a*source + (1-a)*African must lie on one line of slope k."""
    rng = np.random.default_rng(4)
    n = 40000
    r = _refs(n, rng)
    blk = st.assign_blocks(n, NB)
    source, _ = _neanderthal_sources(r)
    afr = np.clip(r["B"] + 0.01 * rng.standard_normal(n), 0, 1)
    tables = {}
    for a in (0.01, 0.02, 0.03, 0.05, 0.10):
        X = a * source + (1 - a) * afr
        tables[f"a{a}"] = sc.build_block_table(f"a{a}", X, r["B"], r["V"],
                                               r["A"], r["C"], blk, NB)
    labels = list(tables)
    k, resid = sc.proportional_fit(tables, labels)
    # one shared source -> all residuals small relative to the spread of D_VA
    spread = np.ptp([tables[l].d_va_theta for l in labels])
    assert np.max(np.abs(resid)) < 0.1 * spread
    assert np.isfinite(k) and k > 0


def test_fit_flags_a_cohort_with_a_different_source():
    """Plant a genuinely different Neanderthal source; the residual must fire."""
    rng = np.random.default_rng(5)
    n = 40000
    r = _refs(n, rng)
    blk = st.assign_blocks(n, NB)
    vindija_like, altai_like = _neanderthal_sources(r)
    afr = r["B"]
    tables = {}
    for a in (0.02, 0.03, 0.05, 0.08):
        tables[f"n{a}"] = sc.build_block_table(
            f"n{a}", a * vindija_like + (1 - a) * afr, r["B"], r["V"], r["A"],
            r["C"], blk, NB)
    tables["odd"] = sc.build_block_table(
        "odd", 0.05 * altai_like + 0.95 * afr, r["B"], r["V"], r["A"], r["C"],
        blk, NB)
    fit, rows = sc.fit_and_residuals(tables, [l for l in tables if l != "odd"],
                                     list(tables), NB)
    by = {r_["label"]: r_ for r_ in rows}
    assert fit["k"] > 0
    assert abs(by["odd"]["residual_z"]) > 3
    for lab in tables:
        if lab != "odd":
            assert abs(by[lab]["residual_z"]) < abs(by["odd"]["residual_z"])


def test_jackknife_over_blocks_matches_cohort_stat():
    rng = np.random.default_rng(6)
    n = 8000
    r = _refs(n, rng)
    blk = st.assign_blocks(n, NB)
    t = sc.build_block_table("x", rng.random(n), r["B"], r["V"], r["A"], r["C"],
                             blk, NB)
    theta, se = sc.jackknife_over_blocks(lambda tb: tb["x"].d_va_theta,
                                         {"x": t}, NB)
    direct = sc.cohort_stat(t, "D_VA")
    assert np.isclose(theta[0], direct["theta"], rtol=1e-12)
    assert np.isclose(se[0], direct["se"], rtol=1e-8)


def test_detection_limit_takes_the_larger_floor():
    # systematic spread dominates
    d = sc.detection_limit([0.01, -0.01, 0.012, -0.011], [0.001] * 4, 0.1)
    assert d["systematic_floor"] > d["statistical_floor"]
    assert np.isclose(d["limit"], d["systematic_floor"])
    assert np.isclose(d["limit_fraction_of_signal"], d["limit"] / 0.1)
    # statistical noise dominates
    d2 = sc.detection_limit([0.0001, -0.0001], [0.05, 0.05], 0.1)
    assert np.isclose(d2["limit"], d2["statistical_floor"])


def test_detection_limit_is_scaled_by_a_D_VA_magnitude():
    """The limit must be expressed against a D_VA magnitude, not a slope.

    Guards a real units bug: dividing a D_VA difference by the fitted slope k
    (which is D_VA per unit D_NEA, ~2.2) instead of by a typical cohort's D_VA
    (~0.07) understates the limit by a factor of ~30.
    """
    d = sc.detection_limit([0.001, -0.001], [0.005, 0.005], 0.0737)
    assert 0.1 < d["limit_fraction_of_signal"] < 0.2      # ~13%, not ~0.4%
    assert d["best_limit"] <= d["limit"]
    assert d["best_fraction_of_signal"] <= d["limit_fraction_of_signal"]


def test_detection_limit_best_case_uses_the_tightest_comparison():
    d = sc.detection_limit([0.0], [0.01, 0.005, 0.001], 0.07)
    assert np.isclose(d["best_case_floor"], 2 * 0.001)
    assert np.isclose(d["statistical_floor"], 2 * 0.005)


def test_technical_covariates_detects_a_planted_trend():
    labels = [f"c{i}" for i in range(12)]
    age = np.arange(12, dtype=float) * 1000
    values = 2.3 - 0.00004 * age              # R falls with age
    noise = np.random.default_rng(0).standard_normal(12)
    out = sc.technical_covariates(labels, values,
                                  {"age": age, "unrelated": noise})
    assert out["age"]["rho"] < -0.95 and out["age"]["p"] < 0.01
    assert abs(out["unrelated"]["rho"]) < 0.9
    # too few points to be meaningful -> reported as nan rather than guessed at
    short = sc.technical_covariates(["a", "b"], [1.0, 2.0], {"age": [1.0, 2.0]})
    assert np.isnan(short["age"]["rho"])


# ------------------------------------------------------- pooled_freq_multi ----
class _FakePG:
    def __init__(self, G):
        self.G = G

    def read(self, rows, cols):
        return self.G[np.ix_(np.asarray(rows), np.asarray(cols))]


class _FakePanel:
    def __init__(self, G):
        self.pg = _FakePG(G)


def test_pooled_freq_multi_matches_pooled_freq():
    rng = np.random.default_rng(7)
    n_snp, n_ind = 500, 40
    G = rng.integers(0, 3, (n_snp, n_ind)).astype(np.int8)
    G[rng.random((n_snp, n_ind)) < 0.2] = -1        # missingness
    panel = _FakePanel(G)
    rows = np.arange(n_snp)
    cohorts = {"a": np.arange(0, 15), "b": np.arange(10, 30),   # overlapping
               "c": np.array([31, 33, 35])}
    fm, cm = pooled_freq_multi(panel, rows, cohorts, chunk=7)
    for lab, cols in cohorts.items():
        p, c = pooled_freq(panel, rows, cols, chunk=7)
        assert np.allclose(fm[lab], p, equal_nan=True)
        assert np.array_equal(cm[lab], c)


def test_pooled_freq_multi_handles_empty_and_all_missing():
    G = np.full((50, 5), -1, dtype=np.int8)
    panel = _FakePanel(G)
    fm, cm = pooled_freq_multi(panel, np.arange(50), {"x": np.arange(5)})
    assert np.all(np.isnan(fm["x"]))
    assert np.all(cm["x"] == 0)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
