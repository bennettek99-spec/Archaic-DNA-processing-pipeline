"""
Unit tests for archaic.stats — validate the core f-statistics math against
analytic properties, independent of any data. Run with: pytest -q
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import stats as st


def _theta(num, den, block, nb=20):
    return st.jackknife_ratio(num, den, block, nb)["theta"]


def test_f4_antisymmetry():
    rng = np.random.default_rng(0); n = 3000
    p = {k: rng.random(n) for k in "ABCD"}
    blk = st.assign_blocks(n, 20); one = np.ones(n)
    abcd = _theta(st.f4_array(p["A"], p["B"], p["C"], p["D"]), one, blk)
    bacd = _theta(st.f4_array(p["B"], p["A"], p["C"], p["D"]), one, blk)
    abdc = _theta(st.f4_array(p["A"], p["B"], p["D"], p["C"]), one, blk)
    assert np.isclose(abcd, -bacd, atol=1e-12)
    assert np.isclose(abcd, -abdc, atol=1e-12)


def test_D_identical_pops_zero():
    rng = np.random.default_rng(1); n = 3000
    X = rng.random(n); Y = rng.random(n); Z = rng.random(n)
    blk = st.assign_blocks(n, 20)
    out = st.dstat({"X": X, "Y": Y, "Z": Z}, "X", "X", "Y", "Z", blk, 20)
    assert np.isclose(out["theta"], 0.0, atol=1e-12)


def test_f4ratio_scale_endpoints():
    """alpha = f4(A,O;X,B)/f4(A,O;Ref,B): X==Ref -> 1, X==B -> 0 (the scale anchors
    behind 'Neanderthals-as-test read ~100%, the African baseline ~0')."""
    rng = np.random.default_rng(2); n = 4000
    p = {"A": rng.random(n), "O": np.zeros(n), "B": rng.random(n), "Ref": rng.random(n)}
    blk = st.assign_blocks(n, 25)
    p1 = dict(p, X=p["Ref"].copy())
    p0 = dict(p, X=p["B"].copy())
    a1 = st.f4_ratio(p1, "A", "O", "X", "B", "Ref", blk, 25)["theta"]
    a0 = st.f4_ratio(p0, "A", "O", "X", "B", "Ref", blk, 25)["theta"]
    assert np.isclose(a1, 1.0, atol=1e-9)
    assert np.isclose(a0, 0.0, atol=1e-9)


def test_D_sign_for_constructed_sharing():
    """If X is pulled toward Y's derived alleles relative to W, D(X,W;Y,Z)>0."""
    rng = np.random.default_rng(3); n = 5000
    Z = np.zeros(n)                     # ancestral outgroup
    Y = (rng.random(n) < 0.5).astype(float)  # 'archaic' derived state
    W = rng.random(n) * 0.5            # baseline
    X = 0.7 * W + 0.3 * Y             # X shares extra derived alleles with Y
    blk = st.assign_blocks(n, 25)
    out = st.dstat({"X": X, "W": W, "Y": Y, "Z": Z}, "X", "W", "Y", "Z", blk, 25)
    assert out["theta"] > 0 and out["z"] > 2


def test_jackknife_matches_batch():
    rng = np.random.default_rng(4); n = 6000
    num = rng.standard_normal(n); den = rng.random(n) + 0.5
    nb = 30
    blk = st.assign_blocks(n, nb); starts = st.block_starts(n, nb)
    a = st.jackknife_ratio(num, den, blk, nb)
    th, se, z, used = st.batch_jackknife_ratio(num[:, None], den[:, None], starts)
    assert np.isclose(a["theta"], th[0], atol=1e-9)
    assert np.isclose(a["se"], se[0], rtol=1e-6)
    assert a["n_used"] == int(used[0])


def test_block_starts():
    starts = st.block_starts(1000, 10)
    assert len(starts) == 10 and starts[0] == 0 and np.all(np.diff(starts) > 0)
