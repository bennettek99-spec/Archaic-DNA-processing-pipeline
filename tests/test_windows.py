"""
Unit tests for archaic.windows — the genomic window aggregator used by the local
archaic-affinity scan. Pure math, no AADR data. Run with: pytest -q
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import windows as w


def test_disjoint_window_means():
    chrom = np.array(["1", "1", "1", "1"])
    pos = np.array([10, 20, 30, 150])
    val = np.array([1.0, 2.0, 3.0, 10.0])
    df = w.window_scan(chrom, pos, val, win=100, min_snp=1)
    # window [0,100): mean of 1,2,3 = 2 (n=3); window [100,200): 10 (n=1)
    d = {int(r.start): r for r in df.itertuples()}
    assert d[0].n_snp == 3 and np.isclose(d[0].mean, 2.0)
    assert d[100].n_snp == 1 and np.isclose(d[100].mean, 10.0)


def test_weighted_mean_and_nan_drop():
    chrom = np.array(["5", "5", "5"])
    pos = np.array([1, 2, 3])
    val = np.array([0.0, 1.0, np.nan])       # NaN must be dropped
    wt = np.array([3.0, 1.0, 99.0])
    df = w.window_scan(chrom, pos, val, weight=wt, win=1000, min_snp=1)
    assert len(df) == 1
    assert df.iloc[0]["n_snp"] == 2
    assert np.isclose(df.iloc[0]["mean"], 0.25)   # (0*3 + 1*1)/4


def test_min_snp_filter():
    chrom = np.array(["1", "1"])
    pos = np.array([10, 20])
    val = np.array([1.0, 2.0])
    assert len(w.window_scan(chrom, pos, val, win=100, min_snp=3)) == 0
    assert len(w.window_scan(chrom, pos, val, win=100, min_snp=2)) == 1


def test_chrom_natural_order():
    chrom = np.array(["10", "2", "1"])
    pos = np.array([5, 5, 5])
    val = np.array([1.0, 1.0, 1.0])
    df = w.window_scan(chrom, pos, val, win=100, min_snp=1)
    assert list(df["chrom"]) == ["1", "2", "10"]


def test_robust_z_and_empirical_p():
    x = np.array([-3.0, -1, 0, 1, 3.0])
    z = w.robust_z(x)
    assert np.isclose(z[2], 0.0)                  # median maps to 0
    assert z[0] < 0 < z[-1]
    p = w.empirical_p(x)
    # the two extremes carry the smallest two-sided p; the centre the largest
    assert p[0] == p.min() and p[-1] == p.min()
    assert p[2] == p.max()
