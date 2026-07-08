"""
Unit tests for archaic.cohort.pooled_freq — the RAM-safe streaming pooled allele
frequency. The invariant that matters: pooling must give *exactly* the same
frequency as panel.frequencies() (mean dosage/2 over non-missing calls), for any
chunk size, on a real on-disk synthetic panel. Run with: pytest -q
"""
import os, sys, shutil, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.synthetic import write_synthetic_panel
from archaic.panel import Panel
from archaic.cohort import pooled_freq


def _panel():
    tmp = tempfile.mkdtemp(prefix="archaic_cohort_")
    prefix, info = write_synthetic_panel(tmp, n_snp=6000, n_test=40,
                                         missing_rate=0.15, seed=3)
    return tmp, Panel(prefix, autosomes_only=True), info


def test_pooled_matches_panel_frequencies():
    tmp, panel, _ = _panel()
    try:
        cols = panel.cols_for(pops=["Test"])
        freq, _ = panel.frequencies({"Test": dict(pops=["Test"])})
        ref = freq["Test"]                                # nanmean(dosage)/2
        p, n = pooled_freq(panel, panel.snp_rows, cols, chunk=8)
        both = np.isfinite(ref) & np.isfinite(p)
        assert both.mean() > 0.99                          # near-full coverage
        assert np.allclose(p[both], ref[both], atol=1e-9)
        # missing everywhere -> NaN with zero count
        assert np.all((n == 0) == ~np.isfinite(p))
    finally:
        del panel; shutil.rmtree(tmp, ignore_errors=True)


def test_pooled_chunk_invariant():
    tmp, panel, _ = _panel()
    try:
        cols = panel.cols_for(pops=["Test"])
        p1, n1 = pooled_freq(panel, panel.snp_rows, cols, chunk=3)
        p2, n2 = pooled_freq(panel, panel.snp_rows, cols, chunk=64)
        assert np.array_equal(n1, n2)
        both = np.isfinite(p1) & np.isfinite(p2)
        assert np.allclose(p1[both], p2[both], atol=1e-12)
    finally:
        del panel; shutil.rmtree(tmp, ignore_errors=True)


def test_empty_cohort():
    tmp, panel, _ = _panel()
    try:
        p, n = pooled_freq(panel, panel.snp_rows, np.array([], dtype=np.int64))
        assert np.all(np.isnan(p)) and np.all(n == 0)
    finally:
        del panel; shutil.rmtree(tmp, ignore_errors=True)
