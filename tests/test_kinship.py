"""
Unit tests for archaic.kinship — READ-style pairwise mismatch classification and
duplicate/relative pruning, exercised on synthetic genotype matrices so no AADR
data is required. Run: pytest -q
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import kinship as kin


def test_haploid_alleles_resolves_homozygotes_and_hets():
    G = np.array([[0, 1, 2, -1],
                  [2, 1, 0, -1]], dtype=np.int8)
    A = kin._haploid_alleles(G, seed=0)
    assert A[0, 0] == 0 and A[0, 2] == 1          # homozygotes deterministic
    assert A[1, 0] == 1 and A[1, 2] == 0
    assert A[0, 3] == -1 and A[1, 3] == -1        # missing preserved
    assert A[0, 1] in (0, 1) and A[1, 1] in (0, 1)  # hets resolved to a haplotype


def test_haploid_alleles_is_reproducible_per_seed():
    G = np.tile(np.array([1], dtype=np.int8), (50, 4))
    a = kin._haploid_alleles(G, seed=7)
    b = kin._haploid_alleles(G, seed=7)
    assert np.array_equal(a, b)


def test_pair_p0_counts_and_mismatch():
    A = np.array([[0, 1],
                  [0, 1],
                  [1, 1],
                  [-1, -1]], dtype=np.int8)
    miss = A < 0
    i, j, p, c = kin._pair_p0(0, 1, A, miss, min_overlap=2)
    assert (i, j) == (0, 1)
    assert c == 3                                 # three non-missing overlapping sites
    assert np.isclose(p, 2.0 / 3.0)               # disagree at two of three sites


def test_pair_p0_below_overlap_is_nan():
    A = np.array([[0, 1], [-1, -1], [-1, -1]], dtype=np.int8)
    miss = A < 0
    _, _, p, c = kin._pair_p0(0, 1, A, miss, min_overlap=3)
    assert c == 0 and np.isnan(p)


def _matrix_with(close_raw_p0, n_extra=4, unrelated_raw_p0=0.9):
    """A 1+(n_extra) matrix: individual 0 is related to 1 (raw P0=close_raw_p0);
    every other pair is unrelated at raw P0 ~ unrelated_raw_p0, so the cohort
    median sits near unrelated_raw_p0 and the close pair normalises down."""
    N = 1 + n_extra
    P0 = np.full((N, N), unrelated_raw_p0)
    np.fill_diagonal(P0, np.nan)
    P0[0, 1] = P0[1, 0] = close_raw_p0
    return P0


def test_classify_identical_is_under_first_degree_threshold():
    # pair (0,1) identical at every shared site -> normalised P0 ~ 0
    P0 = np.array([[np.nan, 0.0],
                   [0.0, np.nan]])
    norm, pairs = kin.classify(P0)
    assert len(pairs) == 1
    _, _, v, deg = pairs[0]
    assert v < kin.READ_1ST and deg == "identical/duplicate"


def test_classify_unrelated_stays_above_threshold():
    # a single unrelated pair has normalised P0 of exactly 1.0 > READ_UNREL
    P0 = np.array([[np.nan, 1.0],
                   [1.0, np.nan]])
    norm, pairs = kin.classify(P0)
    assert pairs == []


def test_classify_second_degree_band():
    # normalised P0 between READ_2ND and READ_UNREL -> second-degree
    P0 = _matrix_with(0.75)                      # 0.75 / 0.9 ~ 0.833
    _, pairs = kin.classify(P0)
    degs = {deg for _, _, v, deg in pairs}
    assert "second-degree" in degs


def test_classify_first_degree_band():
    P0 = _matrix_with(0.60)                      # 0.60 / 0.9 ~ 0.667
    _, pairs = kin.classify(P0)
    degs = {deg for _, _, v, deg in pairs}
    assert "first-degree" in degs


def test_classify_empty_matrix_returns_no_pairs():
    norm, pairs = kin.classify(np.full((0, 0), np.nan))
    assert pairs == []
    assert norm.shape == (0, 0)


def test_classify_all_missing_returns_no_pairs():
    norm, pairs = kin.classify(np.full((3, 3), np.nan))
    assert pairs == []


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} tests passed")
