"""
Unit tests for archaic/lib_eigenstrat.py — the hand-rolled binary TGENO/GENO
packed-genotype reader. This is the most fragile, least-obvious-on-inspection
code in the pipeline (bit-packing/shift arithmetic) and previously had zero
test coverage. Uses archaic/synthetic.py's writer to build small, exactly-known
fixtures so we can assert byte-for-byte round-trip correctness without needing
the real AADR data. Run: pytest -q
"""
import gc
import os
import shutil
import sys
import tempfile
from contextlib import contextmanager

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import lib_eigenstrat as le
from archaic.synthetic import write_packed_geno, write_ind, write_snp


@contextmanager
def tmp_dir():
    """Like tempfile.TemporaryDirectory, but tolerant of Windows keeping a
    memory-mapped .geno file locked until the PackedGeno (and its np.memmap)
    is garbage-collected — a plain TemporaryDirectory raises PermissionError
    on __exit__ in that case."""
    d = tempfile.mkdtemp(prefix="archaic_test_")
    try:
        yield d
    finally:
        gc.collect()
        shutil.rmtree(d, ignore_errors=True)


def _known_dosage(n_ind=6, n_snp=13, seed=0):
    """A small dosage matrix (n_ind, n_snp) with 0/1/2 and some -1 (missing),
    deliberately including a size that is NOT a multiple of 4 (13 snps, 6
    inds) to exercise the packed format's partial-byte padding."""
    rng = np.random.default_rng(seed)
    d = rng.integers(0, 3, size=(n_ind, n_snp)).astype(np.int8)
    miss = rng.random((n_ind, n_snp)) < 0.15
    d[miss] = -1
    return d


@pytest.mark.parametrize("transposed", [True, False])
def test_packed_geno_roundtrip_exact(transposed):
    dosage = _known_dosage()
    n_ind, n_snp = dosage.shape
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=transposed)
        pg = le.PackedGeno(path, n_ind, n_snp)
        assert pg.transposed == transposed
        out = pg.read(np.arange(n_snp), np.arange(n_ind))     # (n_snp, n_ind)
        np.testing.assert_array_equal(out, dosage.T)


@pytest.mark.parametrize("transposed", [True, False])
def test_packed_geno_arbitrary_unsorted_selection(transposed):
    """read() must support an arbitrary (not necessarily sorted) subset/order
    of SNP rows and individual columns, since Panel.cols_for / snp_rows can
    select either."""
    dosage = _known_dosage(n_ind=9, n_snp=27, seed=1)
    n_ind, n_snp = dosage.shape
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=transposed)
        pg = le.PackedGeno(path, n_ind, n_snp)

        snp_rows = np.array([5, 0, 26, 3, 11])
        ind_cols = np.array([8, 2, 0, 7])
        out = pg.read(snp_rows, ind_cols)
        expected = dosage[np.ix_(ind_cols, snp_rows)].T
        np.testing.assert_array_equal(out, expected)


@pytest.mark.parametrize("transposed", [True, False])
def test_packed_geno_chunking_matches_unchunked(transposed):
    """Small ind_chunk/snp_chunk (forcing multiple loop iterations) must give
    the same result as one big chunk."""
    dosage = _known_dosage(n_ind=20, n_snp=200, seed=2)
    n_ind, n_snp = dosage.shape
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=transposed)
        pg = le.PackedGeno(path, n_ind, n_snp)

        snp_rows = np.arange(n_snp)
        ind_cols = np.arange(n_ind)
        full = pg.read(snp_rows, ind_cols, ind_chunk=1000, snp_chunk=1000)
        chunked = pg.read(snp_rows, ind_cols, ind_chunk=3, snp_chunk=17)
        np.testing.assert_array_equal(full, chunked)


def test_packed_geno_missing_code_maps_to_minus_one():
    dosage = np.array([[-1, 0, 1, 2]], dtype=np.int8)  # -1 = missing (packed as code 3)
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=True)
        pg = le.PackedGeno(path, 1, 4)
        out = pg.read(np.arange(4), np.arange(1))
        assert list(out[:, 0]) == [-1, 0, 1, 2]


def test_size_mismatch_raises_tgeno():
    dosage = _known_dosage(n_ind=4, n_snp=10, seed=3)
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=True)
        with pytest.raises(ValueError):
            le.PackedGeno(path, 4, 999)               # wrong nsnp -> size mismatch


def test_size_mismatch_raises_geno():
    dosage = _known_dosage(n_ind=4, n_snp=10, seed=3)
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        write_packed_geno(path, dosage, transposed=False)
        with pytest.raises(ValueError):
            le.PackedGeno(path, 999, 10)               # wrong nind -> size mismatch


def test_unknown_magic_raises():
    with tmp_dir() as d:
        path = os.path.join(d, "test.geno")
        with open(path, "wb") as fh:
            fh.write(b"\x00" * 200)
        with pytest.raises(ValueError):
            le.PackedGeno(path, 4, 10)


def test_read_ind_and_read_snp_roundtrip():
    ids = ["A", "B", "C"]
    pops = ["PopX", "PopY", "PopX"]
    with tmp_dir() as d:
        ind_path = os.path.join(d, "test.ind")
        snp_path = os.path.join(d, "test.snp")
        write_ind(ind_path, ids, pops)
        write_snp(snp_path, n_snp=7)

        ind = le.read_ind(ind_path)
        assert list(ind["id"]) == ids
        assert list(ind["pop"]) == pops

        snp = le.read_snp(snp_path)
        assert len(snp) == 7
        assert list(snp.index) == list(range(7))       # geno_row == file order
        assert snp["chrom"].iloc[0] == "1"
        assert snp["pos"].is_monotonic_increasing
