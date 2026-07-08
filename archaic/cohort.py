"""
cohort.py — RAM-safe pooled allele frequencies for a *set* of individuals.

panel.Panel.frequencies() reads every requested column in one shot, which is
perfect for a handful of reference genomes but impossible for a cohort of
thousands: pooling all 15,443 retained ancients over 1.23M SNPs would allocate
~18 GB of int8 (the machine has ~15.6 GB). This module streams the read in
column chunks and accumulates only two per-SNP vectors (sum of dosages, count of
non-missing calls), so peak memory is O(n_snp * chunk) regardless of cohort size.

The pooled "allele-1 frequency" is sum(dosage)/(2 * n_called), i.e. the mean
dosage/2 over non-missing genotypes — identical to what panel.frequencies()
computes, just built incrementally. As everywhere in this pipeline a single
consistent per-SNP allele coding is used (the f4/D statistics are
polarisation-invariant), so no derived-allele call is needed here.
"""
from __future__ import annotations
import numpy as np


def pooled_freq(panel, rows, cols, chunk: int = 256, log=None):
    """Pooled allele-1 frequency over `cols` individuals at `rows` geno-SNPs.

    rows  : 1-D int array of .geno SNP-row indices.
    cols  : 1-D int array of individual column indices (the cohort).
    Returns (p, n):
      p : float64 (len(rows),) allele-1 frequency, NaN where no individual has a
          call at that SNP.
      n : int64   (len(rows),) number of non-missing calls contributing.
    """
    rows = np.asarray(rows, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    s = np.zeros(len(rows), dtype=np.float64)     # sum of dosages (0/1/2)
    n = np.zeros(len(rows), dtype=np.int64)        # non-missing count
    if len(cols) == 0 or len(rows) == 0:
        return np.full(len(rows), np.nan), n

    for c0 in range(0, len(cols), chunk):
        G = panel.pg.read(rows, cols[c0:c0 + chunk])          # int8 (n_row, k)
        valid = G >= 0
        s += np.where(valid, G, 0).sum(axis=1, dtype=np.float64)
        n += valid.sum(axis=1, dtype=np.int64)
        if log is not None:
            log.info(f"  pooled_freq: {min(c0 + chunk, len(cols)):,}/{len(cols):,} individuals")

    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(n > 0, s / (2.0 * n), np.nan)
    return p, n
