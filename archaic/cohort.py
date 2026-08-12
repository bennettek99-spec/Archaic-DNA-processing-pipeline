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


def pooled_freq_multi(panel, rows, cohorts, chunk: int = 256, log=None):
    """Pooled frequencies for MANY cohorts in a single pass over the genotypes.

    cohorts : mapping label -> 1-D array of individual column indices. Cohorts
              may overlap (a genome can belong to both a named cohort and a
              region-by-period cell); each individual's genotypes are read once
              and added to every cohort that claims it.

    Calling pooled_freq() once per cohort would re-read the .geno file once per
    cohort, which for a few dozen overlapping cohorts spanning most of the AADR
    is the dominant cost of a study. Here the union of requested columns is
    streamed in chunks and accumulated straight into per-cohort dosage/call
    accumulators, so the file is traversed once regardless of cohort count.

    Accumulators are int32: dosages are 0/1/2 and the largest conceivable cohort
    is a few thousand genomes, so the sums stay exact and the whole set of
    accumulators costs ~9 MB per cohort rather than the ~300 MB a genotype block
    would.

    Returns (freqs, counts): label -> float64 frequency array / int64 call count,
    each of length len(rows), matching pooled_freq() exactly.
    """
    rows = np.asarray(rows, dtype=np.int64)
    labels = list(cohorts)
    cols_by = {l: np.asarray(cohorts[l], dtype=np.int64) for l in labels}

    union = np.array(sorted({int(c) for l in labels for c in cols_by[l]}),
                     dtype=np.int64)
    nr = len(rows)
    s = {l: np.zeros(nr, dtype=np.int32) for l in labels}
    n = {l: np.zeros(nr, dtype=np.int32) for l in labels}
    if len(union) == 0 or nr == 0:
        return ({l: np.full(nr, np.nan) for l in labels},
                {l: n[l].astype(np.int64) for l in labels})

    # for each column in the union, which cohorts want it
    pos_in_union = {int(c): i for i, c in enumerate(union)}
    wanted = {l: np.array([pos_in_union[int(c)] for c in cols_by[l]],
                          dtype=np.int64) for l in labels}

    # Gathering a whole cohort's columns at once (Gv[:, take]) would allocate an
    # array the size of the genotype block for every cohort in every chunk, which
    # on a 1.2M-SNP panel is hundreds of MB a time and is what makes the naive
    # version thrash. Gathering in narrow column batches caps the temporary at
    # n_snp * GATHER bytes while keeping the sum vectorised.
    GATHER = 32
    for c0 in range(0, len(union), chunk):
        c1 = min(c0 + chunk, len(union))
        G = panel.pg.read(rows, union[c0:c1])            # int8 (nr, k)
        valid = G >= 0
        Gv = np.where(valid, G, 0)                       # int8, missing -> 0
        for l in labels:
            w = wanted[l]
            take = w[(w >= c0) & (w < c1)] - c0
            if len(take) == 0:
                continue
            for j0 in range(0, len(take), GATHER):
                sub = take[j0:j0 + GATHER]
                s[l] += Gv[:, sub].sum(axis=1, dtype=np.int32)
                n[l] += valid[:, sub].sum(axis=1, dtype=np.int32)
        del G, valid, Gv
        if log is not None:
            log.info(f"  pooled_freq_multi: {c1:,}/{len(union):,} individuals")

    freqs, counts = {}, {}
    for l in labels:
        nn = n[l].astype(np.int64)
        with np.errstate(invalid="ignore", divide="ignore"):
            freqs[l] = np.where(nn > 0, s[l] / (2.0 * np.maximum(nn, 1)), np.nan)
        counts[l] = nn
    return freqs, counts
