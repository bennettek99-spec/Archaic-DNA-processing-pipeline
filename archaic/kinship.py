"""
kinship.py — detect and remove relatives / library duplicates before group-level
estimates (a standard aDNA requirement: related individuals are not independent
observations and bias population frequencies).

Method: READ-style pairwise allele mismatch (Monroy Kuhn et al. 2018). For each
pair we compute P0 = the fraction of SNPs at which their (pseudo-haploidised)
alleles differ, over sites both cover; normalising P0 by the cohort median gives a
relatedness statistic with published thresholds:

    normalised P0   relationship
    > 0.90625       unrelated
    0.8125–0.90625  second degree
    0.625–0.8125    first degree
    < 0.625         identical twin / same individual / duplicate library

`prune` returns the set of individuals to keep (one per related cluster, the
higher-coverage member), so downstream group frequencies use independent samples.
"""
from __future__ import annotations
import numpy as np

READ_UNREL, READ_2ND, READ_1ST = 0.90625, 0.8125, 0.625


def _haploid_alleles(G, seed=0):
    """G (n_snp, n_ind) of 0/1/2/-1 -> alleles 0/1 (het diploid randomly resolved),
    missing kept as -1. Makes pseudo-haploid and diploid samples comparable."""
    rng = np.random.default_rng(seed)
    A = np.full(G.shape, -1, dtype=np.int8)
    A[G == 0] = 0
    A[G == 2] = 1
    het = (G == 1)
    if het.any():
        A[het] = rng.integers(0, 2, size=int(het.sum()), dtype=np.int8)
    return A


def mismatch_matrix(panel, cols, n_snp=40000, min_overlap=500, seed=0):
    """Pairwise P0 (allele-mismatch rate) for the individuals in `cols`."""
    cols = np.asarray(cols, dtype=np.int64)
    sub = panel.snp_rows[np.linspace(0, panel.n_snp - 1, min(n_snp, panel.n_snp)).astype(int)]
    G = panel.pg.read(sub, cols)
    A = _haploid_alleles(G, seed)
    miss = G < 0
    N = len(cols)
    P0 = np.full((N, N), np.nan)
    nov = np.full((N, N), 0)
    for i in range(N):
        for j in range(i + 1, N):
            m = (~miss[:, i]) & (~miss[:, j])
            c = int(m.sum())
            if c >= min_overlap:
                p = float(np.mean(A[m, i] != A[m, j]))
                P0[i, j] = P0[j, i] = p
                nov[i, j] = nov[j, i] = c
    return P0, nov


def classify(P0):
    """Normalise P0 by the cohort median and return (norm_matrix, list of related
    pairs as (i, j, norm, degree))."""
    vals = P0[np.isfinite(P0)]
    if len(vals) == 0:
        return P0, []
    med = np.median(vals)
    norm = P0 / med if med > 0 else P0
    pairs = []
    N = P0.shape[0]
    for i in range(N):
        for j in range(i + 1, N):
            v = norm[i, j]
            if not np.isfinite(v) or v > READ_UNREL:
                continue
            deg = ("identical/duplicate" if v < READ_1ST else
                   "first-degree" if v < READ_2ND else "second-degree")
            pairs.append((i, j, float(v), deg))
    return norm, pairs


def prune(panel, cols, coverage=None, n_snp=120000, seed=0, max_norm=READ_2ND):
    """Return (keep_cols, dropped, pairs): keep one individual per related cluster
    (highest coverage), drop the rest. By default only identical/first-degree pairs
    (normalised P0 < READ_2ND = 0.8125) are merged into clusters — conservative, as
    second-degree calls are noisy on low-overlap pseudo-haploid data. `coverage`
    optional per-col score to choose which member to keep (default = call rate)."""
    cols = np.asarray(cols, dtype=np.int64)
    P0, nov = mismatch_matrix(panel, cols, n_snp=n_snp, seed=seed)
    _, pairs = classify(P0)
    if coverage is None:
        sub = panel.snp_rows[np.linspace(0, panel.n_snp - 1, min(n_snp, panel.n_snp)).astype(int)]
        G = panel.pg.read(sub, cols)
        coverage = (G >= 0).mean(0)
    coverage = np.asarray(coverage, float)
    # union-find only over close kin (<= max_norm); report all pairs separately
    parent = list(range(len(cols)))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i, j, v, _ in pairs:
        if v < max_norm:
            parent[find(i)] = find(j)
    clusters = {}
    for k in range(len(cols)):
        clusters.setdefault(find(k), []).append(k)
    keep_idx, drop_idx = [], []
    for members in clusters.values():
        best = max(members, key=lambda m: coverage[m])
        keep_idx.append(best)
        drop_idx.extend(m for m in members if m != best)
    keep = cols[sorted(keep_idx)]
    dropped = cols[sorted(drop_idx)]
    return keep, dropped, pairs
