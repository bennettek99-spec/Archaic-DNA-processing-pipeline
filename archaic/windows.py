"""
windows.py — aggregate a per-SNP signal into genomic windows.

Used by the local archaic-affinity scan (local_archaic_scan.py) to turn a
per-SNP archaic-allele frequency into a position-resolved landscape whose low
tail is candidate archaic "deserts" and whose high tail is candidate adaptive
introgression. Kept as a small, dependency-light, pure function so it can be
unit-tested without any AADR data.

Windows are physical (base-pair) and tile each chromosome from 0; a step smaller
than the width gives overlapping (smoother) windows, step == width gives disjoint
windows (the default — disjoint windows keep the empirical-p flagging closer to
independent tests). Chromosomes are processed in natural 1,2,...,22 order.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def _chrom_sort_key(c):
    """Natural order for chromosome labels: numerics 1..22 first, then the rest."""
    try:
        return (0, int(c))
    except (TypeError, ValueError):
        return (1, str(c))


def window_scan(chrom, pos, value, weight=None, win=1_000_000, step=None,
                min_snp=10):
    """Weighted-mean of `value` in sliding base-pair windows.

    chrom  : per-SNP chromosome labels (array-like of str/int).
    pos    : per-SNP base-pair positions (int).
    value  : per-SNP signal to average (float; NaN entries are ignored).
    weight : per-SNP weights (default: equal). NaN/<=0 weights are dropped.
    win    : window width in bp. step: window stride (default win == disjoint).
    min_snp: minimum finite-value SNPs for a window to be emitted.

    Returns a DataFrame sorted by (chrom, start) with columns:
      chrom, start, end, mid, n_snp, mean  (weighted mean of value over the
      window's finite SNPs), wsum (sum of weights used).
    """
    chrom = np.asarray(chrom)
    pos = np.asarray(pos, dtype=np.int64)
    value = np.asarray(value, dtype=np.float64)
    weight = np.ones_like(value) if weight is None else np.asarray(weight, float)
    step = win if step is None else int(step)

    ok = np.isfinite(value) & np.isfinite(weight) & (weight > 0)
    chrom, pos, value, weight = chrom[ok], pos[ok], value[ok], weight[ok]

    out = []
    for c in sorted(set(chrom.tolist()), key=_chrom_sort_key):
        m = chrom == c
        cp, cv, cw = pos[m], value[m], weight[m]
        order = np.argsort(cp, kind="stable")
        cp, cv, cw = cp[order], cv[order], cw[order]
        if len(cp) == 0:
            continue
        last = int(cp[-1])
        for start in range(0, last + 1, step):
            end = start + win
            lo = np.searchsorted(cp, start, side="left")
            hi = np.searchsorted(cp, end, side="left")
            if hi - lo < min_snp:
                continue
            w = cw[lo:hi]
            v = cv[lo:hi]
            wsum = float(w.sum())
            mean = float((v * w).sum() / wsum) if wsum > 0 else np.nan
            out.append((str(c), start, end, start + win // 2,
                        int(hi - lo), mean, wsum))

    cols = ["chrom", "start", "end", "mid", "n_snp", "mean", "wsum"]
    df = pd.DataFrame(out, columns=cols)
    if not df.empty:
        df = df.sort_values(
            ["chrom", "start"],
            key=lambda s: s.map(_chrom_sort_key) if s.name == "chrom" else s,
        ).reset_index(drop=True)
    return df


def robust_z(x):
    """Median/MAD standardisation (robust to the heavy tails of a desert scan).

    z = (x - median) / (1.4826 * MAD). Falls back to std if MAD == 0.
    """
    x = np.asarray(x, dtype=np.float64)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale == 0:
        scale = np.nanstd(x)
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(x)
    return (x - med) / scale


def empirical_p(x):
    """Two-sided empirical p-value of each entry within the vector's own ECDF.

    p_i = 2 * min(rank_below, rank_above) / n, clipped to (0, 1]; a value at an
    extreme of the distribution gets a small p. Purely descriptive (windows are
    not independent under overlap), used only to rank candidates.
    """
    x = np.asarray(x, dtype=np.float64)
    n = np.isfinite(x).sum()
    if n == 0:
        return np.full_like(x, np.nan)
    order = np.argsort(np.argsort(x))            # 0-based rank
    below = order + 1
    above = n - order
    p = 2.0 * np.minimum(below, above) / n
    return np.clip(p, 1.0 / n, 1.0)
