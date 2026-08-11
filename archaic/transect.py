"""
transect.py — pooled archaic statistics for time-ordered cohorts.

The rest of the pipeline estimates archaic ancestry one genome at a time. That is
the right unit for outlier detection, but it is the wrong unit for asking how a
*population's* archaic ancestry changed through time: individual ancient genomes
in a single archaeological horizon are often too low-coverage to resolve a
percentage on their own, while the horizon as a whole is well powered.

This module supplies the missing layer: bin individuals by date, pool each bin
into one allele-frequency vector (archaic.cohort.pooled_freq), and run the same
f-statistics the rest of the pipeline uses on the pooled frequencies. Because a
pooled frequency is a ratio-of-sums over individuals and every statistic here is
itself a ratio-of-sums over SNPs, the block jackknife remains valid: it resamples
genomic blocks, which is the correlated axis, and is agnostic to how the
frequency at a block was assembled.

Two facts drive the design:

  * **Effect size, not Z.** Low coverage inflates the standard error of a pooled
    statistic but does not, on its own, move the point estimate. Comparing bins
    by Z would therefore manufacture a trend wherever coverage trends. Every
    routine here returns theta and se separately so a study can compare effect
    sizes and treat power explicitly.

  * **Coverage-matched comparison.** `common_rows` restricts two cohorts to the
    SNPs where both actually have calls, which is the honest way to compare a
    30k-SNP horizon with a 400k-SNP one. `snp_floor_rows` does the same against
    a per-SNP call-count threshold.

Nothing in this module is Oceania-specific; `oceania_transect.py` is one study
built on it.
"""
from __future__ import annotations

from typing import Iterable, NamedTuple, Sequence

import numpy as np

from . import stats as st
from .cohort import pooled_freq


class TimeBin(NamedTuple):
    """A named, half-open date window in calibrated years before present.

    Membership is `hi_bp >= date_bp > lo_bp`, so adjacent bins sharing an edge do
    not double-count an individual.
    """
    label: str
    lo_bp: float
    hi_bp: float

    @property
    def mid_bp(self) -> float:
        return 0.5 * (self.lo_bp + self.hi_bp)


def assign_bins(dates, bins: Sequence[TimeBin]) -> np.ndarray:
    """Label each date with the first matching bin, or "" if none matches.

    Returns an object array of labels the same length as `dates`. Non-finite
    dates never match a bin.
    """
    dates = np.asarray(dates, dtype=np.float64)
    out = np.full(len(dates), "", dtype=object)
    for b in bins:
        hit = (out == "") & np.isfinite(dates) & (dates > b.lo_bp) & (dates <= b.hi_bp)
        out[hit] = b.label
    return out


# ------------------------------------------------------------------ SNP sets --
def common_rows(counts: Iterable[np.ndarray], min_calls: int = 1) -> np.ndarray:
    """Boolean mask of SNPs where *every* cohort has at least `min_calls` calls.

    `counts` are the per-SNP non-missing call counts returned by pooled_freq.
    Use this to compare cohorts of very different coverage on one shared SNP set,
    so a difference cannot be an ascertainment-by-coverage artifact.
    """
    counts = list(counts)
    if not counts:
        return np.zeros(0, dtype=bool)
    mask = np.ones(len(counts[0]), dtype=bool)
    for c in counts:
        mask &= np.asarray(c) >= min_calls
    return mask


def snp_floor_rows(count: np.ndarray, min_calls: int) -> np.ndarray:
    """Boolean mask of SNPs where one cohort has at least `min_calls` calls."""
    return np.asarray(count) >= min_calls


# ----------------------------------------------------------------- statistics --
def pooled_archaic_stats(panel, cols, ref, block, n_blocks: int = 50,
                         chunk: int = 256, rows=None, mask=None, log=None) -> dict:
    """Neanderthal f4-ratio and Denisovan affinity for one pooled cohort.

    cols  : individual column indices making up the cohort.
    ref   : dict of reference frequency vectors over `panel.snp_rows`, with keys
            Altai, Vindija, Denisova, Chimp, Mbuti (as built by panel.frequencies).
    block : per-SNP jackknife block ids over `panel.snp_rows`.
    mask  : optional boolean mask over `panel.snp_rows` restricting the SNP set
            (e.g. from common_rows) — applied to the cohort *and* every reference
            so all statistics share one SNP set.

    Returns alpha (Neanderthal proportion), D_Den (relative Denisovan affinity),
    and D_Nea (relative Neanderthal affinity), each with jackknife se/z/n, plus
    the pooled per-SNP call count so callers can build coverage-matched masks.
    """
    rows = panel.snp_rows if rows is None else rows
    p, n = pooled_freq(panel, rows, np.asarray(cols, dtype=np.int64),
                       chunk=chunk, log=log)

    f = {"POOL": p}
    f.update({k: np.asarray(v, dtype=np.float64) for k, v in ref.items()})
    blk = np.asarray(block)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        f = {k: np.where(mask, v, np.nan) for k, v in f.items()}

    alpha = st.f4_ratio(f, "Altai", "Chimp", "POOL", "Mbuti", "Vindija", blk, n_blocks)
    dden = st.dstat(f, "POOL", "Mbuti", "Denisova", "Chimp", blk, n_blocks)
    dnea = st.dstat(f, "POOL", "Mbuti", "Altai", "Chimp", blk, n_blocks)
    return dict(
        n_ind=int(len(cols)),
        alpha=alpha["theta"], alpha_se=alpha["se"], alpha_z=alpha["z"],
        alpha_nsnp=alpha["n_used"],
        D_Den=dden["theta"], D_Den_se=dden["se"], D_Den_z=dden["z"],
        D_Den_nsnp=dden["n_used"],
        D_Nea=dnea["theta"], D_Nea_se=dnea["se"], D_Nea_z=dnea["z"],
        D_Nea_nsnp=dnea["n_used"],
        _freq=p, _count=n,
    )


def ancestry_fraction(p_pool, ref, block, sahul: str, out: str, source: str,
                      base: str, n_blocks: int = 50, mask=None) -> dict:
    """Fraction of `source`-related ancestry in a pooled cohort, archaic-free.

        frac = f4(sahul, out; POOL, base) / f4(sahul, out; source, base)

    For a two-way mixture of a `source`-related population and a `base`-related
    population, this is the standard f4-ratio admixture proportion (Reich et al.
    2009): the numerator measures how much of the cohort's drift is shared with
    the `sahul` clade relative to `base`, and the denominator rescales that by the
    same quantity measured on an unadmixed `source` reference.

    **No archaic genome enters this statistic.** That is deliberate and load-
    bearing: it is what makes "archaic ancestry tracks source ancestry" a real
    test rather than an identity. `sahul` must be a population related to
    `source` but not a member of the cohort's mixture, and `out` an outgroup to
    both `source` and `base`.
    """
    f = {"POOL": np.asarray(p_pool, dtype=np.float64)}
    f.update({k: np.asarray(v, dtype=np.float64) for k, v in ref.items()})
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        f = {k: np.where(mask, v, np.nan) for k, v in f.items()}
    num = st.f4_array(f[sahul], f[out], f["POOL"], f[base])
    den = st.f4_array(f[sahul], f[out], f[source], f[base])
    r = st.jackknife_ratio(num, den, np.asarray(block), n_blocks)
    r["statistic"] = f"f4({sahul},{out};POOL,{base}) / f4({sahul},{out};{source},{base})"
    return r


def predict_mixture(frac, v_source, v_base):
    """Value expected in a two-way mixture: frac*v_source + (1-frac)*v_base.

    Used to turn an ancestry fraction into a *quantitative* prediction for a
    statistic measured independently, so observed-vs-predicted is falsifiable
    rather than merely correlational.
    """
    return frac * v_source + (1.0 - frac) * v_base


def predict_mixture_se(frac, frac_se, v_source, se_source, v_base, se_base):
    """Standard error of predict_mixture, propagating all three inputs.

    Three terms contribute: uncertainty in the mixture fraction (scaled by how
    far apart the endpoints are) and the sampling error of each endpoint anchor.
    The anchor terms are *common-mode* — the same anchors are used for every
    cohort — so they shift all predictions together rather than adding scatter
    between them. Including them is the conservative choice for a claim that
    observation and prediction agree.
    """
    return float(np.sqrt((abs(v_source - v_base) * frac_se) ** 2
                         + (frac * se_source) ** 2
                         + ((1.0 - frac) * se_base) ** 2))


def mixture_residual_z(observed, obs_se, predicted, pred_se=0.0):
    """Z for observed - predicted, combining both standard errors in quadrature."""
    se = float(np.hypot(obs_se, pred_se))
    if not np.isfinite(se) or se <= 0:
        return np.nan
    return float((observed - predicted) / se)


def cohort_contrast(y_first, se_first, y_last, se_last):
    """Endpoint difference last - first, with combined SE and Z.

    The comparison that works for a control series with only two horizons, where
    a regression slope is undefined. It is also the statistic a reader actually
    wants for "did this population change": a difference between two dated
    cohorts, not a fitted rate.
    """
    diff = float(y_last - y_first)
    se = float(np.hypot(se_first, se_last))
    z = diff / se if (np.isfinite(se) and se > 0) else np.nan
    return dict(diff=diff, se=se, z=float(z) if np.isfinite(z) else np.nan)


def weighted_trend(x, y, se):
    """Inverse-variance weighted least-squares slope of y on x, with its SE.

    Used for "does the statistic trend with time"; weighting by 1/se^2 stops a
    noisy low-coverage bin from dominating the fit. Returns (slope, slope_se, z).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (se > 0)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan
    x, y, w = x[ok], y[ok], 1.0 / se[ok] ** 2
    sw = w.sum()
    mx = (w * x).sum() / sw
    my = (w * y).sum() / sw
    sxx = (w * (x - mx) ** 2).sum()
    if sxx <= 0:
        return np.nan, np.nan, np.nan
    slope = (w * (x - mx) * (y - my)).sum() / sxx
    slope_se = float(np.sqrt(1.0 / sxx))
    return float(slope), slope_se, float(slope / slope_se)
