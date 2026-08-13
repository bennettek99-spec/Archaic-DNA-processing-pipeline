"""
source_contrast.py — *which* Neanderthal population is a cohort descended from?

The rest of the pipeline asks **how much** Neanderthal ancestry a genome carries.
This module asks the orthogonal question: **which** Neanderthal the ancestry came
from. The AADR 1240K release contains two high-coverage Neanderthal genomes from
opposite ends of the range and of the timeline — Altai (Denisova Cave, Siberia,
~120 kya) and Vindija 33.19 (Croatia, ~50 kya) — and Vindija is known to sit
closer to the population that actually introgressed into modern humans (Prufer et
al. 2017, Science 358:655). The contrast between them is therefore a (weak) probe
of the introgressing source itself.

THE STATISTIC

    D_VA = D(X, Yoruba; Vindija, Altai)

is positive when X shares more alleles with Vindija than with Altai, relative to
an African baseline. Three things about it drive every design decision here.

  1. **It scales with how much Neanderthal ancestry X has.** To first order a
     population carrying fraction a of Neanderthal ancestry from source S has
     D_VA(X) ~ a * D_VA(S), so any population with more Neanderthal ancestry is
     *expected* to show a proportionally larger D_VA under a single shared
     source. Reading that as "a different Neanderthal" is the central trap of
     this analysis. It bites hardest on the East Asian question, where published
     estimates put the East-Asian excess over Europeans somewhere around 8-20%
     — a *relative* excess: absolute Neanderthal fractions are ~1.7-2.6% in both
     groups and the disputed gap is a few tenths of a percentage point, not
     8-20 points of genome. Note also that this repository's own f4-ratio on the
     1240K panel puts the excess near 2% relative (2.14% against 2.11%), so the
     size of the confound is itself panel-dependent and should be measured
     rather than assumed. The fix is to divide it out:

         R = D_VA / D_NEA      where     D_NEA = D(X, Yoruba; NeaAvg, Chimp)

     with NeaAvg = (p_Vindija + p_Altai)/2. D_NEA is a relative Neanderthal
     affinity that is **exactly symmetric** in the two archaic genomes, so it
     measures the *quantity* of Neanderthal ancestry without responding to which
     Neanderthal it resembles. R is then the Vindija-over-Altai preference *per
     unit of Neanderthal ancestry* — an estimate of D_VA for the source
     population itself, and the quantity the "same source?" question is about.

  2. **Its absolute value is not interpretable, but differences are.** Vindija is
     released as a pseudo-haploid `.SG` sample and Altai as a diploid `.DG` one;
     Vindija is called at less than half as many sites; and Yoruba is not
     archaic-free (Chen et al. 2020, Cell 180:677). Every one of those offsets is
     *common-mode*: the same two archaic genomes and the same African baseline
     enter every cohort's statistic, so the offsets cancel in a cohort-vs-cohort
     difference and do not cancel in a single cohort's value. Only differences
     are reported as findings.

  3. **The dominant noise is shared between cohorts.** Only ~19k of the 528k
     sites where both archaics are called actually distinguish them, so the
     archaic genomes — not the cohorts — are the limiting sample. That noise is
     identical for every cohort, which means a *paired* block jackknife of the
     difference (resampling the same genomic blocks for both cohorts and
     differencing within each replicate) is far tighter than combining two
     independent standard errors in quadrature; empirically 1.3-2.9x. Every
     comparison in this module is paired for that reason.

DESIGN

Each cohort is reduced to a small **block table**: per-jackknife-block sums of the
four per-SNP quantities (D_VA numerator/denominator, D_NEA numerator/denominator).
That is 4 x n_blocks floats per cohort regardless of cohort size, so once the
genotype pass is done every downstream statistic — point estimates, standard
errors, all pairwise paired differences, the proportionality regression, and
leave-one-block-out versions of all of them — is computed from the tables alone.
`jackknife_over_blocks` generalises this: hand it any function of the full set of
tables and it returns that function's value plus a block-jackknife standard error,
with all shared noise correctly cancelled.

Nothing here is cohort-specific; `neanderthal_source.py` is one study built on it.
"""
from __future__ import annotations

from typing import Callable, Mapping, NamedTuple, Sequence

import numpy as np

from . import stats as st


# ------------------------------------------------------------------ helpers ---
def finite_mask(freqs: Mapping[str, np.ndarray], names: Sequence[str]) -> np.ndarray:
    """Boolean mask of SNPs where every named population has a frequency."""
    names = list(names)
    m = np.isfinite(np.asarray(freqs[names[0]], dtype=np.float64))
    for nm in names[1:]:
        m &= np.isfinite(np.asarray(freqs[nm], dtype=np.float64))
    return m


def neanderthal_average(p_vindija, p_altai) -> np.ndarray:
    """Symmetric Neanderthal frequency (p_V + p_A)/2; NaN unless both are called.

    Averaging the two archaic genomes is what makes the normaliser blind to which
    of them a cohort resembles: swapping Vindija for Altai leaves NeaAvg
    unchanged, so D_NEA responds to the amount of Neanderthal ancestry and not to
    its source. Requiring both calls also forces D_NEA onto exactly the SNP set
    D_VA can use, so the ratio R is a ratio over one shared set of sites.
    """
    v = np.asarray(p_vindija, dtype=np.float64)
    a = np.asarray(p_altai, dtype=np.float64)
    out = 0.5 * (v + a)
    out[~(np.isfinite(v) & np.isfinite(a))] = np.nan
    return out


# --------------------------------------------------------------- block table --
class BlockTable(NamedTuple):
    """Per-block sums of the four per-SNP terms behind D_VA and D_NEA.

    All downstream statistics are functions of these arrays, so a cohort of
    3,000 genomes and a cohort of one cost the same from here on.
    """
    label: str
    n_va: np.ndarray      # (n_blocks,) D_VA numerator sums
    d_va: np.ndarray      # (n_blocks,) D_VA denominator sums
    n_ne: np.ndarray      # (n_blocks,) D_NEA numerator sums
    d_ne: np.ndarray      # (n_blocks,) D_NEA denominator sums
    n_snp: int            # SNPs contributing
    n_ind: int = 0

    @property
    def d_va_theta(self) -> float:
        return float(self.n_va.sum() / self.d_va.sum())

    @property
    def d_ne_theta(self) -> float:
        return float(self.n_ne.sum() / self.d_ne.sum())

    @property
    def ratio_theta(self) -> float:
        """R = D_VA / D_NEA — Vindija preference per unit Neanderthal ancestry."""
        return self.d_va_theta / self.d_ne_theta


def build_block_table(label, p_x, p_base, p_vindija, p_altai, p_chimp, block,
                      n_blocks: int = 50, mask=None, n_ind: int = 0) -> BlockTable:
    """Reduce one cohort to its block table.

    p_x is the cohort's pooled allele frequency; the rest are the fixed
    references. `mask` optionally restricts the SNP set (e.g. a coverage-matched
    shared set). A SNP contributes only where *all* four terms are finite, so
    D_VA and D_NEA — and therefore their ratio — share one SNP set exactly.
    """
    nea = neanderthal_average(p_vindija, p_altai)
    f = {"X": np.asarray(p_x, dtype=np.float64),
         "B": np.asarray(p_base, dtype=np.float64),
         "V": np.asarray(p_vindija, dtype=np.float64),
         "A": np.asarray(p_altai, dtype=np.float64),
         "N": nea,
         "C": np.asarray(p_chimp, dtype=np.float64)}

    n_va = st.d_numerator(f["X"], f["B"], f["V"], f["A"])
    d_va = st.d_denominator(f["X"], f["B"], f["V"], f["A"])
    n_ne = st.d_numerator(f["X"], f["B"], f["N"], f["C"])
    d_ne = st.d_denominator(f["X"], f["B"], f["N"], f["C"])

    ok = (np.isfinite(n_va) & np.isfinite(d_va)
          & np.isfinite(n_ne) & np.isfinite(d_ne))
    if mask is not None:
        ok &= np.asarray(mask, dtype=bool)

    block = np.asarray(block)
    sums = []
    for arr in (n_va, d_va, n_ne, d_ne):
        sums.append(np.bincount(block, weights=np.where(ok, arr, 0.0),
                                minlength=n_blocks))
    return BlockTable(label, sums[0], sums[1], sums[2], sums[3],
                      int(ok.sum()), int(n_ind))


# ------------------------------------------------- generic block jackknife ----
def jackknife_over_blocks(fn: Callable[[Mapping[str, BlockTable]], np.ndarray],
                          tables: Mapping[str, BlockTable],
                          n_blocks: int = 50):
    """Block-jackknife any statistic computed from a set of cohort block tables.

    `fn` receives a dict of BlockTables and returns a float or 1-D array. It is
    called once on the full data and once per deleted block, with every cohort's
    table having that block removed *simultaneously* — which is precisely what
    makes the shared archaic-genome noise cancel in differences between cohorts.

    Returns (theta, se) with the usual (g-1)/g delete-one jackknife variance.
    Blocks where any cohort has no data are skipped.
    """
    theta = np.atleast_1d(np.asarray(fn(tables), dtype=np.float64))
    labels = list(tables)
    usable = np.ones(n_blocks, dtype=bool)
    for lab in labels:
        t = tables[lab]
        usable &= (t.d_va != 0) & (t.d_ne != 0)

    loo = []
    for b in range(n_blocks):
        if not usable[b]:
            continue
        reduced = {}
        for lab in labels:
            t = tables[lab]
            keep = np.ones(n_blocks, dtype=bool)
            keep[b] = False
            reduced[lab] = t._replace(n_va=t.n_va[keep], d_va=t.d_va[keep],
                                      n_ne=t.n_ne[keep], d_ne=t.d_ne[keep])
        loo.append(np.atleast_1d(np.asarray(fn(reduced), dtype=np.float64)))

    g = len(loo)
    if g < 2:
        return theta, np.full_like(theta, np.nan)
    L = np.vstack(loo)
    se = np.sqrt((g - 1) / g * np.sum((L - L.mean(axis=0)) ** 2, axis=0))
    return theta, se


def _loo_values(t: BlockTable, stat: str) -> np.ndarray:
    """Leave-one-block-out values of one cohort statistic, one per block."""
    tn_va, td_va = t.n_va.sum(), t.d_va.sum()
    tn_ne, td_ne = t.n_ne.sum(), t.d_ne.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        va = (tn_va - t.n_va) / (td_va - t.d_va)
        ne = (tn_ne - t.n_ne) / (td_ne - t.d_ne)
    if stat == "D_VA":
        out = va
    elif stat == "D_NEA":
        out = ne
    elif stat == "R":
        out = va / ne
    else:
        raise ValueError(f"unknown statistic {stat!r}")
    return out


def cohort_stat(t: BlockTable, stat: str = "D_VA") -> dict:
    """Point estimate and block-jackknife SE for one cohort statistic."""
    theta = {"D_VA": t.d_va_theta, "D_NEA": t.d_ne_theta,
             "R": t.ratio_theta}[stat]
    loo = _loo_values(t, stat)
    loo = loo[np.isfinite(loo)]
    g = len(loo)
    se = float(np.sqrt((g - 1) / g * np.sum((loo - loo.mean()) ** 2))) if g > 1 else np.nan
    z = theta / se if (np.isfinite(se) and se > 0) else np.nan
    return dict(theta=float(theta), se=se, z=float(z), n_blocks_used=g,
                n_snp=t.n_snp, n_ind=t.n_ind)


def paired_difference(a: BlockTable, b: BlockTable, stat: str = "D_VA") -> dict:
    """Difference a - b with a PAIRED block jackknife.

    The same genomic blocks are deleted from both cohorts in each replicate and
    the statistic differenced *within* the replicate, so noise common to the two
    — above all the sampling noise of the two archaic genomes, which is the
    dominant term here — cancels instead of adding. Compare `se_independent`,
    the naive quadrature combination, to see how much is common-mode.
    """
    th_a = {"D_VA": a.d_va_theta, "D_NEA": a.d_ne_theta, "R": a.ratio_theta}[stat]
    th_b = {"D_VA": b.d_va_theta, "D_NEA": b.d_ne_theta, "R": b.ratio_theta}[stat]
    la, lb = _loo_values(a, stat), _loo_values(b, stat)
    ok = np.isfinite(la) & np.isfinite(lb)
    d_loo = (la - lb)[ok]
    g = len(d_loo)
    diff = float(th_a - th_b)
    if g < 2:
        return dict(diff=diff, se=np.nan, z=np.nan, se_independent=np.nan,
                    n_blocks_used=g)
    se = float(np.sqrt((g - 1) / g * np.sum((d_loo - d_loo.mean()) ** 2)))

    def _se(v):
        v = v[np.isfinite(v)]
        n = len(v)
        return np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2)) if n > 1 else np.nan

    se_ind = float(np.hypot(_se(la), _se(lb)))
    return dict(diff=diff, se=se, z=diff / se if se > 0 else np.nan,
                se_independent=se_ind, n_blocks_used=g)


# ------------------------------------------------------- proportionality fit --
def proportional_fit(tables: Mapping[str, BlockTable], labels: Sequence[str]):
    """Slope k of D_VA = k * D_NEA through the origin, over `labels`.

    Under a single shared Neanderthal source every cohort's Vindija preference is
    the *same* multiple of its Neanderthal quantity, so all cohorts lie on one
    line through the origin and k estimates D_VA of the source population itself.
    The fit is deliberately unweighted and origin-constrained: an intercept would
    absorb exactly the common-mode offsets (Yoruba's own archaic ancestry, the
    `.SG`/`.DG` asymmetry) that the design relies on cancelling, and inverse-
    variance weights would let the highest-coverage cohorts define the line the
    others are then tested against.

    Returns (k, residuals) with residuals in the same order as `labels`; feed
    this to `jackknife_over_blocks` to get standard errors for both.
    """
    x = np.array([tables[l].d_ne_theta for l in labels], dtype=np.float64)
    y = np.array([tables[l].d_va_theta for l in labels], dtype=np.float64)
    k = float((x * y).sum() / (x * x).sum())
    return k, y - k * x


def fit_and_residuals(tables: Mapping[str, BlockTable], fit_labels: Sequence[str],
                      test_labels: Sequence[str] = None, n_blocks: int = 50):
    """Origin-through slope k and per-cohort residuals, all with jackknife SEs.

    `fit_labels` define the line (use the well-powered cohorts assumed to share a
    source); `test_labels` are scored against it and default to `fit_labels`. A
    cohort whose residual is significantly non-zero is one whose Neanderthal
    ancestry has a different Vindija-vs-Altai character than the rest — the
    actual finding this study is looking for.
    """
    fit_labels = list(fit_labels)
    test_labels = list(test_labels) if test_labels is not None else list(fit_labels)
    needed = {l: tables[l] for l in dict.fromkeys(fit_labels + test_labels)}

    def _fn(tb):
        x = np.array([tb[l].d_ne_theta for l in fit_labels])
        y = np.array([tb[l].d_va_theta for l in fit_labels])
        k = (x * y).sum() / (x * x).sum()
        res = [tb[l].d_va_theta - k * tb[l].d_ne_theta for l in test_labels]
        return np.array([k] + res)

    theta, se = jackknife_over_blocks(_fn, needed, n_blocks)
    k, k_se = float(theta[0]), float(se[0])
    rows = []
    for i, lab in enumerate(test_labels, start=1):
        r, r_se = float(theta[i]), float(se[i])
        rows.append(dict(label=lab, residual=r, residual_se=r_se,
                         residual_z=r / r_se if r_se > 0 else np.nan,
                         in_fit=lab in fit_labels))
    return dict(k=k, k_se=k_se, k_z=k / k_se if k_se > 0 else np.nan), rows


# -------------------------------------------------------- detection limit -----
def detection_limit(null_diffs: Sequence[float], comparison_ses: Sequence[float],
                    signal: float, n_sigma: float = 2.0) -> dict:
    """Smallest source difference this data could have resolved.

    Two floors are combined and the larger one wins:

      * a **statistical** floor, `n_sigma` times the standard error of the
        cohort-versus-cohort differences the study actually makes; and
      * a **systematic** floor, the observed spread of differences among control
        splits of single cohorts that share a source by construction, which
        captures whatever coverage-, batch- and ascertainment-linked bias
        survives the design.

    The two inputs are deliberately taken from *different* sources, and getting
    this wrong is easy. Splitting one cohort in half gives two groups that share
    nearly all of their ancestry, so the standard error of that difference is
    several times smaller than the error on a comparison between two genuinely
    different populations. It is the right yardstick for bias and the wrong one
    for power, so `comparison_ses` must be the real comparisons and
    `null_diffs` the splits.

    `signal` must be a **D_VA magnitude** — the typical cohort's Vindija-over-
    Altai displacement — not the fitted slope k, which is a different quantity
    in different units. The limit is expressed as a fraction of it, because a
    cohort that replaced fraction f of its Neanderthal ancestry with ancestry
    from a lineage equidistant between Vindija and Altai would shift its D_VA by
    f times exactly that magnitude.
    """
    d = np.asarray([x for x in null_diffs if np.isfinite(x)], dtype=np.float64)
    s = np.asarray([x for x in comparison_ses if np.isfinite(x)], dtype=np.float64)
    stat_floor = float(n_sigma * np.median(s)) if len(s) else np.nan
    best_floor = float(n_sigma * np.min(s)) if len(s) else np.nan
    sys_floor = float(n_sigma * d.std(ddof=1)) if len(d) > 1 else np.nan
    limit = float(np.nanmax([stat_floor, sys_floor]))
    best = float(np.nanmax([best_floor, sys_floor]))
    sig = abs(float(signal)) if signal else np.nan
    return dict(statistical_floor=stat_floor, systematic_floor=sys_floor,
                best_case_floor=best_floor, limit=limit, best_limit=best,
                signal=sig,
                limit_fraction_of_signal=limit / sig if sig else np.nan,
                best_fraction_of_signal=best / sig if sig else np.nan,
                n_null=int(len(d)), n_comparisons=int(len(s)),
                n_sigma=float(n_sigma),
                max_abs_null=float(np.abs(d).max()) if len(d) else np.nan)


def technical_covariates(labels, values, covariates: Mapping[str, Sequence[float]]):
    """Rank correlation of a cohort statistic against technical covariates.

    An age-correlated signal in an ancient-DNA time series is the one that most
    urgently needs a technical explanation ruled out, because sample age is also
    a proxy for deamination damage, coverage and library preparation. This
    returns Spearman correlations of `values` against each covariate so the
    report can say whether a pattern tracks history or tracks the assay.
    Spearman rather than Pearson because coverage and damage are heavily skewed
    and a single deep genome should not set the slope.
    """
    from scipy.stats import spearmanr

    out = {}
    v = np.asarray(values, dtype=np.float64)
    for name, cov in covariates.items():
        c = np.asarray(cov, dtype=np.float64)
        ok = np.isfinite(v) & np.isfinite(c)
        if ok.sum() < 4:
            out[name] = dict(rho=np.nan, p=np.nan, n=int(ok.sum()))
            continue
        rho, p = spearmanr(v[ok], c[ok])
        out[name] = dict(rho=float(rho), p=float(p), n=int(ok.sum()))
    return out
