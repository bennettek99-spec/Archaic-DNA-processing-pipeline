"""
stats.py — archaic-affinity statistics with block-jackknife standard errors.

All statistics are allele-frequency forms of f-statistics (Patterson et al. 2012,
Genetics 192:1065). Two facts make the implementation simple and robust:

  * f4 and the normalised D-statistic are *invariant* to per-SNP allele
    polarisation: flipping the counted allele at a SNP sends every population's
    frequency p -> 1-p, which leaves (p_A-p_B)(p_C-p_D) and the D denominator
    unchanged. So we never need an explicit "derived allele" call — we just use a
    single consistent allele coding per SNP (which the packed .geno guarantees).

  * A ratio-of-sums estimator theta = sum(num)/sum(den) covers everything we need:
        f4(A,B;C,D)         -> num=(pA-pB)(pC-pD),               den=1
        D(W,X;Y,Z)          -> num=(pW-pX)(pY-pZ),               den=BABA+ABBA
        Neanderthal f4-ratio-> num=f4(Altai,Chimp;X,Mbuti),      den=f4(Altai,Chimp;Vindija,Mbuti)
    so one jackknife routine standard-errors all of them.

Standard errors use a delete-one block jackknife over contiguous, equal-SNP-count
genomic blocks (default 50). Equal block sizes make the unweighted jackknife
variance appropriate (Busing et al. 1999, Stat. Comput. 9:3 reduces to this when
block weights are equal); the blocks are large enough (~10k SNPs each) to absorb
linkage disequilibrium.
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------- blocks ------
def assign_blocks(n_snp: int, n_blocks: int = 50) -> np.ndarray:
    """Contiguous, ~equal-SNP-count block id for each SNP (in .geno/genomic order).

    The AADR .snp file is sorted by (chromosome, position), so contiguous blocks
    are contiguous genomic regions. A block may straddle a chromosome boundary;
    that is harmless for the jackknife (it only ever *removes* a region).
    """
    return (np.arange(n_snp) * n_blocks // max(n_snp, 1)).astype(np.int32)


# -------------------------------------------------- per-SNP statistic builders
def f4_array(pA, pB, pC, pD):
    """Per-SNP f4(A,B;C,D) = (pA-pB)(pC-pD). NaN where any pop has no data."""
    return (pA - pB) * (pC - pD)


def d_numerator(pW, pX, pY, pZ):
    """Per-SNP D numerator = (pW-pX)(pY-pZ) (= BABA - ABBA in freq form)."""
    return (pW - pX) * (pY - pZ)


def d_denominator(pW, pX, pY, pZ):
    """Per-SNP D denominator = (pW+pX-2 pW pX)(pY+pZ-2 pY pZ) (= BABA + ABBA)."""
    return (pW + pX - 2.0 * pW * pX) * (pY + pZ - 2.0 * pY * pZ)


# ------------------------------------------------------- jackknife ratio -------
def jackknife_ratio(num, den, block, n_blocks: int = 50):
    """Block-jackknife a ratio-of-sums estimator theta = sum(num)/sum(den).

    num, den : per-SNP arrays (NaN entries are dropped; a SNP is used only where
               BOTH num and den are finite, so num/den share one SNP set).
    Returns dict: theta, se, z (theta/se), n_used (SNPs), n_blocks_used.
    """
    num = np.asarray(num, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    ok = np.isfinite(num) & np.isfinite(den)
    n_used = int(ok.sum())
    if n_used == 0:
        return dict(theta=np.nan, se=np.nan, z=np.nan, n_used=0, n_blocks_used=0)

    num = np.where(ok, num, 0.0)
    den = np.where(ok, den, 0.0)
    Tn, Td = num.sum(), den.sum()
    theta = Tn / Td if Td != 0 else np.nan

    # per-block sums via bincount (fast, no Python loop over SNPs)
    bn = np.bincount(block, weights=num, minlength=n_blocks)
    bd = np.bincount(block, weights=den, minlength=n_blocks)
    has = np.bincount(block, weights=ok.astype(np.float64), minlength=n_blocks) > 0

    loo = []
    for b in range(n_blocks):
        if not has[b]:
            continue
        d_minus = Td - bd[b]
        if d_minus != 0:
            loo.append((Tn - bn[b]) / d_minus)
    loo = np.asarray(loo, dtype=np.float64)
    g = len(loo)
    if g > 1:
        se = np.sqrt((g - 1) / g * np.sum((loo - loo.mean()) ** 2))
    else:
        se = np.nan
    z = theta / se if (se and np.isfinite(se) and se > 0) else np.nan
    return dict(theta=float(theta), se=float(se), z=float(z),
                n_used=n_used, n_blocks_used=g)


# --------------------------------------- vectorised (batch) jackknife ---------
def block_starts(n_snp: int, n_blocks: int = 50) -> np.ndarray:
    """Row index where each contiguous equal-count block begins (for reduceat)."""
    block = assign_blocks(n_snp, n_blocks)
    return np.unique(block, return_index=True)[1].astype(np.int64)


def batch_jackknife_ratio(num, den, starts):
    """Vectorised block-jackknife of theta = sum(num)/sum(den) for MANY samples.

    num, den : float arrays (n_snp, K); NaN where a SNP is unusable for that
               sample. A SNP is used only where both num and den are finite, so
               num/den share one SNP set per sample (column).
    starts   : block start row indices from block_starts().
    Returns theta, se, z, n_used — each shape (K,). Matches jackknife_ratio()
    column-by-column (only non-empty blocks contribute; (g-1)/g jackknife var).
    """
    ok = np.isfinite(num) & np.isfinite(den)
    num0 = np.where(ok, num, 0.0)
    den0 = np.where(ok, den, 0.0)
    okf = ok.astype(num0.dtype)

    Bn = np.add.reduceat(num0, starts, axis=0)     # (G, K) block sums
    Bd = np.add.reduceat(den0, starts, axis=0)
    Bc = np.add.reduceat(okf, starts, axis=0)      # SNPs used per block
    Tn = Bn.sum(0).astype(np.float64)
    Td = Bd.sum(0).astype(np.float64)
    n_used = Bc.sum(0).astype(np.int64)

    with np.errstate(invalid="ignore", divide="ignore"):
        theta = Tn / Td
        loo = (Tn[None, :] - Bn) / (Td[None, :] - Bd)   # leave-one-block-out (G,K)
        nonempty = Bc > 0
        loo = np.where(nonempty, loo.astype(np.float64), np.nan)
        g = nonempty.sum(0).astype(np.float64)          # blocks with data (K,)
        mean = np.nanmean(loo, axis=0)
        var = (g - 1.0) / g * np.nansum((loo - mean) ** 2, axis=0)
        se = np.sqrt(var)
        z = theta / se
    return theta, se, z, n_used


# ------------------------------------------------- high-level convenience ------
def dstat(p, W, X, Y, Z, block, n_blocks=50):
    """Normalised D(W,X;Y,Z). p is a dict pop-name -> per-SNP freq array."""
    num = d_numerator(p[W], p[X], p[Y], p[Z])
    den = d_denominator(p[W], p[X], p[Y], p[Z])
    out = jackknife_ratio(num, den, block, n_blocks)
    out["statistic"] = f"D({W},{X};{Y},{Z})"
    return out


def f4_ratio(p, A, O, X, B, Ref, block, n_blocks=50):
    """f4-ratio alpha = f4(A,O; X,B) / f4(A,O; Ref,B).

    With A=Altai, O=Chimp, B=Mbuti (African baseline), Ref=Vindija (a second,
    independent high-coverage Neanderthal that scales '100% Neanderthal'),
    alpha estimates the Neanderthal-ancestry fraction of test population X.
    Using two *different* Neanderthals for the statistic (Altai) and the scale
    (Vindija) avoids the bias of using one genome as both source and yardstick
    (cf. Reich et al. 2009; Patterson et al. 2012; Petr et al. 2019 PNAS).
    """
    num = f4_array(p[A], p[O], p[X], p[B])
    den = f4_array(p[A], p[O], p[Ref], p[B])
    out = jackknife_ratio(num, den, block, n_blocks)
    out["statistic"] = f"alpha = f4({A},{O};{X},{B}) / f4({A},{O};{Ref},{B})"
    return out
