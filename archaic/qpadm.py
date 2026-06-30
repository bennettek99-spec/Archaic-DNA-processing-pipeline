"""
qpadm.py — basic qpAdm ancestry modelling (pure Python), to quantify a target
population as a mixture of source populations relative to a set of outgroups.

This is the "rotating outgroup" f4-ratio form of qpAdm (Haak et al. 2015):
  for sources S1..Sn and outgroups R0..Rm, solve for weights w (summing to 1)
        f4(Target, S1; R0, Rj) = sum_{k>1} w_k * f4(Sk, S1; R0, Rj)   for all Rj
  by least squares; weights have block-jackknife SEs and the model fit is a
  GLS chi-square p-value (a plausible model has p > 0.05). It is a simplified
  qpAdm — use it for self-contained estimates and cross-check the authoritative
  ADMIXTOOLS 2 qpadm (tools/qpadm_admixtools.R) where available.
"""
from __future__ import annotations
import math
import numpy as np


def _f4(freq, W, X, Y, Z, sel):
    a = (freq[W] - freq[X]) * (freq[Y] - freq[Z])
    return float(np.nanmean(a[sel]))


def _chi2_sf(x, k):
    if k <= 0 or not np.isfinite(x):
        return float("nan")
    if x <= 0:
        return 1.0
    t = ((x / k) ** (1 / 3) - (1 - 2 / (9 * k))) / math.sqrt(2 / (9 * k))   # Wilson-Hilferty
    return 0.5 * math.erfc(t / math.sqrt(2))


def qpadm(freq, target, sources, outgroups, block, n_blocks=50):
    """Return dict(sources, weights, se, chi2, dof, p, n_snp). freq must contain
    target, every source and every outgroup as per-SNP allele-frequency arrays."""
    S1, others = sources[0], sources[1:]
    R0, Rj = outgroups[0], outgroups[1:]
    need = [target] + list(sources) + list(outgroups)
    mask = np.all([np.isfinite(freq[p]) for p in need], axis=0)

    def solve(sel):
        b = np.array([_f4(freq, target, S1, R0, r, sel) for r in Rj])
        A = np.array([[_f4(freq, s, S1, R0, r, sel) for s in others] for r in Rj])
        a_rest, *_ = np.linalg.lstsq(A, b, rcond=None)
        w = np.concatenate([[1 - a_rest.sum()], a_rest])
        return w, A, b, a_rest

    w, A, b, a_rest = solve(mask)
    # jackknife weights
    loo = []
    for bl in range(n_blocks):
        sel = mask & (block != bl)
        if sel.sum() < 100:
            continue
        try:
            ww, *_ = solve(sel); loo.append(ww)
        except Exception:
            pass
    loo = np.array(loo); B = len(loo)
    se = (np.sqrt((B - 1) / B * np.sum((loo - loo.mean(0)) ** 2, axis=0))
          if B > 1 else np.full(len(w), np.nan))
    # GLS chi-square fit p-value from the jackknife covariance of the residual
    resid = b - A @ a_rest
    lr = []
    for bl in range(n_blocks):
        sel = mask & (block != bl)
        if sel.sum() < 100:
            continue
        bb = np.array([_f4(freq, target, S1, R0, r, sel) for r in Rj])
        AA = np.array([[_f4(freq, s, S1, R0, r, sel) for s in others] for r in Rj])
        lr.append(bb - AA @ a_rest)
    lr = np.array(lr); Bn = len(lr)
    try:
        cov = (Bn - 1) / Bn * np.cov(lr.T, bias=True) * Bn
        chi2 = float(resid @ np.linalg.pinv(np.atleast_2d(cov)) @ resid)
        dof = len(Rj) - len(others)
        p = _chi2_sf(chi2, dof)
    except Exception:
        chi2 = float("nan"); dof = len(Rj) - len(others); p = float("nan")
    return dict(sources=list(sources), weights=w, se=se, chi2=chi2, dof=dof,
                p=p, n_snp=int(mask.sum()))
