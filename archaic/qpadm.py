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

Performance note: every f4(W,X;Y,Z) term needed by the linear system (both at
full data and at each of the 50 leave-one-block-out jackknife replicates) is a
mean over up to ~1.2M SNPs. Recomputing that mean from scratch inside the
jackknife loop costs 50x more genome scans than necessary — sum(mean(A\\B)) can
be decomposed into per-block sums computed ONCE, from which every leave-one-
block-out mean is then a couple of array subtractions (Busing et al. 1999's
usual delete-one-block trick, same idea as archaic.stats.batch_jackknife_ratio).
`_build_system` does this: one O(n_snp) pass per f4 term, not one per block.
"""
from __future__ import annotations
import math
import numpy as np

try:                                             # SciPy is optional at import time
    from scipy.optimize import minimize as _minimize
except Exception:                                # pragma: no cover
    _minimize = None

MIN_SNP_PER_LOO = 100    # minimum SNPs required outside a dropped block


def _chi2_sf(x, k):
    if k <= 0 or not np.isfinite(x):
        return float("nan")
    if x <= 0:
        return 1.0
    t = ((x / k) ** (1 / 3) - (1 - 2 / (9 * k))) / math.sqrt(2 / (9 * k))   # Wilson-Hilferty
    return 0.5 * math.erfc(t / math.sqrt(2))


def _block_sums(freq, quads, mask, block, n_blocks):
    """Per-block sum of (pW-pX)(pY-pZ) over `mask`, for each (W,X,Y,Z) in quads.
    Returns (n_blocks, len(quads)) float64 — one O(n_snp) pass per quad."""
    bsum = np.empty((n_blocks, len(quads)), dtype=np.float64)
    for qi, (W, X, Y, Z) in enumerate(quads):
        a = (freq[W] - freq[X]) * (freq[Y] - freq[Z])
        a = np.where(mask, a, 0.0)
        bsum[:, qi] = np.bincount(block, weights=a, minlength=n_blocks)
    return bsum


def _build_system(freq, target, sources, outgroups, block, n_blocks):
    """Build the qpAdm linear system f4(Target,S1;R0,Rj) = A @ a_rest once, plus
    everything needed for a delete-one-block jackknife of it, using one genome
    scan per f4 term (not one per jackknife block). Returns a dict of arrays."""
    S1, others = sources[0], sources[1:]
    R0, Rj = outgroups[0], outgroups[1:]
    need = [target] + list(sources) + list(outgroups)
    mask = np.all([np.isfinite(freq[p]) for p in need], axis=0)
    total_cnt = float(mask.sum())
    cnt_block = np.bincount(block, weights=mask.astype(np.float64), minlength=n_blocks)

    b_quads = [(target, S1, R0, r) for r in Rj]
    A_quads = [(s, S1, R0, r) for r in Rj for s in others]   # r outer, s inner

    bsum_b = _block_sums(freq, b_quads, mask, block, n_blocks)
    bsum_A = _block_sums(freq, A_quads, mask, block, n_blocks) if others else \
        np.zeros((n_blocks, 0))

    b_full = bsum_b.sum(0) / total_cnt
    A_full = (bsum_A.sum(0) / total_cnt).reshape(len(Rj), len(others))

    denom = total_cnt - cnt_block                              # SNPs left per LOO
    valid = denom >= MIN_SNP_PER_LOO
    with np.errstate(invalid="ignore", divide="ignore"):
        b_loo = (bsum_b.sum(0)[None, :] - bsum_b) / denom[:, None]
        A_loo_flat = (bsum_A.sum(0)[None, :] - bsum_A) / denom[:, None] if others else \
            np.zeros((n_blocks, 0))
    b_loo = np.where(valid[:, None], b_loo, np.nan)
    A_loo = np.where(valid[:, None], A_loo_flat, np.nan).reshape(n_blocks, len(Rj), len(others))

    return dict(S1=S1, others=others, R0=R0, Rj=Rj, mask=mask, valid=valid,
               A=A_full, b=b_full, A_loo=A_loo, b_loo=b_loo, n_snp=int(total_cnt))


def _gls_fit(sys_, a_rest):
    """GLS chi-square fit p-value from the jackknife covariance of the residual
    b - A@a_rest, evaluated at each valid leave-one-block-out replicate."""
    A, b, Rj, others = sys_["A"], sys_["b"], sys_["Rj"], sys_["others"]
    resid = b - A @ a_rest
    lr = []
    for bl in range(len(sys_["valid"])):
        if not sys_["valid"][bl]:
            continue
        Abl, bbl = sys_["A_loo"][bl], sys_["b_loo"][bl]
        lr.append(bbl - Abl @ a_rest)
    lr = np.array(lr); Bn = len(lr)
    try:
        cov = (Bn - 1) / Bn * np.cov(lr.T, bias=True) * Bn
        chi2 = float(resid @ np.linalg.pinv(np.atleast_2d(cov)) @ resid)
        dof = len(Rj) - len(others)
        p = _chi2_sf(chi2, dof)
    except Exception:
        chi2 = float("nan"); dof = len(Rj) - len(others); p = float("nan")
    return chi2, dof, p


def _rank_approx(M, rank):
    """Best rank-r approximation by SVD."""
    if rank <= 0:
        return np.zeros_like(M)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    r = min(rank, len(S))
    return (U[:, :r] * S[:r]) @ Vt[:r, :]


def _f4_matrix_system(freq, lefts, rights, block, n_blocks):
    """Build qpWave's f4 matrix and leave-one-block-out replicates.

    Matrix entries are f4(L_i, L0; R_j, R0), i>0, j>0.  A low-rank matrix means
    the left populations are related to the right outgroups through only a small
    number of ancestry streams.
    """
    L0, Li = lefts[0], lefts[1:]
    R0, Rj = rights[0], rights[1:]
    need = list(lefts) + list(rights)
    mask = np.all([np.isfinite(freq[p]) for p in need], axis=0)
    total_cnt = float(mask.sum())
    cnt_block = np.bincount(block, weights=mask.astype(np.float64), minlength=n_blocks)
    quads = [(l, L0, r, R0) for l in Li for r in Rj]
    bsum = _block_sums(freq, quads, mask, block, n_blocks)
    full = (bsum.sum(0) / total_cnt).reshape(len(Li), len(Rj))
    denom = total_cnt - cnt_block
    valid = denom >= MIN_SNP_PER_LOO
    with np.errstate(invalid="ignore", divide="ignore"):
        loo_flat = (bsum.sum(0)[None, :] - bsum) / denom[:, None]
    loo_flat = np.where(valid[:, None], loo_flat, np.nan)
    loo = loo_flat.reshape(n_blocks, len(Li), len(Rj))
    return dict(matrix=full, loo=loo, valid=valid, n_snp=int(total_cnt),
                lefts=list(lefts), rights=list(rights))


def qpwave(freq, lefts, rights, block, n_blocks=50, max_rank=None):
    """Simplified qpWave rank tests for left/right population sets.

    Returns one row per tested rank with chi-square p-value.  This is a pure
    Python companion diagnostic; authoritative publication runs should still be
    cross-checked against ADMIXTOOLS 2.
    """
    if len(lefts) < 2 or len(rights) < 2:
        return []
    sys_ = _f4_matrix_system(freq, lefts, rights, block, n_blocks)
    M = sys_["matrix"]
    nr, nc = M.shape
    max_rank = min(nr, nc) - 1 if max_rank is None else min(max_rank, min(nr, nc) - 1)
    out = []
    loo = np.array([sys_["loo"][b] for b in range(n_blocks) if sys_["valid"][b]])
    B = len(loo)
    for rank in range(max_rank + 1):
        resid = (M - _rank_approx(M, rank)).reshape(-1)
        lr = np.array([(m - _rank_approx(m, rank)).reshape(-1) for m in loo])
        dof = (nr - rank) * (nc - rank)
        if B > 1 and dof > 0:
            cov = (B - 1) / B * np.cov(lr.T, bias=True) * B
            chi2 = float(resid @ np.linalg.pinv(np.atleast_2d(cov)) @ resid)
            p = _chi2_sf(chi2, dof)
        else:
            chi2 = float("nan"); p = float("nan")
        out.append(dict(rank=rank, chi2=chi2, dof=dof, p=p,
                        n_snp=sys_["n_snp"], n_left=len(lefts), n_right=len(rights)))
    return out


def qpadm(freq, target, sources, outgroups, block, n_blocks=50):
    """Return dict(sources, weights, se, chi2, dof, p, n_snp). freq must contain
    target, every source and every outgroup as per-SNP allele-frequency arrays."""
    sys_ = _build_system(freq, target, sources, outgroups, block, n_blocks)
    A, b = sys_["A"], sys_["b"]
    a_rest, *_ = np.linalg.lstsq(A, b, rcond=None)
    w = np.concatenate([[1 - a_rest.sum()], a_rest])

    loo = []
    for bl in range(n_blocks):
        if not sys_["valid"][bl]:
            continue
        Abl, bbl = sys_["A_loo"][bl], sys_["b_loo"][bl]
        try:
            ww, *_ = np.linalg.lstsq(Abl, bbl, rcond=None)
            loo.append(np.concatenate([[1 - ww.sum()], ww]))
        except Exception:
            pass
    loo = np.array(loo); B = len(loo)
    se = (np.sqrt((B - 1) / B * np.sum((loo - loo.mean(0)) ** 2, axis=0))
          if B > 1 else np.full(len(w), np.nan))

    chi2, dof, p = _gls_fit(sys_, a_rest)
    feasible = bool(np.all(w > -1e-6) and np.all(w < 1 + 1e-6))
    return dict(sources=list(sources), weights=w, se=se, chi2=chi2, dof=dof,
                p=p, n_snp=sys_["n_snp"], feasible=feasible, constrained=False)


def _solve_constrained(A, b, x0=None):
    """min ||A x - b||^2 over x>=0 with sum(x) <= 1 (SLSQP). x are the non-pivot
    source weights a_rest; the pivot weight is 1 - sum(x)."""
    k = A.shape[1]
    if k == 0:
        return np.zeros(0)
    if x0 is None:
        x0, *_ = np.linalg.lstsq(A, b, rcond=None)
        x0 = np.clip(x0, 0, 1)
        if x0.sum() > 1:
            x0 = x0 / x0.sum()
    if _minimize is None:                        # graceful fallback: clip the OLS fit
        return x0
    obj = lambda x: float(np.sum((A @ x - b) ** 2))
    jac = lambda x: 2.0 * A.T @ (A @ x - b)
    cons = ({"type": "ineq", "fun": lambda x: 1.0 - np.sum(x),
             "jac": lambda x: -np.ones_like(x)},)
    res = _minimize(obj, x0, jac=jac, method="SLSQP",
                    bounds=[(0.0, 1.0)] * k, constraints=cons,
                    options=dict(maxiter=200, ftol=1e-12))
    return res.x if res.success else np.clip(res.x, 0, 1)


def qpadm_constrained(freq, target, sources, outgroups, block, n_blocks=50):
    """Constrained qpAdm: identical linear system to qpadm() but the mixture
    weights are forced onto the simplex (each in [0,1], summing to 1) via SLSQP.

    Unconstrained qpAdm can return negative or >1 weights, which are not
    interpretable as ancestry proportions; this variant always returns a valid
    mixture (a "supervised admixture" fit) and is the right thing to report when
    the unconstrained model is close to the simplex boundary. Weights carry
    block-jackknife SEs and the same GLS chi-square fit statistic is reported for
    comparison. Returns the same dict shape as qpadm(), with constrained=True."""
    sys_ = _build_system(freq, target, sources, outgroups, block, n_blocks)
    A, b = sys_["A"], sys_["b"]
    a_rest = _solve_constrained(A, b)
    w = np.concatenate([[1 - a_rest.sum()], a_rest]) if len(a_rest) else np.array([1.0])

    loo = []
    for bl in range(n_blocks):
        if not sys_["valid"][bl]:
            continue
        Abl, bbl = sys_["A_loo"][bl], sys_["b_loo"][bl]
        try:
            ar = _solve_constrained(Abl, bbl, x0=a_rest.copy())
            loo.append(np.concatenate([[1 - ar.sum()], ar]) if len(ar) else np.array([1.0]))
        except Exception:
            pass
    loo = np.array(loo); B = len(loo)
    se = (np.sqrt((B - 1) / B * np.sum((loo - loo.mean(0)) ** 2, axis=0))
          if B > 1 else np.full(len(w), np.nan))

    chi2, dof, p = _gls_fit(sys_, a_rest)
    return dict(sources=list(sources), weights=w, se=se, chi2=chi2, dof=dof,
                p=p, n_snp=sys_["n_snp"], feasible=True, constrained=True)


def compete_models(freq, target, models, outgroups, block, n_blocks=50,
                   constrained=True):
    """Run several candidate source models for one target and rank them.

    models: dict(model_name -> list of source names present in `freq`).
    Ranking key: feasible models first, then by higher fit p-value (a plausible
    model has p > 0.05). Returns a list of result dicts (best first), each with
    model, weights, se, p, feasible added.
    """
    run = qpadm_constrained if constrained else qpadm
    out = []
    for name, srcs in models.items():
        srcs = [s for s in srcs if s in freq]
        if len(srcs) < 1:
            continue
        try:
            r = run(freq, target, srcs, outgroups, block, n_blocks)
        except Exception:
            continue
        r["model"] = name
        out.append(r)
    out.sort(key=lambda r: (not r.get("feasible", False),
                            -(r["p"] if np.isfinite(r["p"]) else -1)))
    return out
