"""
neighborhood.py - local expected-value model for individual outlier scans.

Phase 6 and Phase 9 both ask the same question: given a target genome's ancestry,
geography, and age coordinates, what Neanderthal ancestry would we expect from
nearby high-confidence genomes? This module keeps that residual calculation in
one tested place.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def zscore(x):
    """NaN-aware z-score with a unit fallback for constant vectors."""
    x = np.asarray(x, dtype=np.float64)
    mu, sd = np.nanmean(x), np.nanstd(x)
    return (x - mu) / (sd if sd > 0 else 1.0)


def feature_matrix(df, features, weights):
    """Build the weighted, z-scored feature matrix used for neighbour search."""
    cols = []
    for f in features:
        vals = df[f].fillna(df[f].median()).values
        cols.append(zscore(vals) * weights[f])
    return np.column_stack(cols)


def local_residual_stats(df, X, ref_mask, ids, K, alpha_col="alpha_adj", se_col="alpha_SE"):
    """Expected value and residual z-score from precision-weighted neighbours.

    The expected-value uncertainty uses the *effective* neighbour count

        n_eff = (sum w)^2 / sum(w^2)

    rather than the nominal K, because precision weighting can make a small
    number of high-coverage neighbours dominate the local mean.
    """
    ref_mask = np.asarray(ref_mask, dtype=bool)
    ids = np.asarray(ids, dtype=object)
    n_ref = int(ref_mask.sum())
    if n_ref == 0:
        raise ValueError("No reference samples available for local expectation.")

    Xref = X[ref_mask]
    a_ref = df[alpha_col].to_numpy(dtype=np.float64)[ref_mask]
    se_ref = df[se_col].to_numpy(dtype=np.float64)[ref_mask]
    id_ref = ids[ref_mask]

    tree = cKDTree(Xref)
    k_query = min(int(K) + 1, n_ref)
    _, idx = tree.query(X, k=k_query)
    if k_query == 1:
        idx = idx[:, None]

    A = a_ref[idx]
    SE = se_ref[idx]
    Wn = 1.0 / np.clip(SE, 1e-4, None) ** 2
    Wn = np.where(id_ref[idx] == ids[:, None], 0.0, Wn)

    Wsum = Wn.sum(1)
    W2sum = (Wn ** 2).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        expected = (Wn * A).sum(1) / Wsum
        var_obs = (Wn * (A - expected[:, None]) ** 2).sum(1) / Wsum
        mean_meas_var = (Wn * SE ** 2).sum(1) / Wsum
        n_eff = (Wsum ** 2) / W2sum

    sigma_bio2 = np.clip(var_obs - mean_meas_var, 0.0, None)
    se_expected2 = var_obs / np.clip(n_eff, 1.0, None)
    denom = np.sqrt(df[se_col].to_numpy(dtype=np.float64) ** 2 + sigma_bio2 + se_expected2)
    residual = df[alpha_col].to_numpy(dtype=np.float64) - expected
    z = residual / denom

    bad = ~(np.isfinite(Wsum) & (Wsum > 0))
    expected[bad] = np.nan
    var_obs[bad] = np.nan
    mean_meas_var[bad] = np.nan
    sigma_bio2[bad] = np.nan
    se_expected2[bad] = np.nan
    n_eff[bad] = 0.0
    residual[bad] = np.nan
    z[bad] = np.nan

    return dict(
        expected=expected,
        residual=residual,
        z=z,
        sigma_bio=np.sqrt(sigma_bio2),
        se_expected=np.sqrt(se_expected2),
        n_neighbors_eff=n_eff,
        var_obs=var_obs,
        mean_meas_var=mean_meas_var,
    )
