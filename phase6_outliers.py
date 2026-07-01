#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 5b + Phase 6 — expected introgression and unexpected-individual detection.

The research question: which ancient individuals carry significantly MORE or LESS
Neanderthal ancestry than expected for someone of their ancestry, geography, and
time period? We answer it with a local (nearest-neighbour) expectation plus a
rigorous standardised residual that separates MEASUREMENT noise from BIOLOGICAL
scatter:

  1. data-type adjustment: remove the small additive offset between sequencing
     data types (shotgun SG vs capture AG, ~0.3pp from Phase 4) so a sample is not
     flagged merely for how it was sequenced.
  2. feature space: z-scored ancestry PCs (dominant) + longitude + latitude + age.
  3. expectation: for each sample, the precision-weighted (1/SE^2) mean Neanderthal
     level of its K genetically/geographically/temporally nearest HIGH-CONFIDENCE
     neighbours (self excluded). This is "expected given ancestry+place+time".
  4. variance decomposition at each sample's neighbourhood:
        Var_obs(neighbours) = biological_scatter^2 + mean(measurement SE^2)
     => sigma_bio^2 = max(0, Var_obs - mean SE_n^2);  SE_expected^2 = Var_obs / K.
  5. standardised residual:
        z = (alpha_adj - Expected) / sqrt(SE_self^2 + sigma_bio^2 + SE_expected^2)
     Ranking |z| gives candidate over-/under-introgressed individuals. Only
     HIGH-CONFIDENCE samples are eligible as candidates (others are too noisy);
     Phase 9 will test robustness.

Also fits a weighted ridge (ancestry PCs + geo + time) as an interpretable
"predictive model" and reports its cross-validated R^2.

Outputs:
  results/phase6_<panel>_residuals.csv         all samples: expected, residual, z, neighbourhood stats
  results/phase6_<panel>_top_outliers.txt      ranked +/- candidates (high-confidence)
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.refs import PANELS

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
PCS = [f"PC{i}" for i in range(1, 7)]
FEATURES = PCS + ["lon", "lat", "date_bp"]
FEAT_WEIGHT = {**{p: 1.0 for p in PCS}, "lon": 0.5, "lat": 0.5, "date_bp": 0.7}
K = 80


def zscore(x):
    x = x.astype(np.float64)
    mu, sd = np.nanmean(x), np.nanstd(x)
    return (x - mu) / (sd if sd > 0 else 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--k", type=int, default=K)
    args = ap.parse_args()

    df = pd.read_csv(os.path.join(RESULTS, f"phase4_{args.panel}_analysis.csv"))
    pc = pd.read_csv(os.path.join(RESULTS, f"phase5_{args.panel}_pca.csv"))
    df = df.merge(pc, on="genetic_id")
    print(f"Phase 6 — panel={args.panel}  samples={len(df):,}  K={args.k}")

    # 1. data-type additive adjustment (centre each type to the global weighted mean)
    w = (1.0 / df["alpha_SE"].clip(lower=1e-4) ** 2).values
    gmean = np.average(df["alpha_Nea"], weights=w)
    off = {}
    for dt, g in df.groupby("data_type"):
        wg = (1.0 / g["alpha_SE"].clip(lower=1e-4) ** 2).values
        off[dt] = np.average(g["alpha_Nea"], weights=wg) - gmean if len(g) >= 30 else 0.0
    df["alpha_adj"] = df["alpha_Nea"] - df["data_type"].map(off).fillna(0.0)

    # 2. feature space (z-scored, weighted); impute missing coords to median
    X = np.column_stack([
        zscore(df[f].fillna(df[f].median()).values) * FEAT_WEIGHT[f] for f in FEATURES
    ])

    # 3-4. neighbourhood expectation on the HIGH-CONFIDENCE reference set
    ref = df["high_conf"].values
    ids_np = df["genetic_id"].to_numpy(dtype=object)
    Xref = X[ref]
    a_ref = df["alpha_adj"].values[ref]
    se_ref = df["alpha_SE"].values[ref]
    id_ref = ids_np[ref]
    print(f"  high-confidence reference samples: {ref.sum():,}")

    tree = cKDTree(Xref)
    kq = args.k + 1                                   # +1 to drop self
    dist, idx = tree.query(X, k=kq)                   # (n, kq)

    A = a_ref[idx]                                    # neighbour alpha_adj
    SE = se_ref[idx]
    Wn = 1.0 / np.clip(SE, 1e-4, None) ** 2
    self_mask = id_ref[idx] == ids_np[:, None]
    Wn = np.where(self_mask, 0.0, Wn)

    Wsum = Wn.sum(1)
    expected = (Wn * A).sum(1) / Wsum                 # precision-weighted local mean
    # weighted variance of neighbours around the local mean
    var_obs = (Wn * (A - expected[:, None]) ** 2).sum(1) / Wsum
    mean_meas_var = (Wn * SE ** 2).sum(1) / Wsum
    sigma_bio2 = np.clip(var_obs - mean_meas_var, 0.0, None)
    se_expected2 = var_obs / args.k

    se_self = df["alpha_SE"].values
    denom = np.sqrt(se_self ** 2 + sigma_bio2 + se_expected2)
    residual = df["alpha_adj"].values - expected
    z = residual / denom

    df["expected_Nea"] = expected
    df["residual_Nea"] = residual
    df["z_resid"] = z
    df["sigma_bio_local"] = np.sqrt(sigma_bio2)
    df["n_neighbors_eff"] = (Wsum ** 2) / (Wn ** 2).sum(1)   # effective neighbour count

    # ---- interpretable predictive model (Phase 5): weighted ridge, CV R^2 ------
    Xr = X[ref]
    yr = df["alpha_adj"].values[ref]
    ridge = Ridge(alpha=1.0)
    r2 = cross_val_score(ridge, Xr, yr, cv=5, scoring="r2")
    print(f"  predictive model (ridge on PCs+geo+time) 5-fold R^2 = "
          f"{r2.mean():.3f} +/- {r2.std():.3f}")
    print("    (low R^2 is expected & honest: within-ancestry Neanderthal variation "
          "is small and largely measurement noise)")

    # ---- outputs --------------------------------------------------------------
    rpath = os.path.join(RESULTS, f"phase6_{args.panel}_residuals.csv")
    keep_cols = ["genetic_id", "group_id", "country", "lat", "lon", "date_bp",
                 "data_type", "alpha_Nea", "alpha_adj", "alpha_SE", "alpha_nSNP",
                 "high_conf", "expected_Nea", "residual_Nea", "z_resid",
                 "sigma_bio_local", "n_neighbors_eff", "D_Nea_Z", "D_Den_Z", "flags"]
    df[keep_cols].to_csv(rpath, index=False)

    cand = df[df["high_conf"]].copy()
    cand = cand[np.isfinite(cand["z_resid"])]
    top_pos = cand.sort_values("z_resid", ascending=False).head(25)
    top_neg = cand.sort_values("z_resid").head(25)

    def fmt(s):
        L = []
        for _, r in s.iterrows():
            L.append(f"  z={r['z_resid']:+5.2f}  obs={r['alpha_Nea']*100:5.2f}% "
                     f"exp={r['expected_Nea']*100:5.2f}% (res={r['residual_Nea']*100:+5.2f}pp) "
                     f"nSNP={int(r['alpha_nSNP']):>7,}  {int(r['date_bp']):>6}BP  "
                     f"{str(r['country'])[:14]:14s} {str(r['group_id'])[:42]}")
        return "\n".join(L)

    txt = []
    txt.append(f"Phase 6 — unexpected archaic introgression candidates (panel {args.panel})")
    txt.append(f"high-confidence eligible: {len(cand):,}   K={args.k}   "
               f"data-type offsets: " + ", ".join(f"{k}:{v*100:+.2f}pp" for k, v in off.items() if abs(v) > 1e-9))
    txt.append("z = (observed_adj - expected) / sqrt(SE_self^2 + bio_scatter^2 + SE_expected^2)")
    txt.append("HYPOTHESES ONLY — Phase 7 investigates each, Phase 9 tests robustness.")
    txt.append("")
    txt.append("=== MORE Neanderthal than expected (top 25 positive z) ===")
    txt.append(fmt(top_pos))
    txt.append("")
    txt.append("=== LESS Neanderthal than expected (top 25 negative z) ===")
    txt.append(fmt(top_neg))
    report = "\n".join(txt)
    tpath = os.path.join(RESULTS, f"phase6_{args.panel}_top_outliers.txt")
    with open(tpath, "w", encoding="utf-8") as fh:
        fh.write(report + "\n")

    # brief console view
    print(f"\n  |z| distribution (high-conf): "
          f"95%={np.nanpercentile(np.abs(cand['z_resid']),95):.2f}  "
          f"99%={np.nanpercentile(np.abs(cand['z_resid']),99):.2f}  "
          f"max={np.nanmax(np.abs(cand['z_resid'])):.2f}")
    print("\n  Top 8 MORE-than-expected:")
    print(fmt(top_pos.head(8)))
    print("\n  Top 8 LESS-than-expected:")
    print(fmt(top_neg.head(8)))
    print(f"\nWrote:\n  {rpath}\n  {tpath}")


if __name__ == "__main__":
    main()
