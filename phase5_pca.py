#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 5a — ancestry PCA for the retained samples.

Ancestry is the dominant predictor of archaic ancestry (West- vs East-Eurasian
etc.), so the expected-introgression model needs ancestry coordinates. We compute
a memory-bounded PCA: an even autosomal SNP subset (~30k), Patterson-style
standardisation (centre by 2p, scale 1/sqrt(2p(1-p)), missing -> mean), and a
randomised SVD for the top PCs. PCs are written next to the Phase-4 table.

Sanity checks printed: PC1/PC2 vs longitude/latitude (geographic structure), and
the weighted correlation of alpha_Nea with each PC (the East-Asian excess should
appear as a real ancestry gradient, not as noise).

Output: results/phase5_<panel>_pca.csv   (genetic_id, PC1..PCk)
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic.refs import PANELS

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def randomized_svd(M, k, n_oversamples=10, n_iter=4, seed=0):
    rng = np.random.default_rng(seed)
    n, d = M.shape
    l = k + n_oversamples
    Omega = rng.standard_normal((d, l)).astype(np.float32)
    Y = M @ Omega
    for _ in range(n_iter):
        Y = M @ (M.T @ Y)
    Q, _ = np.linalg.qr(Y)
    B = Q.T @ M
    Ub, S, Vt = np.linalg.svd(B, full_matrices=False)
    U = Q @ Ub
    return U[:, :k], S[:k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--nsnp", type=int, default=30000)
    ap.add_argument("--k", type=int, default=10)
    args = ap.parse_args()
    cfg = PANELS[args.panel]

    df = pd.read_csv(os.path.join(RESULTS, f"phase4_{args.panel}_analysis.csv"))
    ids = df["genetic_id"].tolist()
    print(f"Phase 5a PCA — panel={args.panel}  samples={len(ids):,}  snp_subset={args.nsnp:,}")

    panel = Panel(cfg["prefix"])
    # even autosomal SNP subset
    sub = np.linspace(0, panel.n_snp - 1, args.nsnp).astype(np.int64)
    snp_rows = panel.snp_rows[sub]
    cols = np.array([panel._id_to_col[i] for i in ids], dtype=np.int64)

    print("Reading genotype subset (one pass over the panel)...")
    G = panel.pg.read(snp_rows, cols)                 # (nsnp, nsamp) int8
    nsnp, nsamp = G.shape
    Gf = G.astype(np.float32)
    Gf[G < 0] = np.nan

    print("Standardising (Patterson 2006) ...")
    with np.errstate(invalid="ignore"):
        p = np.nanmean(Gf, axis=1) / 2.0              # allele freq per SNP
    keep = np.isfinite(p) & (p > 0.0) & (p < 1.0)
    Gf = Gf[keep]; p = p[keep]
    scale = np.sqrt(2.0 * p * (1.0 - p)).astype(np.float32)
    # centre by 2p, impute missing -> 0 (the centred mean), scale
    M = (Gf - (2.0 * p)[:, None])
    M[np.isnan(M)] = 0.0
    M /= scale[:, None]
    M = M.T.astype(np.float32)                        # (nsamp, nsnp_kept)
    del Gf, G
    print(f"  PCA matrix: {M.shape[0]:,} samples x {M.shape[1]:,} SNPs")

    U, S = randomized_svd(M, args.k)
    PC = (U * S).astype(np.float64)                   # sample scores
    var = (S ** 2) / (S ** 2).sum()
    print("  variance explained (top PCs): " +
          ", ".join(f"PC{i+1}={var[i]*100:.1f}%" for i in range(min(6, args.k))))

    pcs = pd.DataFrame(PC, columns=[f"PC{i+1}" for i in range(args.k)])
    pcs.insert(0, "genetic_id", ids)
    out = os.path.join(RESULTS, f"phase5_{args.panel}_pca.csv")
    pcs.to_csv(out, index=False)

    # ---- sanity checks -------------------------------------------------------
    j = df.merge(pcs, on="genetic_id")
    w = j["weight"].values
    def wcorr(a, b):
        m = np.isfinite(a) & np.isfinite(b) & np.isfinite(w) & (w > 0)
        if m.sum() < 10:
            return np.nan
        aa, bb, ww = a[m], b[m], w[m]
        ma, mb = np.average(aa, weights=ww), np.average(bb, weights=ww)
        cov = np.average((aa - ma) * (bb - mb), weights=ww)
        va = np.average((aa - ma) ** 2, weights=ww); vb = np.average((bb - mb) ** 2, weights=ww)
        return cov / np.sqrt(va * vb) if va > 0 and vb > 0 else np.nan
    print("\nSanity — geographic structure & archaic gradient (Pearson / weighted):")
    print(f"  {'PC':4s} {'corr(lon)':>10s} {'corr(lat)':>10s} {'wcorr(alpha_Nea)':>17s}")
    for i in range(min(6, args.k)):
        c = f"PC{i+1}"
        rl = np.corrcoef(j[c], j["lon"])[0, 1] if j["lon"].notna().sum() else np.nan
        ra = np.corrcoef(j[c], j["lat"])[0, 1] if j["lat"].notna().sum() else np.nan
        rn = wcorr(j["alpha_Nea"].values, j[c].values)
        print(f"  {c:4s} {rl:10.3f} {ra:10.3f} {rn:17.3f}")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
