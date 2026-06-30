#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — estimate archaic introgression for every retained genome (vectorised).

For each Phase-2 sample:
  alpha_Nea  Neanderthal proportion  = f4(Altai,Chimp; X,Mbuti)/f4(Altai,Chimp; Vindija,Mbuti)
  D_Nea      Neanderthal affinity    = D(X, Mbuti; Altai, Chimp)     (+ jackknife Z)
  D_Den      Denisovan affinity      = D(X, Mbuti; Denisova, Chimp)  (+ jackknife Z)

Reference allele frequencies and all reference-only per-SNP constants are computed
ONCE; test individuals are processed in chunks with a fully vectorised block
jackknife (archaic.stats.batch_jackknife_ratio), which reproduces the per-individual
estimator validated in Phase 1 (verified by --selfcheck). Output:

  results/phase3_<panel>_estimates.csv      (resumable; skips ids already written)

Usage:
  python phase3_estimate.py --panel 1240k
  python phase3_estimate.py --panel 1240k --limit 200 --out results/_tmp.csv   # test
"""
import os, sys, time, argparse, csv
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st
from archaic.refs import PANELS

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
N_BLOCKS = 50
F = np.float32


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=128)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    cfg = PANELS[args.panel]

    meta = pd.read_csv(os.path.join(RESULTS, f"phase2_{args.panel}_metadata.csv"))
    ids = meta["genetic_id"].tolist()
    if args.limit:
        ids = ids[:args.limit]
    print(f"Phase 3 — panel={args.panel}  samples={len(ids):,}  chunk={args.chunk}")

    panel = Panel(cfg["prefix"])
    starts = st.block_starts(panel.n_snp, N_BLOCKS)

    print("Reference allele frequencies (computed once)...")
    rf, ri = panel.frequencies({k: cfg["refs"][k] for k in
                                ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]})
    for k in ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]:
        print(f"  {k:9s} SNPs={ri[k]['n_snp_covered']:,}")

    # reference-only per-SNP constants (float32); NaN propagates where a ref is missing
    pAlt = rf["Altai"].astype(F); pChi = rf["Chimp"].astype(F)
    pVin = rf["Vindija"].astype(F); pDen = rf["Denisova"].astype(F)
    pMb = rf["Mbuti"].astype(F)
    axN = pAlt - pChi                                   # Neanderthal-derived axis
    axD = pDen - pChi                                   # Denisovan-derived axis
    denYN = pAlt + pChi - 2.0 * pAlt * pChi             # D_Nea denom (Y,Z part)
    denYD = pDen + pChi - 2.0 * pDen * pChi             # D_Den denom (Y,Z part)
    a_vmb = axN * (pVin - pMb)                          # f4-ratio denominator (const)

    col_of = panel._id_to_col
    use = [(i, col_of[i]) for i in ids if i in col_of]
    missing = [i for i in ids if i not in col_of]
    if missing:
        print(f"  WARNING: {len(missing)} ids not in .ind (skipped)")

    out_path = args.out or os.path.join(RESULTS, f"phase3_{args.panel}_estimates.csv")
    fields = ["genetic_id", "alpha_Nea", "alpha_SE", "alpha_nSNP",
              "D_Nea", "D_Nea_SE", "D_Nea_Z", "D_Nea_nSNP",
              "D_Den", "D_Den_SE", "D_Den_Z", "D_Den_nSNP"]
    done_ids = set()
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        try:
            done_ids = set(pd.read_csv(out_path, usecols=["genetic_id"])["genetic_id"])
        except Exception:
            done_ids = set()
    if done_ids:
        before = len(use)
        use = [(i, c) for (i, c) in use if i not in done_ids]
        print(f"  resume: {len(done_ids):,} done, {len(use):,}/{before:,} remaining")
    fh = open(out_path, "a", newline="")
    writer = csv.DictWriter(fh, fieldnames=fields)
    if not done_ids:
        writer.writeheader()

    t0 = time.time(); done = 0
    for c0 in range(0, len(use), args.chunk):
        chunk = use[c0:c0 + args.chunk]
        cols = np.array([c for _, c in chunk], dtype=np.int64)
        G = panel.pg.read(panel.snp_rows, cols)                 # (nsnp, K) int8
        pX = G.astype(F); pX[G < 0] = np.nan; pX *= F(0.5)      # allele freq
        diffXM = pX - pMb[:, None]                              # (nsnp, K)
        denX = pX + pMb[:, None] - 2.0 * pX * pMb[:, None]      # D denom (X,M part)

        # alpha_Nea = f4(Altai,Chimp; X,Mbuti) / f4(Altai,Chimp; Vindija,Mbuti)
        num = axN[:, None] * diffXM
        den = np.broadcast_to(a_vmb[:, None], num.shape)
        a_t, a_se, _, a_n = st.batch_jackknife_ratio(num, den, starts)

        # D_Nea = D(X, Mbuti; Altai, Chimp)
        num = diffXM * axN[:, None]
        den = denX * denYN[:, None]
        dn_t, dn_se, dn_z, dn_n = st.batch_jackknife_ratio(num, den, starts)

        # D_Den = D(X, Mbuti; Denisova, Chimp)
        num = diffXM * axD[:, None]
        den = denX * denYD[:, None]
        dd_t, dd_se, dd_z, dd_n = st.batch_jackknife_ratio(num, den, starts)

        for k, (gid, _) in enumerate(chunk):
            writer.writerow(dict(
                genetic_id=gid,
                alpha_Nea=round(float(a_t[k]), 6), alpha_SE=round(float(a_se[k]), 6),
                alpha_nSNP=int(a_n[k]),
                D_Nea=round(float(dn_t[k]), 6), D_Nea_SE=round(float(dn_se[k]), 6),
                D_Nea_Z=round(float(dn_z[k]), 3), D_Nea_nSNP=int(dn_n[k]),
                D_Den=round(float(dd_t[k]), 6), D_Den_SE=round(float(dd_se[k]), 6),
                D_Den_Z=round(float(dd_z[k]), 3), D_Den_nSNP=int(dd_n[k]),
            ))
        done += len(chunk); fh.flush()
        rate = done / (time.time() - t0)
        print(f"  {done:,}/{len(use):,}  ({rate:.0f}/s, ETA {(len(use)-done)/rate/60:.1f} min)")
    fh.close()
    print(f"\nDone in {(time.time()-t0)/60:.1f} min -> {out_path}")


if __name__ == "__main__":
    main()
