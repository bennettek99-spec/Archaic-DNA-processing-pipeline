#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 4 — normalisation / bias control.

Join Phase-3 estimates to Phase-2 metadata, derive technical covariates, and
quantify whether the archaic estimates are confounded by data quality before any
outlier is trusted. We do NOT silently "correct" the estimates; instead we:

  * carry the jackknife SE as each sample's measurement precision (weight 1/SE^2),
    which already absorbs the dominant coverage/SNP-overlap noise;
  * diagnose any *systematic* (directional) bias of alpha_Nea against technical
    covariates — sequencing data type (pseudo-haploid .AG vs diploid .DG/.SG),
    usable-SNP count, DNA damage, contamination — so Phase 5 can include the
    offending ones as nuisance covariates and Phase 6 residuals are net of them;
  * emit an analysis-ready table with a high-confidence flag.

Output:
  results/phase4_<panel>_analysis.csv     estimates + metadata + covariates + weights
  results/phase4_<panel>_bias_report.txt  diagnostics
"""
import os, sys, argparse, re
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.refs import PANELS

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def data_type(gid):
    m = re.search(r"\.([A-Za-z]+)$", str(gid))
    return m.group(1).upper() if m else "NA"


def wmean(x, w):
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    return np.average(x[m], weights=w[m]) if m.sum() else np.nan


def wcorr(x, y, w):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if m.sum() < 10:
        return np.nan
    x, y, w = x[m], y[m], w[m]
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    vx = np.average((x - mx) ** 2, weights=w)
    vy = np.average((y - my) ** 2, weights=w)
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    args = ap.parse_args()
    cfg = PANELS[args.panel]

    est = pd.read_csv(os.path.join(RESULTS, f"phase3_{args.panel}_estimates.csv"))
    meta = pd.read_csv(os.path.join(RESULTS, f"phase2_{args.panel}_metadata.csv"))
    df = meta.merge(est, on="genetic_id", how="inner")
    print(f"Phase 4 — panel={args.panel}  joined {len(df):,} samples")

    # --- derived covariates ---------------------------------------------------
    df["data_type"] = df["genetic_id"].map(data_type)
    df["log_snps"] = np.log10(df["alpha_nSNP"].clip(lower=1))
    df["weight"] = 1.0 / df["alpha_SE"].clip(lower=1e-4) ** 2   # precision weight
    # parse contamination lower bound from ANGSD/hapConX CI strings like "[0.005,0.012]"
    def contam_lb(s):
        m = re.search(r"\[([0-9.]+)", str(s))
        return float(m.group(1)) if m else np.nan
    df["contam_lb"] = df["hapconx_contam"].map(contam_lb)
    df["contam_lb"] = df["contam_lb"].fillna(df["angsd_contam"].map(contam_lb))

    # high-confidence subset for individual-level work (group-level uses all)
    hi_snp = cfg["snp_lowpower"]
    df["high_conf"] = (df["alpha_nSNP"] >= hi_snp) & \
                      (~df["flags"].fillna("").str.contains("questionable"))

    # --- bias diagnostics -----------------------------------------------------
    R = []
    R.append(f"Phase 4 bias / normalisation report — panel {args.panel}")
    R.append(f"samples: {len(df):,}   high-confidence (>= {hi_snp:,} SNP, not questionable): "
             f"{int(df['high_conf'].sum()):,}")
    R.append("")
    R.append("Weighted mean alpha_Nea overall: "
             f"{wmean(df['alpha_Nea'].values, df['weight'].values)*100:.3f}%")
    R.append("")
    R.append("(1) alpha_Nea by sequencing data type (pseudo-haploid AG vs diploid DG/SG):")
    R.append(f"  {'type':6s} {'n':>6s} {'wmean_alpha%':>12s} {'med_SNP':>10s} {'med_SE%':>8s}")
    for dt, g in df.groupby("data_type"):
        if len(g) < 5:
            continue
        R.append(f"  {dt:6s} {len(g):6d} "
                 f"{wmean(g['alpha_Nea'].values, g['weight'].values)*100:12.3f} "
                 f"{g['alpha_nSNP'].median():10,.0f} {g['alpha_SE'].median()*100:8.3f}")
    R.append("")
    R.append("(2) directional bias: weighted corr(alpha_Nea, covariate)")
    R.append("    |r| small => estimate not systematically driven by that technical axis")
    for cov, lab in [("log_snps", "log10 usable SNPs"),
                     ("coverage", "mean coverage"),
                     ("damage", "DNA damage rate"),
                     ("contam_lb", "contamination lower bound"),
                     ("date_bp", "sample age (BP)")]:
        if cov in df:
            r_all = wcorr(df["alpha_Nea"].values, df[cov].values, df["weight"].values)
            hc = df[df["high_conf"]]
            r_hc = wcorr(hc["alpha_Nea"].values, hc[cov].values, hc["weight"].values)
            R.append(f"  {lab:28s} r_all={r_all:+.3f}   r_highconf={r_hc:+.3f}")
    R.append("")
    R.append("Interpretation: a non-trivial r in r_all that SHRINKS in r_highconf indicates a")
    R.append("low-coverage artifact (handled by the SNP floor + precision weighting). A covariate")
    R.append("with persistent |r| in the high-confidence set is a real nuisance axis and is added")
    R.append("to the Phase-5 expected-value model so Phase-6 residuals are net of it.")

    report = "\n".join(R)
    print("\n" + report)

    apath = os.path.join(RESULTS, f"phase4_{args.panel}_analysis.csv")
    rpath = os.path.join(RESULTS, f"phase4_{args.panel}_bias_report.txt")
    df.to_csv(apath, index=False)
    with open(rpath, "w", encoding="utf-8") as fh:
        fh.write(report + "\n")
    print(f"\nWrote:\n  {apath}\n  {rpath}")


if __name__ == "__main__":
    main()
