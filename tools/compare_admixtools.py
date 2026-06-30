#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Concordance: this pipeline vs ADMIXTOOLS 2 (Maier et al. 2023), computed on the
same PLINK export (export_plink.py). Run tools/admixtools_concordance.R first.
Produces results/concordance_summary.csv and figures/fig_concordance.png.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
FIG = os.path.join(RESULTS, "figures")


def main():
    mine = pd.read_csv(os.path.join(RESULTS, "pipeline_fstats_for_concordance.csv"))
    af = os.path.join(RESULTS, "admixtools_fstats.csv")
    if not os.path.exists(af):
        print("Run tools/admixtools_concordance.R first."); return
    at = pd.read_csv(af)
    j = mine.merge(at, on=["stat", "pops"], suffixes=("_mine", "_admix"))

    rows = []
    for s in ["alpha_Nea", "D_Nea"]:
        sub = j[j["stat"] == s]
        r = np.corrcoef(sub["est_mine"], sub["est_admix"])[0, 1]
        rows.append(dict(quantity=s, n=len(sub), pearson_r=round(r, 4),
                         scale_mine_over_admix=round((sub["est_mine"] / sub["est_admix"]).mean(), 2)))
    # qpAdm
    mq = pd.read_csv(os.path.join(RESULTS, "etruscan", "qpadm.csv"))
    aq = pd.read_csv(os.path.join(RESULTS, "admixtools_qpadm.csv"))
    mq["target"] = mq["target"].replace({"Imperial_Roman": "ImperialRoman",
                                          "Italy_BronzeAge": "ItalyBA", "Latin_Italic": "Latin"})
    qj = mq.merge(aq, on="target")
    qdiff = []
    for _, x in qj.iterrows():
        qdiff += [abs(x["Anatolia_N_pct"] - x["Anatolia_N"] * 100),
                  abs(x["Steppe_Yamnaya_pct"] - x["Yamnaya"] * 100),
                  abs(x["WHG_pct"] - x["WHG"] * 100)]
    rows.append(dict(quantity="qpAdm_weights", n=len(qdiff),
                     pearson_r=round(np.corrcoef(
                         qj[["Anatolia_N_pct", "Steppe_Yamnaya_pct", "WHG_pct"]].values.ravel(),
                         qj[["Anatolia_N", "Yamnaya", "WHG"]].values.ravel() * 100)[0, 1], 4),
                     scale_mine_over_admix=f"mean abs diff {np.mean(qdiff):.1f} pp"))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(RESULTS, "concordance_summary.csv"), index=False)
    print(S.to_string(index=False))

    # figure: scatter mine vs admixtools (f-stats z-scored per quantity to overlay) + qpAdm
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for s, c, m in [("alpha_Nea", "#1f77b4", "o"), ("D_Nea", "#2ca02c", "s")]:
        sub = j[j["stat"] == s]
        ax[0].scatter(sub["est_admix"], sub["est_mine"], c=c, marker=m, s=45,
                      label=f"{s} (r={np.corrcoef(sub.est_mine,sub.est_admix)[0,1]:.3f})")
    lims = [min(j.est_admix.min(), j.est_mine.min()), max(j.est_admix.max(), j.est_mine.max())]
    ax[0].plot(lims, lims, "k--", lw=0.8, alpha=0.5)
    ax[0].set_xlabel("ADMIXTOOLS 2"); ax[0].set_ylabel("this pipeline")
    ax[0].set_title("f-statistics concordance"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)
    comps = [("Anatolia_N_pct", "Anatolia_N"), ("Steppe_Yamnaya_pct", "Yamnaya"), ("WHG_pct", "WHG")]
    for mc, ac in comps:
        ax[1].scatter(qj[ac] * 100, qj[mc], s=55, label=ac)
    ax[1].plot([0, 90], [0, 90], "k--", lw=0.8, alpha=0.5)
    ax[1].set_xlabel("ADMIXTOOLS 2 qpAdm %"); ax[1].set_ylabel("this pipeline qpAdm %")
    ax[1].set_title("qpAdm ancestry concordance"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_concordance.png"), dpi=150); plt.close(fig)
    print("\nWrote results/concordance_summary.csv and figures/fig_concordance.png")


if __name__ == "__main__":
    main()
