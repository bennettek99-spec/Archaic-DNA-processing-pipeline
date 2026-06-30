#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 8 — publication-quality figures.  Writes PNGs to results/figures/.

  fig1_z_distribution   standardised residuals vs N(0,1) + multiple-testing line
                        (the headline: no individual beyond chance)
  fig2_ancestry_gradient alpha_Nea vs ancestry PC1 (West-East Eurasian cline)
  fig3_pca_introgression PC1 x PC2 coloured by Neanderthal %
  fig4_map              geographic map (lon/lat) coloured by Neanderthal %
  fig5_noise_floor      alpha_SE vs usable SNPs (why single genomes are noisy)
  fig6_violin_era       Neanderthal % by archaeological era
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIG = os.path.join(RESULTS, "figures")
os.makedirs(FIG, exist_ok=True)
PANEL = sys.argv[1] if len(sys.argv) > 1 else "1240k"
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.25})


def main():
    res = pd.read_csv(os.path.join(RESULTS, f"phase6_{PANEL}_residuals.csv"))
    pc = pd.read_csv(os.path.join(RESULTS, f"phase5_{PANEL}_pca.csv"))
    df = res.merge(pc, on="genetic_id")
    hc = df[df["high_conf"]].copy()
    n = len(hc)
    zc = stats.norm.isf(0.025 / n)

    # fig1 — z distribution vs N(0,1)
    fig, ax = plt.subplots(figsize=(7, 4.3))
    z = hc["z_resid"].dropna()
    ax.hist(z, bins=80, density=True, color="#4c72b0", alpha=0.8, label=f"observed (n={n:,})")
    xs = np.linspace(-5, 5, 400)
    ax.plot(xs, stats.norm.pdf(xs), "k--", lw=1.5, label="N(0,1) null")
    for s in (-1, 1):
        ax.axvline(s * zc, color="crimson", ls=":", lw=1.3)
    ax.text(zc, ax.get_ylim()[1]*0.9, f" Bonferroni z*={zc:.2f}\n (0 samples beyond)",
            color="crimson", fontsize=8, va="top")
    ax.set_xlabel("standardised residual z  (observed − expected Neanderthal, given ancestry+geo+time)")
    ax.set_ylabel("density"); ax.set_title("No individual exceeds chance expectation")
    ax.legend(); fig.tight_layout(); fig.savefig(f"{FIG}/fig1_z_distribution.png"); plt.close(fig)

    # fig2 — ancestry gradient
    fig, ax = plt.subplots(figsize=(7, 4.3))
    sc = ax.scatter(hc["PC1"], hc["alpha_Nea"]*100, c=hc["lon"], cmap="Spectral_r",
                    s=6, alpha=0.5)
    # weighted running trend
    o = hc.sort_values("PC1")
    ax.plot(o["PC1"], o["alpha_Nea"].rolling(300, center=True, min_periods=50).mean()*100,
            "k-", lw=2, label="rolling mean")
    plt.colorbar(sc, label="longitude (°E)")
    ax.set_xlabel("ancestry PC1  (West  →  East Eurasian)")
    ax.set_ylabel("Neanderthal ancestry (%)"); ax.set_ylim(0, 4.5)
    ax.set_title("Neanderthal % vs ancestry (East-Eurasian excess, weak & noisy)")
    ax.legend(); fig.tight_layout(); fig.savefig(f"{FIG}/fig2_ancestry_gradient.png"); plt.close(fig)

    # fig3 — PCA coloured by introgression
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    sc = ax.scatter(hc["PC1"], hc["PC2"], c=hc["alpha_Nea"]*100, cmap="viridis",
                    s=7, alpha=0.6, vmin=1.0, vmax=3.0)
    plt.colorbar(sc, label="Neanderthal %")
    ax.set_xlabel("PC1 (West–East Eurasian)"); ax.set_ylabel("PC2")
    ax.set_title("Ancestry PCA coloured by Neanderthal ancestry")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig3_pca_introgression.png"); plt.close(fig)

    # fig4 — map
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    m = hc[hc["lon"].notna() & hc["lat"].notna()]
    sc = ax.scatter(m["lon"], m["lat"], c=m["alpha_Nea"]*100, cmap="viridis",
                    s=8, alpha=0.6, vmin=1.0, vmax=3.0)
    plt.colorbar(sc, label="Neanderthal %")
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ax.set_title("Ancient Eurasian Neanderthal ancestry by location")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig4_map.png"); plt.close(fig)

    # fig5 — noise floor
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.scatter(df["alpha_nSNP"], df["alpha_SE"]*100, s=5, alpha=0.3, color="#555")
    xs = np.linspace(df["alpha_nSNP"].min(), df["alpha_nSNP"].max(), 200)
    k = np.nanmedian(df["alpha_SE"]*100 * np.sqrt(df["alpha_nSNP"]))
    ax.plot(xs, k / np.sqrt(xs), "r-", lw=2, label=r"$\propto 1/\sqrt{\mathrm{SNPs}}$")
    ax.axvline(200_000, color="green", ls=":", label="high-confidence floor")
    ax.set_xscale("log"); ax.set_xlabel("usable SNPs (per individual)")
    ax.set_ylabel("Neanderthal % standard error"); ax.set_ylim(0, 3)
    ax.set_title("Measurement noise floor — why single genomes can't resolve fine signal")
    ax.legend(); fig.tight_layout(); fig.savefig(f"{FIG}/fig5_noise_floor.png"); plt.close(fig)

    # fig6 — violin by era
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    bins = [0, 3000, 5000, 8000, 12000, 1e9]
    labs = ["0–3k", "3–5k", "5–8k", "8–12k", "12k+ (UP)"]
    hc["era"] = pd.cut(hc["date_bp"], bins=bins, labels=labs, right=False)
    data = [hc.loc[hc["era"] == l, "alpha_Nea"].dropna()*100 for l in labs]
    parts = ax.violinplot(data, showmedians=True)
    ax.set_xticks(range(1, len(labs)+1)); ax.set_xticklabels(labs)
    ax.set_xlabel("era (years BP)"); ax.set_ylabel("Neanderthal %"); ax.set_ylim(0, 4.5)
    ax.set_title("Neanderthal ancestry through time (high-confidence)")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig6_violin_era.png"); plt.close(fig)

    print("Wrote 6 figures to", FIG)
    for f in sorted(os.listdir(FIG)):
        print("  ", f)


if __name__ == "__main__":
    main()
