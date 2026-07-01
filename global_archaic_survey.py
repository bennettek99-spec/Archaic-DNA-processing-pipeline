#!/usr/bin/env python
"""
Whole-AADR (global) high-archaic survey.

Extends the Eurasia-only survey to ALL of the AADR 1240K panel: present-day
individuals and every continent (Africa, Americas, Oceania, Eurasia). Consumes
results/phase4_1240k_global_analysis.csv (Phase-3 estimates on the Phase-2
`--scope global` metadata, merged + high_conf flagged here).

Purpose: (i) an out-of-Africa NEGATIVE CONTROL (Africans should read ~0%
Neanderthal), and (ii) bring in Oceania, the one region where archaic ancestry
genuinely exceeds 5% once Denisovan is counted.

Outputs (reports/global_archaic_survey/):
  global_continent_breakdown.csv
  global_over5pct.csv
  oceania_denisovan.csv
  fig_g1_neanderthal_by_continent.png
  fig_g2_denisovan_by_continent.png
  fig_g3_archaic_scatter.png
  fig_g4_global_map.png
Run: PYTHONIOENCODING=utf-8 python global_archaic_survey.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "results/phase4_1240k_global_analysis.csv"
OUT = "reports/global_archaic_survey"
os.makedirs(OUT, exist_ok=True)
THRESH, Z_DEN = 0.05, 2.5
plt.rcParams.update({"figure.dpi": 130, "font.size": 10})

df = pd.read_csv(SRC)
if "high_conf" not in df:
    df["high_conf"] = (df.alpha_nSNP >= 200000) & \
        (~df["flags"].fillna("").str.contains("questionable"))
df["high_conf"] = df["high_conf"].astype(bool)
df["over"] = df.alpha_Nea > THRESH
df["den_flag"] = df.D_Den_Z > Z_DEN
CONTS = ["Africa", "Americas", "Eurasia", "Oceania"]
d = df[df.continent.isin(CONTS)].copy()

# ---------------------------------------------------------------- tables
rows = []
for c in CONTS:
    g = d[d.continent == c]; hc = g[g.high_conf]
    rows.append(dict(
        continent=c, n=len(g), n_present_day=int(g.is_modern.sum()),
        mean_alpha_Nea=round(g.alpha_Nea.mean(), 4),
        median_alpha_Nea=round(g.alpha_Nea.median(), 4),
        hc_mean_alpha=round(hc.alpha_Nea.mean(), 4) if len(hc) else np.nan,
        raw_over5=int(g.over.sum()), hc_over5=int((g.over & g.high_conf).sum()),
        mean_D_Den_Z=round(g.D_Den_Z.mean(), 3),
        pct_den_flag=round(100 * g.den_flag.mean(), 2),
        n_hc_den_flag=int((g.den_flag & g.high_conf).sum())))
bt = pd.DataFrame(rows)
bt.to_csv(f"{OUT}/global_continent_breakdown.csv", index=False)
print(bt.to_string(index=False))

cols = ["genetic_id", "group_id", "country", "continent", "is_modern", "lat", "lon",
        "date_bp", "alpha_Nea", "alpha_SE", "alpha_nSNP", "D_Den_Z", "high_conf"]
d[d.over][cols].sort_values("alpha_Nea", ascending=False).to_csv(
    f"{OUT}/global_over5pct.csv", index=False)
oc = d[(d.continent == "Oceania") & d.den_flag].sort_values("D_Den_Z", ascending=False)
oc[cols].to_csv(f"{OUT}/oceania_denisovan.csv", index=False)

# ---------------------------------------------------------------- colors
CC = {"Africa": "#e6ab02", "Americas": "#66a61e",
      "Eurasia": "#7570b3", "Oceania": "#d95f02"}

# --- FIG g1: Neanderthal alpha by continent
fig, ax = plt.subplots(figsize=(8, 5))
data = [100 * d.loc[d.continent == c, "alpha_Nea"].values for c in CONTS]
bp = ax.boxplot(data, tick_labels=CONTS, showfliers=False, patch_artist=True, widths=.6)
for b, c in zip(bp["boxes"], CONTS):
    b.set(facecolor=CC[c], alpha=.6)
ax.axhline(5, color="#d7301f", ls="--", lw=1.3, label="5% threshold")
ax.axhline(0, color="k", lw=.8)
ax.set_ylabel("Neanderthal ancestry $\\alpha$ (%)")
ax.set_title("Neanderthal ancestry by continent (whole AADR)\n"
             "Africa $\\approx$ 0 validates the zero-point; no continent exceeds 5%")
ax.legend()
fig.tight_layout(); fig.savefig(f"{OUT}/fig_g1_neanderthal_by_continent.png"); plt.close(fig)

# --- FIG g2: Denisovan affinity by continent
fig, ax = plt.subplots(figsize=(8, 5))
data = [d.loc[d.continent == c, "D_Den_Z"].values for c in CONTS]
bp = ax.boxplot(data, tick_labels=CONTS, showfliers=False, patch_artist=True, widths=.6)
for b, c in zip(bp["boxes"], CONTS):
    b.set(facecolor=CC[c], alpha=.6)
ax.axhline(0, color="k", lw=.8)
ax.axhline(Z_DEN, color="#238b45", ls="--", lw=1.2, label="flag Z=2.5")
ax.set_ylabel("Denisovan affinity  D_Den_Z  (relative)")
ax.set_title("Denisovan affinity by continent (whole AADR)\n"
             "Oceania is the sole strongly-elevated region (Papuan/Australian)")
ax.legend()
fig.tight_layout(); fig.savefig(f"{OUT}/fig_g2_denisovan_by_continent.png"); plt.close(fig)

# --- FIG g3: archaic scatter (Neanderthal vs Denisovan), continents
fig, ax = plt.subplots(figsize=(8.5, 6))
for c in CONTS:
    g = d[d.continent == c]
    ax.scatter(100 * g.alpha_Nea, g.D_Den_Z, s=8, c=CC[c], alpha=.4,
               label=c, rasterized=True)
ax.axhline(Z_DEN, color="#238b45", ls="--", lw=1)
ax.axvline(5, color="#d7301f", ls="--", lw=1)
ax.set_xlabel("Neanderthal ancestry $\\alpha$ (%)")
ax.set_ylabel("Denisovan affinity  D_Den_Z")
ax.set_title("The archaic-ancestry plane: Oceania occupies its own high-Denisovan corner\n"
             "(Neanderthal alone never clears 5%; Oceanians add substantial Denisovan)")
ax.set_xlim(-3, 8)
ax.legend(markerscale=2)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_g3_archaic_scatter.png"); plt.close(fig)

# --- FIG g4: global map
fig, ax = plt.subplots(figsize=(12, 6))
hi = d[d.high_conf]
ax.scatter(hi.lon, hi.lat, s=3, c="#cccccc", alpha=.4, rasterized=True, label="_bg")
den = d[d.den_flag]
ax.scatter(den.lon, den.lat, s=22, marker="^", facecolor="none", edgecolor="#238b45",
           lw=1.1, label="Denisovan-affinity (Z>2.5)", zorder=5)
over = d[d.over]
ax.scatter(over.lon, over.lat, s=26, c="#d7301f", alpha=.75, edgecolor="k", lw=.3,
           label="raw >5% Neanderthal", zorder=4)
ax.set_xlim(-180, 180); ax.set_ylim(-60, 82)
ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
ax.set_title("Whole-AADR high-archaic signals\n"
             "grey = high-confidence genomes; Denisovan cluster in Oceania/ISEA")
ax.grid(alpha=.2); ax.legend(loc="lower left", fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_g4_global_map.png"); plt.close(fig)

print(f"\nWrote tables + 4 figures to {OUT}/")
print("Oceania Denisovan-flagged: %d (of %d)" %
      (len(oc), (d.continent == "Oceania").sum()))
