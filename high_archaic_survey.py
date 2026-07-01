#!/usr/bin/env python
"""
High-archaic-ancestry survey (>5% threshold) over the AADR 1240K panel.

Consumes the completed pipeline output `results/phase4_1240k_analysis.csv`
(per-individual Neanderthal proportion alpha_Nea + SE, Denisovan affinity
D_Den + Z, full geographic/temporal/QC metadata) and answers the specific
question: which ancient/modern Eurasians carry estimated ARCHAIC ancestry
above an ABSOLUTE 5% threshold, where do they cluster, which archaic source
contributed, and how does this change over time?

Distinct from the pipeline's Phase 6 (which asks who is archaic-elevated
*given* ancestry/geography/age -> a near-null). This is the simpler absolute
query, and its honest answer hinges on coverage: raw >5% hits are dominated
by low-SNP noise.

Outputs (all under reports/high_archaic_survey/):
  raw_over5pct.csv          - every sample clearing 5% Neanderthal (raw)
  regional_breakdown.csv    - counts by macro-region (raw vs high-conf vs sig)
  denisovan_candidates.csv  - top Denisovan-affinity samples
  notable_early_up.csv      - biologically meaningful early Upper Palaeolithic set
  fig1_coverage_threshold.png
  fig2_map.png
  fig3_timeline.png
  fig4_denisovan.png
Run: PYTHONIOENCODING=utf-8 python high_archaic_survey.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- config
SRC = "results/phase4_1240k_analysis.csv"
OUT = "reports/high_archaic_survey"
os.makedirs(OUT, exist_ok=True)
THRESH = 0.05          # 5% absolute Neanderthal threshold
HC_SNP = 200_000       # high-confidence SNP floor (matches Phase 4)
Z_DEN = 2.5            # Denisovan relative-affinity flag (D_Den_Z)

plt.rcParams.update({"figure.dpi": 130, "font.size": 10})

# ---------------------------------------------------------------- load
df = pd.read_csv(SRC)
df["high_conf"] = df["high_conf"].astype(bool)
# 95% CI lower bound on the Neanderthal proportion
df["alpha_lcb"] = df["alpha_Nea"] - 1.96 * df["alpha_SE"]
df["alpha_ucb"] = df["alpha_Nea"] + 1.96 * df["alpha_SE"]


# ---------------------------------------------------------------- region map
def macro_region(lat, lon):
    """Coarse lon/lat -> macro-region (reproducible band classifier)."""
    if pd.isna(lat) or pd.isna(lon):
        return "Unknown"
    # Island SE Asia / Oceania (Denisovan-relevant)
    if lon >= 95 and lat < 12:
        return "Island SE Asia / Oceania"
    # Africa (mostly pre-excluded; catch residuals)
    if lat < 12 and -20 <= lon <= 52:
        return "Africa"
    if 55 <= lat and lon >= 55:
        return "Siberia"
    if -12 <= lon <= 45 and 36 <= lat <= 72:
        return "Europe"
    if 45 < lon <= 62 and 33 <= lat <= 60:
        return "Caucasus / W. Central Asia"
    if 20 <= lon < 63 and 12 <= lat < 42:
        return "Near East / Anatolia"
    if 62 <= lon <= 92 and 5 <= lat < 38:
        return "South Asia"
    if 62 <= lon < 100 and 38 <= lat < 55:
        return "Central Asia / Steppe"
    if lon >= 92 and 18 <= lat < 60:
        return "East Asia"
    return "Other Eurasia"


df["region"] = [macro_region(a, b) for a, b in zip(df["lat"], df["lon"])]

# ---------------------------------------------------------------- masks
hc = df["high_conf"]
over = df["alpha_Nea"] > THRESH                    # raw >5%
sig = df["alpha_lcb"] > THRESH                     # 95% CI lower bound >5%
den = df["D_Den_Z"] > Z_DEN                        # Denisovan-affinity flag

print("=" * 64)
print("HIGH-ARCHAIC SURVEY  (n=%d QC-pass Eurasian genomes, 1240K)" % len(df))
print("=" * 64)
print("Neanderthal alpha: mean %.3f  median %.3f  (SE median %.3f)"
      % (df.alpha_Nea.mean(), df.alpha_Nea.median(), df.alpha_SE.median()))
print("raw alpha>5%%              : %d" % over.sum())
print("  ...of which high-conf   : %d" % (over & hc).sum())
print("  ...median nSNP of raw>5%%: %d (vs %d overall)"
      % (df.loc[over, "alpha_nSNP"].median(), df.alpha_nSNP.median()))
print("95%% CI lower bound >5%%    : %d  (high-conf: %d)"
      % (sig.sum(), (sig & hc).sum()))
print("Denisovan D_Den_Z>%.1f     : %d  (high-conf: %d)"
      % (Z_DEN, den.sum(), (den & hc).sum()))

# ---------------------------------------------------------------- CSV 1: raw >5%
cols = ["genetic_id", "group_id", "locality", "country", "region", "lat", "lon",
        "date_bp", "full_date", "alpha_Nea", "alpha_SE", "alpha_lcb", "alpha_ucb",
        "alpha_nSNP", "D_Nea_Z", "D_Den_Z", "coverage", "snps_hit",
        "assessment", "data_type", "high_conf"]
raw = df.loc[over, cols].sort_values("alpha_Nea", ascending=False)
raw.to_csv(f"{OUT}/raw_over5pct.csv", index=False)

# ---------------------------------------------------------------- CSV 2: regional
def region_table(mask):
    return df.loc[mask].groupby("region").size()

reg = pd.DataFrame({
    "n_total": df.groupby("region").size(),
    "n_raw_over5": region_table(over),
    "n_highconf_over5": region_table(over & hc),
    "n_CI_sig_over5": region_table(sig),
    "n_denisovan_flag": region_table(den),
}).fillna(0).astype(int).sort_values("n_raw_over5", ascending=False)
reg["pct_raw_over5"] = (100 * reg.n_raw_over5 / reg.n_total).round(2)
reg.to_csv(f"{OUT}/regional_breakdown.csv")
print("\n--- regional breakdown (raw >5%) ---")
print(reg[["n_total", "n_raw_over5", "n_highconf_over5", "n_denisovan_flag"]].to_string())

# ---------------------------------------------------------------- CSV 3: Denisovan
den_tbl = df.loc[den, ["genetic_id", "group_id", "country", "region", "lat", "lon",
                       "date_bp", "D_Den", "D_Den_Z", "D_Den_nSNP", "alpha_Nea",
                       "alpha_nSNP", "high_conf"]].sort_values("D_Den_Z", ascending=False)
den_tbl.to_csv(f"{OUT}/denisovan_candidates.csv", index=False)

# ---------------------------------------------------------------- CSV 4: notable early-UP
# Individuals >40kBP OR named early-UP sites: the biologically real high-archaic story.
names = ["Oase", "Muierii", "Cioclovina", "Kostenki", "Ust_Ishim", "Tianyuan",
         "Bacho", "BB7", "CC7", "F6-", "BK-1", "ZKU", "Zlaty", "Sunghir",
         "GoyetQ", "Pestera", "Bichon"]
pat = "|".join(names)
early = df[df.genetic_id.str.contains(pat, case=False, na=False)
           | df.group_id.str.contains("IUP|UP|Aurignac|Gravett", case=False, na=False)].copy()
early = early[early.date_bp >= 25000]
early = early[["genetic_id", "group_id", "country", "region", "date_bp",
               "alpha_Nea", "alpha_SE", "alpha_lcb", "alpha_ucb", "alpha_nSNP",
               "D_Nea_Z", "high_conf"]].sort_values("date_bp", ascending=False)
early.to_csv(f"{OUT}/notable_early_up.csv", index=False)

# ================================================================ FIGURES
C_HC = "#2c7fb5"     # high-conf
C_LOW = "#c8c8c8"    # low-conf
C_HIT = "#d7301f"    # >5% hit
C_SIG = "#54278f"    # CI-significant
C_DEN = "#238b45"    # Denisovan

# ---- FIG 1: alpha vs coverage -- the "it's all low-coverage" headline
fig, ax = plt.subplots(figsize=(8.2, 5.4))
lo = df[~hc]; hi = df[hc]
ax.scatter(lo.alpha_nSNP, lo.alpha_Nea, s=6, c=C_LOW, alpha=.45,
           label="Low-confidence (<200k SNP)", rasterized=True)
ax.scatter(hi.alpha_nSNP, hi.alpha_Nea, s=6, c=C_HC, alpha=.45,
           label="High-confidence (>=200k SNP)", rasterized=True)
ax.axhline(THRESH, color=C_HIT, ls="--", lw=1.4, label="5% threshold")
ax.axhline(df.alpha_Nea.median(), color="k", ls=":", lw=1,
           label="panel median (%.1f%%)" % (100 * df.alpha_Nea.median()))
# annotate Oase1
o = df[df.genetic_id.str.startswith("Oase1")]
if len(o):
    r = o.iloc[0]
    ax.annotate("Oase1 (9.8%, recent\nNeanderthal ancestor)",
                (r.alpha_nSNP, r.alpha_Nea), (60000, 0.082),
                arrowprops=dict(arrowstyle="->", color=C_SIG), color=C_SIG, fontsize=9)
ax.set_xlabel("SNPs used in f4-ratio estimate (coverage proxy)")
ax.set_ylabel("Estimated Neanderthal ancestry  $\\alpha$")
ax.set_title("Every raw >5%% Neanderthal call is low-coverage\n"
             "0 of %d high-confidence genomes exceed 5%%" % hc.sum())
ax.set_ylim(-0.06, 0.11)
ax.legend(loc="upper right", framealpha=.9, fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_coverage_threshold.png"); plt.close(fig)

# ---- FIG 2: map
fig, ax = plt.subplots(figsize=(11, 6))
# faint backdrop = all high-conf samples => draws sampled Eurasia
ax.scatter(hi.lon, hi.lat, s=3, c=C_LOW, alpha=.35, label="_bg", rasterized=True)
# raw >5% hits, sized by coverage (small = noisy)
h = df[over]
sizes = 25 + (h.alpha_nSNP / h.alpha_nSNP.max()) * 260
ax.scatter(h.lon, h.lat, s=sizes, c=C_HIT, alpha=.7, edgecolor="k", lw=.3,
           label="Raw >5% Neanderthal (size $\\propto$ coverage)")
# CI-significant (Oase1)
s = df[sig]
ax.scatter(s.lon, s.lat, s=180, marker="*", c=C_SIG, edgecolor="k", lw=.5, zorder=5,
           label="95% CI lower bound >5%")
# Denisovan cluster
d = df[den]
ax.scatter(d.lon, d.lat, s=70, marker="^", facecolor="none",
           edgecolor=C_DEN, lw=1.6, zorder=6, label="Denisovan-affinity (Z>%.1f)" % Z_DEN)
ax.set_xlim(-15, 165); ax.set_ylim(-15, 78)
ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
ax.set_title("Geography of high-archaic signals across sampled Eurasia\n"
             "grey = all high-confidence ancient genomes (backdrop)")
ax.grid(alpha=.2)
ax.legend(loc="lower left", framealpha=.92, fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_map.png"); plt.close(fig)

# ---- FIG 3: timeline (high-conf Neanderthal over time + early-UP annotated)
fig, ax = plt.subplots(figsize=(9.5, 5.6))
hh = hi[(hi.date_bp > 0) & (hi.date_bp < 50000)]
ax.scatter(hh.date_bp, 100 * hh.alpha_Nea, s=7, c=C_LOW, alpha=.4, label="High-conf genomes")
# rolling median in 2.5 kyr bins
bins = np.arange(0, 50001, 2500)
hh = hh.copy(); hh["bin"] = pd.cut(hh.date_bp, bins)
gm = hh.groupby("bin", observed=True).alpha_Nea.agg(["median", "count"])
centers = [iv.mid for iv in gm.index]
ax.plot(centers, 100 * gm["median"], "-o", c=C_HC, lw=2, ms=4,
        label="binned median (2.5 kyr)")
# annotate key early-UP high-conf individuals
for gid, lab in [("Ust_Ishim.DG", "Ust'-Ishim"), ("Tianyuan.AG.BY.AA", "Tianyuan"),
                 ("Kostenki14_d.AG.BY.AA", "Kostenki14"), ("Muierii2.SG", "Muierii2"),
                 ("BB7-240.AG.BY.AA", "Bacho Kiro"), ("F6-620.AG.BY.AA", "Bacho Kiro F6")]:
    r = df[df.genetic_id == gid]
    if len(r):
        r = r.iloc[0]
        ax.annotate(lab, (r.date_bp, 100 * r.alpha_Nea),
                    fontsize=7.5, color=C_SIG,
                    xytext=(4, 4), textcoords="offset points")
        ax.scatter([r.date_bp], [100 * r.alpha_Nea], s=26, c=C_SIG, zorder=5)
ax.axhline(100 * df.alpha_Nea.median(), color="k", ls=":", lw=1,
           label="panel median %.1f%%" % (100 * df.alpha_Nea.median()))
ax.invert_xaxis()
ax.set_xlabel("Age (years before present)")
ax.set_ylabel("Neanderthal ancestry $\\alpha$ (%)")
ax.set_title("Neanderthal ancestry over time (high-confidence genomes)\n"
             "elevated in oldest Upper Palaeolithic, settling to ~2%")
ax.legend(loc="upper left", fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_timeline.png"); plt.close(fig)

# ---- FIG 4: Denisovan affinity by region
fig, ax = plt.subplots(figsize=(9, 5.4))
order = (df.groupby("region").D_Den_Z.median().sort_values().index.tolist())
data = [df.loc[df.region == rn, "D_Den_Z"].dropna().values for rn in order]
bp = ax.boxplot(data, vert=False, tick_labels=order, showfliers=False,
                patch_artist=True, widths=.6)
for b in bp["boxes"]:
    b.set(facecolor="#f0f0f0")
# overlay flagged
for i, rn in enumerate(order, start=1):
    dd = den_tbl[den_tbl.region == rn]
    ax.scatter(dd.D_Den_Z, [i] * len(dd), c=C_DEN, s=28, zorder=5)
ax.axvline(0, color="k", lw=.8)
ax.axvline(Z_DEN, color=C_DEN, ls="--", lw=1, label="flag Z=%.1f" % Z_DEN)
ax.set_xlabel("Denisovan affinity  D_Den_Z  (relative; single Denisova reference)")
ax.set_title("Denisovan affinity by region\n"
             "coherent Island SE Asia elevation; green = flagged samples")
ax.legend(loc="lower right", fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_denisovan.png"); plt.close(fig)

print("\nWrote CSVs + 4 figures to %s/" % OUT)
print("early-UP notable set (n=%d):" % len(early))
print(early.head(12).to_string(index=False))
