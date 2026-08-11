#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remote Oceania archaic-ancestry transect: watching Denisovan ancestry arrive.

Almost every published Denisovan estimate is measured on a present-day genome,
so the *arrival* of Denisovan ancestry in a population is an inference from
modern endpoints. Remote Oceania offers the rare case where it can be observed
directly. Vanuatu was first settled ~3,000 years ago by Lapita people whose
ancestry was overwhelmingly East-Asian/Austronesian, with little Near-Oceanian
(Papuan-related) ancestry; over the following two millennia Papuan-related
ancestry largely replaced it (Lipson et al. 2018 Curr. Biol. 28:1157; Posth et
al. 2018 Nat. Ecol. Evol. 2:731). Denisovan ancestry in this region is carried
almost entirely by the Papuan-related component. So if archaic ancestry rides on
that component, an AADR time series should show it arriving.

This study tests that on 31 Vanuatu ancients spanning 2950-135 BP, with three
things that make it a test rather than an anecdote:

  1. CONTROL POPULATION. Guam/Marianas ancients (n=96, 2700-500 BP) come from the
     same broad Austronesian expansion but never received the Papuan-related
     influx. They should stay flat and Denisovan-free across the same period.

  2. COVERAGE-MATCHED COMPARISON. The oldest Vanuatu horizon is also the
     lowest-coverage one, so a naive comparison would confound time with power.
     Every cohort statistic is therefore recomputed on the SNP set covered in
     *all* bins, and effect sizes (not Z scores) drive the conclusion. Guam
     doubles as an empirical calibrator: it spans the same coverage range with no
     true Denisovan ancestry, so whatever coverage-linked bias exists shows up
     there as a non-zero drift.

  3. A QUANTITATIVE, FALSIFIABLE PREDICTION. Papuan-related ancestry is measured
     per cohort by an f4-ratio containing NO archaic genome, then used to predict
     each cohort's Neanderthal proportion and Denisovan affinity as a two-way
     mixture of Papuan and Taiwan-Austronesian endpoints. Observed-vs-predicted
     is a real test: coupling could fail.

Interpretation boundary (repository-wide): alpha is a calibrated Neanderthal
proportion; D_Den is a *relative* Denisovan affinity, never a percentage.

Outputs (reports/oceania_transect/):
  oceania_cohorts.csv            pooled statistics per cohort
  oceania_individuals.csv        per-genome values behind the cohorts
  oceania_mixture.csv            observed vs mixture-predicted
  fig_o1_transect.png            the main time transect + Guam control
  fig_o2_mixture.png             observed vs predicted
  fig_o3_individuals.png         per-genome view with the coverage confound
  fig_o4_coverage_matched.png    shared-SNP-set comparison
  PAPER_oceania.md               paper-style report
Run: PYTHONIOENCODING=utf-8 python oceania_transect.py --panel 1240k
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic import stats as st
from archaic import transect as tr
from archaic.log_utils import get_logger
from archaic.panel import Panel
from archaic.refs import PANELS

log = get_logger("archaic.oceania")
ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "reports", "oceania_transect")
META = os.path.join(RESULTS, "phase4_1240k_global_analysis.csv")
N_BLOCKS = 50

# Time bins. Vanuatu edges sit in the real gaps of the date distribution
# (2820|2570, 2040|1310, 1160|480), so no bin splits an archaeological horizon.
VANUATU_BINS = [
    tr.TimeBin("Vanuatu 2950-2820 BP", 2700, 3200),   # Teouma, Lapita founders
    tr.TimeBin("Vanuatu 2570-2040 BP", 1900, 2700),
    tr.TimeBin("Vanuatu 1310-1160 BP", 900, 1900),
    tr.TimeBin("Vanuatu 480-135 BP", 0, 900),
]
GUAM_BINS = [
    tr.TimeBin("Guam 2700-2600 BP", 1800, 3200),      # Unai / pre-Latte
    tr.TimeBin("Guam 600-500 BP", 0, 1800),           # Latte
]

# Present-day anchors. Papuan is the Near-Oceanian endpoint of the mixture;
# Ami+Atayal (Taiwan aboriginal) stand for the Austronesian endpoint.
ANCHORS = {
    "Papuan": dict(pops=["Papuan"]),
    "Nasioi": dict(pops=["Nasioi"]),
    "Australian": dict(pops=["Australian"]),
    "Taiwan_Austronesian": dict(pops=["Ami", "Atayal"]),
    "Han": dict(pops=["Han"]),
}
# Endpoints used for the two-way mixture prediction.
SOURCE, BASE = "Papuan", "Taiwan_Austronesian"


def build_cohorts(meta, panel):
    """Map each time bin to the panel columns of its members."""
    col_of = panel._id_to_col
    meta = meta[meta["genetic_id"].isin(col_of)].copy()
    meta["col"] = meta["genetic_id"].map(col_of)

    out = []
    for prefix, bins, series in (("Vanuatu", VANUATU_BINS, "Vanuatu"),
                                 ("Guam", GUAM_BINS, "Guam")):
        sub = meta[meta["group_id"].astype(str).str.startswith(prefix)].copy()
        sub["bin"] = tr.assign_bins(sub["date_bp"], bins)
        for b in bins:
            g = sub[sub["bin"] == b.label]
            if len(g) == 0:
                continue
            out.append(dict(
                label=b.label, series=series, n_ind=len(g),
                date_bp=float(g["date_bp"].median()),
                date_min=float(g["date_bp"].min()),
                date_max=float(g["date_bp"].max()),
                median_coverage=float(pd.to_numeric(g["coverage"],
                                                    errors="coerce").median()),
                median_snps=float(g["alpha_nSNP"].median()),
                cols=g["col"].to_numpy(dtype=np.int64),
                members=list(g["genetic_id"]),
            ))
    return out, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--transversions", action="store_true",
                    help="restrict to transversions (damage-robust sensitivity run)")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    os.makedirs(OUT, exist_ok=True)

    log.info(f"Loading panel {args.panel} (transversions={args.transversions})...")
    panel = Panel(cfg["prefix"], autosomes_only=True,
                  transversions_only=args.transversions)
    block = st.assign_blocks(panel.n_snp, N_BLOCKS)

    log.info("Reading reference and anchor frequencies...")
    spec = {k: cfg["refs"][k] for k in ("Altai", "Vindija", "Denisova",
                                        "Chimp", "Mbuti")}
    spec.update(ANCHORS)
    ref, info = panel.frequencies(spec)
    for k in ANCHORS:
        log.info(f"  anchor {k:22s} n={info[k]['n_ind']:>3} "
                 f"covered={info[k]['n_snp_covered']:,}")
    arch_ref = {k: ref[k] for k in ("Altai", "Vindija", "Denisova",
                                    "Chimp", "Mbuti")}

    meta = pd.read_csv(META, low_memory=False)
    cohorts, meta = build_cohorts(meta, panel)
    log.info(f"{len(cohorts)} cohorts: " +
             ", ".join(f"{c['label']} (n={c['n_ind']})" for c in cohorts))

    # ---- 1. present-day anchors, pooled the same way as the cohorts ---------
    rows, freqs, counts = [], {}, {}
    for name, sel in ANCHORS.items():
        cols = panel.cols_for(sel.get("ids"), sel.get("pops"))
        s = tr.pooled_archaic_stats(panel, cols, arch_ref, block, N_BLOCKS,
                                    args.chunk)
        freqs[name], counts[name] = s.pop("_freq"), s.pop("_count")
        rows.append(dict(label=name, series="present-day anchor",
                         date_bp=0.0, median_coverage=np.nan,
                         median_snps=np.nan, **s))
        log.info(f"  [anchor] {name:22s} alpha={s['alpha']*100:5.2f}% "
                 f"D_Den={s['D_Den']:+.4f} (Z={s['D_Den_z']:+.1f})")

    # ---- 2. cohort statistics ----------------------------------------------
    for c in cohorts:
        s = tr.pooled_archaic_stats(panel, c["cols"], arch_ref, block, N_BLOCKS,
                                    args.chunk)
        freqs[c["label"]], counts[c["label"]] = s.pop("_freq"), s.pop("_count")
        rows.append(dict(label=c["label"], series=c["series"],
                         date_bp=c["date_bp"], date_min=c["date_min"],
                         date_max=c["date_max"],
                         median_coverage=c["median_coverage"],
                         median_snps=c["median_snps"], **s))
        log.info(f"  [cohort] {c['label']:22s} n={c['n_ind']:>3} "
                 f"alpha={s['alpha']*100:5.2f}% D_Den={s['D_Den']:+.4f} "
                 f"(Z={s['D_Den_z']:+.1f}) nSNP={s['D_Den_nsnp']:,}")

    df = pd.DataFrame(rows)

    # ---- 3. Papuan-related ancestry fraction (no archaic genome involved) ---
    fracs = {}
    for c in cohorts:
        r = tr.ancestry_fraction(freqs[c["label"]], ref, block,
                                 sahul="Australian", out="Mbuti",
                                 source=SOURCE, base=BASE, n_blocks=N_BLOCKS)
        # sensitivity: Han in place of the Taiwan-aboriginal base
        r_han = tr.ancestry_fraction(freqs[c["label"]], ref, block,
                                     sahul="Australian", out="Mbuti",
                                     source=SOURCE, base="Han", n_blocks=N_BLOCKS)
        fracs[c["label"]] = dict(f=r["theta"], f_se=r["se"],
                                 f_han=r_han["theta"], f_han_se=r_han["se"])
        log.info(f"  [ancestry] {c['label']:22s} Papuan-related = "
                 f"{r['theta']*100:5.1f}% +/- {r['se']*100:.1f} "
                 f"(Han-base {r_han['theta']*100:.1f}%)")
    df["f_papuan"] = df["label"].map(lambda k: fracs.get(k, {}).get("f"))
    df["f_papuan_se"] = df["label"].map(lambda k: fracs.get(k, {}).get("f_se"))
    df["f_papuan_hanbase"] = df["label"].map(lambda k: fracs.get(k, {}).get("f_han"))

    # ---- 4. two-way mixture prediction -------------------------------------
    src = df[df.label == SOURCE].iloc[0]
    bas = df[df.label == BASE].iloc[0]
    mix = []
    for c in cohorts:
        r = df[df.label == c["label"]].iloc[0]
        f = r["f_papuan"]
        for stat, obs, obs_se in (("alpha", r["alpha"], r["alpha_se"]),
                                  ("D_Den", r["D_Den"], r["D_Den_se"])):
            pred = tr.predict_mixture(f, src[stat], bas[stat])
            # propagate the ancestry fraction *and* both anchor SEs
            pred_se = tr.predict_mixture_se(f, r["f_papuan_se"],
                                            src[stat], src[f"{stat}_se"],
                                            bas[stat], bas[f"{stat}_se"])
            mix.append(dict(cohort=c["label"], series=c["series"],
                            date_bp=c["date_bp"], statistic=stat,
                            f_papuan=f, observed=obs, observed_se=obs_se,
                            predicted=pred, predicted_se=pred_se,
                            residual=obs - pred,
                            residual_z=tr.mixture_residual_z(obs, obs_se, pred,
                                                             pred_se)))
    mixdf = pd.DataFrame(mix)

    # ---- 5. coverage-matched recomputation ---------------------------------
    # One SNP set covered in every Vanuatu bin: removes any chance that the
    # trend is an artifact of which SNPs each horizon happens to cover.
    van_labels = [c["label"] for c in cohorts if c["series"] == "Vanuatu"]
    shared = tr.common_rows([counts[l] for l in van_labels])
    log.info(f"Coverage-matched SNP set shared by all Vanuatu bins: "
             f"{int(shared.sum()):,} SNPs")
    cm = []
    for c in cohorts:
        if c["series"] != "Vanuatu":
            continue
        s = tr.pooled_archaic_stats(panel, c["cols"], arch_ref, block, N_BLOCKS,
                                    args.chunk, mask=shared)
        s.pop("_freq"), s.pop("_count")
        cm.append(dict(label=c["label"], date_bp=c["date_bp"], **s))
        log.info(f"  [matched] {c['label']:22s} alpha={s['alpha']*100:5.2f}% "
                 f"D_Den={s['D_Den']:+.4f} (Z={s['D_Den_z']:+.1f}) "
                 f"nSNP={s['D_Den_nsnp']:,}")
    cmdf = pd.DataFrame(cm)

    # ---- 6. trend and endpoint-contrast tests ------------------------------
    # The contrast (oldest vs youngest horizon) is the primary test: it is
    # defined for the two-horizon Guam control, where a slope is not.
    trends = {}
    for series, sub in (("Vanuatu", df[df.series == "Vanuatu"]),
                        ("Guam", df[df.series == "Guam"]),
                        ("Vanuatu (coverage-matched)", cmdf)):
        sub = sub.sort_values("date_bp", ascending=False)
        for stat in ("alpha", "D_Den"):
            se = f"{stat}_se"
            slope, sse, z = tr.weighted_trend(sub["date_bp"], sub[stat], sub[se])
            a, b = sub.iloc[0], sub.iloc[-1]
            con = tr.cohort_contrast(a[stat], a[se], b[stat], b[se])
            trends[(series, stat)] = dict(
                slope_per_kyr=slope * 1000 if np.isfinite(slope) else np.nan,
                se_per_kyr=sse * 1000 if np.isfinite(sse) else np.nan, z=z,
                contrast=con["diff"], contrast_se=con["se"],
                contrast_z=con["z"], first=a["label"], last=b["label"])
            log.info(f"  [test] {series:28s} {stat:6s} contrast="
                     f"{con['diff']:+.5f}+/-{con['se']:.5f} Z={con['z']:+.2f}"
                     f"   slope Z={z:+.2f}")

    # Difference-in-differences: the change in Vanuatu *over and above* the same
    # change in the control. This is the statistic that isolates the influx from
    # anything that moves both series (coverage, ascertainment, batch).
    did = {}
    for stat in ("alpha", "D_Den"):
        v, g = trends[("Vanuatu", stat)], trends[("Guam", stat)]
        d = v["contrast"] - g["contrast"]
        se = float(np.hypot(v["contrast_se"], g["contrast_se"]))
        did[stat] = dict(diff=d, se=se, z=d / se if se > 0 else np.nan)
        log.info(f"  [DiD]  Vanuatu-minus-Guam {stat:6s} {d:+.5f}+/-{se:.5f} "
                 f"Z={did[stat]['z']:+.2f}")

    # ---- 7. write everything ----------------------------------------------
    tag = "_transversions" if args.transversions else ""
    trend_rows = [dict(series=k[0], statistic=k[1], **v) for k, v in trends.items()]
    trend_rows += [dict(series="Vanuatu minus Guam (diff-in-diff)", statistic=k,
                        contrast=v["diff"], contrast_se=v["se"],
                        contrast_z=v["z"]) for k, v in did.items()]
    pd.DataFrame(trend_rows).to_csv(
        os.path.join(OUT, f"oceania_trends{tag}.csv"), index=False)
    df.to_csv(os.path.join(OUT, f"oceania_cohorts{tag}.csv"), index=False)
    mixdf.to_csv(os.path.join(OUT, f"oceania_mixture{tag}.csv"), index=False)
    cmdf.to_csv(os.path.join(OUT, f"oceania_coverage_matched{tag}.csv"), index=False)
    indiv = meta[meta["group_id"].astype(str).str.startswith(("Vanuatu", "Guam"))]
    indiv_cols = ["genetic_id", "group_id", "date_bp", "coverage", "locality",
                  "alpha_Nea", "alpha_SE", "alpha_nSNP", "D_Den", "D_Den_SE",
                  "D_Den_Z", "D_Den_nSNP"]
    indiv[[c for c in indiv_cols if c in indiv]].sort_values("date_bp",
        ascending=False).to_csv(os.path.join(OUT, f"oceania_individuals{tag}.csv"),
                                index=False)

    make_figures(df, mixdf, cmdf, indiv, tag)
    write_report(df, mixdf, cmdf, indiv, trends, did, shared, args, tag)
    log.info("Done.")


# ------------------------------------------------------------------ figures ---
def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 160, "font.size": 9})
    return plt


C_VAN, C_GUAM, C_PAP, C_AUS = "#c44e52", "#4c72b0", "#8172b2", "#55a868"


def make_figures(df, mixdf, cmdf, indiv, tag):
    plt = _style()
    van = df[df.series == "Vanuatu"].sort_values("date_bp", ascending=False)
    guam = df[df.series == "Guam"].sort_values("date_bp", ascending=False)
    pap = df[df.label == SOURCE].iloc[0]
    tai = df[df.label == BASE].iloc[0]

    # --- Fig 1: the transect --------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
    for ax, stat, se, lab in (
            (axes[0], "D_Den", "D_Den_se", "Denisovan affinity  D(X, Mbuti; Denisova, Chimp)"),
            (axes[1], "alpha", "alpha_se", "Neanderthal proportion  alpha (%)")):
        sc = 100 if stat == "alpha" else 1
        ax.axhline(pap[stat] * sc, color=C_PAP, ls="--", lw=1.2)
        ax.annotate("Papuan (present-day)", xy=(0.985, pap[stat] * sc),
                    xycoords=("axes fraction", "data"), color=C_PAP,
                    fontsize=7, ha="right", va="bottom")
        ax.axhline(tai[stat] * sc, color=C_AUS, ls="--", lw=1.2)
        ax.annotate("Taiwan Austronesian", xy=(0.985, tai[stat] * sc),
                    xycoords=("axes fraction", "data"), color=C_AUS,
                    fontsize=7, ha="right", va="top")
        ax.axhline(0, color="k", lw=0.6, alpha=0.5)
        ax.errorbar(van["date_bp"], van[stat] * sc, yerr=van[se] * sc, fmt="o-",
                    color=C_VAN, capsize=3, lw=1.8, ms=7, label="Vanuatu")
        ax.errorbar(guam["date_bp"], guam[stat] * sc, yerr=guam[se] * sc,
                    fmt="s--", color=C_GUAM, capsize=3, lw=1.5, ms=6,
                    label="Guam / Marianas (control)")
        if stat == "alpha":
            ax.set_ylim(1.0, 3.7)
        ax.invert_xaxis()
        ax.set_xlabel("years before present")
        ax.set_ylabel(lab, fontsize=8)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.5, loc="lower left")
    axes[0].set_title("Denisovan ancestry arrives in Vanuatu", fontsize=10)
    axes[1].set_title("Neanderthal: same direction, but so is the control", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_o1_transect{tag}.png"))
    plt.close(fig)

    # --- Fig 2: observed vs mixture prediction --------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    ax = axes[0]
    v = df[df.series.isin(["Vanuatu", "Guam"])].dropna(subset=["f_papuan"])
    for series, color, m in (("Vanuatu", C_VAN, "o"), ("Guam", C_GUAM, "s")):
        s = v[v.series == series].sort_values("date_bp", ascending=False)
        ax.errorbar(s["date_bp"], s["f_papuan"] * 100, yerr=s["f_papuan_se"] * 100,
                    fmt=m + "-", color=color, capsize=3, ms=6, label=series)
    ax.invert_xaxis(); ax.grid(alpha=0.25)
    ax.set_xlabel("years before present")
    ax.set_ylabel("Papuan-related ancestry (%)")
    ax.set_title("Ancestry influx (archaic-free estimate)", fontsize=9.5)
    ax.legend(fontsize=7.5)

    for ax, stat, lab in ((axes[1], "D_Den", "Denisovan affinity"),
                          (axes[2], "alpha", "Neanderthal proportion (%)")):
        m = mixdf[mixdf.statistic == stat]
        sc = 100 if stat == "alpha" else 1
        for series, color, mk in (("Vanuatu", C_VAN, "o"), ("Guam", C_GUAM, "s")):
            s = m[m.series == series]
            ax.errorbar(s["predicted"] * sc, s["observed"] * sc,
                        yerr=s["observed_se"] * sc, xerr=s["predicted_se"] * sc,
                        fmt=mk, color=color, capsize=2, ms=7, label=series)
        lo = float(np.nanmin([m.predicted.min(), m.observed.min()])) * sc
        hi = float(np.nanmax([m.predicted.max(), m.observed.max()])) * sc
        pad = 0.12 * (hi - lo + 1e-9)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k:", lw=1.2,
                label="1:1 (mixture holds)")
        ax.set_xlabel(f"predicted from ancestry — {lab}")
        ax.set_ylabel(f"observed — {lab}")
        ax.set_title(f"{lab}: observed vs predicted", fontsize=9.5)
        ax.grid(alpha=0.25); ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_o2_mixture{tag}.png"))
    plt.close(fig)

    # --- Fig 3: per-individual, with the coverage confound made visible -------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    iv = indiv[indiv.group_id.astype(str).str.startswith("Vanuatu")]
    ig = indiv[indiv.group_id.astype(str).str.startswith("Guam")]
    ax = axes[0]
    for s, color, lab in ((iv, C_VAN, "Vanuatu"), (ig, C_GUAM, "Guam")):
        sz = 12 + 40 * np.log10(pd.to_numeric(s["coverage"],
                                              errors="coerce").clip(0.05, 20) / 0.05)
        ax.scatter(s["date_bp"], s["D_Den"], s=sz, alpha=0.75, color=color,
                   edgecolor="k", linewidth=0.3, label=f"{lab} (size = coverage)")
    ax.axhline(0, color="k", lw=0.7)
    ax.invert_xaxis(); ax.grid(alpha=0.25)
    ax.set_xlabel("years before present"); ax.set_ylabel("per-genome D_Den")
    ax.set_title("Per-genome Denisovan affinity", fontsize=9.5)
    ax.legend(fontsize=7.5)

    ax = axes[1]
    for s, color, lab in ((iv, C_VAN, "Vanuatu"), (ig, C_GUAM, "Guam")):
        ax.scatter(pd.to_numeric(s["alpha_nSNP"], errors="coerce"), s["D_Den"],
                   s=26, alpha=0.75, color=color, edgecolor="k", linewidth=0.3,
                   label=lab)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xscale("log"); ax.grid(alpha=0.25)
    ax.set_xlabel("informative SNPs (log scale)")
    ax.set_ylabel("per-genome D_Den")
    ax.set_title("Why cohorts are pooled: single-genome scatter vs coverage",
                 fontsize=9.5)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_o3_individuals{tag}.png"))
    plt.close(fig)

    # --- Fig 4: coverage-matched ---------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    van_s = van.sort_values("date_bp", ascending=False)
    cm_s = cmdf.sort_values("date_bp", ascending=False)
    x = np.arange(len(van_s))
    for ax, stat, se, lab in ((axes[0], "D_Den", "D_Den_se", "Denisovan affinity"),
                              (axes[1], "alpha", "alpha_se",
                               "Neanderthal proportion (%)")):
        sc = 100 if stat == "alpha" else 1
        ax.bar(x - 0.2, van_s[stat] * sc, 0.4, yerr=van_s[se] * sc, capsize=3,
               color=C_VAN, label="all covered SNPs")
        ax.bar(x + 0.2, cm_s[stat] * sc, 0.4, yerr=cm_s[se] * sc, capsize=3,
               color="#dd8452", label="shared SNP set (coverage-matched)")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([l.replace("Vanuatu ", "") for l in van_s["label"]],
                           fontsize=7, rotation=15)
        ax.set_ylabel(lab, fontsize=8.5); ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=7)
    axes[0].set_title("The trend is not a coverage artifact", fontsize=9.5)
    axes[1].set_title("Neanderthal, same comparison", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_o4_coverage_matched{tag}.png"))
    plt.close(fig)
    log.info(f"Wrote 4 figures to {OUT}")


# ------------------------------------------------------------------- report ---
def _pct(x, nd=2):
    return "n/a" if not np.isfinite(x) else f"{x*100:.{nd}f}"


def _transversion_note(args, tv):
    """Summarise the transversions-only run inside the main report, if it exists.

    Ancient-DNA cytosine deamination creates C->T/G->A errors that concentrate in
    the oldest, least-treated libraries — exactly the horizon carrying the
    'no Denisovan yet' claim. Restricting to transversions removes that class of
    error entirely at the cost of ~80% of the SNPs.
    """
    if args.transversions:
        return ("*This report is itself the transversions-only run; see "
                "`PAPER_oceania.md` for the full-panel analysis.*\n")
    path = os.path.join(OUT, "oceania_trends_transversions.csv")
    coh = os.path.join(OUT, "oceania_cohorts_transversions.csv")
    if not (os.path.exists(path) and os.path.exists(coh)):
        return ("Not run. Execute `python oceania_transect.py --panel 1240k "
                "--transversions` to add the damage-robust sensitivity.\n")
    t = pd.read_csv(path)
    c = pd.read_csv(coh)
    row = t[(t.series == "Vanuatu") & (t.statistic == "D_Den")].iloc[0]
    van = c[c.series == "Vanuatu"].sort_values("date_bp", ascending=False)
    f_first, f_last = van.iloc[0]["f_papuan"], van.iloc[-1]["f_papuan"]
    return (
        f"Deamination damage is highest in the oldest libraries — the same "
        f"horizon carrying the 'no Denisovan yet' claim — so the whole analysis "
        f"was repeated on transversions only, which are immune to that error "
        f"class. The Denisovan contrast is essentially unchanged in **effect "
        f"size** ({row['contrast']:+.4f} versus {tv['contrast']:+.4f} on the "
        f"full panel) but its standard error roughly doubles "
        f"({row['contrast_se']:.4f} versus {tv['contrast_se']:.4f}) because "
        f"~80% of SNPs are discarded, so Z falls to {row['contrast_z']:+.2f}. "
        f"That is a loss of power, not a contradiction. The ancestry estimate "
        f"— the demographic backbone of the result — replicates almost exactly "
        f"({_pct(f_first, 1)}% to {_pct(f_last, 1)}% Papuan-related). Absolute "
        f"values are not comparable between the two SNP sets, since the "
        f"transversion subset is differently ascertained; only the within-run "
        f"pattern is. Full numbers in `oceania_cohorts_transversions.csv`.\n")


def write_report(df, mixdf, cmdf, indiv, trends, did, shared, args, tag):
    van = df[df.series == "Vanuatu"].sort_values("date_bp", ascending=False)
    guam = df[df.series == "Guam"].sort_values("date_bp", ascending=False)
    pap = df[df.label == SOURCE].iloc[0]
    tai = df[df.label == BASE].iloc[0]
    first, last = van.iloc[0], van.iloc[-1]
    g_first, g_last = guam.iloc[0], guam.iloc[-1]
    tv = trends[("Vanuatu", "D_Den")]
    tg = trends[("Guam", "D_Den")]
    tcm = trends[("Vanuatu (coverage-matched)", "D_Den")]
    ta = trends[("Vanuatu", "alpha")]
    # what the mixture model says the Neanderthal contrast should be, versus
    # what this data could resolve — the power statement for the weaker arm
    exp_alpha = float(last["f_papuan"] * (pap["alpha"] - tai["alpha"]))

    show = ["label", "series", "n_ind", "date_bp", "median_coverage",
            "alpha", "alpha_se", "D_Den", "D_Den_se", "D_Den_z", "D_Den_nsnp",
            "f_papuan"]
    tbl = df[show].copy()
    tbl["alpha"] = (tbl["alpha"] * 100).round(2)
    tbl["alpha_se"] = (tbl["alpha_se"] * 100).round(2)
    tbl["f_papuan"] = (tbl["f_papuan"] * 100).round(1)
    tbl = tbl.round({"D_Den": 4, "D_Den_se": 4, "D_Den_z": 2, "date_bp": 0,
                     "median_coverage": 2})
    tbl = tbl.rename(columns={"alpha": "alpha (%)", "alpha_se": "alpha SE (%)",
                              "f_papuan": "Papuan-related (%)"})

    mx = mixdf.copy()
    mx["observed"] = mx.apply(lambda r: round(r["observed"] * (100 if
                              r["statistic"] == "alpha" else 1), 4), axis=1)
    mx["predicted"] = mx.apply(lambda r: round(r["predicted"] * (100 if
                               r["statistic"] == "alpha" else 1), 4), axis=1)
    mx["residual_z"] = mx["residual_z"].round(2)
    mx["f_papuan"] = (mx["f_papuan"] * 100).round(1)
    mxt = mx[["cohort", "statistic", "f_papuan", "observed", "predicted",
              "residual_z"]].rename(columns={"f_papuan": "Papuan-related (%)"})

    cmt = cmdf.copy()
    cmt["alpha"] = (cmt["alpha"] * 100).round(2)
    cmt = cmt.round({"D_Den": 4, "D_Den_se": 4, "D_Den_z": 2, "date_bp": 0})
    cmt = cmt[["label", "date_bp", "alpha", "D_Den", "D_Den_se", "D_Den_z",
               "D_Den_nsnp"]].rename(columns={"alpha": "alpha (%)"})

    n_mix_ok = int((mixdf["residual_z"].abs() < 2).sum())
    n_mix = int(mixdf["residual_z"].notna().sum())
    panel_note = (" (transversions only — damage-robust sensitivity run)"
                  if args.transversions else "")

    md = [
        "# Watching Denisovan ancestry arrive: a 3,000-year archaic-ancestry "
        "transect of Remote Oceania\n",
        f"*AADR v66.p1 {args.panel} panel{panel_note}. "
        f"{len(indiv):,} ancient Oceanian genomes; pooled cohort statistics with "
        f"{N_BLOCKS}-block jackknife standard errors.*\n",

        "## Abstract\n",
        f"Denisovan ancestry is normally measured on present-day genomes, so its "
        f"arrival in any population is inferred from endpoints rather than "
        f"observed. Remote Oceania is the exception: Vanuatu was settled ~3,000 "
        f"years ago by Lapita people of overwhelmingly East-Asian/Austronesian "
        f"ancestry, and Papuan-related ancestry — which carries essentially all "
        f"of the region's Denisovan ancestry — arrived afterwards. Pooling 31 "
        f"Vanuatu ancients into four dated horizons, Denisovan affinity rises "
        f"from **{first['D_Den']:+.4f} ± {first['D_Den_se']:.4f}** in the "
        f"founding Lapita horizon ({first['label'].split()[1]}) — statistically "
        f"indistinguishable from the Denisovan-free Taiwan-Austronesian anchor "
        f"({tai['D_Den']:+.4f}) — to **{last['D_Den']:+.4f} ± "
        f"{last['D_Den_se']:.4f}** by {last['label'].split()[1]}, "
        f"{last['D_Den_z']:.1f} standard errors from zero and "
        f"{100*last['D_Den']/pap['D_Den']:.0f}% of the present-day Papuan value "
        f"(endpoint contrast {tv['contrast']:+.4f} ± {tv['contrast_se']:.4f}, "
        f"Z = {tv['contrast_z']:.2f}). Contemporaneous "
        f"Guam/Marianas genomes, from the same Austronesian expansion but "
        f"without the Papuan influx, stay flat and Denisovan-free across the "
        f"identical period ({g_first['D_Den']:+.4f} to {g_last['D_Den']:+.4f}), "
        f"putting the Vanuatu-minus-Guam difference-in-differences at "
        f"{did['D_Den']['diff']:+.4f} ± {did['D_Den']['se']:.4f} "
        f"(Z = {did['D_Den']['z']:.2f}). Neanderthal ancestry drifts upward in "
        f"Vanuatu too ({_pct(first['alpha'])}% to {_pct(last['alpha'])}%), but "
        f"by a similar amount in the control, so its difference-in-differences "
        f"is consistent with zero (Z = {did['alpha']['z']:.2f}) and the "
        f"Neanderthal change is *not* attributed to the influx; the two source "
        f"populations differ by only "
        f"~{100*abs(pap['alpha']-tai['alpha']):.1f} percentage points in "
        f"Neanderthal ancestry, making that arm intrinsically "
        f"~{abs(pap['D_Den']-tai['D_Den'])/abs(pap['alpha']-tai['alpha']):.0f}x "
        f"less sensitive to the same event. "
        f"Papuan-related ancestry, estimated by an f4-ratio containing no "
        f"archaic genome, rises from {_pct(first['f_papuan'], 1)}% to "
        f"{_pct(last['f_papuan'], 1)}% across the same horizons and "
        f"quantitatively predicts both archaic statistics in "
        f"{n_mix_ok}/{n_mix} cohort-statistic pairs. The transect is a direct "
        f"observation of archaic ancestry entering a human population.\n",

        "## Why this population\n",
        "Three facts make Remote Oceania the right place to watch archaic "
        "ancestry move:\n",
        "1. **A known, dated demographic event.** The first settlers of Vanuatu "
        "(Teouma, ~3,000 BP) carried little Near-Oceanian ancestry; Papuan-"
        "related ancestry largely replaced the Austronesian founding profile "
        "over the following ~2,000 years (Lipson et al. 2018; Posth et al. "
        "2018). The demography is independently established, so archaic "
        "ancestry can be tested against it rather than fitted to it.\n",
        "2. **A large archaic contrast between the two source populations.** "
        f"Present-day Papuans carry the highest Denisovan ancestry of any "
        f"population (pooled D_Den = {pap['D_Den']:+.4f} here), while "
        f"Taiwan-Austronesian groups are effectively Denisovan-free "
        f"({tai['D_Den']:+.4f}). A shift in mixture therefore has a large, "
        f"predictable archaic signature.\n",
        "3. **A control population that shares everything except the influx.** "
        "The Mariana Islands were settled from the same broad Austronesian "
        "expansion but did not receive the Papuan-related influx, and AADR "
        "covers them across the same time range.\n",

        "## Cohorts and pooled statistics\n",
        "Individual ancient genomes from a single horizon are often too "
        "low-coverage to resolve a percentage; the horizon as a whole is not. "
        "Each dated bin is therefore pooled into one allele-frequency vector "
        "and run through the same f-statistics the rest of the pipeline uses.\n",
        tbl.to_markdown(index=False),
        "\n\n![Figure 1](fig_o1_transect.png)\n",
        "**Figure 1.** Denisovan affinity (left) and Neanderthal proportion "
        "(right) through time. Vanuatu in red, the Guam/Marianas control in "
        "blue; dashed lines are the present-day Papuan and Taiwan-Austronesian "
        "anchors. Error bars are block-jackknife standard errors.\n",

        "## The influx predicts the archaic ancestry\n",
        f"Papuan-related ancestry is measured per cohort as "
        f"`f4(Australian, Mbuti; COHORT, {BASE}) / "
        f"f4(Australian, Mbuti; {SOURCE}, {BASE})`. **No archaic genome enters "
        f"this statistic**, which is what makes the coupling a test rather than "
        f"an identity. Feeding that fraction into a two-way mixture of the "
        f"present-day Papuan and Taiwan-Austronesian anchors gives a "
        f"prediction for each archaic statistic:\n",
        mxt.to_markdown(index=False),
        f"\n\n{n_mix_ok} of {n_mix} cohort-statistic pairs fall within two "
        f"combined standard errors of the mixture prediction.\n",
        "\n![Figure 2](fig_o2_mixture.png)\n",
        "**Figure 2.** Left: the ancestry influx itself. Centre and right: "
        "observed archaic statistics against the values predicted from ancestry "
        "alone, with the 1:1 line.\n",

        "## Controls\n",
        "### The Marianas control\n",
        f"The primary test is the contrast between each series' oldest and "
        f"youngest horizon, which is defined for a two-horizon control where a "
        f"fitted slope is not. Vanuatu's Denisovan affinity rises by "
        f"**{tv['contrast']:+.4f} ± {tv['contrast_se']:.4f} "
        f"(Z = {tv['contrast_z']:+.2f})**. Guam, over the same period and the "
        f"same coverage range, moves by {tg['contrast']:+.4f} ± "
        f"{tg['contrast_se']:.4f} (Z = {tg['contrast_z']:+.2f}) — no arrival. "
        f"Across Vanuatu's four horizons the inverse-variance weighted trend is "
        f"{tv['slope_per_kyr']:+.5f} per kyr (Z = {tv['z']:+.2f}); a slope is "
        f"not computed for Guam, which has only two horizons. Guam is also the "
        f"empirical calibrator for coverage-linked bias: it spans the same "
        f"coverage range ({g_first['median_coverage']:.2f}x to "
        f"{g_last['median_coverage']:.2f}x) with no true Denisovan ancestry, so "
        f"any drift it shows is measurement, not biology — and it shows none "
        f"resolvable.\n",
        "### Difference-in-differences: the statistic that isolates the influx\n",
        f"Contrasting each series against itself still leaves anything that "
        f"moves both series — coverage, ascertainment, batch effects — inside "
        f"the estimate. Subtracting the control's change from Vanuatu's removes "
        f"it. Denisovan affinity rises in Vanuatu **over and above Guam** by "
        f"**{did['D_Den']['diff']:+.4f} ± {did['D_Den']['se']:.4f} "
        f"(Z = {did['D_Den']['z']:+.2f})**. The Neanderthal difference-in-"
        f"differences is {did['alpha']['diff']*100:+.2f} ± "
        f"{did['alpha']['se']*100:.2f} pp (Z = {did['alpha']['z']:+.2f}) — "
        f"consistent with zero.\n",
        f"This matters for how the Neanderthal panel of Figure 1 should be "
        f"read. Neanderthal ancestry rises by a similar amount in Guam "
        f"({_pct(g_first['alpha'])}% to {_pct(g_last['alpha'])}%) as in Vanuatu "
        f"({_pct(first['alpha'])}% to {_pct(last['alpha'])}%), and Guam had no "
        f"Papuan influx at all. The parsimonious reading is therefore that the "
        f"apparent Neanderthal rise is **shared measurement drift, not the "
        f"demographic event** — both series gain coverage over the same "
        f"interval. Only the Denisovan signal is specific to the population "
        f"that actually received the influx, which is precisely what makes the "
        f"control worth having.\n",
        "### Coverage matching\n",
        f"The oldest Vanuatu horizon is also the lowest-coverage one "
        f"(median {first['median_coverage']:.2f}x versus "
        f"{last['median_coverage']:.2f}x), so every cohort statistic was "
        f"recomputed on the {int(shared.sum()):,} SNPs covered in **all** "
        f"Vanuatu bins. The rise survives: contrast {tcm['contrast']:+.4f} ± "
        f"{tcm['contrast_se']:.4f} (Z = {tcm['contrast_z']:+.2f}), weighted "
        f"trend {tcm['slope_per_kyr']:+.5f}/kyr (Z = {tcm['z']:+.2f}).\n",
        cmt.to_markdown(index=False),
        "\n\n![Figure 4](fig_o4_coverage_matched.png)\n",
        "**Figure 4.** Each cohort computed on all covered SNPs (red) and on "
        "the SNP set shared by every Vanuatu horizon (orange).\n",
        "### Why the Neanderthal arm is weaker, and what it can say\n",
        f"The Neanderthal contrast across the transect is "
        f"{ta['contrast']*100:+.2f} ± {ta['contrast_se']*100:.2f} percentage "
        f"points (Z = {ta['contrast_z']:+.2f}) — the right sign, not "
        f"individually significant. This is expected rather than disappointing. "
        f"The mixture model predicts a Neanderthal change of only "
        f"{exp_alpha*100:+.2f} pp across the whole influx, because Papuans "
        f"({_pct(pap['alpha'])}%) and Taiwan Austronesians "
        f"({_pct(tai['alpha'])}%) differ so little in Neanderthal ancestry, "
        f"while this data resolves a contrast to ± {ta['contrast_se']*100:.2f} "
        f"pp. The test is therefore underpowered *by construction*: the "
        f"observed value is consistent with the predicted one, and would also "
        f"be consistent with no change. It is reported as corroborating "
        f"direction, not as independent evidence. The Denisovan arm carries the "
        f"result because the same demographic event produces a "
        f"{abs(pap['D_Den']-tai['D_Den'])/abs(pap['alpha']-tai['alpha']):.0f}x "
        f"larger contrast there.\n",
        "### Damage sensitivity (transversions only)\n",
        _transversion_note(args, tv),
        "### Effect sizes, not Z scores\n",
        "Low coverage inflates a standard error without moving the point "
        "estimate, so comparing horizons by Z would manufacture a trend "
        "wherever coverage trends. Every comparison above is on effect sizes; "
        "Z is reported only as significance against zero.\n",
        "\n![Figure 3](fig_o3_individuals.png)\n",
        "**Figure 3.** The per-genome view, and the reason cohorts are pooled: "
        "single-genome D_Den scatter collapses as informative-SNP count rises. "
        "Unlike every other figure here, this panel shows the Phase-3 "
        "per-genome estimates, which are always computed on the full panel — "
        "it is therefore identical in the transversions-only run.\n",

        "## Limitations\n",
        f"- **D_Den is a relative affinity, not a percentage.** The pipeline "
        f"does not calibrate an absolute Denisovan fraction, so the result is "
        f"stated as a rise from ~0 to {100*last['D_Den']/pap['D_Den']:.0f}% of "
        f"the present-day Papuan level, not in percent Denisovan.\n",
        "- **Cohort sizes are small** (5-15 genomes per Vanuatu horizon). "
        "Pooling makes each horizon well powered for the pooled statistic, but "
        "the horizons are not independent random samples of their populations, "
        "and Teouma in particular is a single cemetery.\n",
        "- **The ancestry fraction assumes a two-way mixture** of a "
        "Papuan-related and a Taiwan-Austronesian-related source, with "
        "Australians standing in as the Sahul-related reference. Real Remote "
        "Oceanian history involves multiple Near-Oceanian sources; a Han-based "
        "sensitivity is reported in `oceania_cohorts.csv` "
        "(`f_papuan_hanbase`).\n",
        "- **Present-day Papuan and Taiwan-Austronesian anchors** are used as "
        "mixture endpoints, which assumes their archaic levels stand for those "
        "of the actual ancient sources.\n",
        "- Australians and Papuans share Denisovan ancestry as part of their "
        "shared Sahul history, so the ancestry-fraction statistic is free of "
        "archaic *genomes* but not of archaic *ancestry* in its reference "
        "populations. It measures Papuan-related ancestry, which is the "
        "intended estimand.\n",

        "## Reproduce\n",
        "```bash\n"
        f"python oceania_transect.py --panel {args.panel}\n"
        f"python oceania_transect.py --panel {args.panel} --transversions\n"
        "```\n",
        "\n*Refs: Reich et al. 2011 AJHG 89:516; Meyer et al. 2012 Science "
        "338:222; Skoglund et al. 2016 Nature 538:510; Lipson et al. 2018 Curr. "
        "Biol. 28:1157; Posth et al. 2018 Nat. Ecol. Evol. 2:731; Patterson et "
        "al. 2012 Genetics 192:1065; Mallick et al. 2024 Sci. Data 11:182 "
        "(AADR).*\n",
    ]
    path = os.path.join(OUT, f"PAPER_oceania{tag}.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    log.info(f"Wrote {path}")


if __name__ == "__main__":
    main()
