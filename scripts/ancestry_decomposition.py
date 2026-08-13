#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ancestry_decomposition.py — West-Eurasian source-population (Steppe/Yamnaya/WHG/
EHG/CHG/Anatolia_N/etc.) admixture survey, spanning early farmers through modern
populations, alongside each cohort's already-validated Neanderthal ancestry.

Complements the archaic (Neanderthal/Denisovan) side of this pipeline with the
human-ancestry side: for a chronological transect of ancient European cohorts
plus a cross-section of modern populations, this fits each as a mixture of
canonical source populations (archaic.ancestry / archaic.qpadm), competing
several candidate models (3-way farmer/steppe/HG through 5-source "deep" models)
and reporting both the classic unconstrained qpAdm fit and a simplex-constrained
("supervised admixture") fit that is always a valid mixture.

Outputs (results/ancestry/, reports/ancestry/):
  ancestry_models.csv     - every model x every target (both fits), long form
  ancestry_west3.csv      - the fixed 3-source model for every target (for a
                             consistent chronological stacked bar) + group-level
                             Neanderthal ancestry (archaic.profiles.group_archaic)
  ancestry_best.csv       - best-fitting model per target (archaic.decompose_best)
  fig_a1_stacked_bar.png  - chronological stacked-bar, west3 model
  fig_a2_steppe_time.png  - Steppe_Yamnaya % vs date, ancient cohorts
  fig_a3_model_fit.png    - qpAdm fit p-value per target x model (feasibility-coded)
  fig_a4_archaic_vs_steppe.png - Neanderthal % vs Steppe % (links both pipeline halves)

Run: PYTHONIOENCODING=utf-8 python scripts/ancestry_decomposition.py
"""
import os
import re
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from archaic.panel import Panel
from archaic.refs import PANELS
from archaic import stats as st, ancestry as anc, profiles as pf

PANEL = "1240k"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results", "ancestry")
FIGDIR = os.path.join(HERE, "reports", "ancestry")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

MODELS_TO_COMPETE = anc.MODELS          # west3, west4, deep4, deep5, hg3
REFERENCE_MODEL = "west3"               # fixed model for cross-target comparison

# ----------------------------------------------------------------- target cohorts
# Ancient chronological transect (verified group_id patterns against the local
# AADR v66.1 1240K release) plus a modern cross-section pulled directly from
# panel population labels. date_bp is filled in after resolution from the meta.
ANCIENT_TARGETS = {
    "Balkans_N": lambda g: ("greece" in g or "bulgaria" in g)
                          and bool(re.search(r"neolithic|_n(_|$)", g)),
    "CentralEuro_EN": lambda g: any(c in g for c in
                          ("germany", "austria", "czech", "hungary", "slovakia"))
                          and bool(re.search(r"_en(_|$)|_lbk|starcevo", g)),
    "CordedWare": lambda g: "corded" in g,
    "BellBeaker": lambda g: "beaker" in g,
    "CentralEuro_BA": lambda g: any(c in g for c in ("germany", "czech", "hungary", "austria"))
                          and ("unetice" in g or "_eba" in g or "_lba" in g
                               or "_mba" in g or "bronze" in g)
                          and not any(x in g for x in ("corded", "beaker", "yamnaya")),
    "Steppe_MLBA": lambda g: any(s in g for s in ("sintashta", "andronovo", "srubnaya")),
    "Italy_BA": lambda g: "italy" in g and any(s in g for s in
                          ("_eba", "_lba", "_mba", "bronze")),
    "Etruscan": lambda g: "etruscan" in g,
    "ImperialRoman": lambda g: "imperialroman" in g,
}
MODERN_TARGETS = {  # name -> panel population label
    "French": "French", "English": "English", "Sardinian": "Sardinian",
    "Spanish": "Spanish", "Italian_North": "Italian_North", "Basque": "Basque",
    "Russian": "Russian", "Orcadian": "Orcadian", "Finnish": "Finnish",
}


def all_source_names(models):
    names = set()
    for srcs in models.values():
        names.update(srcs)
    return names


def build_specs(models):
    """cohort name -> resolve_cohorts() spec, covering every source, every
    modern/near-source outgroup, and every target this run needs."""
    specs = {}
    for name in all_source_names(models):
        specs[name] = anc.SOURCES[name]["pred"]
    for name in anc.EXTRA_RIGHT:                      # near-source outgroups
        specs.setdefault(name, anc.SOURCES[name]["pred"])
    for name, pred in ANCIENT_TARGETS.items():
        specs[name] = pred
    for name, pop in MODERN_TARGETS.items():
        specs[name] = ("pop", pop)
    for name, spec in anc.BASE_RIGHT.items():          # distal outgroups
        specs[name] = spec
    return specs


def add_archaic_refs(specs, refs):
    """Add the archaic reference genomes/pops needed for group_archaic()."""
    specs["Altai"] = ("id", refs["Altai"]["ids"][0])
    specs["Vindija"] = ("id", refs["Vindija"]["ids"][0])
    specs["Denisova"] = ("id", refs["Denisova"]["ids"][0])
    specs["Chimp"] = ("id", refs["Chimp"]["ids"][0])
    specs["Mbuti"] = ("pop", refs["Mbuti"]["pops"][0])
    return specs


def median_date(meta, pred):
    gl = meta["group_id"].str.lower()
    sel = gl.map(pred) if callable(pred) else None
    if sel is None or sel.sum() == 0:
        return np.nan
    return float(meta.loc[sel, "date_bp"].median())


def main():
    t0 = time.time()
    panel = Panel(PANELS[PANEL]["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(HERE, "results", f"phase4_{PANEL}_analysis.csv"))
    print(f"panel loaded: {panel.n_snp} SNP, {len(panel.ind)} individuals "
          f"({time.time()-t0:.1f}s)")

    specs = build_specs(MODELS_TO_COMPETE)
    specs = add_archaic_refs(specs, PANELS[PANEL]["refs"])
    cols = anc.resolve_cohorts(panel, meta, specs, maxn=50, verbose=True)
    freq, info = anc.cohort_freqs(panel, cols)
    print(f"cohort frequencies ready ({time.time()-t0:.1f}s)")

    targets = list(ANCIENT_TARGETS) + list(MODERN_TARGETS)
    dates = {t: (median_date(meta, ANCIENT_TARGETS[t]) if t in ANCIENT_TARGETS else 0.0)
             for t in targets}

    # ---- 1. model competition, every target x every candidate model ----------
    model_rows = []
    best_rows = []
    west3_rows = []
    for tgt in targets:
        if info.get(tgt, {}).get("n", 0) < 2:
            print(f"  skip {tgt}: no individuals resolved")
            continue
        ranked = anc.decompose_best(freq, tgt, MODELS_TO_COMPETE, block, 50)
        for r in ranked:
            if not r["ok"]:
                continue
            for which in ("free", "constrained"):
                rr = r[which]
                row = dict(target=tgt, model=r["model"], fit=which,
                          date_bp=dates[tgt], n=info[tgt]["n"], n_snp=r["n_snp"],
                          p=rr["p"], feasible=rr["feasible"])
                for s, w, se in zip(rr["sources"], rr["weights"], rr["se"]):
                    row[f"{s}_pct"] = w * 100.0
                    row[f"{s}_se"] = se * 100.0
                model_rows.append(row)
        best = ranked[0]
        if best["ok"]:
            rr = best["constrained"]
            brow = dict(target=tgt, model=best["model"], date_bp=dates[tgt],
                       n=info[tgt]["n"], n_snp=best["n_snp"], p=rr["p"])
            for s, w, se in zip(rr["sources"], rr["weights"], rr["se"]):
                brow[f"{s}_pct"] = w * 100.0; brow[f"{s}_se"] = se * 100.0
            best_rows.append(brow)
            print(f"  {tgt:16s} best={best['model']:6s} p={rr['p']:.3f}  " +
                 "  ".join(f"{s}={w*100:.1f}%" for s, w in zip(rr["sources"], rr["weights"])))
        # fixed reference model, for the chronological stacked bar
        w3 = next((r for r in ranked if r["model"] == REFERENCE_MODEL and r["ok"]), None)
        if w3 is not None:
            rr = w3["constrained"]; rf = w3["free"]
            wrow = dict(target=tgt, date_bp=dates[tgt], n=info[tgt]["n"],
                       n_snp=w3["n_snp"], p_free=rf["p"], p_constrained=rr["p"],
                       free_feasible=rf["feasible"])
            for s, w, se in zip(rr["sources"], rr["weights"], rr["se"]):
                wrow[f"{s}_pct"] = w * 100.0; wrow[f"{s}_se"] = se * 100.0
            west3_rows.append(wrow)

    pd.DataFrame(model_rows).to_csv(os.path.join(RESULTS, "ancestry_models.csv"), index=False)
    best_df = pd.DataFrame(best_rows)
    best_df.to_csv(os.path.join(RESULTS, "ancestry_best.csv"), index=False)
    west3_df = pd.DataFrame(west3_rows)

    # ---- 2. attach group-level Neanderthal ancestry (mean-genome, low-noise) --
    arch = pf.group_archaic(freq, targets, block, 50,
                            A="Altai", O="Chimp", B="Mbuti", Ref="Vindija", Den="Denisova")
    arch = arch.rename(columns={"cohort": "target"})
    west3_df = west3_df.merge(arch[["target", "alpha_Nea", "alpha_SE"]], on="target", how="left")
    west3_df = west3_df.sort_values("date_bp", ascending=False).reset_index(drop=True)
    west3_df.to_csv(os.path.join(RESULTS, "ancestry_west3.csv"), index=False)
    print(f"\nWrote {len(model_rows)} model rows, {len(best_rows)} best-model rows, "
          f"{len(west3_df)} reference-model rows to {RESULTS}")

    # ------------------------------------------------------------- figures ----
    _fig_stacked_bar(west3_df)
    _fig_steppe_vs_time(west3_df)
    _fig_model_fit(pd.DataFrame(model_rows))
    _fig_archaic_vs_steppe(west3_df)
    print(f"Figures written to {FIGDIR} ({time.time()-t0:.1f}s total)")


# --------------------------------------------------------------------- figures
def _src_cols(df, sources):
    return [s for s in sources if f"{s}_pct" in df.columns]


def _fig_stacked_bar(df):
    if df.empty:
        return
    sources = anc.MODELS[REFERENCE_MODEL]
    cols = _src_cols(df, sources)
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(df))
    bottom = np.zeros(len(df))
    for s in cols:
        vals = df[f"{s}_pct"].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=s, color=anc.source_color(s), width=0.7)
        bottom += vals
    ax.set_xticks(x)
    labels = [f"{t}\n({'modern' if d == 0 else f'{d/1000:.1f} kya'})"
             for t, d in zip(df["target"], df["date_bp"])]
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"Ancestry proportion (%) — {REFERENCE_MODEL} model")
    ax.set_title("West-Eurasian source-population admixture, early farmers → modern\n"
                 "(constrained qpAdm, block-jackknife)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=9)
    ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_a1_stacked_bar.png"), dpi=150)
    plt.close(fig)


def _fig_steppe_vs_time(df):
    d = df[(df["date_bp"] > 0) & df["Steppe_Yamnaya_pct"].notna()].sort_values("date_bp")
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.errorbar(d["date_bp"] / 1000.0, d["Steppe_Yamnaya_pct"], yerr=d["Steppe_Yamnaya_se"],
               fmt="o-", color="#8e44ad", capsize=3, lw=1.5, ms=6)
    for _, row in d.iterrows():
        ax.annotate(row["target"], (row["date_bp"] / 1000.0, row["Steppe_Yamnaya_pct"]),
                   textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax.invert_xaxis()
    ax.set_xlabel("Date (kya BP)")
    ax.set_ylabel("Steppe/Yamnaya ancestry (%)")
    ax.set_title("Steppe-pastoralist ancestry through time (west3 model)")
    ax.axhline(0, color="gray", lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_a2_steppe_time.png"), dpi=150)
    plt.close(fig)


def _fig_model_fit(df):
    if df.empty:
        return
    d = df[df["fit"] == "constrained"].copy()
    targets = list(dict.fromkeys(d["target"]))
    models = list(MODELS_TO_COMPETE.keys())
    mat = np.full((len(targets), len(models)), np.nan)
    for i, t in enumerate(targets):
        for j, m in enumerate(models):
            sub = d[(d["target"] == t) & (d["model"] == m)]
            if len(sub):
                mat[i, j] = sub["p"].iloc[0]
    fig, ax = plt.subplots(figsize=(7, 0.35 * len(targets) + 2))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=0.2, aspect="auto")
    ax.set_xticks(range(len(models))); ax.set_xticklabels(models, rotation=30, ha="right")
    ax.set_yticks(range(len(targets))); ax.set_yticklabels(targets, fontsize=8)
    for i in range(len(targets)):
        for j in range(len(models)):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="qpAdm fit p-value (>0.05 = plausible)")
    ax.set_title("Model competition: fit p-value per target x model")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_a3_model_fit.png"), dpi=150)
    plt.close(fig)


def _fig_archaic_vs_steppe(df):
    d = df[df["Steppe_Yamnaya_pct"].notna() & df["alpha_Nea"].notna()]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, 6))
    sc = ax.scatter(d["Steppe_Yamnaya_pct"], d["alpha_Nea"] * 100.0,
                    c=d["date_bp"], cmap="viridis_r", s=70, edgecolor="k", linewidth=0.5)
    for _, row in d.iterrows():
        ax.annotate(row["target"], (row["Steppe_Yamnaya_pct"], row["alpha_Nea"] * 100.0),
                   textcoords="offset points", xytext=(5, 3), fontsize=7)
    fig.colorbar(sc, ax=ax, label="date (BP)")
    ax.set_xlabel("Steppe/Yamnaya ancestry (%, west3 model)")
    ax.set_ylabel("Neanderthal ancestry (%, group-level mean-genome)")
    ax.set_title("Human source ancestry vs. archaic ancestry\n(the two halves of this pipeline, side by side)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_a4_archaic_vs_steppe.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
