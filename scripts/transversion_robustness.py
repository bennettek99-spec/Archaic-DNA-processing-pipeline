#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transversion-only robustness sweep for the headline group-level results.

The headline analyses (global AADR survey, Etruscan study) are run on the full
1240K SNP set. A known error class for ancient DNA — cytosine deamination — only
produces C->T / G->A changes, so a rerun restricted to transversions is immune
to it. If the headline effect sizes survive on transversions, they are not a
deamination artifact; if they shrink toward zero, they may be.

This script recomputes the two quantities the headline papers lean on — the
group-level Neanderthal f4-ratio (a calibrated %) and the relative Denisovan
D-statistic — for a small fixed set of key cohorts (Etruscan + the Italian
transect, an African negative control, and Papuan/Oceanian Denisovan-positive
anchors) on the full SNP set and on transversions only, and writes a side-by-side
table so the survival of each effect is explicit.

Reuses the same machinery as the studies: `archaic.panel.Panel` with
`transversions_only=True`, `archaic.profiles` pooled frequencies, and the
validated `archaic.stats` f4-ratio / D-statistic with block jackknife.

Output: results/transversion_robustness/group_stats_full_vs_tv.csv
Run:  PYTHONIOENCODING=utf-8 python scripts/transversion_robustness.py --panel 1240k
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import stats as st
from archaic import profiles as pf
from archaic.panel import Panel
from archaic.refs import PANELS

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "results", "transversion_robustness")
META = os.path.join(HERE, "results", "phase4_1240k_global_analysis.csv")
N_BLOCKS = 50

# Key cohorts: the Etruscan paper's Italian transect, the global survey's
# African negative control, and Oceanian Denisovan-positive anchors.
COHORT_SPECS = {
    "Etruscan": ("Etruscan_context", None),
    "Italy_Neolithic": (None, ("Neolithic/Copper",)),
    "Italy_BronzeAge": (None, ("Bronze Age",)),
    "Italy_Roman": (None, ("Roman",)),
    "Italy_Medieval": (None, ("Late Antique/Medieval",)),
    "African_control": ("all", "Africa"),
    "Papuan": ("pop", "Papuan"),
    "Nasioi": ("pop", "Nasioi"),
}

BINS = [("Neolithic/Copper", 4500, 7000), ("Bronze Age", 3200, 4500),
        ("Roman", 1700, 2300), ("Late Antique/Medieval", 800, 1700)]


def cohort_cols(panel, meta):
    """Resolve each cohort spec to panel column indices (capped, no pruning —
    group means over dozens of genomes are robust to a few relatives)."""
    meta = meta[meta["genetic_id"].isin(panel._id_to_col)].copy()
    out = {}
    rng = np.random.default_rng(0)
    for name, (kind, extra) in COHORT_SPECS.items():
        if kind == "Etruscan_context":
            m = meta["archaeological_cohort"].eq("Etruscan_context") \
                if "archaeological_cohort" in meta.columns else \
                meta["group_id"].astype(str).str.contains("Etruscan", na=False)
        elif kind == "pop":
            out[name] = panel.cols_for(pops=[extra])
            continue
        elif kind == "all":
            m = meta["continent"].eq(extra)
        else:
            # bin labels
            lab = extra[0]
            lo, hi = next((l, h) for lb, l, h in BINS if lb == lab)
            m = (meta["date_bp"] > lo) & (meta["date_bp"] <= hi) \
                & meta["continent"].eq("Eurasia") \
                & meta["country"].eq("Italy")
        cols = np.array([panel._id_to_col[i] for i in meta.loc[m, "genetic_id"]
                         if i in panel._id_to_col], dtype=np.int64)
        if len(cols) > 200:
            cols = np.sort(rng.choice(cols, 200, replace=False))
        out[name] = cols
    return out


def group_stats(panel, refs, cohort_cols):
    """Pooled Neanderthal f4-ratio and D_Den per cohort."""
    ref_cols = {name: panel.cols_for(**spec)
                for name, spec in refs.items()}
    freq, info = pf.cohort_frequencies(panel, {
        name: ref_cols[name]
        for name in ("Altai", "Vindija", "Denisova", "Chimp", "Mbuti")
    })
    ref_freq = {k: freq[k] for k in ("Altai", "Vindija", "Denisova", "Chimp", "Mbuti")}
    block = st.assign_blocks(panel.n_snp, N_BLOCKS)
    rows = []
    for name, cols in cohort_cols.items():
        if len(cols) == 0:
            continue
        cf, cinfo = pf.cohort_frequencies(panel, {name: cols})
        f = dict(ref_freq, POOL=cf[name])
        alpha = st.f4_ratio(f, "Altai", "Chimp", "POOL", "Mbuti", "Vindija", block, N_BLOCKS)
        dden = st.dstat(f, "POOL", "Mbuti", "Denisova", "Chimp", block, N_BLOCKS)
        rows.append(dict(
            cohort=name, n_ind=int(cinfo[name]["n"]),
            alpha_Nea=alpha["theta"], alpha_se=alpha["se"], alpha_nsnp=alpha["n_used"],
            D_Den=dden["theta"], D_Den_se=dden["se"], D_Den_z=dden["z"],
        ))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    args = ap.parse_args()

    cfg = PANELS[args.panel]
    os.makedirs(OUT, exist_ok=True)
    meta = pd.read_csv(META)

    # Cohort columns are resolved once on the full panel; the same individuals
    # are then read from the transversion-restricted panel.
    full_panel = Panel(cfg["prefix"], autosomes_only=True, transversions_only=False)
    cohorts = cohort_cols(full_panel, meta)

    full = group_stats(full_panel, cfg["refs"], cohorts)
    tv_panel = Panel(cfg["prefix"], autosomes_only=True, transversions_only=True)
    tv = group_stats(tv_panel, cfg["refs"], cohorts)

    merged = full.merge(tv, on=["cohort", "n_ind"], suffixes=("_full", "_tv"))
    merged = merged[[
        "cohort", "n_ind",
        "alpha_Nea_full", "alpha_se_full", "alpha_nsnp_full",
        "alpha_Nea_tv", "alpha_se_tv", "alpha_nsnp_tv",
        "D_Den_full", "D_Den_tv", "D_Den_z_full", "D_Den_z_tv",
    ]]
    merged.to_csv(os.path.join(OUT, "group_stats_full_vs_tv.csv"), index=False)

    pd.set_option("display.width", 200)
    print("=== Group-level archaic stats: full vs transversions-only ===")
    print(merged.round(4).to_string(index=False))
    print("\nWrote", os.path.join(OUT, "group_stats_full_vs_tv.csv"))


if __name__ == "__main__":
    main()
