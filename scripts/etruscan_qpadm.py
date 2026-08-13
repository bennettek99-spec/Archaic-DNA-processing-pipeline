#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qpAdm ancestry modelling of the Etruscans and neighbours.

Models each target as a 3-way mixture of the canonical West-Eurasian sources —
Anatolian Neolithic farmers, Steppe pastoralists (Yamnaya), and Western
Hunter-Gatherers (WHG) — relative to a set of distal outgroups, with block-
jackknife SEs and a GLS chi-square fit p-value. Target cohorts are kinship-pruned
first (archaic.kinship) so relatives do not bias the frequencies.

Output: results/etruscan/qpadm.csv  (+ printed table).
Also writes all candidate qpAdm models and qpWave rank tests:
  results/etruscan/qpadm_models.csv
  results/etruscan/qpwave.csv
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.panel import Panel
from archaic import stats as st, qpadm as qp, kinship as kin, cohort_rules
from archaic.refs import PANELS

PANEL = "1240k"
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

# sources (by group_id substrings) and distal outgroups
SOURCES = {
    "Anatolia_N": lambda g: "turkey_n" in g or "anatolia_n" in g,
    "Steppe_Yamnaya": lambda g: "yamnaya" in g,
    "WHG": lambda g: any(s in g for s in (
        "loschbour", "villabruna", "bichon", "labrana", "iberia_mesolithic",
        "france_mesolithic", "england_mesolithic")),
}
MODELS = {
    "west3": ["Anatolia_N", "Steppe_Yamnaya", "WHG"],
    "farmer_steppe": ["Anatolia_N", "Steppe_Yamnaya"],
    "farmer_whg": ["Anatolia_N", "WHG"],
    "steppe_whg": ["Steppe_Yamnaya", "WHG"],
}
OUTGROUPS = {  # pops (modern) or group_id substrings (ancient), distal to the sources
    "Mbuti": ("pop", "Mbuti"), "Han": ("pop", "Han"), "Papuan": ("pop", "Papuan"),
    "Karitiana": ("pop", "Karitiana"), "Onge": ("pop", "Onge"),
    "Iran_N": ("grp", "iran_ganjdareh_n"), "Natufian": ("grp", "israel_natufian"),
    "Ust_Ishim": ("id", "Ust_Ishim.DG"), "MA1": ("id", "MA1.SG"),
}
TARGET_COHORTS = {
    "Etruscan": "Etruscan_context",
    "Latin_Italic": "Latin_context",
    "Imperial_Roman": "Imperial_Roman_context",
    "Italy_BronzeAge": "Preceding_Bronze_Age_Italy",
    "Early_Medieval_Italy": "Early_Medieval_Italy",
}


def main():
    panel = Panel(PANELS[PANEL]["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(RESULTS, f"phase4_{PANEL}_analysis.csv"))
    if "archaeological_cohort" not in meta.columns:
        meta = cohort_rules.apply_cohort_rules(meta)
    if "population_test_keep" not in meta.columns:
        meta = cohort_rules.add_population_test_keep(meta)
    gl = meta["group_id"].str.lower()

    def ids_for(pred):
        return meta.loc[gl.map(pred), "genetic_id"].tolist()

    def ids_for_cohort(cohort):
        m = meta["archaeological_cohort"].eq(cohort) & \
            meta["population_test_keep"].fillna(True).astype(bool)
        return meta.loc[m, "genetic_id"].tolist()

    rng = np.random.default_rng(0)
    MAXN = 60  # cap before O(N^2) kinship pruning; group freqs from ~60 are ample
    cohort_cols = {}
    for name, pred in SOURCES.items():
        ids = [i for i in ids_for(pred) if i in panel._id_to_col]
        cols = np.array([panel._id_to_col[i] for i in ids], dtype=np.int64)
        if len(cols) > MAXN:
            cols = np.sort(rng.choice(cols, MAXN, replace=False))
        if 4 <= len(cols):
            keep, dropped, _ = kin.prune(panel, cols, n_snp=40000)
            cols = keep
        cohort_cols[name] = cols
    for name, cohort in TARGET_COHORTS.items():
        ids = [i for i in ids_for_cohort(cohort) if i in panel._id_to_col]
        cols = np.array([panel._id_to_col[i] for i in ids], dtype=np.int64)
        if len(cols) > MAXN:
            cols = np.sort(rng.choice(cols, MAXN, replace=False))
        if 4 <= len(cols):
            keep, dropped, _ = kin.prune(panel, cols, n_snp=40000)
            cols = keep
        cohort_cols[name] = cols
    # outgroups
    for name, (kind, val) in OUTGROUPS.items():
        if kind == "pop":
            cols = panel.cols_for(pops=[val])
        elif kind == "id":
            cols = panel.cols_for(ids=[val])
        else:
            cols = np.array([panel._id_to_col[i] for i in ids_for(lambda g, v=val: v in g)
                             if i in panel._id_to_col], dtype=np.int64)
        cohort_cols[name] = cols

    # frequencies for everything (one read), via a tiny helper on explicit cols
    from archaic import profiles as pf
    freq, info = pf.cohort_frequencies(panel, {k: v for k, v in cohort_cols.items() if len(v)})
    avail_out = [o for o in OUTGROUPS if info.get(o, {}).get("n", 0) >= 1]
    avail_src = [s for s in SOURCES if info.get(s, {}).get("n", 0) >= 2]
    print("sources:", {s: info[s]["n"] for s in avail_src})
    print("outgroups:", {o: info[o]["n"] for o in avail_out})
    print()

    rows = []
    for tgt in TARGET_COHORTS:
        if info.get(tgt, {}).get("n", 0) < 2:
            continue
        r = qp.qpadm(freq, tgt, avail_src, avail_out, block, 50)
        d = dict(target=tgt, n=info[tgt]["n"], n_snp=r["n_snp"],
                 chi2=r["chi2"], dof=r["dof"], p=r["p"])
        for s, w, se in zip(r["sources"], r["weights"], r["se"]):
            d[f"{s}_pct"] = w * 100; d[f"{s}_se"] = se * 100
        rows.append(d)
        ws = "  ".join(f"{s}={w*100:5.1f}%±{se*100:.1f}"
                       for s, w, se in zip(r["sources"], r["weights"], r["se"]))
        print(f"{tgt:18s} n={info[tgt]['n']:3d}  {ws}   p={r['p']:.3f}  "
              f"({'plausible' if r['p'] > 0.05 else 'rejected' if r['p'] == r['p'] else 'n/a'})")
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "etruscan", "qpadm.csv"), index=False)

    model_rows = []
    wave_rows = []
    for tgt in TARGET_COHORTS:
        if info.get(tgt, {}).get("n", 0) < 2:
            continue
        for model, srcs0 in MODELS.items():
            srcs = [s for s in srcs0 if s in avail_src]
            if len(srcs) < 2 or len(avail_out) <= len(srcs):
                continue
            r = qp.qpadm(freq, tgt, srcs, avail_out, block, 50)
            status = "plausible" if r["p"] > 0.05 and r["feasible"] else "rejected"
            d = dict(target=tgt, model=model, n=info[tgt]["n"], n_snp=r["n_snp"],
                     chi2=r["chi2"], dof=r["dof"], p=r["p"],
                     feasible=r["feasible"], status=status)
            for s, w, se in zip(r["sources"], r["weights"], r["se"]):
                d[f"{s}_pct"] = w * 100
                d[f"{s}_se"] = se * 100
            model_rows.append(d)
            for wr in qp.qpwave(freq, [tgt] + srcs, avail_out, block, 50):
                wave_rows.append({**wr, "target": tgt, "model": model,
                                  "lefts": ";".join([tgt] + srcs),
                                  "rights": ";".join(avail_out)})
    pd.DataFrame(model_rows).to_csv(os.path.join(RESULTS, "etruscan", "qpadm_models.csv"), index=False)
    pd.DataFrame(wave_rows).to_csv(os.path.join(RESULTS, "etruscan", "qpwave.csv"), index=False)
    print("\nWrote results/etruscan/qpadm.csv, qpadm_models.csv, qpwave.csv")


if __name__ == "__main__":
    main()
