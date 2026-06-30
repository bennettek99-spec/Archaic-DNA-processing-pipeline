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
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st, qpadm as qp, kinship as kin
from archaic.refs import PANELS

PANEL = "1240k"
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

# sources (by group_id substrings) and distal outgroups
SOURCES = {
    "Anatolia_N": lambda g: "turkey_n" in g or "anatolia_n" in g,
    "Steppe_Yamnaya": lambda g: "yamnaya" in g,
    "WHG": lambda g: any(s in g for s in (
        "loschbour", "villabruna", "bichon", "labrana", "iberia_mesolithic",
        "france_mesolithic", "england_mesolithic")),
}
OUTGROUPS = {  # pops (modern) or group_id substrings (ancient), distal to the sources
    "Mbuti": ("pop", "Mbuti"), "Han": ("pop", "Han"), "Papuan": ("pop", "Papuan"),
    "Karitiana": ("pop", "Karitiana"), "Onge": ("pop", "Onge"),
    "Iran_N": ("grp", "iran_ganjdareh_n"), "Natufian": ("grp", "israel_natufian"),
    "Ust_Ishim": ("id", "Ust_Ishim.DG"), "MA1": ("id", "MA1.SG"),
}
TARGETS = {
    "Etruscan": lambda g: "etruscan" in g,
    "Etruscan_Tuscany": lambda g: "etruscan" in g and "tuscany" in g,
    "Etruscan_Lazio": lambda g: "etruscan" in g and "lazio" in g,
    "Latin_Italic": lambda g: "latini" in g or ("lazio_ia" in g and "etruscan" not in g),
    "Imperial_Roman": lambda g: "imperialroman" in g,
    "Italy_BronzeAge": lambda g: "italy" in g and any(b in g for b in ("_ba", "_eba", "_mba", "_lba")),
}


def main():
    panel = Panel(PANELS[PANEL]["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(RESULTS, f"phase4_{PANEL}_analysis.csv"))
    gl = meta["group_id"].str.lower()

    def ids_for(pred):
        return meta.loc[gl.map(pred), "genetic_id"].tolist()

    rng = np.random.default_rng(0)
    MAXN = 60  # cap before O(N^2) kinship pruning; group freqs from ~60 are ample
    cohort_cols = {}
    for name, pred in {**SOURCES, **TARGETS}.items():
        ids = [i for i in ids_for(pred) if i in panel._id_to_col]
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
    for tgt in TARGETS:
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
    print("\nWrote results/etruscan/qpadm.csv")


if __name__ == "__main__":
    main()
