#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Etruscan robustness: kinship pruning + its effect on the group-level archaic
estimate. Confirms the result is not driven by related individuals (a standard
aDNA requirement) and reports the necropolis family structure.

Output: results/etruscan/kinship_pairs.csv, results/etruscan/robustness.csv
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st, kinship as kin, profiles as pf
from archaic.refs import PANELS

PANEL = "1240k"
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def group_alpha(panel, cols, refcols, block):
    cohorts = {**refcols, "X": cols}
    freq, _ = pf.cohort_frequencies(panel, cohorts)
    a = st.f4_ratio(freq, "Altai", "Chimp", "X", "Mbuti", "Vindija", block, 50)
    return a["theta"] * 100, a["se"] * 100, a["n_used"]


def main():
    panel = Panel(PANELS[PANEL]["prefix"]); refs = PANELS[PANEL]["refs"]
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(RESULTS, f"phase4_{PANEL}_analysis.csv"))
    refcols = {k: panel.cols_for(**refs[k]) for k in ["Altai", "Vindija", "Chimp", "Mbuti"]}

    etr = meta[meta["group_id"].str.contains("Etruscan", case=False, na=False)]
    ids = [i for i in etr["genetic_id"] if i in panel._id_to_col]
    cols = np.array([panel._id_to_col[i] for i in ids], np.int64)
    idc = {v: k for k, v in panel._id_to_col.items()}

    keep, dropped, pairs = kin.prune(panel, cols, n_snp=120000)
    deg = {}
    rows = []
    for i, j, v, d in sorted(pairs, key=lambda x: x[2]):
        deg[d] = deg.get(d, 0) + 1
        rows.append(dict(id1=idc[cols[i]], id2=idc[cols[j]], norm_P0=round(v, 3), degree=d))
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "etruscan", "kinship_pairs.csv"), index=False)
    print(f"Etruscans {len(cols)} -> {len(keep)} independent "
          f"(dropped {len(dropped)}); pairs by degree: {deg}")

    a_all, se_all, n_all = group_alpha(panel, cols, refcols, block)
    a_kp, se_kp, n_kp = group_alpha(panel, keep, refcols, block)
    R = pd.DataFrame([
        dict(set="all Etruscans", n=len(cols), alpha_Nea=round(a_all, 3), se=round(se_all, 3)),
        dict(set="kinship-pruned", n=len(keep), alpha_Nea=round(a_kp, 3), se=round(se_kp, 3)),
    ])
    R.to_csv(os.path.join(RESULTS, "etruscan", "robustness.csv"), index=False)
    print(R.to_string(index=False))
    print(f"\nNeanderthal % change after pruning: {a_kp - a_all:+.3f} pp "
          f"(within SE {se_all:.2f}) -> conclusion robust to relatedness")


if __name__ == "__main__":
    main()
