#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare this pipeline's estimates to ADMIXTOOLS 2 on the same populations.

Run `tools/admixtools_concordance.R` first (needs R + the admixtools package) to
produce results/admixtools_results.csv, then run this. It computes the same
statistics with the pure-Python estimator and reports the differences — the
concordance check for the methods paper. If the ADMIXTOOLS CSV is absent, it still
computes and saves this pipeline's side and tells you to run the R script.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.panel import Panel
from archaic import stats as st
from archaic.refs import PANELS

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
TESTS = ["French", "Han", "Sardinian", "Papuan", "Yoruba"]


def main():
    cfg = PANELS["1240k"]
    panel = Panel(cfg["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    spec = {k: cfg["refs"][k] for k in ["Altai", "Vindija", "Chimp", "Mbuti", "Yoruba"]}
    for t in TESTS:
        spec[t] = dict(pops=[t])
    freq, info = panel.frequencies(spec)

    mine = []
    for X in TESTS:
        a = st.f4_ratio(freq, "Altai", "Chimp", X, "Mbuti", "Vindija", block, 50)
        dn = st.dstat(freq, X, "Mbuti", "Altai", "Chimp", block, 50)
        mine.append(dict(stat="alpha_Nea", pops=X, est=a["theta"], se=a["se"]))
        mine.append(dict(stat="D_Nea", pops=X, est=dn["theta"], se=dn["se"]))
    dd = st.dstat(freq, "French", "Han", "Altai", "Yoruba", block, 50)
    mine.append(dict(stat="D_diff_FrenchHan", pops="French;Han", est=dd["theta"], se=dd["se"]))
    M = pd.DataFrame(mine)
    M.to_csv(os.path.join(RESULTS, "pipeline_fstats_for_concordance.csv"), index=False)

    apath = os.path.join(RESULTS, "admixtools_results.csv")
    if not os.path.exists(apath):
        print("Saved this pipeline's f-statistics to results/pipeline_fstats_for_concordance.csv")
        print("Run tools/admixtools_concordance.R (needs R + admixtools) to produce")
        print("results/admixtools_results.csv, then re-run this script to see the comparison.")
        return
    A = pd.read_csv(apath)
    j = M.merge(A, on=["stat", "pops"], suffixes=("_mine", "_admixtools"))
    j["abs_diff"] = (j["est_mine"] - j["est_admixtools"]).abs()
    print("Concordance: this pipeline vs ADMIXTOOLS 2")
    print(j[["stat", "pops", "est_mine", "est_admixtools", "abs_diff"]].to_string(index=False))
    print(f"\nmax |difference| = {j['abs_diff'].max():.5f}")
    j.to_csv(os.path.join(RESULTS, "concordance_admixtools.csv"), index=False)


if __name__ == "__main__":
    main()
