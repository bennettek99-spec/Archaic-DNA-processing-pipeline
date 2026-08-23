#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V3 groundwork: continent x period qpAdm/qpWave model-rejection tables.

The global migration paper (ROADMAP_GENETICS.md V3) needs, before any prose:
  * source/target definitions (the same West-Eurasian source library the
    ancestry module already ships), and
  * clear continent/period-specific model-rejection tables.

This script supplies the reusable scaffolding: it defines continent/period
target cohorts (from the phase-4 global analysis table), assembles their pooled
frequencies, and runs the qpAdm model competition plus qpWave rank tests for
each, writing a rejection table per continent/period.

The rejection-table logic itself is pure and unit-tested in
`archaic.ancestry.model_rejection_table` (tests/test_ancestry.py); this script is
the AADR-facing driver that supplies real cohorts and frequencies.

Output (all exploratory, under `results/v3_migration/`):
  rejection_table.csv     one row per (continent, period, model): p, chi2, status
  qpwave.csv              qpWave rank tests per cohort set

Run:
  PYTHONIOENCODING=utf-8 python scripts/v3_migration.py --panel 1240k
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import stats as st
from archaic import qpadm as qp
from archaic import ancestry as anc
from archaic.panel import Panel
from archaic.refs import PANELS

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "results", "v3_migration")

# Periods in years BP (ancient genomes only; present-day `is_modern` excluded so
# the migration analysis compares dated horizons rather than modern proxies).
# Each tuple is (label, older_bp, younger_bp) with membership older_bp <= bp <
# younger_bp, so the bins are contiguous and non-overlapping from oldest to
# youngest (larger BP = older).
PERIODS = [
    ("LGM_and_earlier", 100000, 15000),
    ("Late_Glacial", 15000, 11000),
    ("Early_Holocene", 11000, 7000),
    ("Neolithic", 7000, 4500),
    ("Bronze_Age", 4500, 3200),
    ("Iron_Age", 3200, 2300),
    ("Historic", 2300, 0),
]


def build_period_meta(df):
    """Attach a `period` label to each ancient (non-modern) row."""
    meta = df.copy()
    meta = meta[meta["is_modern"].astype(bool).eq(False)].copy()
    meta["period"] = "other"
    for label, older_bp, younger_bp in PERIODS:
        m = (meta["date_bp"] <= older_bp) & (meta["date_bp"] > younger_bp)
        meta.loc[m, "period"] = label
    return meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", default="1240k")
    parser.add_argument("--maxn", type=int, default=60,
                        help="cap individuals per cohort before kinship pruning")
    args = parser.parse_args()

    os.makedirs(OUT, exist_ok=True)
    panel = Panel(PANELS[args.panel]["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(HERE, "results", f"phase4_{args.panel}_global_analysis.csv"))
    meta = build_period_meta(meta)
    meta["gid_lower"] = meta["group_id"].str.lower()

    # ---- source + outgroup columns (the ancestry module's own library) ----
    source_cols = anc.resolve_cohorts(
        panel, meta,
        {name: spec["pred"] for name, spec in anc.SOURCES.items()},
        maxn=args.maxn, kin_prune=True, verbose=False)
    # distal outgroups
    out_cols = {}
    for name, (kind, val) in anc.BASE_RIGHT.items():
        if kind == "pop":
            out_cols[name] = panel.cols_for(pops=[val])
        else:
            out_cols[name] = panel.cols_for(ids=[val])

    # ---- per-(continent, period) target cohorts ----
    cells = []
    for continent in ("Eurasia", "Africa", "Oceania", "Americas"):
        cmeta = meta[meta["continent"].eq(continent)]
        for label, older_bp, younger_bp in PERIODS:
            pm = cmeta[cmeta["period"].eq(label)]
            if len(pm) < 5:
                continue
            ids = pm["genetic_id"].tolist()
            cols = np.array([panel._id_to_col[i] for i in ids if i in panel._id_to_col],
                            dtype=np.int64)
            if len(cols) > args.maxn:
                rng = np.random.default_rng(0)
                cols = np.sort(rng.choice(cols, args.maxn, replace=False))
            if len(cols) >= 3:
                cells.append((f"{continent}_{label}", cols))
    print(f"target cells (n>=3 individuals): {len(cells)}")

    # ---- one streaming frequency pass over everything ----
    from archaic import profiles as pf
    all_cols = {name: cols for name, cols in source_cols.items() if len(cols)}
    all_cols.update({name: cols for name, cols in out_cols.items() if len(cols)})
    all_cols.update({name: cols for name, cols in cells})
    freq, info = pf.cohort_frequencies(panel, all_cols)

    # ---- rejection table (per cell) ----
    targets = [name for name, _ in cells]
    rows = anc.model_rejection_table(freq, targets, anc.MODELS, block, 50)
    rej = pd.DataFrame(rows)
    rej.to_csv(os.path.join(OUT, "rejection_table.csv"), index=False)
    print(f"rejection rows: {len(rej)}")
    if not rej.empty:
        rejected = (rej["status"] == "rejected").sum()
        plausible = (rej["status"] == "plausible").sum()
        print(f"  plausible models: {plausible}   rejected models: {rejected}")

    # ---- qpWave rank tests per cell (top model source set + target) ----
    wave_rows = []
    for tgt in targets:
        srcs = ["Anatolia_N", "Steppe_Yamnaya", "WHG"]  # canonical west3 sources
        srcs = [s for s in srcs if s in freq and np.isfinite(freq[s]).any()]
        outs = [o for o in anc.BASE_RIGHT if o in freq and np.isfinite(freq[o]).any()]
        if len(srcs) < 2 or len(outs) <= len(srcs):
            continue
        for wr in qp.qpwave(freq, [tgt] + srcs, outs, block, 50):
            wave_rows.append(dict(target=tgt, **wr))
    if wave_rows:
        pd.DataFrame(wave_rows).to_csv(os.path.join(OUT, "qpwave.csv"), index=False)
        print(f"qpWave rows: {len(wave_rows)}")

    print("\nWrote results/v3_migration/")


if __name__ == "__main__":
    main()
