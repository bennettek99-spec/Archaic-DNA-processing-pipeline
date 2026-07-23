#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_dstat_probe.py
=================
A FEASIBILITY PROBE — not the full Tier-1 workflow. It computes a handful of
block-jackknife D-statistics on the AADR 1240K panel to give the feasibility
report real numbers for the power / detectability questions, and to expose the
sign-convention trap that the project spec explicitly warns about.

CONVENTION (matches the validated pipeline, denisovan_survey.py):
    D_Den(X) = D(X, Mbuti; Denisova, Chimp)
    POSITIVE  => X shares MORE derived alleles with the Denisovan than Mbuti does
                 (Papuan is the positive control: ~+3.4%).
    The task spec writes the same quantity as D(African, X; Denisovan, Chimp),
    which is the NEGATIVE of this by the antisymmetry D(W,X;Y,Z) = -D(X,W;Y,Z).
    We use the pipeline convention here for consistency with the already-validated
    positive control, and report the conversion explicitly.

Three families of tests:
  Tier-1 basic affinity:  D(X, Mbuti; Denisova, Chimp)        [positive control scale]
  Tier-2 pairwise split: D(X, Han; Denisova, Chimp)           [does X differ from Han?]
  Tier-3 Neanderthal-cntrl: D(Neanderthal, Denisova; X, Mbuti) [Denisovan-specific axis]
Each is run all-sites and transversions-only. Read-only on AADR; reuses archaic.*.

Outputs (under native_american_denisovan_project/):
  results/tables/table3_basic_denisovan_affinity.tsv   (Table 3 basis)
  results/tables/table4_pairwise_denisovan_contrasts.tsv (Table 4 basis)
  results/logs/dstat_probe.txt
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _c in (os.path.join(os.path.dirname(_REPO), "archaic-introgression"), os.path.dirname(_REPO)):
    if os.path.isdir(os.path.join(_c, "archaic")):
        sys.path.insert(0, _c)
        break

from archaic.panel import Panel
from archaic import stats as st
from archaic.config import panel_prefix
from archaic import profiles as prof

N_BLOCKS = 50


REF_SPEC = {
    "Denisova": {"ids": ["Denisova.SG"]},
    "Altai": {"ids": ["AltaiNeanderthal.DG"]},
    "Vindija": {"ids": ["VindijaG1_final.SG"]},
    "Chimp": {"ids": ["Chimp.REF"]},
    "Mbuti": {"pops": ["Mbuti"]},
    "Yoruba": {"pops": ["Yoruba", "YRI", "YRI-Discovery"]},
}

PRESENT_POOLS = {
    "Papuan": {"pops": ["Papuan"]},
    "Han": {"pops": ["Han"]},
    "Japanese": {"pops": ["Japanese"]},
    "Dai": {"pops": ["Dai"]},
    "Karitiana": {"pops": ["Karitiana"]},
    "French": {"pops": ["French"]},
    "Yoruba": {"pops": ["Yoruba", "YRI", "YRI-Discovery"]},
}

ANCIENT_TARGETS = [
    ("USR1.SG", "Ancient Beringian"),
    ("Kolyma1.SG", "Ancient Paleo-Siberian"),
    ("MA1.SG", "Ancient North Eurasian (Mal'ta)"),
    ("Yana1.SG", "Upper Paleolithic Siberian (Yana)"),
    ("Sumidouro6.SG", "Ancient South America (Lagoa Santa)"),
]


def freqs_for(panel, spec):
    cols = panel.cols_for(spec.get("ids"), spec.get("pops"))
    if len(cols) == 0:
        return None, 0
    f, _ = prof.cohort_frequencies(panel, {"X": cols}, min_ind=1)
    return f["X"], len(cols)


def run_table(panel, ref, block, mode, targets, present):
    rows = []
    freq = dict(ref)
    for name, p in present.items():
        pX, n = freqs_for(panel, p)
        if pX is None:
            continue
        freq[name] = pX
    for sid, role in targets:
        cols = panel.cols_for(ids=[sid])
        if len(cols) == 0:
            cand = [i for i in panel.ind["id"].tolist() if sid.split(".")[0].lower() in i.lower()]
            if cand:
                cols = panel.cols_for(ids=[cand[0]]); sid = cand[0]
            else:
                rows.append(dict(mode=mode, test="D(X,Mbuti;Denisova,Chimp)", X=sid,
                                 role=role, D=np.nan, SE=np.nan, Z=np.nan, nSNP=0, n_ind=0,
                                 note="not in panel"))
                continue
        pX, _ = freqs_for(panel, {"ids": [sid]})
        freq[sid] = pX

    for X in list(present.keys()) + [t[0] for t in targets]:
        if X not in freq or freq[X] is None:
            continue
        if not {"Mbuti", "Denisova", "Chimp"}.issubset({k: v for k, v in freq.items() if v is not None}):
            continue
        o = st.dstat(freq, X, "Mbuti", "Denisova", "Chimp", block, N_BLOCKS)
        rows.append(dict(mode=mode, test="D(X,Mbuti;Denisova,Chimp)", X=X,
                         role=_role(X, present, targets),
                         D=o["theta"], SE=o["se"], Z=o["z"], nSNP=o["n_used"]))
    return pd.DataFrame(rows)


def run_contrasts(panel, ref, block, mode, targets):
    rows = []
    freq = dict(ref)
    pHan, _ = freqs_for(panel, PRESENT_POOLS["Han"])
    freq["Han"] = pHan
    for sid, role in targets:
        cols = panel.cols_for(ids=[sid])
        if len(cols) == 0:
            cand = [i for i in panel.ind["id"].tolist() if sid.split(".")[0].lower() in i.lower()]
            if cand:
                cols = panel.cols_for(ids=[cand[0]]); sid = cand[0]
            else:
                rows.append(dict(mode=mode, test="D(X,Han;Denisova,Chimp)", X=sid,
                                 role=role, D=np.nan, SE=np.nan, Z=np.nan, nSNP=0,
                                 note="not in panel"))
                continue
        pX, _ = freqs_for(panel, {"ids": [sid]})
        freq[sid] = pX
        o = st.dstat(freq, sid, "Han", "Denisova", "Chimp", block, N_BLOCKS)
        rows.append(dict(mode=mode, test="D(X,Han;Denisova,Chimp)", X=sid, role=role,
                         D=o["theta"], SE=o["se"], Z=o["z"], nSNP=o["n_used"]))
    # add Papuan vs Han as a positive contrast anchor
    pPap, _ = freqs_for(panel, PRESENT_POOLS["Papuan"])
    if pPap is not None:
        freq["Papuan"] = pPap
        o = st.dstat(freq, "Papuan", "Han", "Denisova", "Chimp", block, N_BLOCKS)
        rows.append(dict(mode=mode, test="D(X,Han;Denisova,Chimp)", X="Papuan",
                         role="present-day positive contrast",
                         D=o["theta"], SE=o["se"], Z=o["z"], nSNP=o["n_used"]))
    # Karitiana vs Han
    pKar, _ = freqs_for(panel, PRESENT_POOLS["Karitiana"])
    if pKar is not None:
        freq["Karitiana"] = pKar
        o = st.dstat(freq, "Karitiana", "Han", "Denisova", "Chimp", block, N_BLOCKS)
        rows.append(dict(mode=mode, test="D(X,Han;Denisova,Chimp)", X="Karitiana",
                         role="present-day Native American",
                         D=o["theta"], SE=o["se"], Z=o["z"], nSNP=o["n_used"]))
    return pd.DataFrame(rows)


def run_neanderthal_control(panel, ref, block, mode, targets, present):
    """D(Altai, Denisova; X, Mbuti): a Denisovan-vs-Neanderthal axis conditioned on X.
    Expected sign derived in statistic_interpretation.md; reported raw here."""
    rows = []
    freq = dict(ref)
    for name in ("Mbuti",):
        pass
    for X in list(present.keys()) + [t[0] for t in targets]:
        p = present.get(X) or ({"ids": [X]} if not any(X == t[0] for t in targets) else {"ids": [X]})
        pX, _ = freqs_for(panel, p if "pops" in p or "ids" in p else {"ids": [X]})
        if pX is None:
            cand = [i for i in panel.ind["id"].tolist() if X.split(".")[0].lower() in i.lower()]
            if cand:
                pX, _ = freqs_for(panel, {"ids": [cand[0]]}); X = cand[0]
            else:
                continue
        freq[X] = pX
        if "Altai" not in freq or "Denisova" not in freq or "Mbuti" not in freq:
            continue
        o = st.dstat(freq, "Altai", "Denisova", X, "Mbuti", block, N_BLOCKS)
        rows.append(dict(mode=mode, test="D(Altai,Denisova;X,Mbuti)", X=X,
                         role=_role(X, present, targets),
                         D=o["theta"], SE=o["se"], Z=o["z"], nSNP=o["n_used"]))
    return pd.DataFrame(rows)


def _role(X, present, targets):
    if X in present:
        return "present-day"
    for sid, role in targets:
        if sid == X:
            return role
    return "?"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", default="1240k")
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)

    out = os.path.join(_REPO, "results")
    for sub in ("tables", "logs"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    prefix = panel_prefix(args.panel, args.config)
    print(f"Loading panel {args.panel} ({os.path.basename(prefix)})...")
    panel = Panel(prefix, autosomes_only=True)
    block = st.assign_blocks(panel.n_snp, N_BLOCKS)
    print(f"  autosomal SNPs: {panel.n_snp:,}")

    print("Reference frequencies...")
    ref, _ = panel.frequencies(REF_SPEC)

    print("\n[Tier 1] D(X, Mbuti; Denisova, Chimp) — basic Denisovan affinity (all sites)")
    t1_all = run_table(panel, ref, block, "all_snps", ANCIENT_TARGETS, PRESENT_POOLS)
    print(t1_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))

    print("\n[Tier 2] D(X, Han; Denisova, Chimp) — pairwise split vs Han (all sites)")
    t2_all = run_contrasts(panel, ref, block, "all_snps", ANCIENT_TARGETS)
    print(t2_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))

    print("\n[Tier 3] D(Altai, Denisova; X, Mbuti) — Neanderthal-conditioned axis (all sites)")
    t3_all = run_neanderthal_control(panel, ref, block, "all_snps", ANCIENT_TARGETS, PRESENT_POOLS)
    print(t3_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))

    print("\nLoading transversion-only sub-panel...")
    tv = Panel(prefix, autosomes_only=True, transversions_only=True)
    tv_block = st.assign_blocks(tv.n_snp, N_BLOCKS)
    tv_ref, _ = tv.frequencies(REF_SPEC)
    t1_tv = run_table(tv, tv_ref, tv_block, "transversions_only", ANCIENT_TARGETS, PRESENT_POOLS)
    t2_tv = run_contrasts(tv, tv_ref, tv_block, "transversions_only", ANCIENT_TARGETS)
    t3_tv = run_neanderthal_control(tv, tv_ref, tv_block, "transversions_only",
                                    ANCIENT_TARGETS, PRESENT_POOLS)

    t1 = pd.concat([t1_all, t1_tv], ignore_index=True)
    t2 = pd.concat([t2_all, t2_tv], ignore_index=True)
    t3 = pd.concat([t3_all, t3_tv], ignore_index=True)
    t1.to_csv(os.path.join(out, "tables", "table3_basic_denisovan_affinity.tsv"),
              sep="\t", index=False, na_rep="NA")
    t2.to_csv(os.path.join(out, "tables", "table4_pairwise_denisovan_contrasts.tsv"),
              sep="\t", index=False, na_rep="NA")
    t3.to_csv(os.path.join(out, "tables", "table_neanderthal_controlled_axis.tsv"),
              sep="\t", index=False, na_rep="NA")

    lines = []
    def w(s=""): print(s); lines.append(s)
    w("=" * 90)
    w("D-STATISTIC FEASIBILITY PROBE  (convention: D(X, Mbuti; Denisova, Chimp); + = more Denisovan)")
    w("  task-spec D(African,X;Denisovan,Chimp) == -D(X,Mbuti;Denisova,Chimp)  (antisymmetry)")
    w("=" * 90)
    w("\n--- Tier 1: D(X, Mbuti; Denisova, Chimp)  [all sites] ---")
    w(t1_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    w("\n--- Tier 1: transversions only ---")
    w(t1_tv.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    w("\n--- Tier 2: D(X, Han; Denisova, Chimp)  [all sites] ---")
    w(t2_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    w("\n--- Tier 2: transversions only ---")
    w(t2_tv.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    w("\n--- Tier 3: D(Altai, Denisova; X, Mbuti)  [all sites] ---")
    w(t3_all.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    with open(os.path.join(out, "logs", "dstat_probe.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"\nWrote: results/tables/table3*.tsv, table4*.tsv, table_neanderthal_controlled_axis.tsv, results/logs/dstat_probe.txt")


if __name__ == "__main__":
    raise SystemExit(main())
