#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe: which differential-Neanderthal statistic best recovers the East-Asian
excess on the HO array? (methods check, reported honestly — not gate tuning)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st

N_BLOCKS = 50
panel = Panel()
spec = {
    "Altai": dict(ids=["AltaiNeanderthal.DG"]),
    "Vindija": dict(ids=["VindijaG1_final.SG"]),
    "Chimp": dict(ids=["Chimp_HO.HO"]),
    "Mbuti": dict(pops=["Mbuti"]),
    "Yoruba": dict(pops=["Yoruba"]),
    "French": dict(pops=["French"]),
    "Sardinian": dict(pops=["Sardinian"]),
    "Han": dict(pops=["Han"]),
    "Dai": dict(pops=["Dai"]),
    "Papuan": dict(pops=["Papuan"]),
}
freq, info = panel.frequencies(spec)
block = st.assign_blocks(panel.n_snp, N_BLOCKS)
print("coverage:", {k: info[k]["n_snp_covered"] for k in ["Altai","Vindija","Dai"]})

# differential statistics: D(WEur, EAsian; Nea, Outgroup) < 0 => EAsian more Nea
print("\nDifferential-Neanderthal contrasts  (<0 => 2nd pop more Neanderthal):")
combos = [
    ("French","Han","Altai","Chimp"),
    ("French","Han","Altai","Yoruba"),
    ("French","Han","Altai","Mbuti"),
    ("French","Han","Vindija","Yoruba"),
    ("French","Han","Vindija","Chimp"),
    ("Sardinian","Han","Altai","Yoruba"),
    ("Sardinian","Han","Vindija","Yoruba"),
    ("French","Dai","Altai","Yoruba"),
    ("Sardinian","Dai","Vindija","Yoruba"),
]
for W,X,Y,Z in combos:
    if any(k not in freq for k in (W,X,Y,Z)):
        continue
    r = st.dstat(freq, W, X, Y, Z, block, N_BLOCKS)
    print(f"  D({W:9s},{X:4s}; {Y:7s},{Z:6s}) = {r['theta']:+.4f}  "
          f"Z={r['z']:+5.1f}  ({r['n_used']:,} SNPs)")

# rank order of absolute Neanderthal affinity D(X,Mbuti;Altai,Chimp)
print("\nRank order, D(X, Mbuti; Altai, Chimp)  (higher => more Neanderthal):")
res = []
for X in ["Sardinian","French","Han","Dai","Papuan"]:
    if X in freq:
        r = st.dstat(freq, X, "Mbuti", "Altai", "Chimp", block, N_BLOCKS)
        res.append((X, r["theta"], r["z"]))
for X, d, z in sorted(res, key=lambda t: t[1]):
    print(f"  {X:10s} D={d:.4f}  Z={z:.1f}")
