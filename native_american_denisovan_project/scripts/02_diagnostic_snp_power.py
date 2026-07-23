#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_diagnostic_snp_power.py
==========================
Estimate the number and callability of Denisovan-diagnostic SNPs on the AADR
1240K panel — the central power question for this study (Native American
Denisovan ancestry is expected to be small, so the count of informative sites
is the binding constraint).

A "Denisovan-diagnostic" SNP is defined here with the SAME frequency-based rule
already validated in ``archaic.denisovan_genome._marker_sharing`` (no ancestral-
allele call is required, because the AADR .snp file is not reliably polarised to
ancestral/derived; the rule is polarisation-invariant):

    Denisova.SG carries one allele at >= 0.90
    pooled African baseline (Mbuti + Yoruba) carries that allele at <= 0.10
    pooled Neanderthal baseline (Altai + Vindija) carries it at <= 0.50
    the site is finite in all three references

This is the project's **Set B** (Denisovan-shared, Neanderthal-excluded). A
stricter **Set A** additionally requires the Neanderthal mean <= 0.10, and a
shared-archaic **Set F** (negative control) keeps sites where Neanderthal and
Denisovan agree against Africa. We report all three counts so the feasibility
estimate is not tied to one threshold choice.

We then read a handful of representative individuals (ancient American, ANE,
Paleo-Siberian, Jomon, present-day Papuan/Han/Japanese/French) and count how
many diagnostic markers each can actually contribute — the empirical power
number behind every later statistic.

READ-ONLY on AADR genotype data; reuses the validated ``archaic`` reader. Writes
no production files.

Outputs (under ``native_american_denisovan_project/``):
  data/diagnostic_sites/denisovan_diagnostic_sites.tsv  per-site Set A/B/F labels
  results/tables/table_diagnostic_snp_power.tsv         per-individual callability
  results/tables/table_diagnostic_snp_counts.tsv        counts by set + chromosome
  results/logs/diagnostic_snp_power.txt                  summary
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


# Reference sets (AADR 1240K individual IDs / group labels), mirroring refs.py.
REF_SPEC = {
    "Denisova": {"ids": ["Denisova.SG"]},
    "Altai": {"ids": ["AltaiNeanderthal.DG"]},
    "Vindija": {"ids": ["VindijaG1_final.SG"]},
    "Chagyrskaya": {"ids": ["Chagyrskaya8.DG"]},
    "Mbuti": {"pops": ["Mbuti"]},
    "Yoruba": {"pops": ["Yoruba", "YRI", "YRI-Discovery"]},
    "Chimp": {"ids": ["Chimp.REF"]},
}

# Representative individuals for the callability / power table.
REPRESENTATIVE = [
    ("USR1.SG", "Ancient Beringian"),
    ("Kolyma1.SG", "Ancient Paleo-Siberian"),
    ("MA1.SG", "Ancient North Eurasian (Mal'ta)"),
    ("Yana1.SG", "Upper Paleolithic Siberian (Yana)"),
    ("Sumidouro6.SG", "Ancient South America (Lagoa Santa)"),
    ("Anzick-1.SG", "Anzick-1 (Clovis)"),
    ("SpiritCave.SG", "Spirit Cave"),
]


def build_diagnostic(freq):
    p_den = freq["Denisova"]
    with np.errstate(invalid="ignore"):
        p_nea = np.nanmean(np.vstack([freq["Altai"], freq["Vindija"]]), axis=0)
        p_afr = np.nanmean(np.vstack([freq["Mbuti"], freq["Yoruba"]]), axis=0)
    den_is_a1 = p_den > 0.5
    den_extreme = np.where(den_is_a1, p_den, 1.0 - p_den)
    afr_of_den = np.where(den_is_a1, p_afr, 1.0 - p_afr)
    nea_of_den = np.where(den_is_a1, p_nea, 1.0 - p_nea)
    finite = (np.isfinite(p_den) & np.isfinite(p_afr) & np.isfinite(p_nea))
    set_B = finite & (den_extreme >= 0.90) & (afr_of_den <= 0.10) & (nea_of_den <= 0.50)
    set_A = set_B & (nea_of_den <= 0.10)
    set_F = finite & (den_extreme >= 0.90) & (afr_of_den <= 0.10) & (nea_of_den >= 0.90)
    return den_is_a1, set_A, set_B, set_F, {"p_den": p_den, "p_nea": p_nea, "p_afr": p_afr}


def oriented_freq(panel, cols, den_is_a1, mask):
    """Mean frequency of the Denisovan-high allele across masked diagnostic SNPs."""
    cols = np.asarray(cols, dtype=np.int64)
    if len(cols) == 0:
        return np.nan, 0
    G = panel.pg.read(panel.snp_rows[mask], cols)
    Gf = G.astype(np.float32); Gf[G < 0] = np.nan
    with np.errstate(invalid="ignore"):
        p = np.nanmean(Gf, axis=1) / 2.0
    sub_den_a1 = den_is_a1[mask]
    oriented = np.where(sub_den_a1, p, 1.0 - p)
    callable_mask = np.isfinite(oriented)
    n = int(callable_mask.sum())
    return (float(np.nanmean(oriented[callable_mask])) if n else np.nan), n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", default="1240k")
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)

    out = os.path.join(_REPO, "results")
    for sub in ("tables", "logs"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    dsite = os.path.join(_REPO, "data", "diagnostic_sites")
    os.makedirs(dsite, exist_ok=True)

    prefix = panel_prefix(args.panel, args.config)
    print(f"Loading panel {args.panel} ({os.path.basename(prefix)})...")
    panel = Panel(prefix, autosomes_only=True)
    print(f"  autosomal SNPs: {panel.n_snp:,}")

    print("Reading reference frequencies (Denisova/Altai/Vindija/Mbuti/Yoruba/Chimp)...")
    freq, info = panel.frequencies(REF_SPEC)
    for k, v in info.items():
        print(f"    {k}: n_ind={v['n_ind']}, snp_covered={v['n_snp_covered']:,}")

    den_is_a1, set_A, set_B, set_F, raw = build_diagnostic(freq)
    n_A, n_B, n_F = int(set_A.sum()), int(set_B.sum()), int(set_F.sum())
    print(f"\nDenisovan-diagnostic SNP counts (autosomal 1240K):")
    print(f"  Set A (strict, Neanderthal <= 0.10):  {n_A:,}")
    print(f"  Set B (shared, Neanderthal <= 0.50):  {n_B:,}")
    print(f"  Set F (shared-archaic negative ctrl): {n_F:,}")

    chrom = panel.snp.loc[panel.snp_rows, "chrom"].to_numpy()
    counts_by_chrom = []
    for c in [str(i) for i in range(1, 23)]:
        isc = chrom == c
        counts_by_chrom.append({
            "chromosome": c, "n_SNP_panel": int(isc.sum()),
            "n_set_A": int((set_A & isc).sum()),
            "n_set_B": int((set_B & isc).sum()),
            "n_set_F": int((set_F & isc).sum()),
        })
    cdf = pd.DataFrame(counts_by_chrom)
    cdf.to_csv(os.path.join(out, "tables", "table_diagnostic_snp_counts.tsv"),
               sep="\t", index=False)

    snp = panel.snp.loc[panel.snp_rows, ["name", "chrom", "pos", "a1", "a2"]].copy()
    snp = snp.reset_index(drop=True)
    snp["set_A"] = set_A
    snp["set_B"] = set_B
    snp["set_F"] = set_F
    snp["denisova_high_allele_is_a1"] = den_is_a1
    snp["p_denisova"] = raw["p_den"]
    snp["p_neanderthal_mean"] = raw["p_nea"]
    snp["p_african_mean"] = raw["p_afr"]
    diag = snp[snp["set_B"]].copy()
    diag.to_csv(os.path.join(dsite, "denisovan_diagnostic_sites.tsv"),
                sep="\t", index=False)

    print("\nPer-individual callability of Set B diagnostic markers:")
    rows = []
    # present-day pools for context
    present_pools = {
        "Papuan": ("pops", ["Papuan"]),
        "Han": ("pops", ["Han"]),
        "Japanese": ("pops", ["Japanese"]),
        "Dai": ("pops", ["Dai"]),
        "French": ("pops", ["French"]),
        "Mbuti": ("pops", ["Mbuti"]),
        "Yoruba": ("pops", ["Yoruba"]),
    }
    for name, (kind, vals) in present_pools.items():
        if kind == "pops":
            cols = panel.cols_for(pops=vals)
        else:
            cols = panel.cols_for(ids=vals)
        if len(cols) == 0:
            print(f"    {name}: not present in panel")
            continue
        mean, n = oriented_freq(panel, cols, den_is_a1, set_B)
        rows.append({"sample": f"[pool] {name}", "role": "present-day control",
                     "n_callable_diag_markers": n,
                     "mean_denisovan_marker_allele": mean,
                     "pct": 100 * (mean if np.isfinite(mean) else np.nan)})
        print(f"    [pool] {name:<10} n_diag={n:>6}  marker_allele={mean:.4f}  ({100*mean:.2f}%)")

    for sid, role in REPRESENTATIVE:
        cols = panel.cols_for(ids=[sid])
        if len(cols) == 0:
            # try a fuzzy lookup
            cand = [i for i in panel.ind["id"].tolist() if sid.split(".")[0].lower() in i.lower()]
            if cand:
                cols = panel.cols_for(ids=[cand[0]]); sid = cand[0]
            else:
                print(f"    {sid}: not present in panel")
                rows.append({"sample": sid, "role": role,
                             "n_callable_diag_markers": 0,
                             "mean_denisovan_marker_allele": np.nan, "pct": np.nan})
                continue
        g = panel.pg.read(panel.snp_rows[set_B], cols)
        g = g.astype(np.float32); g[g < 0] = np.nan
        sub_den_a1 = den_is_a1[set_B]
        with np.errstate(invalid="ignore"):
            p = g[:, 0] / 2.0
        oriented = np.where(sub_den_a1, p, 1.0 - p)
        n = int(np.isfinite(oriented).sum())
        mean = float(np.nanmean(oriented)) if n else np.nan
        rows.append({"sample": sid, "role": role,
                     "n_callable_diag_markers": n,
                     "mean_denisovan_marker_allele": mean, "pct": 100 * mean})
        print(f"    {sid:<22} {role:<42} n_diag={n:>5}  marker_allele={mean if np.isfinite(mean) else float('nan'):.4f}")

    pow = pd.DataFrame(rows)
    pow.to_csv(os.path.join(out, "tables", "table_diagnostic_snp_power.tsv"),
               sep="\t", index=False)

    # transversion-only diagnostic count (damage-robust subset)
    tv = Panel(prefix, autosomes_only=True, transversions_only=True)
    print(f"\nTransversion-only sub-panel: {tv.n_snp:,} autosomal SNPs")
    tv_freq, _ = tv.frequencies(REF_SPEC)
    _, tvA, tvB, tvF, _ = build_diagnostic(tv_freq)
    print(f"  Set A (TV only): {int(tvA.sum()):,}   Set B (TV only): {int(tvB.sum()):,}   "
          f"Set F (TV only): {int(tvF.sum()):,}")

    lines = []
    def w(s=""): print(s); lines.append(s)
    w("=" * 78)
    w("DENISOVAN-DIAGNOSTIC SNP POWER ESTIMATE (1240K autosomal panel)")
    w("=" * 78)
    w(f"Autosomal panel SNPs: {panel.n_snp:,}")
    w(f"Set A (strict):  {n_A:,}   Set B (shared): {n_B:,}   Set F (neg ctrl): {n_F:,}")
    w(f"TV-only sub-panel SNPs: {tv.n_snp:,}")
    w(f"TV-only Set A: {int(tvA.sum()):,}   Set B: {int(tvB.sum()):,}   Set F: {int(tvF.sum()):,}")
    w("")
    w("Per-chromosome Set B counts:")
    w(cdf.to_string(index=False))
    w("")
    w("Per-individual callability (Set B):")
    w(pow.to_string(index=False))
    with open(os.path.join(out, "logs", "diagnostic_snp_power.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"\nWrote: data/diagnostic_sites/denisovan_diagnostic_sites.tsv, "
          f"results/tables/table_diagnostic_snp_*.tsv, results/logs/diagnostic_snp_power.txt")


if __name__ == "__main__":
    raise SystemExit(main())
