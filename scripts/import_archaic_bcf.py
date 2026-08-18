#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Put Chagyrskaya on the 1240K panel, and count what each contrast axis is worth.

Every power result in this study so far ends in the same place: the binding
constraint is the number of sites that distinguish the two archaic genomes being
contrasted. Genome-wide that is 18,926 of 528,283 for Vindija against Altai.
The three-way subsample on the tract statistic (f3917f0) found the same axis
binding there. So the only lever left is a better contrast, and the question
this script answers is a counting question:

    how many 1240K sites separate each available pair of archaic genomes?

Zenodo 7246376 supplies Altai, Vindija, Chagyrskaya and Denisova on hg19, the
AADR's own build, so the join is on chromosome and position with no liftover.

WHAT MAKES THIS FIDDLY

  * A BCF record with no ALT is not missing data. It means every sample is
    homozygous reference, which is a call, and on chr22 that case covers nearly
    half the panel's sites. An earlier concordance check discarded them and this
    one must not.
  * The panel's two alleles and the BCF's REF/ALT have to describe the same
    pair before a genotype can be turned into a dosage of the panel's counted
    allele. Sites where they do not are counted and dropped rather than
    force-matched; a strand flip guessed wrong would invent a difference
    between two archaic genomes, which is precisely the signal being measured.
  * Indels and the '-' allele that appears in these files are not SNPs and are
    excluded.

SELF-CHECK, NOT SELF-REPORT

AADR already contains Altai and Vindija, so those two are read twice by
independent routes and compared genome-wide. The comparison is stratified by
zygosity because AADR ships Vindija as pseudo-haploid `.SG`: it records one
randomly drawn allele, so it *cannot* agree at a true heterozygote, and an
unstratified concordance number would look like a parser bug. Altai is diploid
`.DG` and should agree essentially everywhere, including at heterozygotes.

Outputs (results/archaic_panel/):
  archaic_1240k_freqs.npz    per-archaic allele frequency on the panel's SNP order
  archaic_pair_counts.csv    for every pair: sites both called, sites differing
  archaic_import_qc.csv      concordance against AADR, and what was dropped

Run: PYTHONIOENCODING=utf-8 python scripts/import_archaic_bcf.py
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from archaic.bcf_reader import BCFReader
from archaic.log_utils import get_logger
from archaic.panel import Panel
from archaic.refs import PANELS

log = get_logger("archaic.import_archaic_bcf")
BCF_DIR = os.path.join(ROOT, "data", "archaic_hg19")
OUT = os.path.join(ROOT, "results", "archaic_panel")
BASES = frozenset("ACGT")

# BCF sample name -> the short label used throughout this pipeline
SAMPLES = {"AltaiNeandertal": "Altai", "Vindija33.19": "Vindija",
           "Denisova": "Denisova", "Chagyrskaya-Phalanx": "Chagyrskaya"}
# AADR's own copies of two of them, for the ground-truth check
AADR_EQUIV = {"Altai": "AltaiNeanderthal.DG", "Vindija": "VindijaG1_final.SG"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chroms", default="1-22")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    chroms = []
    for part in args.chroms.split(","):
        if "-" in part:
            a, b = part.split("-")
            chroms.extend(range(int(a), int(b) + 1))
        else:
            chroms.append(int(part))

    cfg = PANELS[args.panel]
    log.info(f"Loading panel {args.panel}...")
    panel = Panel(cfg["prefix"], autosomes_only=True)
    snp = panel.snp
    n_snp = len(snp)
    log.info(f"  {n_snp:,} autosomal panel sites")

    chrom_s = snp["chrom"].astype(str).to_numpy()
    pos = snp["pos"].to_numpy(np.int64)
    a1 = snp["a1"].astype(str).to_numpy()
    a2 = snp["a2"].astype(str).to_numpy()
    # position -> row, per chromosome
    index = {}
    for c in np.unique(chrom_s):
        m = np.flatnonzero(chrom_s == c)
        index[c] = dict(zip(pos[m].tolist(), m.tolist()))

    labels = list(SAMPLES.values())
    dose = {l: np.full(n_snp, -1, dtype=np.int8) for l in labels}
    stats = dict(records=0, matched=0, no_alt_used=0, allele_mismatch=0,
                 non_snp=0, multiallelic_ok=0)

    t0 = time.time()
    for c in chroms:
        path = os.path.join(BCF_DIR, f"highcov_ind_{c}.bcf")
        if not os.path.exists(path):
            log.warning(f"  chr{c}: {os.path.basename(path)} absent, skipped")
            continue
        want = {str(c): set(index.get(str(c), {}).keys())}
        n_before = stats["matched"]
        with BCFReader(path) as r:
            order = [r.samples.index(s) for s in SAMPLES]
            names = [SAMPLES[s] for s in SAMPLES]
            for rec in r.iter_at(want):
                stats["records"] += 1
                row = index[str(c)].get(rec.pos)
                if row is None:
                    continue
                ok = _assign(rec, row, a1, a2, order, names, dose, stats)
                if ok:
                    stats["matched"] += 1
        log.info(f"  chr{c:<3} {stats['matched']-n_before:>7,} panel sites "
                 f"assigned  ({time.time()-t0:6.1f}s elapsed)")

    freqs = {l: np.where(dose[l] >= 0, dose[l] / 2.0, np.nan) for l in labels}
    np.savez(os.path.join(OUT, "archaic_1240k_freqs.npz"),
             **{l: freqs[l].astype(np.float32) for l in labels},
             pos=pos, chrom=chrom_s.astype("U2"))
    for l in labels:
        n = int((dose[l] >= 0).sum())
        log.info(f"  {l:12s} called at {n:>9,} panel sites "
                 f"({100*n/n_snp:.1f}%)")
    log.info(f"  dropped: {stats['allele_mismatch']:,} allele-set mismatches, "
             f"{stats['non_snp']:,} non-SNP; used "
             f"{stats['no_alt_used']:,} no-ALT (homozygous-reference) records")

    # ---- ground truth: AADR's own Altai and Vindija -------------------------
    qc = []
    for lab, aadr_id in AADR_EQUIV.items():
        col = panel._id_to_col.get(aadr_id)
        if col is None:
            continue
        g = panel.pg.read(np.arange(n_snp), np.array([col]))[:, 0]
        mine = dose[lab]
        both = (g >= 0) & (mine >= 0)
        het_bcf = mine == 1
        for kind, m in (("het", both & het_bcf), ("hom", both & ~het_bcf)):
            n = int(m.sum())
            agree = int((g[m] == mine[m]).sum())
            qc.append(dict(sample=lab, zygosity=kind, n=n, concordant=agree,
                           pct=100 * agree / n if n else np.nan))
            log.info(f"  [QC] {lab:8s} {kind:3s} {agree:>8,}/{n:<8,} = "
                     f"{100*agree/max(n,1):7.3f}%")
    pd.DataFrame(qc).to_csv(os.path.join(OUT, "archaic_import_qc.csv"),
                            index=False)

    # ---- the counting question ---------------------------------------------
    rows = []
    for i, x in enumerate(labels):
        for y in labels[i + 1:]:
            both = (dose[x] >= 0) & (dose[y] >= 0)
            diff = both & (dose[x] != dose[y])
            # "fixed" differences are the cleanest contrast: both homozygous
            # and opposite, so no heterozygote ambiguity enters the statistic.
            fixed = both & (((dose[x] == 0) & (dose[y] == 2))
                            | ((dose[x] == 2) & (dose[y] == 0)))
            rows.append(dict(a=x, b=y, both_called=int(both.sum()),
                             differing=int(diff.sum()),
                             fixed_differing=int(fixed.sum()),
                             pct_differing=100 * diff.sum() / max(both.sum(), 1)))
    cdf = pd.DataFrame(rows).sort_values("differing", ascending=False)
    cdf.to_csv(os.path.join(OUT, "archaic_pair_counts.csv"), index=False)
    log.info("Sites distinguishing each pair on the 1240K panel:")
    for _, r in cdf.iterrows():
        log.info(f"  {r['a']:12s} vs {r['b']:12s}  both called "
                 f"{r['both_called']:>9,}  differing {r['differing']:>8,} "
                 f"({r['pct_differing']:.2f}%)  fixed "
                 f"{r['fixed_differing']:>8,}")
    base = cdf[((cdf.a == "Vindija") & (cdf.b == "Altai"))
               | ((cdf.a == "Altai") & (cdf.b == "Vindija"))]
    if len(base):
        b = int(base["differing"].iloc[0])
        log.info(f"  baseline Vindija-vs-Altai here = {b:,} "
                 f"(the study's AADR-derived figure is 18,926)")
        for _, r in cdf.iterrows():
            if {r["a"], r["b"]} != {"Vindija", "Altai"}:
                log.info(f"    {r['a']}-vs-{r['b']}: "
                         f"{r['differing']/b:.2f}x the baseline")
    log.info(f"Wrote archaic_1240k_freqs.npz, archaic_pair_counts.csv and "
             f"archaic_import_qc.csv to {OUT}")


def _assign(rec, row, a1, a2, order, names, dose, stats):
    """Turn one BCF record into a dosage of the panel's counted allele."""
    pa, pb = a1[row], a2[row]
    alleles = [rec.ref] + list(rec.alts)
    if any((len(x) != 1 or x not in BASES) for x in alleles):
        stats["non_snp"] += 1
        return False
    if not rec.alts:
        # No ALT: every sample is homozygous reference. This is a call, not a
        # gap, and on this panel it covers a large share of the sites.
        if rec.ref == pa:
            d = 2
        elif rec.ref == pb:
            d = 0
        else:
            stats["allele_mismatch"] += 1
            return False
        for idx, nm in zip(order, names):
            g = rec.genotypes[idx]
            dose[nm][row] = d if not any(x < 0 for x in g) else -1
        stats["no_alt_used"] += 1
        return True
    if not set(alleles) <= {pa, pb}:
        stats["allele_mismatch"] += 1
        return False
    if len(alleles) > 2:
        stats["multiallelic_ok"] += 1
    for idx, nm in zip(order, names):
        g = rec.genotypes[idx]
        if any(x < 0 or x >= len(alleles) for x in g):
            continue
        dose[nm][row] = sum(1 for x in g if alleles[x] == pa)
    return True


if __name__ == "__main__":
    main()
