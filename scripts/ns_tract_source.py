#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1: the Vindija-vs-Altai contrast measured on introgressed tracts.

The genome-wide study (`PAPER_neanderthal_source.md`) asks which Neanderthal a
population descends from by comparing whole genomes, in which the Neanderthal
component is ~2% of the sequence and 98% is noise for this purpose. Its
calibrated detection limit is f50 ~ 37%: a population would have to have
re-sourced more than a third of its Neanderthal ancestry before the panel would
notice. `ns_phase0_checks.py` projected that restricting the statistic to called
introgressed tracts — where local archaic ancestry is ~0.5 rather than ~0.021 —
should buy roughly 18x.

This runs that experiment on real tracts.

THE DATA, AND WHY IT WAS HIDING

Skov et al. (2018, PLoS Genetics 14:e1007641) table S5 is already cached in this
repository, downloaded for the Papuan Denisovan project. That project's importer
filters it to `pop == "papuans"`, so what the repository has been using is 89
Papuan individuals. The raw table holds **271 individuals from 94 populations**
across West Eurasia, East Asia, South Asia, Central Asia/Siberia and Melanesia,
and — decisively — it carries per-tract counts of the tract's variants shared
with Altai, with Vindija and with Denisova. The Vindija-versus-Altai contrast is
therefore already tabulated, per tract, on material that is essentially pure
archaic ancestry.

THE STATISTIC

    D_tract = (sum Shared_with_Vindija - sum Shared_with_Altai)
              / (sum Shared_with_Vindija + sum Shared_with_Altai)

summed over a group's Neanderthal-affinity tracts. It is the same question the
genome-wide D_VA asks, on the same two archaic genomes, but posed to the
introgressed sequence directly instead of to a 2% dilution of it.

As in the genome-wide study, only *differences between groups* are interpreted.
The absolute value carries the same common-mode offsets — Vindija and Altai are
not called at the same sites or to the same depth — and those cancel between
groups and do not cancel within one.

WHAT IS MEASURED RATHER THAN ASSUMED

  * kappa, again. The lesson of `POWER_mixture_calibration.md` is that a
    detection limit is arithmetic until the response of the statistic to a real
    mixture has been measured. Re-sourcing a tract to a lineage equidistant
    between Vindija and Altai makes its two sharing counts equal, so the mixture
    is injected by replacing (V, A) with ((V+A)/2, (V+A)/2) for a random
    fraction f of tracts, and kappa is read off the result.
  * The error bar, by block jackknife over genomic blocks, paired between groups
    exactly as the genome-wide study pairs it. A delete-one-individual jackknife
    is reported alongside, because with two or three genomes per population the
    between-individual variance is the term most likely to be understated.

Outputs (reports/neanderthal_source/):
  ns_tract_groups.csv      per-population and per-region D_tract with SEs
  ns_tract_pairs.csv       paired differences between the groups that matter
  ns_tract_limit.csv       kappa, the detection limit, and the genome-wide one
  fig_n7_tract.png         group values and the limit against the genome-wide one
  POWER_tract_restricted.md

Run: PYTHONIOENCODING=utf-8 python scripts/ns_tract_source.py
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from archaic import source_contrast as sc
from archaic.log_utils import get_logger

log = get_logger("archaic.ns_tract_source")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")
S5 = os.path.join(ROOT, "archaic_admixture_dating", "outputs",
                  "papuan_skov_real_v1", "downloads",
                  "skov_2018_s5_segments.tsv")
N_BLOCKS = 50
N_SIGMA = 2.0
N_BOOT = 2000
SEED = 20260817

# Regions compared. Melanesia is held out of the headline comparisons and used
# as a positive control instead: Papuans carry Denisovan ancestry, the Altai
# genome carries Denisovan ancestry too, and the genome-wide study already
# showed that combination pulls a population towards Altai. If this statistic
# does not reproduce that pull it is not measuring what it claims to.
HEADLINE = [("WestEurasia", "EastAsia"),
            ("WestEurasia", "SouthAsia"),
            ("EastAsia", "CentralAsiaSiberia"),
            ("WestEurasia", "CentralAsiaSiberia")]
POP_PAIRS = [("French", "Han"), ("Sardinian", "Japanese"),
             ("English", "Dai"), ("French", "Sardinian")]


# ------------------------------------------------------------------- loading --
def load_tracts(path=S5, min_snps=1):
    """Skov S5, restricted to the HMM caller and classified by archaic affinity.

    The published table mixes three callers (HMM, CRF, S*). Only the HMM rows
    are used, because that is the caller this repository has already calibrated
    (`skov_hmm.py`) and mixing callers would mix their biases.

    A tract is called Neanderthal-affinity when it shares more variants with
    Vindija than with Denisova, which is the same rule `tract_import.py` applies
    in the other direction for the Papuan project. Tracts where the two tie are
    reported separately rather than split arbitrarily: they are overwhelmingly
    low-count tracts carrying little information either way.
    """
    d = pd.read_csv(path, sep="\t")
    d = d[d["method"].astype(str).str.upper() == "HMM"].copy()
    for c in ("Shared_with_Altai", "Shared_with_Vindija",
              "Shared_with_Denisova", "snps", "start", "end"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["Shared_with_Altai", "Shared_with_Vindija",
                         "Shared_with_Denisova", "start", "end"])
    d = d[d["snps"] >= min_snps]
    v, a, dn = (d["Shared_with_Vindija"], d["Shared_with_Altai"],
                d["Shared_with_Denisova"])
    d["cls"] = np.where(v > dn, "neanderthal",
                        np.where(dn > v, "denisovan", "unresolved"))
    d["chrom"] = pd.to_numeric(d["chrom"], errors="coerce")
    d = d.dropna(subset=["chrom"])
    d["block"] = assign_blocks(d["chrom"].to_numpy(int),
                               d["start"].to_numpy(float), N_BLOCKS)
    return d


def assign_blocks(chrom, start, n_blocks):
    """Contiguous genomic blocks over a concatenated genome.

    Tracts are megabase-scale and strongly autocorrelated within an individual,
    so the resampling unit has to be a genomic region rather than a tract. Bins
    are equal-width in concatenated coordinates, which keeps a block from
    spanning a chromosome boundary in any way that matters and mirrors the
    fixed-width block scheme the genome-wide study uses on SNP indices.
    """
    order = np.argsort(chrom, kind="stable")
    offs, run = {}, 0.0
    for c in np.unique(chrom):
        offs[int(c)] = run
        run += float(start[chrom == c].max()) + 1e6
    pos = np.array([offs[int(c)] + s for c, s in zip(chrom, start)])
    edges = np.linspace(pos.min(), pos.max() + 1.0, n_blocks + 1)
    return np.clip(np.digitize(pos, edges) - 1, 0, n_blocks - 1)


# ----------------------------------------------------------------- statistic --
def block_sums(sub, n_blocks=N_BLOCKS):
    """Per-block (Vindija, Altai) sharing sums for one group of tracts."""
    v = np.bincount(sub["block"], weights=sub["Shared_with_Vindija"],
                    minlength=n_blocks)
    a = np.bincount(sub["block"], weights=sub["Shared_with_Altai"],
                    minlength=n_blocks)
    return v, a


def theta(v, a):
    tot = v.sum() + a.sum()
    return float((v.sum() - a.sum()) / tot) if tot > 0 else np.nan


def jackknife(v, a):
    """Delete-one-block jackknife of D_tract."""
    loo = []
    for b in range(len(v)):
        keep = np.ones(len(v), dtype=bool)
        keep[b] = False
        if (v[keep].sum() + a[keep].sum()) > 0:
            loo.append(theta(v[keep], a[keep]))
    loo = np.array([x for x in loo if np.isfinite(x)])
    g = len(loo)
    if g < 2:
        return np.nan
    return float(np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum()))


def paired_difference(va, aa, vb, ab):
    """a - b with the same blocks deleted from both, as the study does."""
    th = theta(va, aa) - theta(vb, ab)
    loo = []
    for b in range(len(va)):
        keep = np.ones(len(va), dtype=bool)
        keep[b] = False
        if (va[keep].sum() + aa[keep].sum()) > 0 and \
           (vb[keep].sum() + ab[keep].sum()) > 0:
            loo.append(theta(va[keep], aa[keep]) - theta(vb[keep], ab[keep]))
    loo = np.array([x for x in loo if np.isfinite(x)])
    g = len(loo)
    if g < 2:
        return dict(diff=th, se=np.nan, z=np.nan)
    se = float(np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum()))
    return dict(diff=float(th), se=se, z=float(th / se) if se > 0 else np.nan)


def jackknife_individuals(sub):
    """Delete-one-individual jackknife, as a check on the block SE.

    With two or three genomes in some populations, the between-individual term
    is the one a genomic block jackknife is most likely to understate. If the
    two SEs disagree badly the group is too small to quote.
    """
    names = sub["name"].unique()
    if len(names) < 3:
        return np.nan
    loo = []
    for nm in names:
        s = sub[sub["name"] != nm]
        v, a = block_sums(s)
        loo.append(theta(v, a))
    loo = np.array([x for x in loo if np.isfinite(x)])
    g = len(loo)
    return float(np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum()))


# --------------------------------------------------------- mixture calibration
def measure_kappa(sub, fractions, rng, n_rep=8):
    """kappa for this statistic, by re-sourcing a fraction of tracts.

    A tract drawn from a lineage equidistant between Vindija and Altai shares
    its variants equally with the two, so re-sourcing is injected by replacing
    that tract's (V, A) with ((V+A)/2, (V+A)/2). Tracts are the unit because
    tracts are what a source change would actually replace.

    The published conversion assumes kappa = 1. On the genome-wide statistic it
    was 0.33, and the reason was that D_VA's absolute value is inflated by
    offsets no mixture can move. Whether the same is true here has to be
    measured, not inherited.
    """
    v0 = sub["Shared_with_Vindija"].to_numpy(float)
    a0 = sub["Shared_with_Altai"].to_numpy(float)
    blk = sub["block"].to_numpy(int)
    base = theta(*_bs(v0, a0, blk))
    rows = []
    for f in fractions:
        if f == 0:
            continue
        for _ in range(n_rep):
            pick = rng.random(len(v0)) < f
            mid = 0.5 * (v0 + a0)
            v = np.where(pick, mid, v0)
            a = np.where(pick, mid, a0)
            th = theta(*_bs(v, a, blk))
            rows.append(dict(fraction=f, realised=th - base,
                             predicted=-f * base,
                             kappa=(th - base) / (-f * base)))
    return base, pd.DataFrame(rows)


def _bs(v, a, blk, n_blocks=N_BLOCKS):
    return (np.bincount(blk, weights=v, minlength=n_blocks),
            np.bincount(blk, weights=a, minlength=n_blocks))


def detection_curve(va, aa, vb, ab, kappa, base, rng, n_boot=N_BOOT):
    """f50/f80 for this statistic, by paired block bootstrap under injection.

    The same scheme as `ns_mixture_power.py`: resample the blocks with
    replacement, identically for both groups, and within each replicate
    recompute the difference and its jackknife SE, so the simulated analyst
    tests with the error bar they would really have had. The injected effect is
    a shift of kappa * f * base applied to group a.
    """
    nb = len(va)
    w = rng.multinomial(nb, np.full(nb, 1.0 / nb), size=n_boot).astype(float)

    def boot(v, a, shift=0.0):
        tv, ta = w @ v, w @ a
        with np.errstate(invalid="ignore", divide="ignore"):
            th = (tv - ta) / (tv + ta) + shift
            lv = tv[:, None] - w * v
            la = ta[:, None] - w * a
            loo = (lv - la) / (lv + la) + shift
        return th, loo

    out = []
    for f in (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
        sa, la_ = boot(va, aa, shift=-kappa * f * base)
        sb, lb_ = boot(vb, ab)
        diff = sa - sb
        dl = la_ - lb_
        ok = np.isfinite(dl)
        g = ok.sum(axis=1)
        dl0 = np.where(ok, dl, 0.0)
        m = dl0.sum(axis=1) / np.maximum(g, 1)
        var = (((dl0 - m[:, None]) * ok) ** 2).sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se = np.sqrt((g - 1) / np.maximum(g, 1) * var)
            z = diff / se
        good = np.isfinite(z) & (se > 0)
        out.append(dict(fraction=f,
                        detect_rate=float(np.mean(np.abs(z[good]) > N_SIGMA))
                        if good.any() else np.nan))
    return pd.DataFrame(out)


# --------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-snps", type=int, default=1)
    ap.add_argument("--boot", type=int, default=N_BOOT)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    log.info(f"Loading Skov S5 from {os.path.relpath(S5, ROOT)}")
    d = load_tracts(min_snps=args.min_snps)
    log.info(f"  {len(d):,} HMM tracts, {d['name'].nunique()} individuals, "
             f"{d['pop'].nunique()} populations")
    log.info("  classes: " + ", ".join(f"{k}={v:,}" for k, v in
                                        d["cls"].value_counts().items()))
    nea = d[d["cls"] == "neanderthal"]
    log.info(f"  {len(nea):,} Neanderthal-affinity tracts used")

    # ---- per-group statistics ----------------------------------------------
    rows = []
    for key, sub in list(nea.groupby("region")) + \
            [(p, s) for p, s in nea.groupby("pop") if s["name"].nunique() >= 2]:
        v, a = block_sums(sub)
        rows.append(dict(group=key,
                         level="region" if key in set(nea["region"]) else "pop",
                         n_ind=int(sub["name"].nunique()), n_tracts=len(sub),
                         V=float(v.sum()), A=float(a.sum()),
                         D_tract=theta(v, a), se_block=jackknife(v, a),
                         se_individual=jackknife_individuals(sub)))
    gdf = pd.DataFrame(rows).sort_values(["level", "group"])
    gdf.to_csv(os.path.join(OUT, "ns_tract_groups.csv"), index=False)
    for _, r in gdf[gdf.level == "region"].iterrows():
        log.info(f"  [{r['group']:20s}] n={r['n_ind']:>3} "
                 f"tracts={r['n_tracts']:>6,} D_tract={r['D_tract']:+.5f} "
                 f"+/- {r['se_block']:.5f} (block) / "
                 f"{r['se_individual']:.5f} (individual)")

    # Positive control: Papuan Denisovan-affinity tracts must pull toward Altai,
    # because Altai carries Denisovan ancestry. Same logic as the genome-wide
    # study's Denisova control, which read -0.147 at Z = -6.7.
    den = d[(d["cls"] == "denisovan") & (d["region"] == "Melanesia")]
    if len(den):
        v, a = block_sums(den)
        log.info(f"  [CONTROL Papuan Denisovan tracts] n={den['name'].nunique()}"
                 f" tracts={len(den):,} D_tract={theta(v, a):+.5f} "
                 f"+/- {jackknife(v, a):.5f}  (must be below the Neanderthal "
                 f"groups: Altai carries Denisovan ancestry)")

    # ---- kappa for this statistic ------------------------------------------
    pool = nea[nea["region"] != "Melanesia"]
    base, kdf = measure_kappa(pool, [0.05, 0.10, 0.20, 0.30, 0.50], rng)
    kappa = float(kdf["kappa"].median())
    log.info(f"  kappa = {kappa:.4f} (median over "
             f"{len(kdf)} injections; flat in f: " +
             ", ".join(f"{f:.2f}->{g['kappa'].median():.3f}"
                       for f, g in kdf.groupby('fraction')) + ")")

    # ---- pairs and the limit ------------------------------------------------
    prows = []
    for lvl, pairs in (("region", HEADLINE), ("pop", POP_PAIRS)):
        col = "region" if lvl == "region" else "pop"
        for a_, b_ in pairs:
            sa, sb = nea[nea[col] == a_], nea[nea[col] == b_]
            if not len(sa) or not len(sb):
                continue
            va, aa = block_sums(sa)
            vb, ab = block_sums(sb)
            r = paired_difference(va, aa, vb, ab)
            f50 = N_SIGMA * r["se"] / (kappa * abs(base)) if r["se"] > 0 else np.nan
            prows.append(dict(level=lvl, a=a_, b=b_, **r,
                              f50_analytic=f50))
            log.info(f"  [{lvl}] {a_:20s} - {b_:20s} "
                     f"{r['diff']:+.5f} +/- {r['se']:.5f} (Z={r['z']:+.2f})  "
                     f"f50~{100*f50:.1f}%")
    pdf = pd.DataFrame(prows)
    pdf.to_csv(os.path.join(OUT, "ns_tract_pairs.csv"), index=False)

    # Empirical curve on the headline regional comparison.
    sa, sb = nea[nea.region == "WestEurasia"], nea[nea.region == "EastAsia"]
    va, aa = block_sums(sa)
    vb, ab = block_sums(sb)
    curve = detection_curve(va, aa, vb, ab, kappa, base, rng, args.boot)
    f50 = sc.power_crossing(curve["fraction"], curve["detect_rate"], 0.50)
    f80 = sc.power_crossing(curve["fraction"], curve["detect_rate"], 0.80)
    log.info("  detection curve (WestEurasia vs EastAsia): " +
             ", ".join(f"{100*r.fraction:.0f}%->{100*r.detect_rate:.0f}%"
                       for r in curve.itertuples()))
    log.info(f"  f50 = {100*f50:.1f}%   f80 = {100*f80:.1f}%   "
             f"(genome-wide: 37% and 64%)")

    ldf = pd.DataFrame([dict(statistic="D_tract", kappa=kappa, base=base,
                             f50=f50, f80=f80,
                             genome_wide_f50=0.374, genome_wide_f80=0.641,
                             gain_f50=0.374 / f50 if f50 else np.nan)])
    ldf.to_csv(os.path.join(OUT, "ns_tract_limit.csv"), index=False)
    curve.to_csv(os.path.join(OUT, "ns_tract_curve.csv"), index=False)

    # ---- would more genomes help? -------------------------------------------
    # The decision that follows this run is whether to go and get tract calls
    # for thousands more individuals (1000G via IBDmix). That is only worth
    # doing if the error bar is set by how many genomes there are. The
    # genome-wide study found it was not - sites bound, genomes did not - and
    # the same question has to be asked again of this statistic rather than
    # assumed to have the same answer. Individuals are thinned within both
    # groups at once, mirroring the two-way subsample.
    trows = []
    for q in (1.0, 0.5, 0.25, 0.125):
        reps = 1 if q == 1.0 else 6
        for rep in range(reps):
            sub = []
            for reg in ("WestEurasia", "EastAsia"):
                names = nea[nea.region == reg]["name"].unique()
                k = max(2, int(round(q * len(names))))
                pick = names if q == 1.0 else rng.choice(names, k, replace=False)
                sub.append(nea[(nea.region == reg) & nea["name"].isin(pick)])
            va_, aa_ = block_sums(sub[0])
            vb_, ab_ = block_sums(sub[1])
            r = paired_difference(va_, aa_, vb_, ab_)
            trows.append(dict(fraction=q, rep=rep, se=r["se"],
                              n_ind=sub[0]["name"].nunique()
                              + sub[1]["name"].nunique(),
                              n_tracts=len(sub[0]) + len(sub[1])))
    tdf = pd.DataFrame(trows)
    tdf.to_csv(os.path.join(OUT, "ns_tract_genome_scaling.csv"), index=False)
    fit = sc.subsample_exponent(tdf["fraction"], tdf["se"])
    g = tdf.groupby("fraction")["se"].mean()
    log.info("  genome-count scaling of the paired SE: " +
             ", ".join(f"q={q:g}->{v:.5f}" for q, v in g.items()))
    log.info(f"  b = {fit['b']:+.3f} +/- {fit['b_se']:.3f}  "
             f"(0.5 = genomes bind, 0 = saturated)")
    # A cutoff on b would turn a continuous quantity into a binary claim, and at
    # b = 0.145 +/- 0.026 - weakly binding, but many sigma from zero - that claim
    # would come out the wrong way round. Report the requirement instead, and
    # only call the axis dead if it is genuinely consistent with zero.
    n_ind_now = int(tdf[tdf.fraction == 1.0]["n_ind"].iloc[0])
    if fit["b"] < 2 * fit["b_se"]:
        log.info("  => genome count is consistent with saturated; more "
                 "individuals would not tighten this limit.")
    else:
        mult = (f50 / 0.10) ** (1.0 / fit["b"])
        log.info(f"  => genomes bind weakly. Reaching f50 = 10% needs about "
                 f"{mult:.1f}x the current {n_ind_now} individuals "
                 f"(~{mult*n_ind_now:.0f}), which 1000G-scale tract calls "
                 f"would supply.")
        for n_target, name in ((2504, "1000G"), (279, "SGDP full")):
            proj = f50 * (n_target / n_ind_now) ** (-fit["b"])
            log.info(f"     projected f50 at {name} scale "
                     f"({n_target} individuals): {100*proj:.1f}%")

    log.info(f"Wrote ns_tract_*.csv to {OUT}")


if __name__ == "__main__":
    main()
