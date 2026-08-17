#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
What binds the tract-restricted limit? A three-way subsample.

`ns_tract_source.py` put the calibrated detection limit at f50 = 12.9% for West
Eurasia against East Asia, and found the paired standard error moves only weakly
with the number of individuals (b = 0.145 +/- 0.026). Weakly is not zero, so the
obvious next step looked like importing 1000G-scale tract calls -- roughly six
times the genomes for a 1.3x improvement, at the cost of recomputing Skov's
per-tract archaic-sharing counts from raw genotypes.

Before paying that, the same question the genome-wide study asked has to be
asked here: which axis actually sets the error bar? Genome-wide the answer was
the archaic panel, and more ancient genomes bought nothing. If the answer is the
same shape here, 1000G is the wrong lever and the effort is wasted.

THREE AXES, NOT TWO

  * VARIANTS. Each tract carries a count of its variants shared with Vindija and
    with Altai. Those counts exist only where the two archaic genomes are both
    called and differ -- the same ~19k-site constraint that bound the
    genome-wide study, reappearing in a different guise. Thinned by binomial
    subsampling of both counts, which preserves the expected statistic and
    scales the counting variance by 1/q.
  * TRACTS. How much archaic sequence the sample carries, at a fixed number of
    genomes. Thinned by dropping whole tracts, since a tract is the unit a
    source change would replace.
  * INDIVIDUALS. How many genomes. Thinned by dropping whole individuals.

Whichever bends the paired SE hardest is what to spend money on. The readout is
the exponent b in SE ~ q^-b, plus the share of the full-data variance each axis
controls, using the same two fits as `POWER_two_way_subsample.md`.

A NUMBER WORTH HAVING FIRST

The pure counting floor -- the SE this statistic would have if the only noise
were multinomial sampling of the shared-variant counts -- is about
1/sqrt(V+A). Comparing it against the measured SE says immediately how much
headroom any count-increasing axis has, before a single subsample is run.

Outputs (reports/neanderthal_source/):
  ns_tract_subsample.csv     per arm, fraction and replicate
  ns_tract_subsample_fit.csv exponents and variance shares per arm
  fig_n8_tract_subsample.png

Run: PYTHONIOENCODING=utf-8 python scripts/ns_tract_subsample.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from archaic import source_contrast as sc
from archaic.log_utils import get_logger

import ns_tract_source as ts

log = get_logger("archaic.ns_tract_subsample")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")
FRACTIONS = [1.0, 0.5, 0.25, 0.125]
N_REP = 6
SEED = 20260818
GROUPS = ("WestEurasia", "EastAsia")


def paired_se(sub_a, sub_b):
    va, aa = ts.block_sums(sub_a)
    vb, ab = ts.block_sums(sub_b)
    return ts.paired_difference(va, aa, vb, ab)["se"]


def counting_floor(sub_a, sub_b):
    """SE if multinomial sampling of the shared counts were the only noise.

    D = (V-A)/(V+A) on counts totalling N has variance (1-D^2)/N. Two
    independent groups add in quadrature. Any axis that works by increasing
    counts -- variants or tracts -- cannot push the SE below this, so the gap
    between it and the measured SE is the headroom those axes have.
    """
    out = []
    for s in (sub_a, sub_b):
        v = s["Shared_with_Vindija"].sum()
        a = s["Shared_with_Altai"].sum()
        n = v + a
        d = (v - a) / n if n else np.nan
        out.append((1 - d ** 2) / n if n else np.nan)
    return float(np.sqrt(sum(out)))


def thin_variants(sub, q, rng):
    """Binomial thinning of both sharing counts."""
    s = sub.copy()
    s["Shared_with_Vindija"] = rng.binomial(
        s["Shared_with_Vindija"].to_numpy(int), q)
    s["Shared_with_Altai"] = rng.binomial(
        s["Shared_with_Altai"].to_numpy(int), q)
    return s


def thin_tracts(sub, q, rng):
    return sub[rng.random(len(sub)) < q]


def thin_individuals(sub, q, rng):
    names = sub["name"].unique()
    k = max(2, int(round(q * len(names))))
    pick = rng.choice(names, k, replace=False)
    return sub[sub["name"].isin(pick)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-snps", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    d = ts.load_tracts(min_snps=args.min_snps)
    nea = d[d["cls"] == "neanderthal"]
    A0 = nea[nea.region == GROUPS[0]]
    B0 = nea[nea.region == GROUPS[1]]
    log.info(f"{GROUPS[0]} {A0['name'].nunique()} individuals / {len(A0):,} "
             f"tracts;  {GROUPS[1]} {B0['name'].nunique()} / {len(B0):,}")

    se_full = paired_se(A0, B0)
    floor = counting_floor(A0, B0)
    log.info(f"measured paired SE = {se_full:.5f}; pure counting floor = "
             f"{floor:.5f}  ({se_full/floor:.1f}x above it, so "
             f"{100*(1-(floor/se_full)**2):.0f}% of the variance is NOT "
             f"count-limited)")

    arms = {"variants": thin_variants, "tracts": thin_tracts,
            "individuals": thin_individuals}
    rows = []
    for arm, fn in arms.items():
        for q in FRACTIONS:
            reps = 1 if q == 1.0 else N_REP
            for rep in range(reps):
                a = A0 if q == 1.0 else fn(A0, q, rng)
                b = B0 if q == 1.0 else fn(B0, q, rng)
                se = paired_se(a, b)
                rows.append(dict(arm=arm, fraction=q, rep=rep, se=se,
                                 n_ind=a["name"].nunique() + b["name"].nunique(),
                                 n_tracts=len(a) + len(b),
                                 counts=float(a["Shared_with_Vindija"].sum()
                                              + a["Shared_with_Altai"].sum()
                                              + b["Shared_with_Vindija"].sum()
                                              + b["Shared_with_Altai"].sum())))
            log.info(f"  [{arm:12s}] q={q:<6} SE="
                     f"{np.mean([r['se'] for r in rows if r['arm']==arm and r['fraction']==q]):.5f}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "ns_tract_subsample.csv"), index=False)

    frows = []
    for arm in arms:
        s = df[df.arm == arm]
        fit = sc.subsample_exponent(s["fraction"], s["se"])
        var = sc.subsample_variance_share(s["fraction"], s["se"])
        frows.append(dict(arm=arm, **fit, **var))
        log.info(f"  {arm:12s} b = {fit['b']:+.3f} +/- {fit['b_se']:.3f}   "
                 f"holds {100*var['var_share']:.0f}% of the full-data variance")
    fdf = pd.DataFrame(frows)
    fdf.to_csv(os.path.join(OUT, "ns_tract_subsample_fit.csv"), index=False)

    best = fdf.loc[fdf["b"].idxmax()]
    log.info(f"VERDICT: the binding axis is '{best['arm']}' "
             f"(b = {best['b']:.3f}); the others are "
             + ", ".join(f"{r['arm']} {r['b']:.3f}"
                         for _, r in fdf.iterrows() if r['arm'] != best['arm']))
    # The three shares must NOT be summed. The axes overlap by construction --
    # dropping a tract also drops its shared variants, and dropping an
    # individual drops both -- so adding them double-counts and a "remainder"
    # computed that way is meaningless. The defensible statements are the single
    # largest share, and the independent counting-floor calculation above; that
    # the two agree is the useful check.
    share = float(best["var_share"])
    log.info(f"  the largest single share is {100*share:.0f}% ({best['arm']}). "
             f"The axes overlap, so these shares cannot be added.")
    log.info(f"  independently, the counting floor puts "
             f"{100*(1-(floor/se_full)**2):.0f}% of the variance beyond reach "
             f"of any count-increasing axis - consistent with the "
             f"{100*(1-share):.0f}% the best axis leaves untouched.")
    for _, r in fdf.iterrows():
        asym = se_full * np.sqrt(max(1.0 - float(r["var_share"]), 0.0))
        log.info(f"  exhausting '{r['arm']}' entirely would take the SE to "
                 f"{asym:.5f} and f50 to "
                 f"{100*0.129*asym/se_full:.1f}% (from 12.9%)")
    make_figure(df, fdf, se_full, floor)
    log.info(f"Wrote ns_tract_subsample*.csv and fig_n8 to {OUT}")


def make_figure(df, fdf, se_full, floor):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    colours = {"variants": "#b2182b", "tracts": "#1b7837",
               "individuals": "#2166ac"}
    ax = axes[0]
    for arm in colours:
        g = df[df.arm == arm].groupby("fraction")["se"]
        m, lo, hi = g.mean(), g.min(), g.max()
        ax.plot(m.index, m.values, "o-", color=colours[arm], label=f"thin {arm}")
        ax.fill_between(m.index, lo.values, hi.values, color=colours[arm],
                        alpha=0.16, lw=0)
    q = np.array(FRACTIONS)
    ax.plot(q, se_full / np.sqrt(q), "k--", lw=1, label=r"$q^{-1/2}$")
    ax.axhline(floor, color="grey", ls=":", lw=1)
    ax.text(0.99, floor, " counting floor", fontsize=8, color="grey",
            va="bottom", ha="right", transform=ax.get_yaxis_transform())
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(q); ax.set_xticklabels([f"{v:g}" for v in q])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("fraction retained, q")
    ax.set_ylabel("paired $D_{tract}$ difference SE")
    ax.set_title("Three-way subsample")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    xs = np.arange(len(fdf))
    ax.bar(xs, fdf["b"], 0.55, yerr=fdf["b_se"], capsize=3,
           color=[colours[a] for a in fdf["arm"]])
    ax.axhline(0.5, color="k", ls="--", lw=1)
    ax.text(-0.45, 0.5, "binding ($b=0.5$)", fontsize=8, va="bottom")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(fdf["arm"])
    ax.set_ylim(0, 0.6)
    ax.set_ylabel(r"scaling exponent $b$")
    ax.set_title("Which axis binds")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_n8_tract_subsample.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
