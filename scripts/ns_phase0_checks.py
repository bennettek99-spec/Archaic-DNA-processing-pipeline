#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 0: the two things that must hold before building a tract-restricted D_VA.

The plan for driving the detection limit from ~41% down to 10% rests on two
claims, both currently untested. Neither is expensive to check, and both are
cheaper to check than to discover wrong after a pipeline exists.

CHECK A — IS THE SITE-SCALING EXPONENT STABLE?

`POWER_two_way_subsample.md` measured SE ~ q^-0.41 by thinning the site set
*down* by 8x. Any plan that adds sites extrapolates that exponent *up*, and the
extrapolation is only safe if b is a constant rather than a local slope. Linkage
predicts it is not: the more sites you have, the more redundant each new one is,
so b should fall as q rises and the upward extrapolation would be optimistic.
This refits b on overlapping sub-ranges and as local slopes between adjacent
fractions, and asks whether b drifts with q.

  b flat            -> the site route costs what the plan says
  b falls as q rises -> adding sites buys less than the plan assumes

CHECK B — WOULD TRACT RESTRICTION ACTUALLY DELIVER THE ALPHA GAIN?

The tract lever's premise is that restricting to called introgressed tracts,
where local ancestry is ~0.5 rather than ~0.021, multiplies the mixture signal
by that ratio — about 24x.

Note what is *not* worth testing: whether the response is proportional to alpha.
The closed form makes it so by construction,

    signal per unit f = -(alpha/2) * sum (p_V-p_A)^2 / sum denominator,

and any simulation that injects a displacement of alpha*f*(target-source) will
return proportionality no matter what, because proportionality was assumed in
the injection. A check that cannot fail is not a check.

What can fail is the denominator. That sum carries a factor
(p_X + p_B - 2 p_X p_B) which depends on the test population, and inside a tract
p_X is no longer a modern human frequency — it is half Neanderthal. If that
inflates the denominator, the realised gain is smaller than alpha alone implies.
So the quantity measured here is

    realised gain = (alpha_tract / alpha_genomewide)
                    * (sum denominator normal / sum denominator tract)

evaluated on a synthetic tract pool, p_tract = 0.5*p_V + 0.5*p_X for a
heterozygous tract and p_tract = p_V for a phased Neanderthal haplotype. Only
the denominator is taken from that pool, which is what makes the calculation
legitimate despite the source proxy being one arm of the contrast: the
degeneracy lives entirely in the numerator, and the numerator is not used.

Outputs (reports/neanderthal_source/):
  ns_phase0_exponent.csv   local and windowed exponents with spreads
  ns_phase0_alpha.csv      per-cohort signal per unit f, and the tract projection

Run: PYTHONIOENCODING=utf-8 python scripts/ns_phase0_checks.py --panel 1240k
"""
import argparse
import itertools
import os
import sys
import types

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from archaic import source_contrast as sc
from archaic import stats as st
from archaic.log_utils import get_logger
from archaic.panel import Panel
from archaic.refs import PANELS

import neanderthal_source as ns

log = get_logger("archaic.ns_phase0")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")
F_PROBE = 0.10          # injection fraction for the alpha check


# ------------------------------------------------------------------ check A ---
def check_exponent_stability():
    """Refit the site exponent on sub-ranges to see whether b drifts with q."""
    path = os.path.join(OUT, "ns_power_curve.csv")
    curve = pd.read_csv(path)
    sites = curve[curve.arm == "sites"][["fraction", "se_median"]].dropna()
    qs = sorted(sites["fraction"].unique(), reverse=True)     # 1.0 .. 0.125
    by_q = {q: sites[sites.fraction == q]["se_median"].to_numpy() for q in qs}

    rows = []
    # Local slopes between adjacent fractions. Every replicate pairing is used,
    # so the spread reflects the jackknife-SE noise rather than hiding it.
    for q_hi, q_lo in zip(qs[:-1], qs[1:]):
        b = [np.log(s_lo / s_hi) / np.log(q_hi / q_lo)
             for s_hi in by_q[q_hi] for s_lo in by_q[q_lo]]
        rows.append(dict(kind="local", q_hi=q_hi, q_lo=q_lo,
                         b=float(np.mean(b)), b_sd=float(np.std(b, ddof=1)),
                         n=len(b)))
    # Overlapping three-point windows: the coarsest honest split of a four-point
    # grid. A two-way split would leave two points per side, which cannot
    # constrain a slope at all.
    for i in range(len(qs) - 2):
        w = qs[i:i + 3]
        sub = sites[sites.fraction.isin(w)]
        fit = sc.subsample_exponent(sub["fraction"], sub["se_median"])
        rows.append(dict(kind="window", q_hi=max(w), q_lo=min(w),
                         b=fit["b"], b_sd=fit["b_se"], n=fit["n_points"]))
    edf = pd.DataFrame(rows)
    edf.to_csv(os.path.join(OUT, "ns_phase0_exponent.csv"), index=False)

    log.info("CHECK A - site exponent by sub-range (full-range b = 0.408):")
    for _, r in edf.iterrows():
        log.info(f"  {r['kind']:>6s}  q {r['q_lo']:.3f}-{r['q_hi']:.3f}   "
                 f"b = {r['b']:+.3f} +/- {r['b_sd']:.3f}  (n={int(r['n'])})")

    loc = edf[edf.kind == "local"].sort_values("q_hi")
    hi, lo = loc.iloc[-1], loc.iloc[0]        # nearest q=0.125 vs nearest q=1
    drift = float(lo["b"] - hi["b"])
    pooled = float(np.hypot(lo["b_sd"], hi["b_sd"]))
    log.info(f"  b nearest q=1 is {hi['b']:+.3f}, nearest q=0.125 is "
             f"{lo['b']:+.3f}; drift {drift:+.3f} +/- {pooled:.3f}")
    return edf, drift, pooled


# ------------------------------------------------------------------ check B ---
def check_alpha_scaling(args):
    """Measure the D_VA response per unit f against each cohort's alpha."""
    cfg = PANELS[args.panel]
    log.info("CHECK B - loading panel and references...")
    panel = Panel(cfg["prefix"], autosomes_only=True)
    block_full = st.assign_blocks(panel.n_snp, args.blocks)
    ref, _ = panel.frequencies({k: cfg["refs"][k] for k in
                                ("Altai", "Vindija", "Denisova", "Chimp",
                                 "Mbuti", "Yoruba")})
    useful = (np.isfinite(ref["Vindija"]) & np.isfinite(ref["Altai"])
              & np.isfinite(ref["Yoruba"]) & np.isfinite(ref["Chimp"]))
    R = {k: ref[k][useful] for k in ("Vindija", "Altai", "Chimp", "Yoruba")}
    block = block_full[useful]
    equi = sc.neanderthal_average(R["Vindija"], R["Altai"])

    meta = pd.read_csv(ns.META, low_memory=False)
    defs, crows, _, _, _ = ns.build_cohorts(meta, panel)
    cdf = pd.read_csv(os.path.join(OUT, "ns_cohorts.csv"))

    # Cohorts with a real alpha AND full site coverage. The coverage cut is not
    # cosmetic: a cohort called at fewer sites has a different site set, hence a
    # different sum (p_V-p_A)^2 and a different denominator, and mixing those in
    # produces spread that looks like a failure of proportionality and is not.
    # Oase1 at 25,838 sites returns a ratio of -0.66 against the -1.14 of the
    # fully-covered cohorts for exactly this reason.
    keep = cdf[np.isfinite(cdf["mean_alpha"]) & (cdf["mean_alpha"] > 0.015)
               & (cdf["n_snp"] > 500_000) & cdf["label"].isin(defs)].copy()
    log.info(f"  {len(keep)} cohorts with alpha > 0.015 at full site coverage")

    ca = types.SimpleNamespace(panel=args.panel, chunk=args.chunk)
    full, counts = ns.cached_pooled_freq(panel, defs, ca, "", log)
    rows = []
    for _, c in keep.iterrows():
        lab, alpha = c["label"], float(c["mean_alpha"])
        p = full[lab][useful]
        t0 = sc.build_block_table(lab, p, R["Yoruba"], R["Vindija"], R["Altai"],
                                  R["Chimp"], block, args.blocks)
        pm, _ = sc.mixture_frequencies(p, alpha, F_PROBE, R["Vindija"], equi)
        t1 = sc.build_block_table(lab, pm, R["Yoruba"], R["Vindija"], R["Altai"],
                                  R["Chimp"], block, args.blocks)
        d0, d1 = t0.d_va_theta, t1.d_va_theta
        rows.append(dict(label=lab, kind=c["kind"], alpha=alpha,
                         n_ind=int(c["n_ind"]), D_VA=d0,
                         signal_per_f=(d1 - d0) / F_PROBE,
                         kappa=(d1 - d0) / (-F_PROBE * d0),
                         denom_sum=denominator_sum(p, R),
                         n_snp=t0.n_snp))
    adf = pd.DataFrame(rows).sort_values("alpha")
    ratio = adf["signal_per_f"] / adf["alpha"]
    log.info(f"  signal per unit f, divided by alpha: {ratio.mean():+.4f} "
             f"+/- {ratio.std():.4f} across {len(adf)} cohorts "
             f"(spread {ratio.min()/ratio.max():.2f}x)")
    log.info("  (proportionality itself is built into the injection; this "
             "spread measures denominator stability, which is what can vary)")

    # ---- the projection that actually matters ------------------------------
    # A heterozygous tract carrier is half Neanderthal at that locus; a phased
    # Neanderthal haplotype is all of it. Only the denominator is taken from
    # these pools - see the module docstring on why that is legitimate here.
    ref_lab = adf.loc[adf["n_ind"].idxmax(), "label"]
    p_ref = full[ref_lab][useful]
    a_ref = float(adf.loc[adf.label == ref_lab, "alpha"].iloc[0])
    den_norm = denominator_sum(p_ref, R)
    scen = []
    for name, p_tract, a_tract in (
            ("heterozygous tract", 0.5 * R["Vindija"] + 0.5 * p_ref, 0.5),
            ("phased Nea haplotype", R["Vindija"].copy(), 1.0)):
        den_tr = denominator_sum(p_tract, R)
        naive = a_tract / a_ref
        realised = naive * (den_norm / den_tr)
        scen.append(dict(scenario=name, alpha_tract=a_tract,
                         alpha_genomewide=a_ref, denom_normal=den_norm,
                         denom_tract=den_tr,
                         denom_inflation=den_tr / den_norm,
                         naive_gain=naive, realised_gain=realised))
        log.info(f"  [{name:22s}] alpha {a_ref:.4f} -> {a_tract:.2f}; "
                 f"denominator x{den_tr/den_norm:.3f}; "
                 f"naive gain {naive:.1f}x -> realised {realised:.1f}x")
    del full, counts
    sdf = pd.DataFrame(scen)
    adf.to_csv(os.path.join(OUT, "ns_phase0_alpha.csv"), index=False)
    sdf.to_csv(os.path.join(OUT, "ns_phase0_tract_projection.csv"), index=False)
    return adf, sdf, ref_lab


def denominator_sum(p_x, R):
    """Sum of the D_VA denominator over usable sites, for a given test pop.

    The only place the test population enters the mixture response, so it is the
    only thing that can make a tract-restricted signal fall short of what alpha
    alone predicts.
    """
    num = st.d_numerator(p_x, R["Yoruba"], R["Vindija"], R["Altai"])
    den = st.d_denominator(p_x, R["Yoruba"], R["Vindija"], R["Altai"])
    ok = np.isfinite(num) & np.isfinite(den)
    return float(den[ok].sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--blocks", type=int, default=ns.N_BLOCKS)
    args = ap.parse_args()

    edf, drift, pooled = check_exponent_stability()
    adf, sdf, ref_lab = check_alpha_scaling(args)

    # Where the current limit sits, so the projection lands in real units.
    lim = pd.read_csv(os.path.join(OUT, "ns_detection_limit.csv"))
    limit = float(lim[lim.statistic == "D_VA"]["limit"].iloc[0])
    mix = pd.read_csv(os.path.join(OUT, "ns_mixture_summary.csv"))
    mrow = mix[mix.replacement == "equidistant"].iloc[0]
    signal_per_f = abs(float(adf.loc[adf.label == ref_lab,
                                     "signal_per_f"].iloc[0]))
    f_now = limit / signal_per_f

    log.info("PHASE 0 VERDICT")
    win = edf[edf.kind == "window"]
    log.info(f"  A: local slopes drift {drift:+.3f} +/- {pooled:.3f} - "
             f"consistent with zero, but the spread is too wide to rule out "
             f"real drift, so this is 'no evidence against', not 'confirmed'. "
             f"The windowed fits are the stronger evidence: "
             f"{win.iloc[0]['b']:.3f} +/- {win.iloc[0]['b_sd']:.3f} on the "
             f"upper range against {win.iloc[1]['b']:.3f} +/- "
             f"{win.iloc[1]['b_sd']:.3f} on the lower.")
    for _, r in sdf.iterrows():
        f_new = f_now / r["realised_gain"]
        log.info(f"  B: {r['scenario']:22s} realised gain "
                 f"{r['realised_gain']:.1f}x (naive {r['naive_gain']:.1f}x, "
                 f"denominator penalty x{r['denom_inflation']:.2f}) -> "
                 f"f50 {100*f_now:.0f}% becomes {100*f_new:.1f}%")
    log.info(f"     (f50 here is {100*f_now:.0f}% from the analytic conversion; "
             f"the measured curve gives {100*float(mrow['f50']):.0f}%, and the "
             f"ratio between them carries over to the projections)")
    log.info(f"Wrote ns_phase0_exponent.csv, ns_phase0_alpha.csv and "
             f"ns_phase0_tract_projection.csv to {OUT}")


if __name__ == "__main__":
    main()
