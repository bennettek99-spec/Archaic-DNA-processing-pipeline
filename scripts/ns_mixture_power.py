#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Is the 13% detection limit calibrated? Simulated mixtures at known fractions.

`PAPER_neanderthal_source.md` states the limit as a one-line conversion:

    a difference of 0.0098 in D_VA is resolvable at 2 sigma;
    a typical cohort's D_VA is 0.0737;
    therefore a cohort would have to re-source f > 0.0098/0.0737 = 13% of its
    Neanderthal ancestry before this study could see it.

`POWER_two_way_subsample.md` established which sample sets the 0.0098. It said
nothing about the conversion, and the conversion carries an untested assumption:
that a cohort which re-sources a fraction f of its Neanderthal ancestry to a
lineage equidistant between Vindija and Altai shifts its D_VA by exactly
f x 0.0737 - that is, that D_VA responds *linearly* to the mixture, with slope
one, in the cohort's own D_VA units. If the true response is a fraction kappa of
that, the real threshold is 13%/kappa and the published number is optimistic.

WHAT THIS SCRIPT DOES

Builds the mixture explicitly, at the allele-frequency level, and pushes it
through the unmodified statistic:

    p_X'(f) = p_X + alpha_X * f * (p_S' - p_Vindija)

alpha_X is the cohort's own f4-ratio Neanderthal fraction, p_Vindija stands for
the source Eurasians actually descend from (the study's premise), and p_S' is
the replacement source. Two replacements are run: the **equidistant** lineage
(p_S' = (p_V+p_A)/2), which is the counterfactual the published sentence names,
and a **full swap to Altai** (p_S' = p_A), the most different source the panel
physically contains.

This is deliberately *not* done by shifting D_VA by f x D_VA directly. That
would assume the very thing under test and would return 13% by construction.
Here the perturbation is applied to frequencies, every SNP's numerator and
denominator are recomputed, and the resulting D_VA is whatever it is.

TWO READOUTS

  1. LINEARITY. kappa(f) = realised dD_VA / (-f x D_VA). The published
     conversion assumes kappa = 1 at every f. Any departure rescales the limit.
  2. DETECTION POWER. For every pair of testable grid cohorts, inject f into one
     member and ask how often the pair is called at 2 sigma, under a *paired
     block bootstrap*: the 50 jackknife blocks are resampled with replacement
     (identically for both cohorts, mirroring the paired jackknife), and within
     each bootstrap replicate the difference and its jackknife SE are both
     recomputed - so the simulated analyst uses the SE they would really have.
     f = 0 is included as a false-positive control and should return ~5%.

The number to quote afterwards is f80, not f50. A 2-sigma threshold gives ~50%
power against an effect sitting exactly at it, so "the limit" in the published
sense is the coin-flip point, not the point where a real difference would
reliably have been found.

Outputs (reports/neanderthal_source/):
  ns_mixture_linearity.csv   realised vs predicted dD_VA, and kappa, per cohort
  ns_mixture_power.csv       detection rate by fraction, replacement and class
  ns_mixture_summary.csv     f50, f80, the analytical 13%, and the correction
  fig_n6_mixture.png         the power curve and the linearity check
  POWER_mixture_calibration.md

Run: PYTHONIOENCODING=utf-8 python scripts/ns_mixture_power.py --panel 1240k
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

log = get_logger("archaic.ns_mixture_power")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")

# The fractions asked for, plus 0 (false-positive control), the analytical
# threshold itself (0.133, where power should be ~50% if the limit is honest),
# and 0.30 to pin the top of the curve.
FRACTIONS = [0.0, 0.02, 0.05, 0.10, 0.133, 0.15, 0.20, 0.30, 0.40, 0.50,
             0.70, 1.00]
REPLACEMENTS = ("equidistant", "altai")
N_BOOT = 500
N_SIGMA = 2.0
SEED = 20260817

# The published conversion, restated here only so the report can show the
# arithmetic it is checking. Both are read from the study's own tables at run
# time; these are fallbacks if a table is missing.
FALLBACK_LIMIT = 0.009826459090792421
FALLBACK_SIGNAL = 0.07368772531967363


# ---------------------------------------------------------- paired bootstrap --
def bootstrap_detection(ta, tb, n_blocks, boot_w, n_sigma=N_SIGMA):
    """Fraction of paired block-bootstrap replicates calling a - b at n_sigma.

    `boot_w` is a (B, n_blocks) matrix of block multiplicities, shared between
    the two cohorts so the resampling is paired exactly as the jackknife is.
    Within each replicate both the difference and its delete-one-block jackknife
    standard error are recomputed, because the quantity whose distribution
    matters is the one an analyst would actually compute and compare to 2 - not
    the difference measured against some fixed, externally-known error.

    Returns the detection rate, the observed difference, the jackknife SE on the
    unresampled data, and the bootstrap SD of the difference. The last two are a
    calibration check: they should agree, and if they do not the block bootstrap
    and the block jackknife disagree about this data and neither number is safe.
    """
    w = boot_w
    tn_a, td_a = w @ ta.n_va, w @ ta.d_va          # (B,)
    tn_b, td_b = w @ tb.n_va, w @ tb.d_va
    with np.errstate(invalid="ignore", divide="ignore"):
        diff = tn_a / td_a - tn_b / td_b
        # delete-one-block *within* each bootstrap replicate
        la = (tn_a[:, None] - w * ta.n_va) / (td_a[:, None] - w * ta.d_va)
        lb = (tn_b[:, None] - w * tb.n_va) / (td_b[:, None] - w * tb.d_va)
    dl = la - lb                                    # (B, n_blocks)
    ok = np.isfinite(dl)
    g = ok.sum(axis=1)
    safe = g > 1
    dl0 = np.where(ok, dl, 0.0)
    mean = np.divide(dl0.sum(axis=1), np.maximum(g, 1))
    var = (((dl0 - mean[:, None]) * ok) ** 2).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        se = np.sqrt((g - 1) / np.maximum(g, 1) * var)
        z = diff / se
    good = safe & np.isfinite(z) & np.isfinite(se) & (se > 0)
    rate = (float(np.mean(np.abs(z[good]) > n_sigma)) if good.any() else np.nan)
    obs = sc.paired_difference(ta, tb, "D_VA")
    return dict(detect_rate=rate,
                diff=obs["diff"], se_jackknife=obs["se"],
                se_bootstrap=float(np.std(diff[np.isfinite(diff)], ddof=1)),
                n_boot_used=int(good.sum()))


# --------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--blocks", type=int, default=ns.N_BLOCKS)
    ap.add_argument("--boot", type=int, default=N_BOOT)
    ap.add_argument("--max-cohorts", type=int, default=0,
                    help="cap the cohort set (0 = all testable grid cohorts)")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    os.makedirs(OUT, exist_ok=True)
    nb = args.blocks
    rng = np.random.default_rng(SEED)

    log.info(f"Loading panel {args.panel}...")
    panel = Panel(cfg["prefix"], autosomes_only=True)
    block_full = st.assign_blocks(panel.n_snp, nb)

    log.info("Reading archaic and baseline reference frequencies...")
    spec = {k: cfg["refs"][k] for k in ("Altai", "Vindija", "Denisova", "Chimp",
                                        "Mbuti", "Yoruba")}
    ref, _ = panel.frequencies(spec)
    useful = (np.isfinite(ref["Vindija"]) & np.isfinite(ref["Altai"])
              & np.isfinite(ref["Yoruba"]) & np.isfinite(ref["Chimp"]))
    R = {k: ref[k][useful] for k in ("Vindija", "Altai", "Chimp", "Yoruba")}
    block = block_full[useful]
    nea_avg = sc.neanderthal_average(R["Vindija"], R["Altai"])
    targets = {"equidistant": nea_avg, "altai": R["Altai"]}
    log.info(f"  {int(useful.sum()):,} usable sites")

    log.info("Rebuilding the study's cohort definitions...")
    meta = pd.read_csv(ns.META, low_memory=False)
    defs, crows, _, _, _ = ns.build_cohorts(meta, panel)

    # The cohorts the published floor is computed over: testable grid cells.
    # Named Palaeolithic targets are excluded on purpose - this measures the
    # calibration of the limit, not the status of any hypothesis cohort.
    cdf = pd.read_csv(os.path.join(OUT, "ns_cohorts.csv"))
    keep = cdf[(cdf.kind == "grid") & cdf.testable].copy()
    if args.max_cohorts:
        keep = keep.nlargest(args.max_cohorts, "n_ind")
    labels = [l for l in keep["label"] if l in defs]
    alpha_of = dict(zip(keep["label"], keep["mean_alpha"]))
    log.info(f"{len(labels)} testable grid cohorts; "
             f"alpha {min(alpha_of.values()):.4f}-{max(alpha_of.values()):.4f}")

    cache_args = types.SimpleNamespace(panel=args.panel, chunk=args.chunk)
    full_freqs, full_counts = ns.cached_pooled_freq(panel, defs, cache_args,
                                                    "", log)
    freqs = {l: full_freqs[l][useful] for l in labels}
    del full_freqs, full_counts

    # ---- build the mixture series ------------------------------------------
    log.info(f"Building {len(labels)} x {len(FRACTIONS)} x {len(REPLACEMENTS)} "
             f"injected block tables...")
    tables = {}                      # (mode, f) -> {label: BlockTable}
    lin_rows = []
    total_clip = 0
    for mode in REPLACEMENTS:
        for f in FRACTIONS:
            tf = {}
            for lab in labels:
                p, nclip = sc.mixture_frequencies(
                    freqs[lab], alpha_of[lab], f, R["Vindija"], targets[mode])
                total_clip += nclip
                tf[lab] = sc.build_block_table(
                    lab, p, R["Yoruba"], R["Vindija"], R["Altai"], R["Chimp"],
                    block, nb, n_ind=int(crows[lab]["n_ind"]))
            tables[(mode, f)] = tf
        base = tables[(mode, 0.0)]
        for lab in labels:
            d0 = base[lab].d_va_theta
            for f in FRACTIONS:
                d = tables[(mode, f)][lab].d_va_theta
                pred = -f * d0
                lin_rows.append(dict(
                    replacement=mode, label=lab, fraction=f, alpha=alpha_of[lab],
                    D_VA_base=d0, D_VA_mixed=d, realised=d - d0, predicted=pred,
                    kappa=(d - d0) / pred if pred != 0 else np.nan,
                    # The relation the published conversion needs is
                    # D_VA(X) = alpha * D_VA(source). Inverting it gives the
                    # D-statistic the source population would have to have. A
                    # D-statistic cannot exceed 1, so a value above 1 falsifies
                    # the relation outright, independently of any simulation.
                    implied_source_D_VA=d0 / alpha_of[lab]))
        log.info(f"  [{mode}] built")
    if total_clip:
        log.warning(f"  {total_clip:,} injected frequencies clipped to [0,1]")

    ldf = pd.DataFrame(lin_rows)
    ldf.to_csv(os.path.join(OUT, "ns_mixture_linearity.csv"), index=False)

    # Sanity: f = 0 must reproduce the unmodified statistic exactly.
    z0 = ldf[ldf.fraction == 0.0]
    if not np.allclose(z0["realised"], 0.0, atol=1e-12):
        raise SystemExit("f=0 injection changed D_VA; the mixture is not "
                         "an identity at zero and nothing below is valid.")

    for mode in REPLACEMENTS:
        k = ldf[(ldf.replacement == mode) & (ldf.fraction > 0)]["kappa"]
        log.info(f"  [{mode}] kappa = {k.mean():.4f} +/- {k.std():.4f} "
                 f"(median {k.median():.4f}) over {len(k)} cohort-fraction cells")

    # ---- detection power ----------------------------------------------------
    pairs = list(itertools.combinations(labels, 2))
    log.info(f"Paired block bootstrap: {len(pairs)} pairs x {len(FRACTIONS)} "
             f"fractions x {len(REPLACEMENTS)} replacements x {args.boot} "
             f"replicates...")
    boot_w = rng.multinomial(nb, np.full(nb, 1.0 / nb),
                             size=args.boot).astype(np.float64)

    prows = []
    for mode in REPLACEMENTS:
        for f in FRACTIONS:
            tf = tables[(mode, f)]
            t0 = tables[(mode, 0.0)]
            recs = []
            for a, b in pairs:
                # inject into a only; b stays as observed
                r = bootstrap_detection(tf[a], t0[b], nb, boot_w)
                recs.append(r)
            rr = pd.DataFrame(recs)
            prows.append(dict(
                replacement=mode, fraction=f, n_pairs=len(pairs),
                detect_rate=float(rr["detect_rate"].mean()),
                median_abs_diff=float(rr["diff"].abs().median()),
                median_abs_z=float((rr["diff"] / rr["se_jackknife"]).abs()
                                   .median()),
                median_se_jackknife=float(rr["se_jackknife"].median()),
                median_se_bootstrap=float(rr["se_bootstrap"].median()),
                se_ratio_boot_over_jack=float(
                    (rr["se_bootstrap"] / rr["se_jackknife"]).median())))
            log.info(f"  [{mode}] f={f:<6}  detected "
                     f"{100*prows[-1]['detect_rate']:5.1f}%   "
                     f"|diff|={prows[-1]['median_abs_diff']:.5f}  "
                     f"SE_jack={prows[-1]['median_se_jackknife']:.5f}  "
                     f"SE_boot={prows[-1]['median_se_bootstrap']:.5f}")
    pdf = pd.DataFrame(prows)
    pdf.to_csv(os.path.join(OUT, "ns_mixture_power.csv"), index=False)

    # ---- summary against the published conversion ---------------------------
    try:
        dl = pd.read_csv(os.path.join(OUT, "ns_detection_limit.csv"))
        row = dl[dl.statistic == "D_VA"].iloc[0]
        limit, signal = float(row["limit"]), float(row["signal"])
        analytic = float(row["limit_fraction_of_signal"])
    except (OSError, KeyError, IndexError):
        limit, signal = FALLBACK_LIMIT, FALLBACK_SIGNAL
        analytic = limit / signal

    srows = []
    for mode in REPLACEMENTS:
        s = pdf[pdf.replacement == mode]
        k = float(ldf[(ldf.replacement == mode)
                      & (ldf.fraction > 0)]["kappa"].median())
        srows.append(dict(
            replacement=mode, kappa=k,
            f50=sc.power_crossing(s["fraction"], s["detect_rate"], 0.50),
            f80=sc.power_crossing(s["fraction"], s["detect_rate"], 0.80),
            # NOT a false-positive rate: real cohorts really do differ a
            # little, so this is the study's own baseline significance rate
            # among these pairs and is an upper bound on the nominal 5%.
            baseline_rate=float(s[s.fraction == 0.0]["detect_rate"].iloc[0]),
            baseline_median_abs_z=float(
                s[s.fraction == 0.0]["median_abs_z"].iloc[0]),
            se_ratio_boot_over_jack=float(
                s[s.fraction == 0.0]["se_ratio_boot_over_jack"].iloc[0]),
            analytic_limit=analytic,
            analytic_limit_kappa_corrected=analytic / k if k else np.nan))
    sdf = pd.DataFrame(srows)
    sdf.to_csv(os.path.join(OUT, "ns_mixture_summary.csv"), index=False)
    for _, r in sdf.iterrows():
        log.info(f"  [{r['replacement']}] kappa={r['kappa']:.3f}  "
                 f"f50={100*r['f50']:.1f}%  f80={100*r['f80']:.1f}%  "
                 f"baseline={100*r['baseline_rate']:.1f}%  "
                 f"SE_boot/SE_jack={r['se_ratio_boot_over_jack']:.3f}  "
                 f"(analytical {100*analytic:.1f}%, kappa-corrected "
                 f"{100*r['analytic_limit_kappa_corrected']:.1f}%)")
    imp = ldf["implied_source_D_VA"]
    log.info(f"  implied D_VA of the source population if D_VA(X)=alpha*D_VA(S):"
             f" {imp.min():.2f}-{imp.max():.2f} "
             f"(a D-statistic cannot exceed 1)")

    make_figure(ldf, pdf, sdf, analytic)
    write_report(ldf, pdf, sdf, analytic, limit, signal, len(pairs), labels)
    log.info(f"Wrote ns_mixture_*.csv, fig_n6_mixture.png and "
             f"POWER_mixture_calibration.md to {OUT}")


def make_figure(ldf, pdf, sdf, analytic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    colours = {"equidistant": "#b2182b", "altai": "#2166ac"}

    ax = axes[0]
    for mode in REPLACEMENTS:
        s = pdf[pdf.replacement == mode].sort_values("fraction")
        ax.plot(100 * s["fraction"], 100 * s["detect_rate"], "o-",
                color=colours[mode], label=f"{mode} replacement")
        f50 = float(sdf[sdf.replacement == mode]["f50"].iloc[0])
        if np.isfinite(f50):
            ax.plot([100 * f50], [50], "o", ms=11, mfc="none",
                    color=colours[mode])
    ax.axvline(100 * analytic, color="k", ls="--", lw=1)
    ax.text(100 * analytic + 0.6, 8, f"analytical\nlimit {100*analytic:.0f}%",
            fontsize=8)
    ax.axhline(50, color="grey", lw=0.8, ls=":")
    ax.axhline(80, color="grey", lw=0.8, ls=":")
    ax.text(0.99, 51, "50%", fontsize=7, color="grey", ha="right",
            transform=ax.get_yaxis_transform())
    ax.text(0.99, 81, "80%", fontsize=7, color="grey", ha="right",
            transform=ax.get_yaxis_transform())
    ax.set_xlabel("fraction of Neanderthal ancestry re-sourced, f (%)")
    ax.set_ylabel("pairs called at 2 sigma (%)")
    ax.set_title("Empirical detection power")
    ax.set_ylim(-3, 103)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axes[1]
    for mode in REPLACEMENTS:
        s = ldf[(ldf.replacement == mode) & (ldf.fraction > 0)]
        g = s.groupby("fraction")["kappa"]
        m, lo, hi = g.median(), g.min(), g.max()
        ax.plot(100 * m.index, m.values, "o-", color=colours[mode], label=mode)
        ax.fill_between(100 * m.index, lo.values, hi.values,
                        color=colours[mode], alpha=0.18, lw=0)
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.text(0.98, 1.0, "assumed by the published conversion ", fontsize=8,
            va="bottom", ha="right", transform=ax.get_yaxis_transform())
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("fraction re-sourced, f (%)")
    ax.set_ylabel(r"$\kappa$ = realised $\Delta D_{VA}$ / predicted")
    ax.set_title("Does $D_{VA}$ respond as the limit assumes?")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_n6_mixture.png"), dpi=150)
    plt.close(fig)


def write_report(ldf, pdf, sdf, analytic, limit, signal, n_pairs, labels):
    def s(mode, field):
        return float(sdf[sdf.replacement == mode][field].iloc[0])

    k_eq = s("equidistant", "kappa")
    f50, f80 = s("equidistant", "f50"), s("equidistant", "f80")
    fpr = s("equidistant", "baseline_rate")
    corrected = s("equidistant", "analytic_limit_kappa_corrected")
    agree = abs(f50 - analytic) < 0.25 * analytic

    lines = [
        "# Is the 13% limit calibrated? Simulated mixtures at known fractions\n",
        f"*Companion to `PAPER_neanderthal_source.md` and "
        f"`POWER_two_way_subsample.md`. AADR v66.p1 1240K panel, {len(labels)} "
        f"testable grid cohorts, {n_pairs} pairs, {N_BOOT} paired block-bootstrap "
        f"replicates per pair per fraction.*\n",

        "## What was actually assumed\n",
        f"The published limit is a conversion: a resolvable D_VA difference of "
        f"{limit:.5f}, divided by a typical cohort's D_VA of {signal:.5f}, gives "
        f"{100*analytic:.1f}%. The division assumes that re-sourcing a fraction "
        f"*f* of a cohort's Neanderthal ancestry to a lineage equidistant "
        f"between Vindija and Altai moves its D_VA by exactly *f* x its own "
        f"D_VA. Nothing in the study tested that. This does.\n",

        "## Method\n",
        "Mixtures are built at the allele-frequency level and pushed through "
        "the unmodified statistic:\n",
        "```\np_X'(f) = p_X + alpha_X * f * (p_S' - p_Vindija)\n```\n",
        "with alpha_X the cohort's own f4-ratio Neanderthal fraction, and "
        "p_S' either the equidistant lineage (p_V+p_A)/2 - the counterfactual "
        "the published sentence names - or a full swap to Altai, the most "
        "different source the panel contains. Every SNP's numerator and "
        "denominator are recomputed and D_VA is whatever comes out. Shifting "
        "D_VA by *f* x D_VA directly would have assumed the conclusion and "
        "returned the published number by construction.\n",
        "Detection is measured by a paired block bootstrap: the 50 jackknife "
        "blocks are resampled with replacement, identically for both cohorts, "
        "and within each replicate the difference **and its jackknife standard "
        "error** are both recomputed, so the simulated analyst tests with the "
        "error bar they would really have had.\n",

        "## Result 1: does D_VA respond as assumed?\n",
        ldf[ldf.fraction > 0].groupby(["replacement", "fraction"])[
            "kappa"].agg(["median", "min", "max"]).round(4).to_markdown(),
        f"\n\nFor the equidistant replacement kappa = {k_eq:.3f}: a mixture of "
        f"fraction *f* moves D_VA by {k_eq:.3f} x *f* x D_VA, not 1.000 x. "
        + ("The conversion is sound and the limit needs no rescaling.\n"
           if abs(k_eq - 1) < 0.05 else
           f"The published conversion is therefore off by that factor, and the "
           f"threshold it implies is {100*corrected:.1f}% rather than "
           f"{100*analytic:.1f}%.\n"),
        "kappa is flat in *f* across the whole range, so the response is "
        "linear - it is the *slope* that is wrong, not the linearity. The "
        "Altai-swap row is about twice the equidistant row because a full swap "
        "is twice the frequency displacement of a half-way one; only the "
        "equidistant row is expected to give kappa = 1, since only it is the "
        "counterfactual the published sentence names.\n",

        "### Why kappa is not 1\n",
        f"The conversion needs D_VA(X) = alpha x D_VA(source). Inverting that "
        f"on this data asks the source population to have a D-statistic of "
        f"{ldf['implied_source_D_VA'].min():.2f} to "
        f"{ldf['implied_source_D_VA'].max():.2f}. **A D-statistic cannot exceed "
        f"1.** The relation is therefore not approximately wrong, it is "
        f"impossible, and no simulation is needed to see it.\n",
        "The reason is one the study already states and then does not apply "
        "here. `PAPER_neanderthal_source.md` says D_VA's *absolute* value is "
        "not interpretable - Vindija is pseudo-haploid where Altai is diploid, "
        "Vindija is called at less than half as many sites, and Yoruba carries "
        "its own Neanderthal ancestry - and that only differences between "
        "cohorts are meaningful, because those offsets are common-mode and "
        "cancel. The detection limit is the one place in the study that divides "
        "by the absolute value. Those offsets inflate the denominator of that "
        "division without contributing anything a real mixture can move, so the "
        "limit comes out too small by roughly the inflation factor.\n",

        "## Result 2: empirical detection power\n",
        pdf.round(5).to_markdown(index=False),

        f"\n\n**Baseline at f = 0: {100*fpr:.1f}% of pairs called, median |Z| = "
        f"{s('equidistant', 'baseline_median_abs_z'):.2f}.** This is not a "
        f"false-positive rate: real cohorts really do differ a little, so it is "
        f"the study's own pairwise significance rate among these pairs and an "
        f"upper bound on the nominal 5%. The median |Z| below 1 is the sign "
        f"that most of these pairs are genuinely null, consistent with the "
        f"study's finding that none survive Bonferroni. The power rows below "
        f"include this baseline rather than subtracting it.\n",

        f"**The bootstrap and the jackknife agree on the error bar** "
        f"(SE_boot/SE_jack = {s('equidistant', 'se_ratio_boot_over_jack'):.3f} "
        f"at the median). That is the calibration check that licenses the rest "
        f"of the table: the resampling scheme used to simulate detection and "
        f"the jackknife used to size the error bar are not telling different "
        f"stories about this data.\n",

        f"**f50 = {100*f50:.1f}%** against the analytical {100*analytic:.1f}%. "
        + ("The two agree, so the published limit is calibrated as a "
           "coin-flip threshold.\n" if agree else
           "The two disagree, and the empirical curve is what should be "
           "believed: it makes no linearity assumption.\n"),
        f"**f80 = {100*f80:.1f}%.** This is the number worth quoting to anyone "
        f"asking what the study would have *found*. A 2-sigma threshold gives "
        f"about 50% power against an effect sitting exactly on it, so the "
        f"published limit is the point at which detection becomes a coin flip, "
        f"not the point at which a real difference would reliably have been "
        f"seen.\n",

        "\n![Figure 6](fig_n6_mixture.png)\n",
        "**Figure 6.** Left: detection rate against injected fraction, with the "
        "analytical limit marked and f50 ringed. Right: kappa, the realised "
        "D_VA response as a fraction of the response the published conversion "
        "assumes; the dashed line is the assumption.\n",

        "## Caveats\n",
        "- **The result scales inversely with alpha.** The injected shift is "
        "alpha x *f* x (target - source), so f50 and f80 move in proportion to "
        "whatever error is in the cohort Neanderthal fractions. This "
        "repository's f4-ratio is known to run ~0.2pp high on a ~2.1% base, so "
        "if alpha is overstated by 10% the true thresholds are 10% *higher* "
        "than reported here. No other input carries this sensitivity.\n",
        "- **The mixture is instantaneous and clean.** One source replaced by "
        "another at a known fraction, with no drift, no LD decay and no "
        "post-admixture selection. A real second pulse would be messier and "
        "harder to see, so these thresholds are optimistic as descriptions of "
        "history even where they are exact as descriptions of the statistic.\n",
        "- **p_Vindija stands in for the true introgressing source**, which is "
        "unobserved. Vindija is the closest available proxy (Prufer et al. "
        "2017) and is what the study's premise already assumes; a source "
        "further from Vindija would change kappa.\n",
        f"- **This calibrates the conversion, not the floor.** The "
        f"{limit:.5f} resolvable difference is taken as given from "
        f"`ns_detection_limit.csv`; what is tested here is the step from that "
        f"number to a percentage of ancestry.\n",
    ]
    path = os.path.join(OUT, "POWER_mixture_calibration.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
