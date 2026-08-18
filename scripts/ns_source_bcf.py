#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The source contrast recomputed on full-quality archaic genomes.

The AADR ships Vindija as `VindijaG1_final.SG`, pseudo-haploid: one randomly
drawn allele per site, called at under half the panel. Every detection limit in
this study rests on the sites where Vindija and Altai are both called and
differ, and on the AADR's own copies that is 18,926 sites. The published BCFs of
the same two genomes are diploid and call 79.5% of the panel, giving 32,674 --
1.73x on the one axis that three separate experiments have found binding.

This recomputes the contrast with those frequencies substituted in.

WHY BOTH ARMS ARE RUN HERE

The obvious shortcut is to compute the new limit and compare it against the 37%
already in `POWER_mixture_calibration.md`. That would confound the change of
data with every difference between this script and that one -- cohort
selection, block assignment, which pairs enter the median. So the AADR
frequencies are run through this same code as a baseline, and only the two
numbers produced here are compared with each other. The published figure is
reported alongside as a cross-check on the baseline, not as the comparator.

WHAT IS SUBSTITUTED, AND WHAT IS NOT

Only the two contrast genomes change. Cohort allele frequencies are the AADR
ancient genomes exactly as before, reused from the study's own cache; Yoruba
stays as the baseline and Chimp as the outgroup for the symmetric normaliser.
The question, the cohorts and the statistic are untouched.

KAPPA IS RE-MEASURED, NOT CARRIED OVER

`POWER_mixture_calibration.md` found the published conversion from a resolvable
D_VA difference to a percentage of ancestry optimistic by a factor of three,
because it divides by an absolute D_VA that is inflated by offsets no real
mixture can move. Those offsets are properties of the archaic genomes, and this
run changes the archaic genomes -- so kappa is measured again on both arms
rather than assumed to transfer. A limit quoted before that is arithmetic.

Outputs (reports/neanderthal_source/):
  ns_bcf_comparison.csv   both arms: sites, limit, kappa, f50, f80
  ns_bcf_cohorts.csv      per-cohort D_VA under each arm
  fig_n9_bcf.png

Run: PYTHONIOENCODING=utf-8 python scripts/ns_source_bcf.py
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
import ns_mixture_power as nmp

log = get_logger("archaic.ns_source_bcf")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")
FREQS = os.path.join(ROOT, "results", "archaic_panel", "archaic_1240k_freqs.npz")
N_BOOT = 400
SEED = 20260818
MIX_F = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 1.00]


def arm(name, V, A, freqs, labels, alpha_of, block_full, nb, yor, chimp, rng,
        split_pairs, log):
    """One complete pass: limit, kappa and detection curve for one V/A pair."""
    useful = (np.isfinite(V) & np.isfinite(A) & np.isfinite(yor)
              & np.isfinite(chimp))
    R = dict(Vindija=V[useful], Altai=A[useful], Yoruba=yor[useful],
             Chimp=chimp[useful])
    block = block_full[useful]
    both = int(useful.sum())
    informative = int((R["Vindija"] != R["Altai"]).sum())
    log.info(f"[{name}] {both:,} sites with all four called; "
             f"{informative:,} distinguish Vindija from Altai "
             f"({100*informative/max(both,1):.2f}%)")

    def table_for(lab, p):
        return sc.build_block_table(lab, p, R["Yoruba"], R["Vindija"],
                                    R["Altai"], R["Chimp"], block, nb)

    tables = {lab: table_for(lab, freqs[lab][useful]) for lab in labels}
    # The control-split cohorts are needed too, or the systematic floor cannot
    # be computed and `detection_limit` silently falls back to the statistical
    # one -- which is the larger of the two here, so the error would not show
    # up in the answer.
    for lab in {x for pair in split_pairs for x in pair}:
        if lab not in tables and lab in freqs:
            tables[lab] = table_for(lab, freqs[lab][useful])
    d_va = {l: tables[l].d_va_theta for l in labels}
    signal = float(np.median([d_va[l] for l in labels]))

    # statistical floor: real cohort-vs-cohort comparisons
    ses = []
    for a, b in itertools.combinations(labels, 2):
        ses.append(sc.paired_difference(tables[a], tables[b], "D_VA")["se"])
    # systematic floor: the study's own same-cohort control splits
    nulls = []
    for a, b in split_pairs:
        if a in tables and b in tables:
            nulls.append(sc.paired_difference(tables[a], tables[b],
                                              "D_VA")["diff"])
    lim = sc.detection_limit(nulls, ses, signal)
    log.info(f"[{name}] typical D_VA {signal:.5f}; limit {lim['limit']:.5f} "
             f"(stat {lim['statistical_floor']:.5f}, sys "
             f"{lim['systematic_floor']:.5f}); median paired SE "
             f"{float(np.median(ses)):.5f}")

    # Injected tables are built once per (cohort, fraction) and reused by both
    # the kappa estimate and the detection curve. Building them inside the pair
    # loop instead costs 595 x 9 builds per arm rather than 35 x 9 -- a 17x
    # difference that turns a few minutes into hours.
    equi = sc.neanderthal_average(R["Vindija"], R["Altai"])
    inj, ks = {}, []
    for lab in labels:
        p = freqs[lab][useful]
        d0 = tables[lab].d_va_theta
        for f in MIX_F:
            if f == 0:
                inj[(lab, f)] = tables[lab]
                continue
            pm, _ = sc.mixture_frequencies(p, alpha_of[lab], f, R["Vindija"],
                                           equi)
            t = table_for(lab, pm)
            inj[(lab, f)] = t
            ks.append((t.d_va_theta - d0) / (-f * d0))
    kappa = float(np.median(ks))
    log.info(f"[{name}] kappa = {kappa:.4f} +/- {float(np.std(ks)):.4f}")

    pairs = list(itertools.combinations(labels, 2))
    boot_w = rng.multinomial(nb, np.full(nb, 1.0 / nb),
                             size=N_BOOT).astype(np.float64)
    curve = []
    for f in MIX_F:
        rates = [nmp.bootstrap_detection(inj[(a, f)], tables[b], nb,
                                         boot_w)["detect_rate"]
                 for a, b in pairs]
        curve.append(dict(fraction=f, detect_rate=float(np.mean(rates))))
    cdf = pd.DataFrame(curve)
    f50 = sc.power_crossing(cdf["fraction"], cdf["detect_rate"], 0.50)
    f80 = sc.power_crossing(cdf["fraction"], cdf["detect_rate"], 0.80)
    log.info(f"[{name}] detection: " +
             ", ".join(f"{100*r.fraction:.0f}%->{100*r.detect_rate:.0f}%"
                       for r in cdf.itertuples()))
    log.info(f"[{name}] f50 = {100*f50:.1f}%   f80 = {100*f80:.1f}%")
    return dict(arm=name, sites_called=both, informative=informative,
                typical_D_VA=signal, limit=lim["limit"],
                statistical_floor=lim["statistical_floor"],
                systematic_floor=lim["systematic_floor"],
                median_paired_se=float(np.median(ses)), kappa=kappa,
                f50=f50, f80=f80,
                analytic_f50=lim["limit"] / (kappa * abs(signal))), cdf, d_va


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--blocks", type=int, default=ns.N_BLOCKS)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    nb = args.blocks
    rng = np.random.default_rng(SEED)
    cfg = PANELS[args.panel]

    log.info("Loading panel and references...")
    panel = Panel(cfg["prefix"], autosomes_only=True)
    block_full = st.assign_blocks(panel.n_snp, nb)
    ref, _ = panel.frequencies({k: cfg["refs"][k] for k in
                                ("Altai", "Vindija", "Denisova", "Chimp",
                                 "Mbuti", "Yoruba")})

    # The importer indexes the *full* .snp table (1,233,013 rows); everything
    # downstream here -- the cached cohort frequencies, the reference arrays and
    # the block assignment -- is on the autosomal subset that `panel.snp_rows`
    # selects. Subset rather than assume, then verify against the stored
    # positions: a silent misalignment of two arrays of plausible length is the
    # kind of error that produces a confident wrong answer.
    z = np.load(FREQS, allow_pickle=False)
    rows = np.asarray(panel.snp_rows)
    stored_pos = z["pos"].astype(np.int64)
    if len(stored_pos) == panel.n_snp:
        take = slice(None)
    elif len(stored_pos) == len(panel.snp):
        take = rows
    else:
        raise SystemExit(f"frequency array length {len(stored_pos)} matches "
                         f"neither the full .snp table ({len(panel.snp)}) nor "
                         f"the autosomal subset ({panel.n_snp})")
    bcf = {k: z[k][take].astype(np.float64)
           for k in ("Altai", "Vindija", "Chagyrskaya", "Denisova")}
    want_pos = panel.snp["pos"].to_numpy(np.int64)[rows]
    if not np.array_equal(stored_pos[take], want_pos):
        raise SystemExit("stored positions do not match the panel's SNP order")
    log.info(f"  BCF frequencies aligned to {len(bcf['Altai']):,} autosomal "
             f"panel sites (positions verified)")

    meta = pd.read_csv(ns.META, low_memory=False)
    defs, crows, split_meta, _, _ = ns.build_cohorts(meta, panel)
    cdf0 = pd.read_csv(os.path.join(OUT, "ns_cohorts.csv"))
    keep = cdf0[(cdf0.kind == "grid") & cdf0.testable
                & cdf0.label.isin(defs)]
    labels = list(keep["label"])
    alpha_of = dict(zip(keep["label"], keep["mean_alpha"]))
    split_pairs = [(a, b) for _, _, a, b in split_meta]
    log.info(f"  {len(labels)} testable grid cohorts, "
             f"{len(split_pairs)} control splits")

    ca = types.SimpleNamespace(panel=args.panel, chunk=args.chunk)
    full, counts = ns.cached_pooled_freq(panel, defs, ca, "", log)
    # The control-split cohorts must be carried through as well, or the
    # systematic floor cannot be computed at all. It does not bind here -- the
    # statistical floor is several times larger -- but a limit reported as
    # max(stat, sys) with sys silently NaN is a limit that has not been checked.
    need = set(labels) | {x for pair in split_pairs for x in pair}
    freqs = {l: full[l] for l in need if l in full}
    missing = sorted(need - set(freqs))
    if missing:
        log.warning(f"  {len(missing)} cohorts absent from the cache: "
                    f"{missing[:3]}...")
    del full, counts

    rows, curves, dvas = [], {}, {}
    for name, V, A in (("AADR", ref["Vindija"], ref["Altai"]),
                       ("BCF", bcf["Vindija"], bcf["Altai"])):
        r, c, d = arm(name, V, A, freqs, labels, alpha_of, block_full, nb,
                      ref["Yoruba"], ref["Chimp"], rng, split_pairs, log)
        rows.append(r)
        curves[name] = c
        dvas[name] = d
    cmp = pd.DataFrame(rows)
    cmp.to_csv(os.path.join(OUT, "ns_bcf_comparison.csv"), index=False)
    pd.DataFrame({"label": labels,
                  "D_VA_AADR": [dvas["AADR"][l] for l in labels],
                  "D_VA_BCF": [dvas["BCF"][l] for l in labels]}
                 ).to_csv(os.path.join(OUT, "ns_bcf_cohorts.csv"), index=False)

    a, b = cmp.iloc[0], cmp.iloc[1]
    log.info("SUMMARY")
    log.info(f"  informative sites  {a['informative']:>9,} -> "
             f"{b['informative']:>9,}  ({b['informative']/a['informative']:.2f}x)")
    log.info(f"  median paired SE   {a['median_paired_se']:.5f} -> "
             f"{b['median_paired_se']:.5f}  "
             f"({b['median_paired_se']/a['median_paired_se']:.2f}x)")
    log.info(f"  kappa              {a['kappa']:.3f} -> {b['kappa']:.3f}")
    log.info(f"  f50                {100*a['f50']:.1f}% -> {100*b['f50']:.1f}%")
    log.info(f"  f80                {100*a['f80']:.1f}% -> {100*b['f80']:.1f}%")
    log.info(f"  (published AADR-based figures, different code: 37% / 64%)")
    make_figure(curves, cmp)
    log.info(f"Wrote ns_bcf_comparison.csv, ns_bcf_cohorts.csv and fig_n9 to {OUT}")


def make_figure(curves, cmp):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    col = {"AADR": "#777777", "BCF": "#b2182b"}
    for name, c in curves.items():
        r = cmp[cmp.arm == name].iloc[0]
        ax.plot(100 * c["fraction"], 100 * c["detect_rate"], "o-",
                color=col[name],
                label=f"{name}: {r['informative']:,} sites, f50={100*r['f50']:.0f}%")
    ax.axhline(50, color="grey", lw=0.8, ls=":")
    ax.axhline(80, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("fraction of Neanderthal ancestry re-sourced, f (%)")
    ax.set_ylabel("pairs called at 2 sigma (%)")
    ax.set_title("Source contrast on AADR vs full-quality archaic genomes")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_n9_bcf.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
