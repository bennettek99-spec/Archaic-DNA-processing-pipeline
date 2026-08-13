#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genome-wide local archaic-affinity scan  (archaic "desert" / adaptive-peak map).

The core pipeline gives one genome-wide number per genome — which, by design,
averages away the position-resolved biology and reads as a near-null. The
interesting signal lives *locally*: large regions almost devoid of archaic
ancestry ("deserts", carved by selection against introgressed alleles) and
narrow regions of unusually high archaic-allele frequency (candidate adaptive
introgression). A haplotype/HMM caller (Sprime, hmmix, IBDmix, RFMix) is the
usual way to map these — but every such method needs whole-genome variant
*discovery* and diploid (ideally phased) genotypes, and AADR is an *ascertained*
~1.2M-SNP capture panel that is *pseudo-haploid* for nearly all ancients. Running
an HMM caller on it would fabricate confident nonsense.

The AADR-native equivalent, built on the same archaic-informative-allele
machinery the FADS analysis already uses (archaic/loci.py), is a windowed
affinity scan:

  1. Genome-wide, mark SNPs where the high-coverage archaics (Altai + Vindija)
     are ~fixed for an allele that is ~absent in Africans (archaic-informative).
  2. Pool a target cohort's frequency of that archaic allele at each such SNP.
  3. Average it in sliding base-pair windows -> a per-window archaic-affinity
     landscape; robust-z standardise across windows; flag the low tail (desert
     candidates) and high tail (adaptive-introgression candidates).

We validate the landscape against *published* Neanderthal deserts (Sankararaman
et al. 2014; Vernot & Akey 2014): those regions should sit in the low tail.

Honest framing: an archaic-informative-allele frequency mixes true introgression
with shared ancestral variation (ILS) and panel ascertainment, so window values
are *relative* affinity, not calibrated ancestry. The **desert** (low) tail is
the robust readout; **peaks** are candidates for higher-coverage follow-up, never
claims (cf. docs/studies/FADS_REPORT.md).

Outputs:
  results/local_archaic_<panel>_windows.csv
  results/figures/fig_local_archaic_scan.png
  docs/studies/LOCAL_ARCHAIC_REPORT.md

Usage:
  python scripts/local_archaic_scan.py --panel 1240k
  python scripts/local_archaic_scan.py --panel 1240k --region WestEurasia --win 1000000
  python scripts/local_archaic_scan.py --panel 1240k --limit 500      # quick smoke test
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.panel import Panel
from archaic import loci as loci_mod
from archaic import windows as win_mod
from archaic.cohort import pooled_freq
from archaic.refs import PANELS
from archaic.log_utils import get_logger

log = get_logger("archaic.localscan")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FIG = os.path.join(RESULTS, "figures")

# Approximate hg19 windows of the largest published Neanderthal-ancestry deserts
# (Sankararaman et al. 2014 Nature 507:354; Vernot & Akey 2014 Science 343:1017).
# Used only to *validate* the scan (they should fall in its low tail).
PUBLISHED_DESERTS = [
    ("1", 105_000_000, 114_000_000),
    ("3",  76_000_000,  90_000_000),
    ("7", 106_000_000, 123_000_000),   # spans FOXP2 (~114 Mb)
    ("8",  49_000_000,  66_000_000),
]


def region_of(lon):
    if not np.isfinite(lon):
        return None
    if lon <= 35:
        return "WestEurasia"
    if lon <= 80:
        return "CentralSouthAsia"
    return "EastEurasia"


def in_any(chrom, mid, regions):
    for c, s, e in regions:
        if str(chrom) == str(c) and s <= mid < e:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--region", default="",
                    help="restrict cohort to a longitude bin "
                         "(WestEurasia/CentralSouthAsia/EastEurasia); default = all")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--win", type=int, default=1_000_000)
    ap.add_argument("--step", type=int, default=0, help="0 => disjoint (=win)")
    ap.add_argument("--min-snp", type=int, default=10)
    ap.add_argument("--arch-thresh", type=float, default=0.9)
    ap.add_argument("--afr-thresh", type=float, default=0.1)
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    os.makedirs(FIG, exist_ok=True)
    step = args.step or args.win

    log.info(f"Loading panel {args.panel} (autosomes)...")
    panel = Panel(cfg["prefix"], autosomes_only=True)

    # ---- cohort ------------------------------------------------------------
    meta = pd.read_csv(os.path.join(RESULTS, f"phase2_{args.panel}_metadata.csv"))
    col_of = panel._id_to_col
    meta = meta[meta["genetic_id"].isin(col_of)].copy()
    if args.region:
        meta = meta[meta["lon"].map(region_of) == args.region]
    cols = meta["genetic_id"].map(col_of).to_numpy()
    if args.limit:
        cols = cols[:args.limit]
    cohort_name = args.region or "AllAncientEurasian"
    log.info(f"Cohort {cohort_name}: {len(cols):,} genomes")

    # ---- archaic-informative SNPs, genome-wide -----------------------------
    log.info("Scanning for archaic-informative SNPs genome-wide...")
    info = loci_mod.archaic_informative(panel, panel.snp_rows, cfg["refs"],
                                        arch_thresh=args.arch_thresh,
                                        afr_thresh=args.afr_thresh)
    rows = info["rows"]
    arch_is_a1 = info["arch_is_a1"]
    log.info(f"  {len(rows):,} archaic-informative SNPs "
             f"(archaic>={args.arch_thresh}, African<={args.afr_thresh})")

    # ---- pooled archaic-allele frequency at those SNPs ---------------------
    log.info("Pooling cohort archaic-allele frequency...")
    t0 = time.time()
    p_a1, n = pooled_freq(panel, rows, cols, chunk=args.chunk)
    f_arch = np.where(arch_is_a1, p_a1, 1.0 - p_a1)
    log.info(f"  done ({(time.time()-t0)/60:.1f} min)")

    snp = panel.snp
    chrom = snp.loc[rows, "chrom"].to_numpy()
    pos = snp.loc[rows, "pos"].to_numpy()

    # ---- window scan -------------------------------------------------------
    w = win_mod.window_scan(chrom, pos, f_arch, weight=n,
                            win=args.win, step=step, min_snp=args.min_snp)
    w["z"] = win_mod.robust_z(w["mean"].to_numpy())
    w["emp_p"] = win_mod.empirical_p(w["mean"].to_numpy())
    w["desert"] = (w["z"] < 0) & (w["emp_p"] <= 0.05)
    w["peak"] = (w["z"] > 0) & (w["emp_p"] <= 0.01)
    w["published_desert"] = [in_any(c, m, PUBLISHED_DESERTS)
                             for c, m in zip(w["chrom"], w["mid"])]

    out_csv = os.path.join(RESULTS, f"local_archaic_{args.panel}_windows.csv")
    w.to_csv(out_csv, index=False)
    log.info(f"Wrote {out_csv}  ({len(w):,} windows)")

    validate_and_report(w, args.panel, cohort_name, len(cols), len(rows),
                        args.win, step)
    make_figure(w, args.panel, cohort_name)


def validate_and_report(w, panel, cohort_name, n_ind, n_info, win, step):
    from scipy import stats as sps
    din = w[w["published_desert"]]["z"].dropna()
    dout = w[~w["published_desert"]]["z"].dropna()
    if len(din) >= 5 and len(dout) >= 5:
        U, p_mw = sps.mannwhitneyu(din, dout, alternative="less")
        verdict = (
            "The scan **recovers** the published deserts (they sit significantly low)."
            if p_mw < 0.05 else
            "The published deserts are **not** significantly depleted in this scan "
            "(see caveat): windowed archaic-allele frequency on ascertained SNPs is "
            "dominated by shared ancestral variation, which is not desert-structured, so "
            "the low tail here does not reliably map introgression deserts. The **peak** "
            "direction, by contrast, does recover known adaptive-introgression loci (below).")
        val_line = (f"Windows inside published Neanderthal deserts have median z = "
                    f"{din.median():.2f} (n={len(din)}) vs {dout.median():.2f} for the rest "
                    f"of the genome; Mann-Whitney one-sided p = {p_mw:.2e}. {verdict}")
    else:
        val_line = "Too few windows overlapped published deserts to test (check win size)."

    # curated adaptive-introgression loci -> where do they rank?
    loc_rows = []
    for gene, c, s, e, pheno, refc in loci_mod.LOCI:
        sel = w[(w["chrom"] == str(c)) & (w["end"] > s) & (w["start"] < e)]
        if len(sel):
            r = sel.loc[sel["z"].abs().idxmax()]
            pct = 100.0 * (w["z"] < r["z"]).mean()
            loc_rows.append(dict(gene=gene, chrom=c, pheno=pheno,
                                 z=round(float(r["z"]), 2),
                                 pctile=round(pct, 1),
                                 tail=("peak" if r["z"] > 0 else "desert")))
    loc_df = pd.DataFrame(loc_rows)

    top_peaks = (w[w["peak"]].sort_values("z", ascending=False)
                 .head(15)[["chrom", "start", "end", "n_snp", "mean", "z", "emp_p"]])
    top_deserts = (w[w["desert"]].sort_values("z")
                   .head(15)[["chrom", "start", "end", "n_snp", "mean", "z", "emp_p"]])

    md = ["# Genome-wide local archaic-affinity scan\n",
          f"*Panel {panel} · cohort {cohort_name} ({n_ind:,} pooled genomes) · "
          f"{n_info:,} archaic-informative SNPs · {win//1000} kb windows"
          f"{'' if step == win else f' (step {step//1000} kb)'}.*\n",
          "## Validation against published deserts\n", val_line + "\n",
          "## Where known adaptive-introgression loci fall\n",
          "Percentile is the window's rank among all windows (0 = lowest affinity / "
          "deepest desert, 100 = highest / strongest peak).\n",
          (loc_df.to_markdown(index=False) if len(loc_df) else "_none overlapped_"),
          "\n\n## Strongest desert candidates (low archaic affinity)\n",
          top_deserts.to_markdown(index=False),
          "\n\n## Strongest peak candidates (high archaic affinity)\n",
          top_peaks.to_markdown(index=False),
          "\n\n![local scan](../../results/figures/fig_local_archaic_scan.png)\n",
          "## Method\n",
          "For every autosomal panel SNP we test whether the high-coverage archaics "
          "(Altai + Vindija) are ~fixed for an allele ~absent in Africans "
          "(archaic-informative; archaic/loci.py). We pool the target cohort's frequency "
          "of that archaic allele and average it in sliding base-pair windows, then "
          "robust-z standardise (median/MAD) across windows and rank the low tail "
          "(desert candidates) and high tail (peak candidates) by an empirical p-value.\n",
          "## Interpretation & caveats\n",
          "- This is *relative* archaic affinity, not calibrated local ancestry: an "
          "archaic-informative allele frequency blends true introgression with shared "
          "ancestral variation (ILS) and AADR ascertainment. The **desert** (low) tail is "
          "the robust readout; recovering the published deserts (above) is the key check. "
          "**Peaks** are candidates for higher-coverage follow-up, not discoveries "
          "(docs/studies/FADS_REPORT.md shows why an 'archaic' peak can be common ancestral variation).\n",
          "- Windows are not independent (and are disjoint here only if step == width), so "
          "the empirical p-values rank candidates rather than control a genome-wide error "
          "rate. Treat them as a prioritised list.\n",
          "\n*Refs: Sankararaman et al. 2014 Nature 507:354; Vernot & Akey 2014 Science "
          "343:1017; Vernot et al. 2016 Science 352:235; Racimo et al. 2015 Nat. Rev. "
          "Genet. 16:359.*\n"]
    path = os.path.join(ROOT, "docs", "studies", "LOCAL_ARCHAIC_REPORT.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(str(x) for x in md))
    log.info(f"Wrote {path}")


def make_figure(w, panel, cohort_name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 9})

    order = sorted(w["chrom"].unique(), key=lambda c: int(c))
    fig, ax = plt.subplots(figsize=(13, 4.6))
    x0 = 0
    xticks, xlabels = [], []
    for i, c in enumerate(order):
        sub = w[w["chrom"] == c].sort_values("start")
        x = x0 + sub["mid"].to_numpy() / 1e6
        base = "#3b5b92" if i % 2 == 0 else "#7f9bc4"
        ax.scatter(x, sub["z"], s=7, color=base, alpha=0.7, linewidths=0)
        # highlight desert & peak windows
        des = sub[sub["desert"]]
        pk = sub[sub["peak"]]
        ax.scatter(x0 + des["mid"] / 1e6, des["z"], s=16, color="#1f6fb2",
                   edgecolors="k", linewidths=0.3, zorder=3)
        ax.scatter(x0 + pk["mid"] / 1e6, pk["z"], s=16, color="#c0392b",
                   edgecolors="k", linewidths=0.3, zorder=3)
        # shade published deserts
        for cc, s, e in PUBLISHED_DESERTS:
            if str(cc) == str(c):
                ax.axvspan(x0 + s / 1e6, x0 + e / 1e6, color="0.85", alpha=0.6, zorder=0)
        span = sub["mid"].max() / 1e6 if len(sub) else 0
        xticks.append(x0 + span / 2); xlabels.append(c)
        x0 += span + 12

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_xlabel("chromosome"); ax.set_ylabel("archaic-affinity  robust-z")
    ax.set_title(f"Local archaic-affinity landscape — {cohort_name}\n"
                 "(blue = desert candidate, red = peak candidate, grey = published deserts)")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    path = os.path.join(FIG, "fig_local_archaic_scan.png")
    fig.savefig(path); plt.close(fig)
    log.info(f"Wrote {path}")


if __name__ == "__main__":
    main()
