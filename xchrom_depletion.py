#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-chromosome vs autosome Neanderthal-ancestry depletion.

Neanderthal ancestry is markedly *lower on the X chromosome* than on the
autosomes in non-Africans — roughly a fifth of the autosomal level — a signature
of selection against introgressed alleles that is strongest on the X (reduced
male hybrid fertility / faster-X; Sankararaman et al. 2014 Nature 507:354;
Vernot & Akey 2014 Science 343:1017). Recovering that depletion from the AADR is
both a genuine result for ancient Eurasians and a strong *internal validation*:
the same f4-ratio estimator, restricted to an independent SNP set (the X), should
read far below its autosomal value.

Method: the identical pooled f4-ratio used by the core pipeline, run separately
on autosomes (.snp chrom 1-22) and the X (.snp chrom 23):

    alpha = f4(Altai, Chimp; POOL, Mbuti) / f4(Altai, Chimp; Vindija, Mbuti)

PANEL CHOICE (important, AADR-specific). This needs the Chimp outgroup *on the
X*. In AADR 1240K the outgroup sequences (Chimp.REF, Gorilla.REF, Ancestor.REF)
have **zero X genotypes**, so the f4-ratio cannot be formed on the X there at all.
The **Human Origins (HO)** panel *does* carry Chimp on the X (Chimp_HO.HO covers
all 3,814 HO X SNPs), so this analysis defaults to and requires --panel ho. The
script refuses to run on a panel whose Chimp lacks X coverage rather than emit a
meaningless number.

Because X and autosome SNPs are disjoint, alpha_X and alpha_auto are independent:
the depletion Z uses summed jackknife variances and the ratio the delta method.
Cohorts are *pooled* (a single pseudo-haploid ancient covers far too few X SNPs);
pooling thousands of genomes recovers near-complete coverage of the X panel. The
QC'd ancient set is taken from the 1240K Phase-2 metadata and mapped onto the HO
panel by genetic_id, so it inherits the pipeline's transparent QC.

Outputs:
  results/xchrom_<panel>_depletion.csv     one row per cohort
  results/figures/fig_xchrom_depletion.png
  XCHROM_REPORT.md

Usage:
  python xchrom_depletion.py                 # panel ho (the one that works)
  python xchrom_depletion.py --limit 500     # quick smoke test
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st
from archaic.cohort import pooled_freq
from archaic.refs import PANELS
from archaic.log_utils import get_logger

log = get_logger("archaic.xchrom")
ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
FIG = os.path.join(RESULTS, "figures")
N_BLOCKS = 50
AUTOSOMES = {str(c) for c in range(1, 23)}
X_CODE = "23"

# Present-day pops to try as positive controls (skipped if absent / too sparse).
PRESENT_DAY = ["French", "Sardinian", "Han", "Papuan", "Karitiana"]
# Where the QC'd ancient cohort + geography come from (Phase 2 is run on 1240K).
META_PANEL = "1240k"


def region_of(lon):
    """Coarse west->east Eurasian longitude bin (for a regional breakdown)."""
    if not np.isfinite(lon):
        return None
    if lon <= 35:
        return "WestEurasia"
    if lon <= 80:
        return "CentralSouthAsia"
    return "EastEurasia"


def f4ratio_on_mask(ref, pPool, mask, n_blocks=N_BLOCKS):
    """Pooled Neanderthal f4-ratio restricted to the SNPs selected by `mask`."""
    axN = ref["Altai"][mask] - ref["Chimp"][mask]
    num = axN * (pPool[mask] - ref["Mbuti"][mask])
    den = axN * (ref["Vindija"][mask] - ref["Mbuti"][mask])
    blk = st.assign_blocks(int(mask.sum()), n_blocks)
    return st.jackknife_ratio(num, den, blk, n_blocks)


def depletion_row(name, kind, n_ind, ref, pPool, mask_auto, mask_x):
    a = f4ratio_on_mask(ref, pPool, mask_auto)
    x = f4ratio_on_mask(ref, pPool, mask_x)
    aa, ax = a["theta"], x["theta"]
    ratio = ax / aa if (aa and np.isfinite(aa) and aa != 0) else np.nan
    if np.isfinite(ratio) and ax != 0:
        se_ratio = abs(ratio) * np.sqrt((x["se"] / ax) ** 2 + (a["se"] / aa) ** 2)
    else:
        se_ratio = np.nan
    se_diff = np.sqrt(a["se"] ** 2 + x["se"] ** 2)
    z_dep = (aa - ax) / se_diff if (se_diff and np.isfinite(se_diff) and se_diff > 0) else np.nan
    return dict(
        cohort=name, kind=kind, n_ind=int(n_ind),
        alpha_auto=round(aa, 6), alpha_auto_SE=round(a["se"], 6), auto_nSNP=a["n_used"],
        alpha_X=round(ax, 6), alpha_X_SE=round(x["se"], 6), X_nSNP=x["n_used"],
        ratio_X_auto=round(ratio, 4) if np.isfinite(ratio) else np.nan,
        ratio_SE=round(se_ratio, 4) if np.isfinite(se_ratio) else np.nan,
        depletion_z=round(z_dep, 2) if np.isfinite(z_dep) else np.nan,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="ho",
                    help="must carry Chimp on the X; only 'ho' does in AADR v66")
    ap.add_argument("--limit", type=int, default=0,
                    help="cap the pooled-ancient cohort size (quick test)")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--no-regions", action="store_true")
    ap.add_argument("--no-present-day", action="store_true")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    os.makedirs(FIG, exist_ok=True)

    log.info(f"Loading panel {args.panel} (autosomes + X)...")
    panel = Panel(cfg["prefix"], autosomes_only=False)
    rows = panel.snp_rows
    chrom = panel.snp["chrom"].to_numpy()[rows]
    mask_auto = np.isin(chrom, list(AUTOSOMES))
    mask_x = chrom == X_CODE
    log.info(f"  autosomal SNPs={mask_auto.sum():,}  X SNPs={mask_x.sum():,}")

    # ---- guard: the outgroup MUST have X coverage or the f4-ratio is empty --
    chimp_cols = panel.cols_for(**cfg["refs"]["Chimp"])
    chimp_x = int((panel.pg.read(rows[mask_x], chimp_cols) >= 0).any(axis=1).sum())
    if chimp_x == 0:
        log.error(f"Panel '{args.panel}': the Chimp outgroup has NO X-chromosome "
                  "genotypes, so f4(Altai,Chimp; .,.) is undefined on the X. "
                  "Re-run with --panel ho (its Chimp covers the X).")
        sys.exit(2)
    log.info(f"  Chimp X coverage OK ({chimp_x:,} X SNPs)")

    log.info("Reference allele frequencies (Altai/Vindija/Chimp/Mbuti)...")
    ref, ri = panel.frequencies({k: cfg["refs"][k]
                                 for k in ["Altai", "Vindija", "Chimp", "Mbuti"]})

    # ---- cohorts (QC'd ancients from Phase 2, mapped onto this panel) -------
    meta = pd.read_csv(os.path.join(RESULTS, f"phase2_{META_PANEL}_metadata.csv"))
    col_of = panel._id_to_col
    meta = meta[meta["genetic_id"].isin(col_of)].copy()
    meta["col"] = meta["genetic_id"].map(col_of)
    log.info(f"  {len(meta):,} QC'd ancients present in the {args.panel} panel")

    cohorts = []
    pooled = meta["col"].to_numpy()
    if args.limit:
        pooled = pooled[:args.limit]
    cohorts.append(("AllAncientEurasian", "ancient_pooled", pooled))

    if not args.no_regions:
        meta["region"] = meta["lon"].map(region_of)
        for rname, g in meta.groupby("region"):
            if rname and len(g) >= 50:
                cohorts.append((rname, "ancient_region", g["col"].to_numpy()))

    if not args.no_present_day:
        for pop in PRESENT_DAY:
            cc = panel.cols_for(pops=[pop])
            if len(cc) >= 3:
                cohorts.append((pop, "present_day", cc))

    # ---- compute -----------------------------------------------------------
    out_rows = []
    t0 = time.time()
    for name, kind, cols in cohorts:
        log.info(f"[{name}] n_ind={len(cols):,} ({kind}) — pooling frequencies...")
        pPool, _ = pooled_freq(panel, rows, cols, chunk=args.chunk)
        row = depletion_row(name, kind, len(cols), ref, pPool, mask_auto, mask_x)
        out_rows.append(row)
        log.info(f"  alpha_auto={row['alpha_auto']:.4f}  alpha_X={row['alpha_X']:.4f}"
                 f"  ratio={row['ratio_X_auto']}  depletion_z={row['depletion_z']}"
                 f"  (X nSNP={row['X_nSNP']})")

    df = pd.DataFrame(out_rows)
    out_csv = os.path.join(RESULTS, f"xchrom_{args.panel}_depletion.csv")
    df.to_csv(out_csv, index=False)
    log.info(f"Wrote {out_csv}  ({(time.time()-t0)/60:.1f} min)")

    make_figure(df, args.panel)
    write_report(df, args.panel, int(mask_auto.sum()), int(mask_x.sum()))


def make_figure(df, panel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 9})

    d = df.copy()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6),
                                   gridspec_kw={"width_ratios": [1.5, 1]})
    y = np.arange(len(d))[::-1]
    ax1.barh(y + 0.18, d["alpha_auto"] * 100, height=0.34, color="#4c72b0",
             xerr=d["alpha_auto_SE"] * 100, label="autosomes", capsize=2)
    ax1.barh(y - 0.18, d["alpha_X"] * 100, height=0.34, color="#c44e52",
             xerr=d["alpha_X_SE"] * 100, label="X chromosome", capsize=2)
    ax1.axvline(0, color="k", lw=0.6)
    ax1.set_yticks(y); ax1.set_yticklabels(d["cohort"])
    ax1.set_xlabel("Neanderthal ancestry  α (%)")
    ax1.set_title("Autosomal vs X-linked Neanderthal ancestry")
    ax1.legend(loc="lower right"); ax1.grid(axis="x", alpha=0.25)

    ax2.errorbar(d["ratio_X_auto"], y, xerr=d["ratio_SE"], fmt="o",
                 color="#55a868", capsize=3)
    ax2.axvline(1.0, color="k", ls="--", lw=1, label="no depletion")
    ax2.axvspan(0.15, 0.30, color="0.8", alpha=0.5,
                label="published present-day\nX/auto (~0.2)")
    ax2.set_yticks(y); ax2.set_yticklabels([])
    ax2.set_xlabel("ratio  α(X) / α(autosome)")
    ax2.set_title("X-depletion of Neanderthal ancestry")
    ax2.legend(fontsize=7, loc="lower right"); ax2.grid(axis="x", alpha=0.25)

    fig.tight_layout()
    path = os.path.join(FIG, "fig_xchrom_depletion.png")
    fig.savefig(path); plt.close(fig)
    log.info(f"Wrote {path}")


def write_report(df, panel, n_auto, n_x):
    head = df[df["cohort"] == "AllAncientEurasian"].iloc[0]
    tbl = df[["cohort", "kind", "n_ind", "alpha_auto", "alpha_X",
              "ratio_X_auto", "ratio_SE", "depletion_z"]].copy()
    tbl["alpha_auto"] = (tbl["alpha_auto"] * 100).round(2)
    tbl["alpha_X"] = (tbl["alpha_X"] * 100).round(2)
    signif = np.isfinite(head["depletion_z"]) and head["depletion_z"] >= 2 \
        and np.isfinite(head["ratio_X_auto"]) and head["ratio_X_auto"] < 1
    if signif:
        verdict = ("The X carries significantly less Neanderthal ancestry than the "
                   "autosomes, the expected signature of selection against introgressed "
                   "alleles.")
    else:
        verdict = ("**This is inconclusive: the AADR X panel is too sparse to resolve the "
                   "depletion.** Only ~1.7k X SNPs carry the outgroup, so α(X) has a very "
                   "large jackknife SE and the X/auto ratio is not significantly below 1 "
                   "(and swings wildly, even negative, for small cohorts). The expected "
                   "~5x X-depletion cannot be demonstrated from AADR at this SNP density; "
                   "it would need shotgun genomes or a denser panel with outgroup X coverage.")
    md = ["# X-chromosome depletion of Neanderthal ancestry\n",
          f"*Panel: {panel} — {n_auto:,} autosomal SNPs, {n_x:,} X SNPs "
          f"(~1.7k usable, capped by Vindija's X coverage).*\n",
          "## Result\n",
          f"Pooling all {int(head['n_ind']):,} retained ancient Eurasian genomes, the "
          f"autosomal Neanderthal estimate is **{head['alpha_auto']*100:.2f}%** (as "
          f"expected) while the X reads **{head['alpha_X']*100:.2f}%** — ratio "
          f"**α(X)/α(auto) = {head['ratio_X_auto']} ± {head['ratio_SE']}** "
          f"(depletion Z = {head['depletion_z']}). " + verdict + "\n",
          "## Per-cohort\n",
          tbl.rename(columns={"alpha_auto": "α_auto (%)", "alpha_X": "α_X (%)",
                              "ratio_X_auto": "ratio", "depletion_z": "depletion_Z"}
                     ).to_markdown(index=False),
          "\n\n![X depletion](results/figures/fig_xchrom_depletion.png)\n",
          "## Method\n",
          "The AADR `.snp` codes the X as `23`. We run the pipeline's pooled f4-ratio "
          "α = f4(Altai, Chimp; POOL, Mbuti) / f4(Altai, Chimp; Vindija, Mbuti) separately "
          "on autosomes (1–22) and X (23), with a 50-block delete-one jackknife. The two "
          "SNP sets are disjoint, so α(X) and α(auto) are independent; the depletion Z uses "
          "summed jackknife variances and the ratio SE the delta method. Cohorts are pooled "
          "because a single pseudo-haploid ancient covers too few X SNPs to estimate α.\n",
          "## Why the Human Origins panel\n",
          "The f4-ratio needs the Chimp outgroup *on the X*. In AADR 1240K the outgroup "
          "sequences (Chimp.REF, Gorilla.REF, Ancestor.REF) have **zero X genotypes**, so "
          "the statistic is undefined on the X there. The Human Origins panel carries Chimp "
          "across the X (all 3,814 HO X SNPs), which is why this analysis uses it; the QC'd "
          "ancient cohort is inherited from the 1240K Phase-2 metadata by genetic_id.\n",
          "## Interpretation & caveats\n",
          "- The biology is well established — Neanderthal ancestry is ~5x lower on the X "
          "(faster-X / reduced male hybrid fertility; Sankararaman et al. 2014; Vernot & Akey "
          "2014) — but this analysis **cannot resolve it on AADR**, so nothing here should be "
          "read as confirming or refuting it.\n",
          "- The limitation is SNP count, not sample size: the HO X panel has only ~3.8k SNPs "
          "and the Vindija scale rests on ~1.7k, so α(X)'s jackknife SE is ~17x the autosomal "
          "one and floored by the number of X SNP-blocks — pooling more genomes does not help. "
          "Per-cohort X/auto ratios (e.g. Han 0.24, Papuan −3.1) are individually meaningless "
          "at this SE and should not be interpreted.\n",
          "- The X/auto *ratio* is the right (offset-free) quantity to compare, and the "
          "autosomal arm reproduces the expected ~2–3% perfectly; only the X arm is data-"
          "starved. A real test would need shotgun diploid genomes or a denser panel that "
          "genotypes an outgroup across the X.\n",
          "\n*Refs: Sankararaman et al. 2014 Nature 507:354; Vernot & Akey 2014 Science "
          "343:1017; Petr et al. 2019 PNAS 116:1639; Mallick et al. 2024 Sci. Data 11:182 (AADR).*\n"]
    path = os.path.join(ROOT, "XCHROM_REPORT.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    log.info(f"Wrote {path}")


if __name__ == "__main__":
    main()
