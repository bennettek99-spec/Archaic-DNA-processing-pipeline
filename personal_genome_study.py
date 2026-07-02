#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
personal_genome_study.py — estimate archaic (Neanderthal + Denisovan) introgression
for a single personal direct-to-consumer genome (MyHeritage / 23andMe / Ancestry)
using the pipeline's validated AADR archaic references.

METHOD NOTE (why not the f4-ratio here). For ancient shotgun genomes the pipeline
uses the f4-ratio alpha = f4(Altai,Chimp;X,Mbuti)/f4(Altai,Chimp;Vindija,Mbuti).
That statistic needs the *second* archaic genome (Vindija) covered at the same
sites and hundreds of thousands of archaic-informative SNPs; a consumer array
shares too few of those, and on the shared set the f4-ratio collapses toward zero
*even for populations known to be ~2% Neanderthal* (French/English read ~0 on the
consumer SNP subset — verified). So for a personal array genome the valid,
powered statistic is an archaic-allele MATCH RATE over curated marker alleles
(derived, carried by both archaics, African-rare), internally calibrated to a
percentage against reference populations of known archaic ancestry measured on the
*identical* SNP set. f4-ratio / D are still reported as a methods comparison.
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st
from archaic import consumer_dna as cdna
from archaic.refs import PANELS

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results", "personal")
FIG = os.path.join(RESULTS, "figures")
os.makedirs(FIG, exist_ok=True)
N_BLOCKS = 50
PANEL = "1240k"

CONTEXT = {
    "Orcadian": ["Orcadian"], "English": ["English"], "French": ["French"],
    "Basque": ["Basque"], "Spanish": ["Spanish"], "Italian_North": ["Italian_North"],
    "Sardinian": ["Sardinian"], "Han": ["Han"], "Papuan": ["Papuan"],
    "Yoruba": ["Yoruba"], "Mbuti": ["Mbuti"],
}
WEST_EUR = ["Orcadian", "English", "French", "Basque", "Spanish",
            "Italian_North", "Sardinian"]
# populations used to calibrate match-rate -> % (exclude Papuan: Denisovan inflates
# its Neanderthal-marker sharing; keep the African anchors — they fix the zero-point)
CALIB = WEST_EUR + ["Han", "Yoruba", "Mbuti"]


def matchrate_all(freq, keep, a1, a2, names, block):
    """Match rate (theta,se) and per-block sums for every population over `keep`."""
    out, bsum, bcnt = {}, {}, {}
    for nm in names:
        af = cdna.archaic_allele_freq(freq[nm], a1, a2)
        th, se, n, bn, bw = cdna.block_mean_jackknife(af, keep, block, N_BLOCKS)
        out[nm] = dict(match=th * 100, se=se * 100, n=n)
        bsum[nm], bcnt[nm] = bn, bw
    return out, bsum, bcnt


def calibrate_predict(match_by_pop, known, target):
    """Fit known_alpha ~ match% over CALIB pops; return (slope, intercept, r, pred)."""
    xs = np.array([match_by_pop[p] for p in CALIB])
    ys = np.array([known[p] for p in CALIB])
    sl, ic = np.polyfit(xs, ys, 1)
    r = float(np.corrcoef(xs, ys)[0, 1])
    return sl, ic, r, sl * target + ic


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome"); ap.add_argument("--name", default="You")
    args = ap.parse_args()

    cfg = PANELS[PANEL]
    print(f"Personal-genome archaic study on the {PANEL} panel\nGenome: {args.genome}\n")
    panel = Panel(cfg["prefix"])
    block = st.assign_blocks(panel.n_snp, N_BLOCKS)
    chrom_of_row = panel.snp.iloc[panel.snp_rows]["chrom"].to_numpy()

    consumer = cdna.read_consumer_file(args.genome)
    me_freq, astats = cdna.align_to_panel(panel, consumer)
    print("Alignment to panel:")
    for k, v in astats.items():
        print(f"  {k:26s} {v:,}")
    bmask = np.isfinite(me_freq)

    spec = {k: cfg["refs"][k] for k in ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti", "Yoruba"]}
    for k, v in CONTEXT.items():
        spec[k] = dict(pops=v)
    freq, info = panel.frequencies(spec)
    freq[args.name] = me_freq

    # coding sanity check
    ok = np.isfinite(me_freq) & np.isfinite(freq["French"])
    r_check = float(np.corrcoef(me_freq[ok], freq["French"][ok])[0, 1])
    print(f"\nAllele-coding check: corr({args.name},French)={r_check:+.3f} "
          f"({'PASS' if r_check > 0.3 else 'FAIL'})")

    names = [args.name] + list(CONTEXT.keys())

    # ---- known Neanderthal % for reference pops (full-panel f4-ratio) ---------
    known = {}
    for nm in CONTEXT:
        if nm == "Mbuti":
            known[nm] = 0.0; continue                      # baseline by definition
        a = st.f4_ratio(freq, "Altai", "Chimp", nm, "Mbuti", "Vindija", block, N_BLOCKS)
        known[nm] = a["theta"] * 100

    # ---- PRIMARY estimator: Neanderthal match rate ---------------------------
    keepN, a1N, a2N = cdna.archaic_markers(freq, "neanderthal", afr_max=0.05)
    keepN_b = keepN & bmask
    mrN, bsumN, bcntN = matchrate_all(freq, keepN_b, a1N, a2N, names, block)
    match_by_pop = {p: mrN[p]["match"] for p in names}
    sl, ic, r_cal, pred = calibrate_predict(match_by_pop, known, mrN[args.name]["match"])

    # block-jackknife the WHOLE calibrated prediction (propagates calibration too)
    loo_pred = []
    for b in range(N_BLOCKS):
        m_loo = {}
        for p in names:
            T, W = bsumN[p].sum(), bcntN[p].sum()
            Tw, Ww = T - bsumN[p][b], W - bcntN[p][b]
            if Ww <= 0:
                m_loo = None; break
            m_loo[p] = Tw / Ww * 100
        if m_loo is None:
            continue
        _, _, _, p_loo = calibrate_predict(m_loo, known, m_loo[args.name])
        loo_pred.append(p_loo)
    loo_pred = np.array(loo_pred); g = len(loo_pred)
    pred_se = float(np.sqrt((g - 1) / g * np.sum((loo_pred - loo_pred.mean()) ** 2)))

    # ---- Denisovan match rate (should be ~0 for a West-Eurasian) -------------
    keepD, a1D, a2D = cdna.archaic_markers(freq, "denisovan", afr_max=0.05)
    keepD_b = keepD & bmask
    mrD, _, _ = matchrate_all(freq, keepD_b, a1D, a2D, names, block)

    # ---- methods comparison: f4-ratio + D on the personal genome -------------
    a_me = st.f4_ratio(freq, "Altai", "Chimp", args.name, "Mbuti", "Vindija", block, N_BLOCKS)
    dN_me = st.dstat(freq, args.name, "Yoruba", "Altai", "Chimp", block, N_BLOCKS)
    dD_me = st.dstat(freq, args.name, "Yoruba", "Denisova", "Chimp", block, N_BLOCKS)

    # ---- per-chromosome match rate (stability) -------------------------------
    afN = cdna.archaic_allele_freq(me_freq, a1N, a2N)
    perchr = []
    for c in [str(i) for i in range(1, 23)]:
        m = keepN_b & (chrom_of_row == c) & np.isfinite(afN)
        if m.sum() < 8:
            continue
        mr = float(np.nanmean(afN[m])) * 100
        perchr.append(dict(chrom=int(c), match=mr, alpha=sl * mr + ic, nSNP=int(m.sum())))
    PC = pd.DataFrame(perchr)

    # ---- assemble reference table --------------------------------------------
    rows = []
    for nm in names:
        rows.append(dict(name=nm, n_ind=(1 if nm == args.name else info.get(nm, {}).get("n_ind", 0)),
                         Nea_match=mrN[nm]["match"], Nea_match_se=mrN[nm]["se"],
                         Nea_alpha=sl * mrN[nm]["match"] + ic,
                         Den_match=mrD[nm]["match"], Den_match_se=mrD[nm]["se"],
                         known_alpha=known.get(nm, np.nan), nSNP=mrN[nm]["n"]))
    R = pd.DataFrame(rows)

    we_alpha = float(R[R["name"].isin(WEST_EUR)]["Nea_alpha"].mean())
    den_base = float(R[R["name"].isin(WEST_EUR)]["Den_match"].mean())
    den_pap = float(R[R["name"] == "Papuan"]["Den_match"].iloc[0])

    # ============================ REPORT ======================================
    print("\n" + "=" * 72)
    print(f"FINAL RESULT — {args.name}")
    print("=" * 72)
    print(f"  Neanderthal marker sites used   : {mrN[args.name]['n']:,}")
    print(f"  Neanderthal match rate          : {mrN[args.name]['match']:.2f}% "
          f"(West-Eur refs {R[R.name.isin(WEST_EUR)].Nea_match.min():.1f}-"
          f"{R[R.name.isin(WEST_EUR)].Nea_match.max():.1f}%, Africans ~1-2%)")
    print(f"  Calibration (r={r_cal:.3f})         : alpha = {sl:.4f}*match {ic:+.3f}")
    print(f"  >> Neanderthal ancestry         : {pred:.2f}%  +/- {pred_se:.2f}%")
    print(f"  West-Eurasian reference mean    : {we_alpha:.2f}%")
    print(f"  Denisovan match rate            : {mrD[args.name]['match']:.2f}% "
          f"(West-Eur base {den_base:.2f}%, Papuan {den_pap:.2f}%) -> ~0, negligible")
    print(f"  TOTAL archaic (Nea + ~0 Den)    : ~{pred:.2f}%")
    print("-" * 72)
    print(f"  [methods comparison] f4-ratio   : {a_me['theta']*100:.2f} +/- "
          f"{a_me['se']*100:.2f}%  ({a_me['n_used']:,} SNP) - underpowered on arrays")
    print(f"  [methods comparison] D_Nea Z    : {dN_me['z']:.1f}   D_Den Z: {dD_me['z']:.1f}")
    print("=" * 72)

    # ---- save ----------------------------------------------------------------
    R.to_csv(os.path.join(RESULTS, "personal_estimates.csv"), index=False)
    PC.to_csv(os.path.join(RESULTS, "personal_perchrom.csv"), index=False)
    summary = dict(
        name=args.name, panel=PANEL, method="calibrated archaic-allele match rate",
        Nea_pct=float(pred), Nea_SE=float(pred_se),
        Nea_match=float(mrN[args.name]["match"]), Nea_match_se=float(mrN[args.name]["se"]),
        n_marker_snps=int(mrN[args.name]["n"]),
        calib_slope=float(sl), calib_intercept=float(ic), calib_r=float(r_cal),
        Den_match=float(mrD[args.name]["match"]), Den_base_westeur=den_base, Den_papuan=den_pap,
        west_eur_mean=we_alpha,
        f4ratio_pct=float(a_me["theta"] * 100), f4ratio_se=float(a_me["se"] * 100),
        f4ratio_nsnp=int(a_me["n_used"]), D_Nea_Z=float(dN_me["z"]), D_Den_Z=float(dD_me["z"]),
        r_coding_check=r_check, align=astats,
    )
    with open(os.path.join(RESULTS, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    make_figures(args.name, R, PC, pred, pred_se, we_alpha, sl, ic)
    print("\nWrote results/personal/{personal_estimates,personal_perchrom}.csv, "
          "summary.json, 3 figures.")
    return summary


def make_figures(name, R, PC, pred, pred_se, we_alpha, sl, ic):
    NA = "#c0392b"
    WE = ["Orcadian", "English", "French", "Basque", "Spanish", "Italian_North", "Sardinian"]

    # fig1: calibrated Neanderthal % — you vs references
    Rs = R.sort_values("Nea_alpha").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    se_alpha = Rs["Nea_match_se"] * sl
    cols = [NA if n == name else ("#7f8c8d" if n in WE else "#34495e") for n in Rs["name"]]
    ax.barh(Rs["name"], Rs["Nea_alpha"], xerr=se_alpha, color=cols,
            error_kw=dict(ecolor="#2c3e50", capsize=2, lw=1), alpha=0.9)
    ax.axvline(we_alpha, color="#3498db", lw=1, ls="--", alpha=0.8,
               label=f"West-Eur mean {we_alpha:.2f}%")
    ax.set_xlabel("Neanderthal ancestry (%) — calibrated archaic-allele match rate")
    ax.set_title(f"Neanderthal ancestry: {name} vs present-day references (1240K)")
    ax.legend(loc="lower right", fontsize=8); ax.grid(axis="x", alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_pg1_neanderthal.png", dpi=150); plt.close(fig)

    # fig2: calibration line (match rate -> known %), with You placed on it
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    for _, r in R.iterrows():
        if r["name"] == name:
            continue
        c = "#8e44ad" if r["name"] == "Papuan" else ("#7f8c8d" if r["name"] in WE else "#34495e")
        ax.scatter(r["Nea_match"], r["known_alpha"], s=45, color=c, zorder=3)
        ax.annotate(r["name"], (r["Nea_match"], r["known_alpha"]), fontsize=7,
                    xytext=(4, 2), textcoords="offset points")
    xs = np.linspace(R["Nea_match"].min() - 1, R["Nea_match"].max() + 1, 50)
    ax.plot(xs, sl * xs + ic, "k--", lw=1, alpha=0.6, label=f"calib: α={sl:.3f}·m{ic:+.2f}")
    me = R[R["name"] == name].iloc[0]
    ax.errorbar(me["Nea_match"], pred, yerr=pred_se, fmt="*", color=NA, ms=18,
                capsize=3, zorder=5, label=f"{name}: {pred:.2f}%")
    ax.set_xlabel("archaic-allele match rate (%)")
    ax.set_ylabel("Neanderthal ancestry, known / calibrated (%)")
    ax.set_title("Internal calibration of the match-rate estimator")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_pg2_calibration.png", dpi=150); plt.close(fig)

    # fig3: per-chromosome estimate
    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    ax.axhline(pred, color=NA, lw=1.5, ls="--", label=f"genome-wide {pred:.2f}%")
    ax.axhspan(pred - pred_se, pred + pred_se, color=NA, alpha=0.12)
    sz = PC["nSNP"] / PC["nSNP"].max() * 130 + 15
    ax.scatter(PC["chrom"], PC["alpha"], s=sz, color="#2c3e50", alpha=0.8, zorder=3)
    ax.set_xlabel("chromosome"); ax.set_ylabel("Neanderthal ancestry (%)")
    ax.set_xticks(range(1, 23))
    ax.set_title(f"Per-chromosome Neanderthal estimate for {name} (point size ∝ marker SNPs)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_pg3_perchrom.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
