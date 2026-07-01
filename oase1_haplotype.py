#!/usr/bin/env python
"""
Oase1 archaic-haplotype segment analysis (1240K array resolution).

Reproduces, at capture-array resolution, the *spatial segment* method that
Fu et al. (2015, Nature 524:216) used on the shotgun Oase1 genome: plot the
positions at which the individual carries the Neanderthal-derived allele along
each chromosome, find long contiguous Neanderthal-ancestry segments, and use
their genetic lengths to date the admixture (expected introgressed segment
length after g generations of recombination ~ 100/g cM).

We use the AADR 1240K genotypes we already have for Oase1 plus the in-panel
archaic (Altai, Vindija, Denisova) and African-outgroup (Mbuti+Yoruba) samples,
and the genetic map (Morgans) stored in the .snp file. This is *array
resolution* (~1 informative SNP per Mb) -- it can localise the long blocks that
distinguish a recent admixture from background introgression, but cannot measure
short (<~5 cM) segments. A read-level (BAM) version using hmmix is provided
separately in oase1_bam_pipeline/ for full resolution.

Output (reports/oase1_haplotype/):
  oase1_segments.csv          detected Neanderthal segments (chrom, cM span, ...)
  oase1_informative_sites.csv per-site carrier track (Oase1 vs baselines)
  fig_oase1_karyogram.png     genome-wide segment map (Oase1 vs baseline)
  fig_oase1_segment_lengths.png  segment-length distribution + generations estimate
Run: PYTHONIOENCODING=utf-8 python oase1_haplotype.py
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic.refs import PANELS

OUT = "reports/oase1_haplotype"
os.makedirs(OUT, exist_ok=True)

# ---- parameters -------------------------------------------------------------
AFR_MAX = 0.05        # African-outgroup derived-freq ceiling for an informative site
NEA_MIN = 0.5         # Altai & Vindija must each carry the derived allele (>=1 copy)
# HMM emissions are set from the DATA: E0 = the observed genome-wide carrier rate
# at informative sites in ~2% baselines (ILS + pseudo-haploid error background);
# E1 = the within-block carrier rate (limited by the same noise + single-read
# pseudo-haploid sampling, so ~2x background, not ~1). Sensitivity reported below.
E0 = 0.085            # P(carry archaic allele | background)
E1 = 0.50             # P(carry archaic allele | archaic segment)
TAU_CM = 25.0         # HMM prior mean segment length (cM)
PI1 = 0.08            # prior fraction of genome in archaic state (~Oase1 alpha)
MIN_SITES = 4         # a reported segment must span >= this many informative SNPs

# Oase1's two AADR libraries are merged (union of covered sites, damage-restricted
# call preferred where both exist) to maximise the sparse informative-site overlap.
OASE_LIBS = ["Oase1_d.AG.BY.AA", "Oase1.AG.BY.AA"]
COMPARE = {           # baselines for contrast (penecontemporaneous + Holocene)
    "Ust_Ishim.DG":     "Ust'-Ishim (45 kya, ~2%)",
    "Kostenki14.SG":    "Kostenki14 (38 kya, ~2%)",
    "Loschbour.AG":     "Loschbour (Mesolithic, ~2%)",
}


def viterbi(carry, gpos_cm):
    """2-state Viterbi over ordered informative sites.
    carry: 0/1 array (carries archaic-derived allele). gpos_cm: cM positions.
    Returns state path (0=background, 1=archaic)."""
    n = len(carry)
    logE = np.array([[np.log(1 - E0), np.log(E0)],
                     [np.log(1 - E1), np.log(E1)]])          # [state][obs]
    d = np.diff(gpos_cm)
    d = np.clip(d, 0, None)
    V = np.full((n, 2), -np.inf); B = np.zeros((n, 2), dtype=int)
    V[0, 0] = np.log(1 - PI1) + logE[0, carry[0]]
    V[0, 1] = np.log(PI1) + logE[1, carry[0]]
    for i in range(1, n):
        pstay = np.exp(-d[i - 1] / TAU_CM)                   # stay in same state
        pswap = 1 - pstay
        # transition matrix row-normalised toward stationary (PI1)
        T = np.array([[pstay + pswap * (1 - PI1), pswap * PI1],
                      [pswap * (1 - PI1),        pstay + pswap * PI1]])
        logT = np.log(np.clip(T, 1e-12, None))
        for s in range(2):
            cand = V[i - 1] + logT[:, s]
            B[i, s] = int(np.argmax(cand))
            V[i, s] = cand[B[i, s]] + logE[s, carry[i]]
    path = np.zeros(n, dtype=int); path[-1] = int(np.argmax(V[-1]))
    for i in range(n - 1, 0, -1):
        path[i - 1] = B[i, path[i]]
    return path


def main():
    cfg = PANELS["1240k"]
    print("Loading 1240K panel ...")
    panel = Panel(cfg["prefix"])
    snp = panel.snp.iloc[panel.snp_rows].reset_index(drop=True)
    chrom = snp["chrom"].to_numpy()
    gpos_cm = snp["gpos"].to_numpy() * 100.0                 # Morgans -> cM
    pos = snp["pos"].to_numpy()

    # reference frequencies (allele-1 dosage/2)
    rf, ri = panel.frequencies({
        "Altai": cfg["refs"]["Altai"], "Vindija": cfg["refs"]["Vindija"],
        "Denisova": cfg["refs"]["Denisova"], "Chimp": cfg["refs"]["Chimp"],
        "Afr": {"pops": ["Mbuti", "Yoruba", "YRI", "YRI-Discovery"]},
    })
    for k in rf:
        print(f"  {k:9s} covered SNPs = {ri[k]['n_snp_covered']:,}")

    # polarise by chimp: derived-allele freq = allele-1 freq if chimp~0 else 1-freq
    pChi = rf["Chimp"]
    anc_is_a1 = pChi < 0.5                                   # chimp carries a1 -> a1 ancestral
    def derived(p):
        return np.where(anc_is_a1, p, 1.0 - p)
    dAlt, dVin, dDen, dAfr = map(derived, (rf["Altai"], rf["Vindija"],
                                           rf["Denisova"], rf["Afr"]))

    # Neanderthal-informative sites (NIS): Altai (high-coverage, 1.15M SNPs)
    # carries the derived allele, African outgroup lacks it, chimp state defined.
    # A high-confidence subset additionally requires Vindija to confirm.
    ok = np.isfinite(pChi) & (np.abs(pChi - 0.5) > 0.4) \
        & np.isfinite(dAlt) & np.isfinite(dAfr)
    nis = ok & (dAlt >= NEA_MIN) & (dAfr < AFR_MAX)
    nis_hc = nis & np.isfinite(dVin) & (dVin >= NEA_MIN)    # Altai+Vindija confirmed
    dis = ok & np.isfinite(dDen) & (dDen >= NEA_MIN) & (dAfr < AFR_MAX) \
        & (dAlt < NEA_MIN)                                  # Denisovan-specific
    print(f"\nNeanderthal-informative sites (Altai-based): {nis.sum():,}"
          f"  (Vindija-confirmed subset: {nis_hc.sum():,})")
    print(f"Denisovan-specific sites:                    {dis.sum():,}")

    # read Oase1 libraries + baselines
    want = [i for i in OASE_LIBS + list(COMPARE) if i in panel._id_to_col]
    cols = np.array([panel._id_to_col[i] for i in want])
    G = panel.pg.read(panel.snp_rows, cols)                 # (n_snp, k) 0/1/2/-1
    dose_all = np.where(anc_is_a1[:, None], G, np.where(G < 0, -1, 2 - G))
    colof = {sid: k for k, sid in enumerate(want)}

    # merge Oase1 libraries: prefer damage-restricted call, fall back to standard
    libs = [colof[l] for l in OASE_LIBS if l in colof]
    oase = np.full(len(chrom), -1, dtype=np.int8)
    for k in libs:                                          # later libs fill gaps
        oase = np.where((oase < 0) & (dose_all[:, k] >= 0), dose_all[:, k], oase)

    # working set: merged Oase1 first, then baselines
    ids = ["Oase1"] + [c for c in COMPARE if c in colof]
    LABELS = {"Oase1": "Oase1 (merged libraries)", **COMPARE}
    derived_dose = np.column_stack(
        [oase] + [dose_all[:, colof[c]] for c in COMPARE if c in colof])

    # ------------------------------------------------------------------ diag
    print("\nGenome-wide Neanderthal-allele carrier rate at informative sites:")
    site_track = {"chrom": chrom[nis], "pos": pos[nis], "gpos_cm": gpos_cm[nis]}
    for j, sid in enumerate(ids):
        dd = derived_dose[nis, j]
        cov = dd >= 0
        carry = (dd >= 1).astype(int)
        rate = carry[cov].mean() if cov.sum() else np.nan
        print(f"  {LABELS[sid]:32s} nSNP={cov.sum():6d}  carrier_rate={rate:.3f}")
        site_track[f"carry__{sid}"] = np.where(cov, carry, -1)

    pd.DataFrame(site_track).to_csv(f"{OUT}/oase1_informative_sites.csv", index=False)

    # ------------------------------------------------- segment calling (HMM)
    all_segs = []
    for j, sid in enumerate(ids):
        dd = derived_dose[nis, j]
        cov = dd >= 0
        for ch in [str(c) for c in range(1, 23)]:
            m = (chrom[nis] == ch) & cov
            if m.sum() < MIN_SITES:
                continue
            g = gpos_cm[nis][m]
            order = np.argsort(g)
            g = g[order]; c = (dd[m][order] >= 1).astype(int)
            path = viterbi(c, g)
            # extract maximal state-1 runs
            i = 0
            while i < len(path):
                if path[i] == 1:
                    k = i
                    while k + 1 < len(path) and path[k + 1] == 1:
                        k += 1
                    span = g[k] - g[i]
                    nsit = k - i + 1
                    if nsit >= MIN_SITES:
                        all_segs.append(dict(
                            sample=sid, label=LABELS[sid], chrom=ch,
                            start_cM=round(float(g[i]), 2),
                            end_cM=round(float(g[k]), 2),
                            span_cM=round(float(span), 2), n_sites=int(nsit),
                            carrier_frac=round(float(c[i:k + 1].mean()), 3)))
                    i = k + 1
                else:
                    i += 1

    segs = pd.DataFrame(all_segs).sort_values(["sample", "span_cM"],
                                              ascending=[True, False])
    segs.to_csv(f"{OUT}/oase1_segments.csv", index=False)

    print("\nDetected Neanderthal segments (>= %d informative sites):" % MIN_SITES)
    for sid in ids:
        s = segs[segs["sample"] == sid]
        if not len(s):
            print(f"  {COMPARE[sid]:32s} none"); continue
        tot = s["span_cM"].sum()
        print(f"  {LABELS[sid]:32s} n={len(s):3d}  total={tot:7.1f} cM  "
              f"longest={s['span_cM'].max():5.1f} cM  "
              f">50cM={int((s.span_cM>50).sum())}  >30cM={int((s.span_cM>30).sum())}")

    # ---- generations estimate from Oase1 segment lengths (mean L ~ 100/g)
    key = "Oase1"
    so = segs[(segs["sample"] == key) & (segs.span_cM > 0)]
    Lmax = so.span_cM.max() if len(so) else np.nan
    Ltop = so.span_cM.nlargest(3).mean() if len(so) else np.nan
    print(f"\nOase1 longest segment = {Lmax:.1f} cM  ->  g ~= 100/L = {100/Lmax:.1f} gen")
    print(f"Oase1 mean of top-3   = {Ltop:.1f} cM  ->  g ~= {100/Ltop:.1f} gen  "
          f"(array resolution fragments blocks -> lengths are LOWER bounds, "
          f"so g is an UPPER bound)")

    # ============================================================== FIGURES
    CH = [str(c) for c in range(1, 23)]
    chrom_nis = chrom[nis]; g_nis = gpos_cm[nis]
    C_OASE = "#d7301f"; C_TICK = "#333333"; C_BASE = "#2c7fb5"

    # --- FIG 1: two-column karyogram, Oase1 vs a ~2% baseline ---------------
    base = "Loschbour.AG" if "Loschbour.AG" in ids else ids[-1]
    panels = [("Oase1", "Oase1 (merged) — Neanderthal segments"),
              (base, f"{LABELS[base]} — Neanderthal segments")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
    for ax, (sid, ttl) in zip(axes, panels):
        j = ids.index(sid)
        dd = derived_dose[nis, j]
        for yi, ch in enumerate(CH):
            y = len(CH) - yi
            m = chrom_nis == ch
            gmax = g_nis[m].max() if m.any() else 0
            ax.plot([0, gmax], [y, y], color="#dddddd", lw=6, solid_capstyle="round", zorder=1)
            car = m & (dd >= 1)
            ax.scatter(g_nis[car], np.full(car.sum(), y), s=5, c=C_TICK, zorder=3)
            s = segs[(segs["sample"] == sid) & (segs.chrom == ch)]
            for _, r in s.iterrows():
                ax.plot([r.start_cM, r.end_cM], [y, y],
                        color=C_OASE if sid == "Oase1" else C_BASE,
                        lw=7, solid_capstyle="butt", zorder=2, alpha=.9)
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("genetic position (cM)")
        ax.set_xlim(-2, None)
    axes[0].set_yticks(range(1, len(CH) + 1))
    axes[0].set_yticklabels(CH[::-1]); axes[0].set_ylabel("chromosome")
    fig.suptitle("Oase1 carries Neanderthal ancestry in a few LONG blocks; a ~2% "
                 "baseline carries it in short scattered ones\n"
                 "(black ticks = Neanderthal-allele sites; bars = HMM-detected "
                 "segments; 1240K array resolution)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUT}/fig_oase1_karyogram.png", dpi=130); plt.close(fig)

    # --- FIG 2: segment-length distribution + generations ------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    order = ["Oase1"] + [c for c in COMPARE if c in ids]
    colors = [C_OASE] + ["#7fb0d3", "#9ecae1", "#c6dbef"][:len(order) - 1]
    for sid, col in zip(order, colors):
        s = segs[segs["sample"] == sid].span_cM.values
        if len(s):
            a1.hist(s, bins=np.arange(0, 45, 3), alpha=.55, label=LABELS[sid],
                    color=col, edgecolor="k", lw=.4)
    a1.axvline(30, color="k", ls=":", lw=1)
    a1.set_xlabel("segment length (cM)"); a1.set_ylabel("number of segments")
    a1.set_title("Detected Neanderthal segment lengths"); a1.legend(fontsize=8)
    # longest-segment bar
    lm = [(LABELS[s], segs[segs["sample"] == s].span_cM.max() if
           (segs["sample"] == s).any() else 0) for s in order]
    a2.barh([x[0] for x in lm][::-1], [x[1] for x in lm][::-1],
            color=[C_OASE if "Oase1" in x[0] else C_BASE for x in lm][::-1])
    a2.set_xlabel("longest detected segment (cM)")
    a2.set_title("Longest Neanderthal segment (Oase1 vs ~2%% baselines)\n"
                 "%.1f cM $\\rightarrow$ ancestor within ~%.0f generations"
                 % (Lmax, round(100 / Lmax)), fontsize=10)
    a2.margins(x=0.15)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_oase1_segment_lengths.png", dpi=130)
    plt.close(fig)

    print(f"\nWrote figures + CSVs to {OUT}/")
    return dict(nis=int(nis.sum()), nis_hc=int(nis_hc.sum()), segs=segs,
                Lmax=float(Lmax), key=key, labels=LABELS)


if __name__ == "__main__":
    res = main()
    print("\n(diagnostics complete; figures written by oase1_haplotype_figs)")
