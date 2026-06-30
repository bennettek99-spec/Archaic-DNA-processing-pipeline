#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADMIXTURE-style ancestry components via LEA::snmf (sparse non-negative matrix
factorisation, an ADMIXTURE equivalent that runs natively on Windows R).

  python admixture_lea.py export    # write results/lea/etr.geno + samples.csv
  Rscript tools/run_lea.R           # snmf, picks K, writes results/lea/Q.csv
  python admixture_lea.py plot      # -> results/figures/fig_admixture.png

The clustering set is West-Eurasian ancestry sources + Etruscan/Italian targets
(archaics excluded so they don't dominate K); SNPs are thinned and call-rate
filtered.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results"); FIG = os.path.join(RESULTS, "figures")
LEA = os.path.join(RESULTS, "lea"); os.makedirs(LEA, exist_ok=True)

# curated West-Eurasian ancestry set: label -> (kind, value)
_whg = lambda g: any(s in g for s in ("loschbour", "villabruna", "bichon", "labrana",
                                      "france_mesolithic", "england_mesolithic", "iberia_mesolithic"))
_latin = lambda g: ("latini" in g) or (("lazio_ia" in g) and ("etruscan" not in g))
_itba = lambda g: "italy" in g and any(b in g for b in ("_ba", "_eba", "_mba", "_lba"))
SET = [
    ("Anatolia_N", "grp", "turkey_n"), ("WHG", "grpfn", _whg), ("Yamnaya", "grp", "yamnaya"),
    ("Iran_N", "grp", "iran_ganjdareh_n"), ("Natufian", "grp", "israel_natufian"),
    ("Etruscan", "grp", "etruscan"), ("Latin", "grpfn", _latin),
    ("ImperialRoman", "grp", "imperialroman"), ("ItalyBA", "grpfn", _itba),
    ("Sardinian", "pop", "Sardinian"),
]
ORDER = ["WHG", "Anatolia_N", "Iran_N", "Natufian", "Yamnaya", "ItalyBA",
         "Etruscan", "Latin", "ImperialRoman", "Sardinian"]
MAXN = 40
SNP_STEP = 12     # thin ~1.15M -> ~96k
MIN_CALLRATE = 0.5


def export():
    from archaic.panel import Panel
    from archaic.refs import PANELS
    panel = Panel(PANELS["1240k"]["prefix"])
    gl = panel.ind["pop"].str.lower(); pops = panel.ind["pop"].values
    rng = np.random.default_rng(0)
    cols, labels = [], []
    for label, kind, val in SET:
        if kind == "pop":
            idx = np.where(pops == val)[0]
        elif kind == "grpfn":
            idx = np.where(gl.map(val).fillna(False).values)[0]
        else:
            idx = np.where(gl.str.contains(val.lower(), na=False).values)[0]
        idx = [int(i) for i in idx if int(i) not in cols]
        if len(idx) > MAXN:
            idx = sorted(rng.choice(idx, MAXN, replace=False).tolist())
        for i in idx:
            cols.append(i); labels.append(label)
    cols = np.array(cols, dtype=np.int64)
    snp_rows = panel.snp_rows[::SNP_STEP]
    print(f"reading {len(cols)} individuals x {len(snp_rows):,} SNPs...")
    G = panel.pg.read(snp_rows, cols)                      # (nsnp, nind) int8
    callrate = (G >= 0).mean(1)
    keep = callrate >= MIN_CALLRATE
    G = G[keep]
    print(f"SNPs after call-rate >= {MIN_CALLRATE}: {keep.sum():,}")
    # LEA .geno: one row per SNP, char per individual (0/1/2/9)
    chars = np.where(G < 0, 9, G).astype(np.int8).astype(str)
    with open(os.path.join(LEA, "etr.geno"), "w") as fh:
        for row in chars:
            fh.write("".join(row) + "\n")
    pd.DataFrame({"genetic_id": panel.ind["id"].values[cols], "pop": labels}).to_csv(
        os.path.join(LEA, "samples.csv"), index=False)
    print(f"Wrote {LEA}/etr.geno ({keep.sum():,} SNPs x {len(cols)} ind) + samples.csv")


def compute(K=4):
    """ADMIXTURE-style ancestry coefficients via sparse NMF (the algorithm behind
    LEA::snmf; Frichot et al. 2014), since Bioconductor/LEA will not install on the
    bleeding-edge R. Reads the exported .geno, imputes missing to the per-SNP mean,
    runs non-negative matrix factorisation, and writes per-individual coefficients."""
    from sklearn.decomposition import NMF
    rows = [list(map(int, line.strip())) for line in open(os.path.join(LEA, "etr.geno"))]
    G = np.array(rows, dtype=np.float32).T            # (n_ind, n_snp); 9 = missing
    miss = G == 9
    G[miss] = np.nan
    with np.errstate(invalid="ignore"):
        colmean = np.nanmean(G, axis=0)
    G = np.where(np.isnan(G), colmean[None, :], G)    # impute to per-SNP mean
    var = G.var(0) > 0
    G = G[:, var]
    print(f"NMF on {G.shape[0]} individuals x {G.shape[1]:,} SNPs, K={K}")
    W = NMF(n_components=K, init="nndsvda", max_iter=500, random_state=0).fit_transform(G)
    Q = W / W.sum(1, keepdims=True)                   # rows sum to 1 = ancestry coefficients
    pd.DataFrame(Q).to_csv(os.path.join(LEA, "Q.csv"), index=False, header=False)
    print(f"Wrote {LEA}/Q.csv  (K={K})")


def plot():
    samp = pd.read_csv(os.path.join(LEA, "samples.csv"))
    Q = pd.read_csv(os.path.join(LEA, "Q.csv"), header=None).values
    K = Q.shape[1]
    samp = samp.assign(**{f"Q{k}": Q[:, k] for k in range(K)})
    # order individuals: by population (ORDER), then by dominant component
    samp["popord"] = samp["pop"].map({p: i for i, p in enumerate(ORDER)}).fillna(99)
    samp["dom"] = Q.argmax(1)
    samp = samp.sort_values(["popord", "dom"] + [f"Q{k}" for k in range(K)]).reset_index(drop=True)
    Qs = samp[[f"Q{k}" for k in range(K)]].values
    palette = ["#e0a458", "#7b9e89", "#5c7a99", "#b5651d", "#8c6d8c", "#9aa14a"][:K]
    # label each component by the source population where it is most abundant
    srcrefs = ["WHG", "Anatolia_N", "Iran_N", "Yamnaya", "Natufian"]
    means = samp.groupby("pop", sort=False)[[f"Q{k}" for k in range(K)]].mean()
    klab = {}
    for k in range(K):
        m = means[f"Q{k}"]
        top = m[[p for p in srcrefs if p in m.index]].idxmax()
        klab[k] = {"WHG": "WHG-like", "Anatolia_N": "Anatolian farmer",
                   "Iran_N": "Iran/CHG", "Yamnaya": "Steppe", "Natufian": "Levantine"}.get(top, top)
    import matplotlib.patches as mp
    fig, ax = plt.subplots(figsize=(11, 3.8))
    bottom = np.zeros(len(samp))
    for k in range(K):
        ax.bar(range(len(samp)), Qs[:, k], bottom=bottom, width=1.0,
               color=palette[k], linewidth=0)
        bottom += Qs[:, k]
    # population dividers + labels
    bounds = samp.groupby("pop", sort=False).size().cumsum().values
    starts = np.concatenate([[0], bounds[:-1]])
    for b in bounds[:-1]:
        ax.axvline(b, color="white", lw=1.2)
    ax.set_xticks((starts + bounds) / 2)
    ax.set_xticklabels(samp.groupby("pop", sort=False).size().index, rotation=40, ha="right", fontsize=8)
    ax.set_xlim(0, len(samp)); ax.set_ylim(0, 1); ax.set_ylabel(f"ancestry (K={K})")
    ax.legend(handles=[mp.Patch(color=palette[k], label=klab[k]) for k in range(K)],
              ncol=K, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.16), frameon=False)
    ax.set_title(f"Sparse-NMF ancestry components (K={K}, snmf-style) — descriptive", pad=24)
    fig.tight_layout(); p = os.path.join(FIG, "fig_admixture.png"); fig.savefig(p, dpi=150); plt.close(fig)
    print(f"Wrote {p}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "export"
    {"export": export, "compute": compute, "plot": plot}.get(cmd, export)()
