#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 7 — automated investigation reports for the top residual candidates.

For each of the strongest +/- standardised-residual individuals (high-confidence
only) we emit a structured report that deliberately separates EVIDENCE from
SPECULATION and lists both biological and technical explanations, plus the
genetically nearest neighbours used to form the expectation. Because Phase 6/9
found nothing significant after multiple-testing, every report states up front
that the candidate is a HYPOTHESIS, not a discovery.

Output:
  results/phase7_reports/<rank>_<id>.md   one per candidate
  results/phase7_reports/INDEX.md
"""
import os, sys, re
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import stats

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
OUT = os.path.join(RESULTS, "phase7_reports")
os.makedirs(OUT, exist_ok=True)
PANEL = sys.argv[1] if len(sys.argv) > 1 else "1240k"
PCS = [f"PC{i}" for i in range(1, 7)]
FEATURES = PCS + ["lon", "lat", "date_bp"]
FEAT_WEIGHT = {**{p: 1.0 for p in PCS}, "lon": 0.5, "lat": 0.5, "date_bp": 0.7}
N_EACH = 10
K = 80

ANCESTRY_FLAGS = ["-o", "african", "levant", "steppe", "siberia", "han", "east",
                  "natufian", "iran", "aasi", "andamanese", "jomon"]


def zscore(x):
    x = x.astype(np.float64); mu, sd = np.nanmean(x), np.nanstd(x)
    return (x - mu) / (sd if sd > 0 else 1.0)


def main():
    df = pd.read_csv(os.path.join(RESULTS, f"phase6_{PANEL}_residuals.csv"))
    pc = pd.read_csv(os.path.join(RESULTS, f"phase5_{PANEL}_pca.csv"))
    an = pd.read_csv(os.path.join(RESULTS, f"phase4_{PANEL}_analysis.csv"))[
        ["genetic_id", "locality", "hapconx_contam", "angsd_contam", "coverage"]]
    df = df.merge(pc, on="genetic_id").merge(an, on="genetic_id").reset_index(drop=True)
    n_hc = int(df["high_conf"].sum())
    zc = stats.norm.isf(0.025 / n_hc)

    # alpha_adj already computed in Phase 6; build the same feature space
    X = np.column_stack([zscore(df[f].fillna(df[f].median()).values) * FEAT_WEIGHT[f]
                         for f in FEATURES])

    hc = df["high_conf"].values
    Xref, idref = X[hc], np.where(hc)[0]
    tree = cKDTree(Xref)

    cand = df[df["high_conf"] & df["z_resid"].notna()]
    top = pd.concat([cand.sort_values("z_resid", ascending=False).head(N_EACH),
                     cand.sort_values("z_resid").head(N_EACH)])

    index = ["# Phase 7 — candidate investigation reports",
             "", f"Panel {PANEL}. **None of these is statistically significant** "
             f"(all |z| < Bonferroni z*={zc:.2f}; 0 pass FDR). They are the most "
             "deviant individuals and the only defensible targets for future deeper "
             "sequencing — hypotheses, not discoveries.", "",
             "| rank | id | z | obs% | exp% | group | likely explanation |",
             "|---|---|---|---|---|---|---|"]

    for rank, (_, r) in enumerate(top.iterrows(), 1):
        gid = r["genetic_id"]
        ri = df.index[df["genetic_id"] == gid][0]
        _, nbr = tree.query(X[ri], k=K + 1)
        nbr_rows = df.iloc[idref[nbr]]
        nbr_rows = nbr_rows[nbr_rows["genetic_id"] != gid].head(8)

        # heuristic explanations
        gl = str(r["group_id"]).lower()
        tech = []
        if r["alpha_nSNP"] < 300_000: tech.append(f"modest coverage ({int(r['alpha_nSNP']):,} SNP)")
        if "AG" == r["data_type"]: tech.append("pseudo-haploid capture (.AG)")
        cl = re.search(r"\[([0-9.]+)", str(r.get("hapconx_contam", "")))
        if cl and float(cl.group(1)) > 0.01: tech.append(f"contamination LB {cl.group(1)} (deflates Neanderthal)")
        bio = []
        if "-o" in gl: bio.append("AADR-flagged ancestry OUTLIER ('-o') — atypical ancestry for its context")
        if r["z_resid"] < 0 and any(a in gl for a in ["african", "levant"]):
            bio.append("recent African/Levantine admixture → LESS Neanderthal (expected, not anomalous)")
        if r["z_resid"] > 0 and any(a in gl for a in ["han", "east", "siberia", "steppe"]):
            bio.append("East-Eurasian-related ancestry carries slightly MORE Neanderthal")
        if not bio:
            bio.append("no obvious ancestry driver; most consistent with statistical noise "
                       "(extreme order statistic of ~9,000 tests)")

        likely = (bio[0] if bio else "noise")[:60]
        index.append(f"| {rank} | {gid} | {r['z_resid']:+.2f} | {r['alpha_Nea']*100:.2f} "
                     f"| {r['expected_Nea']*100:.2f} | {str(r['group_id'])[:34]} | {likely} |")

        md = []
        md.append(f"# Candidate {rank}: {gid}")
        md.append(f"\n> **Status: HYPOTHESIS, not a discovery.** |z|={abs(r['z_resid']):.2f} "
                  f"< Bonferroni z*={zc:.2f}; does not survive multiple-testing correction.\n")
        md.append("## Evidence (measured)")
        md.append(f"- **Group / culture:** {r['group_id']}")
        md.append(f"- **Site / locality:** {r.get('locality','?')}")
        md.append(f"- **Country:** {r['country']}    **Coordinates:** {r['lat']}, {r['lon']}")
        md.append(f"- **Date:** {int(r['date_bp']) if pd.notna(r['date_bp']) else '?'} BP")
        md.append(f"- **Observed Neanderthal:** {r['alpha_Nea']*100:.2f}% "
                  f"(adj {r['alpha_adj']*100:.2f}%, SE {r['alpha_SE']*100:.2f}%)")
        md.append(f"- **Expected (ancestry+geo+time):** {r['expected_Nea']*100:.2f}%  "
                  f"→ residual {r['residual_Nea']*100:+.2f} pp,  **z = {r['z_resid']:+.2f}**")
        md.append(f"- **Quality:** {int(r['alpha_nSNP']):,} usable SNP, data type {r['data_type']}, "
                  f"contamination {r.get('hapconx_contam','?')}, flags: {r.get('flags','') or 'none'}")
        md.append(f"- **Neanderthal D Z:** {r['D_Nea_Z']:.1f}   **Denisovan D Z:** {r['D_Den_Z']:.1f}")
        md.append("\n## Genetically nearest neighbours (the expectation)")
        md.append("| neighbour | group | Neanderthal% | SNP |")
        md.append("|---|---|---|---|")
        for _, nr in nbr_rows.iterrows():
            md.append(f"| {nr['genetic_id']} | {str(nr['group_id'])[:36]} | "
                      f"{nr['alpha_Nea']*100:.2f} | {int(nr['alpha_nSNP']):,} |")
        md.append("\n## Possible biological explanations")
        for b in bio: md.append(f"- {b}")
        md.append("\n## Possible technical explanations")
        for t in (tech or ["coverage/contamination unremarkable for a high-confidence sample"]):
            md.append(f"- {t}")
        md.append("\n## Confidence assessment")
        md.append("- **Not significant** after multiple-testing (Phase 6) and **not robust** "
                  "to perturbation as a *significant* hit (Phase 9).")
        md.append("- Recommended only as a target for higher-coverage shotgun re-sequencing "
                  "and explicit local-ancestry / admixture analysis before any claim.")
        with open(os.path.join(OUT, f"{rank:02d}_{re.sub(r'[^A-Za-z0-9_.-]','_',gid)}.md"),
                  "w", encoding="utf-8") as fh:
            fh.write("\n".join(md) + "\n")

    with open(os.path.join(OUT, "INDEX.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(index) + "\n")
    print(f"Wrote {2*N_EACH} candidate reports + INDEX to {OUT}")


if __name__ == "__main__":
    main()
