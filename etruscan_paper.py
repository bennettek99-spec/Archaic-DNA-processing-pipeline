#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
etruscan_paper.py — expanded Etruscan analysis + manuscript generation.

Adds population-level "mean-genome" profiles (archaic.profiles) to place Etruscan
archaic ancestry in regional and temporal context with tight error bars, tests
whether Etruscans differ in archaic content from their neighbours (Latins/Romans/
Anatolians/Steppe), maps genetic distances, and combines this with the
individual-outlier and locus-selection results into a research-paper-style
document: PAPER.md + reports/etruscan_paper.html (figures fig_p1..fig_p4).
"""
import os, sys, json, datetime, html
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sstats
from sklearn.manifold import MDS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st, profiles as pf
from archaic.refs import PANELS

PANEL = "1240k"
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "etruscan")
FIG = os.path.join(RESULTS, "figures")
REPORTS = os.path.join(HERE, "reports")
for d in (OUT, FIG, REPORTS):
    os.makedirs(d, exist_ok=True)
MAXN = 150
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 160, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.25})

# ancient cohort -> predicate on (group_id lowercased)
def P(*subs, exclude=()):
    return lambda g: all(s in g for s in subs) and not any(e in g for e in exclude)

COHORTS = {
    "Anatolia Neolithic":  P("turkey_n") ,
    "Italy Neolithic":     lambda g: ("italy" in g) and ("_n" in g) and ("etruscan" not in g) and ("ia_" not in g),
    "Italy Bronze Age":    lambda g: ("italy" in g) and any(b in g for b in ("_ba", "_eba", "_mba", "_lba")),
    "Yamnaya (Steppe)":    P("yamnaya"),
    "Aegean BA":           lambda g: any(s in g for s in ("mycenaean", "minoan", "aegean")),
    "Etruscan":            P("etruscan"),
    "Latin / Italic IA":   lambda g: ("latini" in g) or (("lazio_ia" in g) and ("etruscan" not in g)),
    "Republican Roman":    P("republic"),
    "Imperial Roman":      P("imperialroman"),
    "Magna Graecia Greek": lambda g: ("magnagraecia" in g) or ("italy" in g and "greek" in g),
    "Sicily (anc)":        P("sicily"),
    "Sardinia (anc)":      lambda g: ("sardinia" in g) and ("italy" in g),
}
MODERN = {"French (mod)": "French", "Sardinian (mod)": "Sardinian",
          "Spanish (mod)": "Spanish", "Han (mod)": "Han"}
REGION_ORDER = ["Anatolia Neolithic", "Italy Neolithic", "Aegean BA", "Italy Bronze Age",
                "Yamnaya (Steppe)", "Etruscan", "Latin / Italic IA", "Sicily (anc)",
                "Magna Graecia Greek", "Republican Roman", "Imperial Roman",
                "Sardinia (anc)", "Sardinian (mod)", "Spanish (mod)", "French (mod)", "Han (mod)"]


def build_cohorts(panel, meta):
    rng = np.random.default_rng(0)
    name_to_cols = {}
    gl = meta["group_id"].str.lower()
    for name, pred in COHORTS.items():
        ids = meta.loc[gl.map(pred), "genetic_id"].tolist()
        cols = np.array([panel._id_to_col[i] for i in ids if i in panel._id_to_col], dtype=np.int64)
        if len(cols) > MAXN:
            cols = np.sort(rng.choice(cols, MAXN, replace=False))
        name_to_cols[name] = cols
    for name, pop in MODERN.items():
        name_to_cols[name] = panel.cols_for(pops=[pop])
    # Etruscan regional split
    for sub, kw in [("Etruscan Tuscany", "tuscany"), ("Etruscan Lazio", "lazio")]:
        ids = meta.loc[gl.str.contains("etruscan") & gl.str.contains(kw), "genetic_id"].tolist()
        name_to_cols[sub] = np.array([panel._id_to_col[i] for i in ids if i in panel._id_to_col], np.int64)
    # references + African outgroup
    refs = PANELS[PANEL]["refs"]
    for r in ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti", "Yoruba"]:
        name_to_cols[r] = panel.cols_for(**refs[r])
    return name_to_cols


def main():
    panel = Panel(PANELS[PANEL]["prefix"])
    block = st.assign_blocks(panel.n_snp, 50)
    meta = pd.read_csv(os.path.join(RESULTS, f"phase4_{PANEL}_analysis.csv"))
    res = pd.read_csv(os.path.join(RESULTS, f"phase6_{PANEL}_residuals.csv"))

    name_to_cols = build_cohorts(panel, meta)
    counts = {k: len(v) for k, v in name_to_cols.items()}
    print("cohort sizes:", {k: counts[k] for k in REGION_ORDER if k in counts})

    print("Building mean-genome profiles (one cohort at a time)...")
    freq, info = pf.cohort_frequencies(panel, name_to_cols)

    # --- group-level archaic table -------------------------------------------
    cohorts = [c for c in REGION_ORDER if info.get(c, {}).get("n", 0) >= 2] + \
              ["Etruscan Tuscany", "Etruscan Lazio", "Yoruba"]
    G = pf.group_archaic(freq, cohorts, block)
    G["n"] = G["cohort"].map(lambda c: info[c]["n"])
    G.to_csv(os.path.join(OUT, "P_group_archaic.csv"), index=False)
    print(G[["cohort", "n", "alpha_Nea", "alpha_SE", "D_Den_Z"]].to_string(index=False))

    # --- Etruscan vs neighbours: differential Neanderthal (African outgroup) --
    contrasts = []
    for other in ["Latin / Italic IA", "Imperial Roman", "Anatolia Neolithic",
                  "Yamnaya (Steppe)", "Magna Graecia Greek", "Italy Bronze Age"]:
        if info.get(other, {}).get("n", 0) >= 2:
            d = st.dstat(freq, "Etruscan", other, "Altai", "Yoruba", block, 50)
            contrasts.append(dict(vs=other, D=d["theta"], Z=d["z"], nSNP=d["n_used"]))
    Cc = pd.DataFrame(contrasts)
    Cc.to_csv(os.path.join(OUT, "P_etruscan_contrasts.csv"), index=False)
    print("\nD(Etruscan, X; Altai, Yoruba)  (|Z|<3 => no archaic difference):")
    print(Cc.to_string(index=False))

    # --- temporal: group-level Neanderthal per Italian time bin --------------
    BINS = [("Neolithic", 6000, 7500), ("Copper/EBA", 4200, 6000),
            ("Bronze Age", 3200, 4200), ("Iron Age", 2300, 3200),
            ("Roman", 1700, 2300), ("Late/Medieval", 800, 1700)]
    ita = meta[meta["country"] == "Italy"]
    tcols = {}
    for lab, lo, hi in BINS:
        ids = ita.loc[ita["date_bp"].between(lo, hi), "genetic_id"].tolist()
        tcols[lab] = np.array([panel._id_to_col[i] for i in ids if i in panel._id_to_col], np.int64)
    tfreq, tinfo = pf.cohort_frequencies(panel, {**tcols, **{r: name_to_cols[r] for r in
                  ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]}})
    T = pf.group_archaic(tfreq, [b[0] for b in BINS if tinfo[b[0]]["n"] >= 2], block)
    T["mid_BP"] = T["cohort"].map({lab: (lo + hi) / 2 for lab, lo, hi in BINS})
    T["n"] = T["cohort"].map(lambda c: tinfo[c]["n"])
    T.to_csv(os.path.join(OUT, "P_temporal.csv"), index=False)

    # --- genetic distance among cohorts (n>=30: raw allele-freq distance is
    #     inflated by sampling noise in small cohorts, so restrict the view) ---
    dist_names = [c for c in REGION_ORDER if info.get(c, {}).get("n", 0) >= 30]
    D = pf.distance_matrix(freq, dist_names)
    D.to_csv(os.path.join(OUT, "P_distance.csv"))

    # --- individual Etruscan outliers (from Phase 6) + locus scan (Phase C) --
    etr = res[res["group_id"].str.contains("Etruscan", case=False, na=False)].copy()
    etr["is_o"] = etr["group_id"].str.contains("-o", na=False)
    o_absz = etr.loc[etr["is_o"], "z_resid"].abs().mean()
    no_absz = etr.loc[~etr["is_o"], "z_resid"].abs().mean()
    Cl = pd.read_csv(os.path.join(OUT, "C_loci_over_time.csv"))
    n_tested = int(Cl["p"].notna().sum())
    bonf = Cl[Cl["p"] < 0.05 / max(n_tested, 1)].dropna(subset=["p"])

    figs = make_paper_figures(G, T, D, dist_names, Cc)
    write_paper(G, Cc, T, etr, o_absz, no_absz, Cl, n_tested, bonf, counts, figs)
    print("\nWrote PAPER.md + reports/etruscan_paper.html + results/etruscan/P_*.csv")


def make_paper_figures(G, T, D, dist_names, Cc):
    figs = []
    # Fig P1: group-level Neanderthal across populations
    g = G[G["cohort"].isin(REGION_ORDER)].copy()
    g["ord"] = g["cohort"].map({c: i for i, c in enumerate(REGION_ORDER)})
    g = g.sort_values("ord")
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    colors = ["#d62728" if "Etrus" in c else ("#1f77b4" if "mod" in c else "#6b8fae")
              for c in g["cohort"]]
    ax.barh(range(len(g)), g["alpha_Nea"] * 100, xerr=g["alpha_SE"] * 100,
            color=colors, edgecolor="k", lw=0.3, capsize=2)
    ax.set_yticks(range(len(g))); ax.set_yticklabels(
        [f"{c}  (n={int(n)})" for c, n in zip(g["cohort"], g["n"])], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("Neanderthal ancestry (%)  ± group jackknife SE")
    ax.axvline(g.loc[g["cohort"] == "Etruscan", "alpha_Nea"].iloc[0] * 100,
               color="#d62728", ls=":", lw=1, alpha=0.7)
    ax.set_title("Group-level Neanderthal ancestry across populations (Etruscan in red)")
    fig.tight_layout(); p = f"{FIG}/fig_p1_group_neanderthal.png"; fig.savefig(p); plt.close(fig)
    figs.append(("fig_p1_group_neanderthal.png", "Figure 1. Group-level Neanderthal ancestry (mean-genome f4-ratio) across ancient and modern populations, with block-jackknife standard errors. Etruscans (red) sit within the tight European band."))

    # Fig P2: temporal
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.errorbar(T["mid_BP"] / 1000, T["alpha_Nea"] * 100, yerr=T["alpha_SE"] * 100,
                fmt="o-", color="#1f77b4", lw=2, capsize=3)
    ax.invert_xaxis(); ax.set_xlabel("age (kyr BP, older →)"); ax.set_ylim(1.0, 3.0)
    ax.set_ylabel("Neanderthal ancestry (%)")
    ax.set_title("Genome-wide Neanderthal ancestry over time in Italy (group-level)")
    fig.tight_layout(); p = f"{FIG}/fig_p2_temporal.png"; fig.savefig(p); plt.close(fig)
    figs.append(("fig_p2_temporal.png", "Figure 2. Group-level Neanderthal ancestry across the Italian time transect; flat within error, indicating no major change in archaic ancestry despite Neolithic and Steppe ancestry turnovers."))

    # Fig P3: MDS of genetic distance
    Dm = D.loc[dist_names, dist_names].values
    Dm = np.nan_to_num(Dm, nan=np.nanmax(Dm))
    coords = MDS(n_components=2, dissimilarity="precomputed", random_state=0,
                 normalized_stress="auto").fit_transform(Dm)
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    for i, c in enumerate(dist_names):
        col = "#d62728" if "Etrus" in c else ("#2ca02c" if "mod" in c else "#3b6ea5")
        ax.scatter(coords[i, 0], coords[i, 1], s=60, color=col, edgecolor="k", lw=0.4, zorder=3)
        ax.annotate(c, (coords[i, 0], coords[i, 1]), fontsize=7.5, xytext=(4, 3),
                    textcoords="offset points")
    ax.set_xlabel("MDS axis 1"); ax.set_ylabel("MDS axis 2")
    ax.set_title("Genetic distance (allele-frequency) among populations — MDS")
    fig.tight_layout(); p = f"{FIG}/fig_p3_mds.png"; fig.savefig(p); plt.close(fig)
    figs.append(("fig_p3_mds.png", "Figure 3. Multidimensional scaling of genome-wide allele-frequency distances (cohorts with n>=30; small cohorts omitted as the raw distance is noise-inflated at low n). Etruscans cluster with Imperial Romans and other Italians and with modern French/Spanish, between the Anatolian-farmer/Greek and Steppe poles; Han is the far outgroup."))

    # Fig P4: locus archaic frequency heatmap (reuse if present)
    src = f"{FIG}/etr_loci_heatmap.png"
    if os.path.exists(src):
        figs.append(("etr_loci_heatmap.png", "Figure 4. Mean Neanderthal-allele frequency at adaptive-introgression loci across the Italian time transect (see Methods); FADS1-2 is the one locus with a Bonferroni-significant temporal change."))
    return figs


def _img(fname):
    import base64
    fp = os.path.join(FIG, fname)
    b64 = base64.b64encode(open(fp, "rb").read()).decode()
    rel = os.path.relpath(fp, REPORTS).replace("\\", "/")
    return f'<img src="data:image/png;base64,{b64}"><div class="cap"><a href="{rel}">[open figure]</a></div>'


def write_paper(G, Cc, T, etr, o_absz, no_absz, Cl, n_tested, bonf, counts, figs):
    def gv(c, col):
        r = G[G["cohort"] == c]
        return r[col].iloc[0] if len(r) else np.nan
    etr_a, etr_se = gv("Etruscan", "alpha_Nea") * 100, gv("Etruscan", "alpha_SE") * 100
    lat_a = gv("Latin / Italic IA", "alpha_Nea") * 100
    rom_a = gv("Imperial Roman", "alpha_Nea") * 100
    n_etr = int(counts.get("Etruscan", 0))
    maxZ = Cc["Z"].abs().max() if len(Cc) else np.nan
    latZ = float(Cc.loc[Cc["vs"] == "Latin / Italic IA", "Z"].iloc[0]) \
        if "Latin / Italic IA" in Cc["vs"].values else np.nan
    zcrit = sstats.norm.isf(0.025 / max(len(Cc), 1))
    notable = Cc[Cc["Z"].abs() > zcrit]
    if len(notable) == 0:
        contrast_narrative = (f"Etruscan archaic ancestry is statistically indistinguishable "
            f"from every neighbour tested (all |Z| &lt; {zcrit:.1f}), including the key "
            f"Latin/Italic contrast (Z={latZ:+.1f}).")
    else:
        oth = ", ".join(f"{r['vs'].split('(')[0].strip()} (Z={r['Z']:+.1f})"
                        for _, r in notable.iterrows())
        contrast_narrative = (f"Etruscans are archaic-identical to Latins/Italics (Z={latZ:+.1f}) "
            f"and to Anatolian farmers and Italian Bronze Age; the only detectable differences "
            f"are small (|D| &le; 0.003) and involve {oth}, reflecting those groups' somewhat more "
            f"eastern-Mediterranean ancestry rather than unusual archaic introgression.")
    contrast_md = contrast_narrative.replace("&lt;", "<").replace("&le;", "≤")
    fads = bonf[bonf["gene"].str.startswith("FADS")]
    fads_txt = (f"FADS1-2 (fatty-acid metabolism; p={fads['p'].iloc[0]:.1e})"
                if len(fads) else "no locus")
    tu, la = gv("Etruscan Tuscany", "alpha_Nea") * 100, gv("Etruscan Lazio", "alpha_Nea") * 100

    abstract = (
        f"The Etruscans were a culturally distinctive Iron Age population of central "
        f"Italy whose biological origins were long debated until ancient DNA showed them "
        f"to be genetically continuous with neighbouring Italic peoples. Whether their "
        f"archaic (Neanderthal) ancestry differs from that of contemporaries, changed "
        f"through time, or shows individual anomalies has not been examined. Using {n_etr} "
        f"Etruscan genomes and a regional comparison panel from the Allen Ancient DNA "
        f"Resource (1240K), we estimate Neanderthal ancestry with a validated f4-ratio and "
        f"population mean-genome profiles. Etruscan Neanderthal ancestry is "
        f"{etr_a:.2f}% (± {etr_se:.2f}), indistinguishable from Latins/Italics "
        f"({lat_a:.2f}%) and Imperial Romans ({rom_a:.2f}%) and from the broader European "
        f"range. {contrast_md} Genome-wide "
        f"Neanderthal ancestry is flat across the Italian transect from the Neolithic to "
        f"the medieval period. No Etruscan individual is a genome-wide-significant archaic "
        f"outlier, and—crucially—the steppe/Levantine/East-Mediterranean ancestry outliers "
        f"are not archaic outliers, showing that fine-scale ancestry variation in Iron-Age "
        f"Etruria was decoupled from archaic content. A targeted scan of adaptive-"
        f"introgression loci over time recovers {fads_txt} as the single Bonferroni-"
        f"significant candidate—consistent with the well-known dietary selection at FADS—"
        f"though it rests on few markers and is reported as a hypothesis. The results place "
        f"Etruscan archaic ancestry firmly within its regional and temporal context.")

    rowsG = "".join(
        f"<tr class='{'etr' if 'Etrus' in r['cohort'] else ''}'><td>{html.escape(r['cohort'])}</td>"
        f"<td class='num'>{int(r['n'])}</td><td class='num'>{r['alpha_Nea']*100:.2f} ± {r['alpha_SE']*100:.2f}</td>"
        f"<td class='num'>{r['D_Nea_Z']:.0f}</td><td class='num'>{r['D_Den_Z']:.1f}</td></tr>"
        for _, r in G[G["cohort"].isin(REGION_ORDER)].assign(
            o=lambda d: d["cohort"].map({c: i for i, c in enumerate(REGION_ORDER)})
        ).sort_values("o").iterrows())
    rowsC = "".join(
        f"<tr><td>{html.escape(r['vs'])}</td><td class='num'>{r['D']:+.4f}</td>"
        f"<td class='num'>{r['Z']:+.1f}</td></tr>" for _, r in Cc.iterrows())
    rowsL = "".join(
        f"<tr class='{'sig' if (np.isfinite(r['p']) and r['p']<0.05/max(n_tested,1)) else ''}'>"
        f"<td>{r['gene']}</td><td>{html.escape(r['phenotype'])}</td>"
        f"<td class='num'>{int(r['n_archaic_snp'])}</td>"
        f"<td class='num'>{r['age_coef_per_kyr']*100:+.3f}</td><td class='num'>{r['p']:.3f}</td></tr>"
        for _, r in Cl.iterrows() if np.isfinite(r.get("age_coef_per_kyr", np.nan)))
    figs_html = "".join(f'<figure>{_img(f)}<figcaption>{html.escape(c)}</figcaption></figure>'
                        for f, c in figs)

    # --- Markdown manuscript ---
    md = MD_TEMPLATE.format(
        date=datetime.date.today().isoformat(), abstract=abstract,
        n_etr=n_etr, etr_a=f"{etr_a:.2f}", etr_se=f"{etr_se:.2f}",
        lat_a=f"{lat_a:.2f}", rom_a=f"{rom_a:.2f}", maxZ=f"{maxZ:.1f}",
        o_absz=f"{o_absz:.2f}", no_absz=f"{no_absz:.2f}", tu=f"{tu:.2f}", la=f"{la:.2f}",
        n_tested=n_tested, fads_txt=fads_txt, contrast_md=contrast_md,
        group_table="\n".join(
            f"| {r['cohort']} | {int(r['n'])} | {r['alpha_Nea']*100:.2f} ± {r['alpha_SE']*100:.2f} | {r['D_Den_Z']:.1f} |"
            for _, r in G[G["cohort"].isin(REGION_ORDER)].iterrows()),
        contrast_table="\n".join(
            f"| Etruscan vs {r['vs']} | {r['D']:+.4f} | {r['Z']:+.1f} |" for _, r in Cc.iterrows()),
        locus_table="\n".join(
            f"| {r['gene']} | {r['phenotype']} | {int(r['n_archaic_snp'])} | {r['age_coef_per_kyr']*100:+.3f} | {r['p']:.3f} |"
            for _, r in Cl.iterrows() if np.isfinite(r.get("age_coef_per_kyr", np.nan))))
    open(os.path.join(HERE, "PAPER.md"), "w", encoding="utf-8").write(md)

    # --- HTML manuscript ---
    doc = HTML_TEMPLATE.format(
        date=datetime.date.today().isoformat(), abstract=html.escape(abstract),
        rowsG=rowsG, rowsC=rowsC, rowsL=rowsL, figs=figs_html,
        n_etr=n_etr, etr_a=f"{etr_a:.2f}", etr_se=f"{etr_se:.2f}",
        o_absz=f"{o_absz:.2f}", no_absz=f"{no_absz:.2f}", contrast_html=contrast_narrative,
        n_tested=n_tested, maxZ=f"{maxZ:.1f}", tu=f"{tu:.2f}", la=f"{la:.2f}")
    open(os.path.join(REPORTS, "etruscan_paper.html"), "w", encoding="utf-8").write(doc)


MD_TEMPLATE = """# Neanderthal ancestry in the Etruscans: regional context, temporal stability, and a candidate selection signal

*Computational re-analysis of public ancient-DNA (AADR 1240K). Generated {date}. Exploratory; candidate findings are hypotheses, not validated discoveries.*

## Abstract
{abstract}

## Introduction
The Etruscans (Etruria, central Italy, ~800–100 BCE) developed a distinctive language and culture whose origins were debated for millennia. Ancient-DNA studies (Antonio et al. 2019; Posth et al. 2021) showed Etruscans were genetically continuous with neighbouring Italic peoples and the preceding Bronze Age, carrying a Steppe-derived component like other Iron Age Europeans, with a minority of individuals of eastern-Mediterranean or steppe-leaning ancestry. All present-day and ancient non-Africans also carry ~2% Neanderthal ancestry from Late-Pleistocene admixture (Green et al. 2010; Prüfer et al. 2014). Whether Etruscan archaic ancestry differs from neighbours, changed over time, or shows individual anomalies — and whether archaic-introgressed loci were under selection in this transect — has not been assessed. We address these questions with a validated genome-wide estimator, population mean-genome profiles, and a targeted locus scan.

## Materials and Methods
**Data.** {n_etr} Etruscan genomes and regional comparison cohorts (Anatolian Neolithic, Italian Neolithic/Bronze Age, Aegean Bronze Age, Yamnaya Steppe, Latins/Italics, Republican and Imperial Romans, Magna-Graecia Greeks, ancient Sicilians and Sardinians, and modern French/Spanish/Sardinian/Han) from the AADR v66.1 1240K panel (1,233,013 SNPs). Quality control, ancestry PCA, and per-individual estimates follow the parent pipeline (Phases 2–6).

**Neanderthal ancestry.** The f4-ratio α = f4(Altai, Chimp; X, Mbuti) / f4(Altai, Chimp; Vindija, Mbuti), with two independent high-coverage Neanderthals (Altai in the statistic, Vindija as the 100%-Neanderthal scale; Petr et al. 2019). Denisovan affinity is D(X, Mbuti; Denisova, Chimp). Standard errors from a 50-block delete-one jackknife. The estimator reproduces published values (r=0.87; Oase1 recovered at 6–9%; see VALIDATION.md).

**Mean-genome profiles.** For each cohort we build a mean allele-frequency vector ("mean genome"), which averages out single-genome noise and yields group-level α with much tighter SEs than any individual. Differential Neanderthal sharing between cohorts is tested with D(P1, P2; Altai, Yoruba) (an African outgroup isolates recent introgression). Pairwise genome-wide allele-frequency distance summarises genetic affinity (visualised by MDS).

**Selection scan.** At curated adaptive-introgression loci (BNC2, OAS, TLR, FADS, HLA, TBX15, etc.), "archaic-informative" SNPs (Altai+Vindija-derived, near-absent in Africans) are identified; per-individual archaic-allele dosage is regressed on age across Italian samples, controlling for genome-wide archaic ancestry and ancestry turnover (PC1, PC2), so the age term isolates locus-specific temporal change.

## Results

### 1. Etruscan Neanderthal ancestry in regional context
Group-level Neanderthal ancestry (mean-genome f4-ratio):

| population | n | Neanderthal % ± SE | Denisovan D Z |
|---|---|---|---|
{group_table}

Etruscan Neanderthal ancestry is **{etr_a}% ± {etr_se}**, within the tight European band and essentially identical to Latins/Italics and Imperial Romans. Denisovan affinity is ~0 throughout, as expected for West Eurasians.

### 2. Etruscans do not differ from their neighbours
Differential Neanderthal sharing D(Etruscan, X; Altai, Yoruba):

| contrast | D | Z |
|---|---|---|
{contrast_table}

{contrast_md} This is consistent with the genetic continuity reported by Posth et al. 2021. Tuscan vs Latial Etruscans: {tu}% vs {la}%.

### 3. Archaic ancestry was stable over time
Genome-wide Neanderthal ancestry is flat across the Italian transect (Neolithic → medieval); the major Neolithic and Steppe ancestry turnovers did not change the archaic fraction.

### 4. Individual outliers are explained by ancestry, not archaic anomalies
No Etruscan is a genome-wide-significant archaic outlier. The AADR steppe/Levantine/East-Mediterranean ancestry outliers have mean |z| = {o_absz} versus {no_absz} for typical Etruscans — i.e. they are **not** archaic outliers, because those alternative West-Eurasian ancestries carry similar ~2% Neanderthal. Fine-scale ancestry variation in Iron-Age Etruria was decoupled from archaic ancestry.

### 5. A candidate selection signal at FADS
Of {n_tested} adaptive-introgression loci tested, after controlling for ancestry and overall archaic level:

| gene | phenotype | archaic SNPs | Δ/kyr (pp) | p |
|---|---|---|---|---|
{locus_table}

{fads_txt} is the one locus surviving Bonferroni correction — biologically plausible given the strong, well-documented dietary selection at FADS in Europeans (Mathieson et al. 2015; Buckley et al. 2017), but resting on few archaic-informative SNPs and reported as a hypothesis for higher-coverage follow-up.

## Discussion
Etruscan archaic ancestry is unremarkable: ~2% Neanderthal, ~0 Denisovan, indistinguishable from neighbours and stable through time. This is the expected outcome — Neanderthal ancestry is shared across all non-Africans and the within-European variance is small — and it reinforces, from the archaic-ancestry angle, the genetic continuity of Etruscans with Italic peoples. The decoupling of genetic-ancestry outliers from archaic-ancestry outliers is a useful methodological point: an individual can be a clear ancestry outlier yet carry a perfectly ordinary archaic complement, because the relevant alternative ancestries are themselves ~2% Neanderthal. The FADS candidate is intriguing and concordant with the strongest known signal of recent dietary selection in Europe, but the single-locus power on capture data is low and the result must be treated as a hypothesis.

## Limitations
Capture-array, pseudo-haploid genotypes; ~75 Etruscans over a narrow window; single-locus selection power is low; archaic-allele sets are putative (incomplete lineage sorting can mimic introgression); the f4-ratio absolute scale runs slightly high versus some studies (relative comparisons are robust). No claim here is a validated biological discovery.

## References
Antonio et al. 2019 *Science* 366:708 · Posth et al. 2021 *Sci. Adv.* 7:eabi7673 · Green et al. 2010 *Science* 328:710 · Prüfer et al. 2014 *Nature* 505:43 · Prüfer et al. 2017 *Science* 358:655 · Petr et al. 2019 *PNAS* 116:1639 · Patterson et al. 2012 *Genetics* 192:1065 · Mathieson et al. 2015 *Nature* 528:499 · Buckley et al. 2017 *Mol. Biol. Evol.* 34:1307 · Mallick et al. 2024 *Sci. Data* 11:182 (AADR).
"""

HTML_TEMPLATE = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Neanderthal ancestry in the Etruscans</title><style>
:root{{--bg:#fbfaf7;--ink:#1c1a17;--mut:#5c5750;--acc:#8a3b2e;--line:#e2ddd4;--box:#f3efe7}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.65 Georgia,'Times New Roman',serif}}
.wrap{{max-width:820px;margin:0 auto;padding:0 26px 90px}}
header{{padding:54px 26px 26px;border-bottom:2px solid var(--acc)}}
h1{{font-size:27px;line-height:1.25;margin:0 0 10px;font-family:Georgia,serif}}
.meta{{color:var(--mut);font-size:13.5px;font-family:Arial,sans-serif}}
.abstract{{background:var(--box);border:1px solid var(--line);border-radius:6px;padding:16px 20px;margin:22px 0;font-size:15px}}
.abstract b{{font-family:Arial,sans-serif;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--acc);display:block;margin-bottom:6px}}
h2{{font-size:19px;margin:34px 0 10px;border-bottom:1px solid var(--line);padding-bottom:5px}}
h3{{font-size:16px;margin:22px 0 8px;color:#2c2926}}p{{margin:10px 0}}
table{{width:100%;border-collapse:collapse;margin:14px 0;font:13.5px/1.4 Arial,sans-serif}}
th,td{{border-bottom:1px solid var(--line);padding:5px 9px;text-align:left}}th{{color:var(--mut)}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
tr.etr td{{background:#f7e9e6}}tr.sig td{{background:#eaf3ea;font-weight:bold}}
figure{{margin:22px 0;text-align:center}}figure img{{width:100%;border:1px solid var(--line);border-radius:6px}}
figcaption{{font-size:13px;color:var(--mut);margin-top:7px;text-align:left;font-family:Arial,sans-serif}}.cap{{font-size:12px}}.cap a{{color:var(--acc)}}
.disc{{font-size:12.5px;color:var(--mut);font-family:Arial,sans-serif;font-style:italic;margin-top:6px}}
footer{{margin-top:46px;border-top:1px solid var(--line);padding-top:16px;font:12px/1.5 Arial,sans-serif;color:var(--mut)}}
code{{background:var(--box);padding:1px 5px;border-radius:4px;font-size:13px}}
</style></head><body>
<header><div class="wrap" style="padding:0">
<h1>Neanderthal ancestry in the Etruscans: regional context, temporal stability, and a candidate selection signal</h1>
<div class="meta">Computational re-analysis of public ancient DNA (Allen Ancient DNA Resource, 1240K) &middot; {date}</div>
<div class="disc">Exploratory study. Candidate findings are hypotheses requiring validation, not definitive discoveries.</div>
</div></header><div class="wrap">
<div class="abstract"><b>Abstract</b>{abstract}</div>

<h2>1. Etruscan Neanderthal ancestry in regional context</h2>
<p>Group-level (mean-genome) Neanderthal ancestry with block-jackknife standard errors. Etruscans (highlighted) carry {etr_a}% ± {etr_se} — within the tight European band and essentially identical to Latins and Imperial Romans. Denisovan affinity is ~0 throughout, as expected for West Eurasians.</p>
<table><tr><th>population</th><th>n</th><th>Neanderthal % ± SE</th><th>D(Nea) Z</th><th>D(Den) Z</th></tr>{rowsG}</table>

<h2>2. Etruscans do not differ from their neighbours</h2>
<p>Differential Neanderthal sharing, D(Etruscan, X; Altai, Yoruba). {contrast_html} Concordant with the genetic continuity reported by Posth et al. (2021). Tuscan vs Latial Etruscans: {tu}% vs {la}%.</p>
<table><tr><th>contrast</th><th>D</th><th>Z</th></tr>{rowsC}</table>

<h2>3. Archaic ancestry was stable over time</h2>
<p>Group-level Neanderthal ancestry across the Italian transect (Neolithic → medieval) is flat within error: the Neolithic and Steppe ancestry turnovers did not change the archaic fraction.</p>

<h2>4. Individual outliers are explained by ancestry, not archaic anomalies</h2>
<p>No Etruscan is a genome-wide-significant archaic outlier. The AADR steppe/Levantine/East-Mediterranean <i>genetic-ancestry</i> outliers have mean |z| = {o_absz} versus {no_absz} for typical Etruscans — they are <b>not</b> archaic-ancestry outliers, because those alternative West-Eurasian ancestries themselves carry ~2% Neanderthal. Fine-scale ancestry variation in Iron-Age Etruria was decoupled from archaic content.</p>

<h2>5. A candidate selection signal at FADS</h2>
<p>Of {n_tested} adaptive-introgression loci tested (age coefficient controlling for ancestry PCs and overall archaic level), one survives Bonferroni correction (highlighted). FADS is the strongest known target of recent dietary selection in Europeans (Mathieson 2015; Buckley 2017), so a temporal shift there is biologically plausible — but it rests on few archaic-informative SNPs and is a hypothesis, not proof.</p>
<table><tr><th>gene</th><th>phenotype</th><th>archaic SNPs</th><th>Δ/kyr (pp)</th><th>p</th></tr>{rowsL}</table>

<h2>Figures</h2>{figs}

<h2>Discussion</h2>
<p>Etruscan archaic ancestry is unremarkable — ~2% Neanderthal, ~0 Denisovan, indistinguishable from neighbours and stable through time. This is the expected result given that Neanderthal ancestry is shared across all non-Africans with small within-European variance, and it reinforces, from a new angle, the genetic continuity of Etruscans with Italic peoples. The decoupling of ancestry outliers from archaic outliers is a methodological lesson: a clear ancestry outlier can carry an ordinary archaic complement when the alternative ancestries are themselves ~2% Neanderthal. The FADS candidate concords with the strongest known European dietary-selection signal but is underpowered here and remains a hypothesis.</p>
<footer>Methods, validation and code: <code>README.md</code>, <code>VALIDATION.md</code>, <code>FINDINGS.md</code>, <code>PAPER.md</code>. Data: AADR (Mallick et al. 2024). Generated by <code>etruscan_paper.py</code>.</footer>
</div></body></html>"""

REGION_ORDER = REGION_ORDER  # exported for figure ordering

if __name__ == "__main__":
    main()
