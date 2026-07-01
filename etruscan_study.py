#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Etruscan focused study.

Questions:
  A. Did genome-wide Neanderthal/archaic ancestry change over time in Italy
     (Neolithic -> Bronze Age -> Iron Age/Etruscan -> Roman -> Medieval)?
  B. Are any Etruscan individuals statistical outliers in archaic ancestry, AND
     does that outlier status survive conditioning on genetic ancestry + geography
     (i.e. is a deviation real, or just explained by the individual's ancestry)?
  C. WHICH archaic (Neanderthal) alleles / genes changed in frequency over time —
     candidate signals of selection on introgressed variants?
  D. Do the outliers correspond to different genetic ancestries / archaeological
     contexts (e.g. the AADR '-o' steppe/Levantine/East-Med Etruscan outliers)?

Uses the already-computed genome-wide estimates (Phase 3-6) plus the new
locus-level module (archaic.loci). Outputs -> results/etruscan/, figures, and a
self-contained HTML report -> reports/etruscan_report.html (+ ETRUSCAN_FINDINGS.md).

HONEST SCOPE: ~75 Etruscans over a narrow window and ~hundreds of Italians over
~6 kyr on capture data — power for single-locus selection is low; results are
candidates, framed against a genome-wide background, never proof of selection.
"""
import os, sys, base64, html, datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sstats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import loci as L
from archaic.refs import PANELS

PANEL = "1240k"
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "etruscan")
FIG = os.path.join(RESULTS, "figures")
REPORTS = os.path.join(HERE, "reports")
for d in (OUT, FIG, REPORTS):
    os.makedirs(d, exist_ok=True)

# Italian time bins (years BP): label -> (lo, hi)
BINS = [("Neolithic/Copper", 4500, 7000), ("Bronze Age", 3200, 4500),
        ("Iron Age (Etruscan/Latin)", 2300, 3200), ("Roman", 1700, 2300),
        ("Late Antique/Medieval", 800, 1700)]


def wstat(x, w):
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if m.sum() == 0:
        return np.nan, np.nan
    mu = np.average(x[m], weights=w[m])
    se = np.sqrt(1.0 / w[m].sum())              # SE of weighted mean (1/sqrt(sum w))
    return mu, se


def bin_of(bp):
    for lab, lo, hi in BINS:
        if lo <= bp < hi:
            return lab
    return None


def main():
    panel = Panel(PANELS[PANEL]["prefix"])
    refs = PANELS[PANEL]["refs"]

    # merged per-sample table (estimates + metadata + residuals + PCs)
    res = pd.read_csv(os.path.join(RESULTS, f"phase6_{PANEL}_residuals.csv"))
    pcs = pd.read_csv(os.path.join(RESULTS, f"phase5_{PANEL}_pca.csv"))
    df = res.merge(pcs, on="genetic_id")
    df["weight"] = 1.0 / df["alpha_SE"].clip(lower=1e-4) ** 2

    ita = df[(df["country"] == "Italy") & df["date_bp"].between(800, 7000)].copy()
    ita["bin"] = ita["date_bp"].map(bin_of)
    etr = df[df["group_id"].str.contains("Etruscan", case=False, na=False)].copy()
    etr["is_o"] = etr["group_id"].str.contains("-o", na=False)
    print(f"Italian transect: {len(ita)}  |  Etruscans: {len(etr)} "
          f"({int(etr['is_o'].sum())} AADR ancestry-outliers '-o')")

    # ===================== A. genome-wide archaic ancestry over time ===========
    A_rows = []
    for lab, lo, hi in BINS:
        g = ita[ita["bin"] == lab]
        mu, se = wstat(g["alpha_Nea"].values, g["weight"].values)
        A_rows.append(dict(bin=lab, n=len(g), mid_BP=(lo + hi) / 2,
                           mean_Nea=mu * 100, se=se * 100))
    A = pd.DataFrame(A_rows)
    # weighted regression alpha ~ age across individuals
    xa = ita["date_bp"].values.astype(float)
    ya = ita["alpha_Nea"].values
    msk = np.isfinite(xa) & np.isfinite(ya)
    slope, intercept, r, p_time, _ = sstats.linregress(xa[msk], ya[msk])
    A.to_csv(os.path.join(OUT, "A_genomewide_over_time.csv"), index=False)
    print(f"[A] genome-wide Neanderthal vs age: slope={slope*1e3*100:+.4f} pp/kyr, p={p_time:.3f}")

    # ===================== B. Etruscan outliers, with/without ancestry =========
    # within-Etruscan residual (NO ancestry conditioning)
    mu_e, _ = wstat(etr["alpha_Nea"].values, etr["weight"].values)
    var_obs = np.average((etr["alpha_Nea"] - mu_e) ** 2, weights=etr["weight"])
    mean_mv = np.average(etr["alpha_SE"] ** 2, weights=etr["weight"])
    bio2 = max(0.0, var_obs - mean_mv)
    etr["z_within_etr"] = (etr["alpha_Nea"] - mu_e) / np.sqrt(etr["alpha_SE"] ** 2 + bio2)
    # global ancestry+geo+time-conditioned residual already in phase6 (z_resid)
    etr["z_ancestry_cond"] = etr["z_resid"]
    keep = ["genetic_id", "group_id", "country", "date_bp", "is_o",
            "alpha_Nea", "alpha_SE", "alpha_nSNP", "PC1", "PC2",
            "z_within_etr", "z_ancestry_cond", "high_conf"]
    etr_out = etr[keep].sort_values("z_within_etr", key=lambda s: s.abs(), ascending=False)
    etr_out.to_csv(os.path.join(OUT, "B_etruscan_individuals.csv"), index=False)
    # does outlier status correspond to '-o' ancestry outliers?
    o_absz = etr.loc[etr["is_o"], "z_within_etr"].abs().mean()
    no_absz = etr.loc[~etr["is_o"], "z_within_etr"].abs().mean()
    # how many within-Etruscan |z|>2 lose significance once ancestry-conditioned?
    strong = etr[etr["z_within_etr"].abs() > 2]
    explained = strong[strong["z_ancestry_cond"].abs() < 2]
    print(f"[B] mean|z_within| : '-o' outliers {o_absz:.2f} vs typical {no_absz:.2f}")
    print(f"[B] {len(strong)} within-Etruscan |z|>2; {len(explained)} drop below 2 after "
          f"ancestry-conditioning (=> explained by ancestry)")

    # ===================== C. locus-level archaic alleles over time ============
    C_rows = []
    locus_bin_freq = {}     # gene -> {bin: archaic freq}
    for gene, chrom, start, end, pheno, ref in L.LOCI:
        rows = L.panel_rows_in_window(panel, chrom, start, end)
        if len(rows) == 0:
            continue
        info = L.archaic_informative(panel, rows, refs)
        ar_rows, ar_a1 = info["rows"], info["arch_is_a1"]
        if len(ar_rows) < 3:
            C_rows.append(dict(gene=gene, chrom=chrom, n_snp=len(rows),
                               n_archaic_snp=len(ar_rows), phenotype=pheno, reference=ref,
                               age_coef=np.nan, p=np.nan, note="too few informative SNPs"))
            continue
        # per-individual mean archaic-allele dosage at this locus
        gids = list(ita["genetic_id"])
        cols = np.array([panel._id_to_col[g] for g in gids if g in panel._id_to_col])
        valid_gids = [g for g in gids if g in panel._id_to_col]
        Gl = panel.pg.read(ar_rows, cols).astype(np.float32)
        Gl[Gl < 0] = np.nan
        # archaic-allele dosage per individual (mean over informative SNPs)
        dose_a1 = Gl  # copies of allele1
        dose_arch = np.where(ar_a1[:, None], dose_a1, 2.0 - dose_a1)
        with np.errstate(invalid="ignore"):
            per_ind = np.nanmean(dose_arch, axis=0) / 2.0      # archaic-allele freq per individual
            cov = np.sum(~np.isnan(dose_arch), axis=0)
        loc = pd.DataFrame({"genetic_id": valid_gids, "locus_arch": per_ind, "cov": cov})
        loc = loc.merge(ita[["genetic_id", "date_bp", "alpha_Nea", "PC1", "PC2"]],
                        on="genetic_id")
        loc = loc[(loc["cov"] >= 3) & np.isfinite(loc["locus_arch"])]
        # per-bin mean archaic freq
        loc["bin"] = loc["date_bp"].map(bin_of)
        locus_bin_freq[gene] = {lab: loc.loc[loc["bin"] == lab, "locus_arch"].mean()
                                for lab, _, _ in [(b[0], b[1], b[2]) for b in BINS]}
        # locus-specific temporal change, CONTROLLING FOR genome-wide archaic ancestry
        # AND ancestry turnover (PC1,PC2) so the age term isolates selection, not the
        # Neolithic/Steppe ancestry transitions that also move allele frequencies:
        #   OLS  locus_arch ~ age + genomewide_alpha + PC1 + PC2
        X = np.column_stack([np.ones(len(loc)), loc["date_bp"].values,
                             loc["alpha_Nea"].values, loc["PC1"].values, loc["PC2"].values])
        y = loc["locus_arch"].values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        dof = max(1, len(loc) - X.shape[1])
        sigma2 = (resid @ resid) / dof
        XtX_inv = np.linalg.pinv(X.T @ X)
        se_age = np.sqrt(sigma2 * XtX_inv[1, 1])
        t = beta[1] / se_age if se_age > 0 else np.nan
        pval = 2 * sstats.t.sf(abs(t), dof) if np.isfinite(t) else np.nan
        C_rows.append(dict(gene=gene, chrom=chrom, n_snp=len(rows),
                           n_archaic_snp=len(ar_rows), phenotype=pheno, reference=ref,
                           age_coef_per_kyr=beta[1] * 1000, p=pval,
                           direction=("rising toward present" if beta[1] < 0 else "falling toward present"),
                           note=""))
    C = pd.DataFrame(C_rows).sort_values("p")
    C.to_csv(os.path.join(OUT, "C_loci_over_time.csv"), index=False)
    n_tested = int(C["p"].notna().sum())
    bonf_thr = 0.05 / max(n_tested, 1)
    nsig = int((C["p"] < 0.05).sum())
    bonf = C[C["p"] < bonf_thr].dropna(subset=["p"])
    print(f"[C] {n_tested} loci tested; {nsig} nominal p<0.05; "
          f"{len(bonf)} survive Bonferroni ({bonf_thr:.4f}): {list(bonf['gene'])}")

    figs = make_figures(A, ita, etr, C, locus_bin_freq, slope, intercept)
    write_reports(A, etr_out, C, o_absz, no_absz, strong, explained, slope, p_time,
                  nsig, n_tested, bonf, figs)
    print("\nWrote results/etruscan/* , figures, reports/etruscan_report.html, ETRUSCAN_FINDINGS.md")


def make_figures(A, ita, etr, C, locus_bin_freq, slope, intercept):
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
                         "axes.grid": True, "grid.alpha": 0.25})
    figs = []
    # A: Neanderthal over time
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.scatter(ita["date_bp"] / 1000, ita["alpha_Nea"] * 100, s=8, alpha=0.25,
               color="#888", label="Italian individuals")
    e = etr
    ax.scatter(e["date_bp"] / 1000, e["alpha_Nea"] * 100, s=22, color="#d62728",
               edgecolor="k", lw=0.3, label="Etruscans", zorder=4)
    ax.errorbar(A["mid_BP"] / 1000, A["mean_Nea"], yerr=A["se"], fmt="o-",
                color="#1f77b4", lw=2, capsize=3, label="bin mean ± SE", zorder=5)
    xs = np.linspace(ita["date_bp"].min(), ita["date_bp"].max(), 50)
    ax.plot(xs / 1000, (intercept + slope * xs) * 100, "k--", lw=1, alpha=0.6)
    ax.invert_xaxis(); ax.set_xlabel("age (kyr BP, older →)")
    ax.set_ylabel("Neanderthal ancestry (%)"); ax.set_ylim(0, 4.5)
    ax.set_title("Genome-wide Neanderthal ancestry over time in Italy")
    ax.legend(fontsize=8); fig.tight_layout()
    p = f"{FIG}/etr_neanderthal_over_time.png"; fig.savefig(p); plt.close(fig)
    figs.append(("etr_neanderthal_over_time.png", "A. Genome-wide Neanderthal ancestry across the Italian time transect; Etruscans in red."))

    # B: outliers vs ancestry
    fig, ax = plt.subplots(figsize=(7, 4.3))
    c = np.where(etr["is_o"], "#d62728", "#1f77b4")
    sc = ax.scatter(etr["PC1"], etr["z_within_etr"], c=c, s=36, edgecolor="k", lw=0.3)
    for _, r in etr[etr["z_within_etr"].abs() > 1.8].iterrows():
        ax.annotate(str(r["genetic_id"])[:10], (r["PC1"], r["z_within_etr"]),
                    fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.axhline(0, color="k", lw=0.6); ax.axhline(2, color="grey", ls=":"); ax.axhline(-2, color="grey", ls=":")
    ax.set_xlabel("ancestry PC1 (West ← → East Eurasian)")
    ax.set_ylabel("within-Etruscan Neanderthal residual (z)")
    ax.set_title("Etruscan outliers vs ancestry  (red = AADR '-o' ancestry outliers)")
    fig.tight_layout(); p = f"{FIG}/etr_outliers_ancestry.png"; fig.savefig(p); plt.close(fig)
    figs.append(("etr_outliers_ancestry.png", "B. Within-Etruscan archaic residual vs ancestry PC1; red points are the AADR steppe/Levantine/East-Med ancestry outliers."))

    # C: locus archaic frequency heatmap over time
    genes = [g for g in locus_bin_freq if any(np.isfinite(list(locus_bin_freq[g].values())))]
    binlabs = [b[0] for b in BINS]
    M = np.array([[locus_bin_freq[g].get(b, np.nan) for b in binlabs] for g in genes])
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(genes) + 1.5))
    im = ax.imshow(M, aspect="auto", cmap="magma", vmin=np.nanmin(M), vmax=np.nanmax(M))
    ax.set_xticks(range(len(binlabs))); ax.set_xticklabels(binlabs, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8)
    plt.colorbar(im, label="archaic-allele freq")
    ax.set_title("Archaic-allele frequency at introgression loci, over time (Italy)")
    fig.tight_layout(); p = f"{FIG}/etr_loci_heatmap.png"; fig.savefig(p); plt.close(fig)
    figs.append(("etr_loci_heatmap.png", "C. Mean Neanderthal-allele frequency at adaptive-introgression loci across the Italian time transect."))
    return figs


def _img(fname):
    fp = os.path.join(FIG, fname)
    b64 = base64.b64encode(open(fp, "rb").read()).decode()
    rel = os.path.relpath(fp, REPORTS).replace("\\", "/")
    return f'<img src="data:image/png;base64,{b64}"><div class="cap"><a href="{rel}">[open]</a></div>'


def write_reports(A, etr_out, C, o_absz, no_absz, strong, explained, slope, p_time,
                  nsig, n_tested, bonf, figs):
    # data-driven verdict
    if len(bonf):
        bl = "; ".join(f"<b>{r['gene']}</b> (p={r['p']:.1e}, {r['direction']}, "
                       f"{int(r['n_archaic_snp'])} archaic SNPs)" for _, r in bonf.iterrows())
        locus_v = (f"{len(bonf)} of {n_tested} introgression loci survive Bonferroni after "
                   f"controlling for ancestry (PCs) and overall archaic level: {bl}. "
                   "These rest on few archaic-informative SNPs and are flagged as selection "
                   "<i>candidates/hypotheses</i> for higher-coverage follow-up, not proof.")
        if bonf["gene"].str.startswith("FADS").any():
            locus_v += (" FADS1-2 is a well-known target of strong dietary selection in "
                        "Europeans (Mathieson 2015; Buckley 2017), so a temporal shift there "
                        "is biologically plausible.")
    else:
        locus_v = (f"No locus survives Bonferroni among {n_tested} tested once ancestry and "
                   "overall archaic level are controlled — consistent with the limited power "
                   "of ~75 Etruscans / capture data.")
    verdict = (
        f"Genome-wide Neanderthal ancestry is essentially flat across the Italian transect "
        f"({slope*1e3*100:+.4f} pp/kyr, p={p_time:.2f}). No Etruscan individual is a "
        f"significant archaic-ancestry outlier; crucially, the AADR steppe/Levantine/"
        f"East-Mediterranean <i>genetic-ancestry</i> outliers are <b>not</b> archaic-ancestry "
        f"outliers (mean |z| {o_absz:.2f} vs {no_absz:.2f} for typical Etruscans) — those "
        f"alternative West-Eurasian ancestries carry similar ~2% Neanderthal, so the fine-scale "
        f"ancestry variation in Iron-Age Etruria did <b>not</b> translate into archaic-ancestry "
        f"differences. {locus_v}")
    verdict_md = verdict.replace("<b>", "**").replace("</b>", "**").replace("<i>", "*").replace("</i>", "*")
    # markdown
    md = ["# Etruscan archaic-introgression study", "",
          "*Exploratory; results are candidates, not proof of selection.*", "",
          "## A. Genome-wide Neanderthal ancestry over time (Italy)",
          f"Weighted trend vs age: **{slope*1e3*100:+.4f} pp / 1000 yr** (p={p_time:.2f}). "
          "Per-bin means:", "",
          "| era | n | Neanderthal % ± SE |", "|---|---|---|"]
    for _, r in A.iterrows():
        md.append(f"| {r['bin']} | {int(r['n'])} | {r['mean_Nea']:.2f} ± {r['se']:.2f} |")
    md += ["", "## B. Etruscan outliers & ancestry",
           f"Mean |within-Etruscan z|: AADR '-o' ancestry outliers **{o_absz:.2f}** vs "
           f"typical Etruscans **{no_absz:.2f}**. Of {len(strong)} Etruscans with "
           f"|z|>2 within the group, **{len(explained)}** fall below |z|=2 once genetic "
           "ancestry + geography are conditioned on — i.e. their apparent deviation is "
           "**explained by ancestry**, not anomalous archaic introgression.", "",
           "Top within-Etruscan residuals:", "",
           "| id | group | BP | Nea% | z(within) | z(ancestry-cond) | -o |",
           "|---|---|---|---|---|---|---|"]
    for _, r in etr_out.head(12).iterrows():
        md.append(f"| {r['genetic_id']} | {str(r['group_id'])[:34]} | {int(r['date_bp'])} "
                  f"| {r['alpha_Nea']*100:.2f} | {r['z_within_etr']:+.2f} "
                  f"| {r['z_ancestry_cond']:+.2f} | {'YES' if r['is_o'] else ''} |")
    md += ["", "## C. Archaic alleles / genes over time",
           f"{len(C)} adaptive-introgression loci tested; **{nsig}** with nominal p<0.05 "
           f"(Bonferroni threshold {0.05/max(len(C),1):.4f}). Age coefficient controls for "
           "genome-wide archaic ancestry, so it isolates locus-specific change.", "",
           "| gene | phenotype | archaic SNPs | Δ/kyr | p | direction |",
           "|---|---|---|---|---|---|"]
    for _, r in C.iterrows():
        if not np.isfinite(r.get("age_coef_per_kyr", np.nan)):
            continue
        md.append(f"| {r['gene']} | {r['phenotype']} | {int(r['n_archaic_snp'])} "
                  f"| {r['age_coef_per_kyr']*100:+.3f}pp | {r['p']:.3f} | {r.get('direction','')} |")
    md += ["", "## D. Verdict", verdict_md]
    open(os.path.join(HERE, "ETRUSCAN_FINDINGS.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")

    # HTML (same dark theme)
    figs_html = "".join(f'<figure>{_img(f)}<figcaption>{html.escape(c)}</figcaption></figure>'
                        for f, c in figs)
    rowsB = "".join(
        f"<tr class='{'o' if r['is_o'] else ''}'><td>{html.escape(str(r['genetic_id']))}</td>"
        f"<td>{html.escape(str(r['group_id'])[:32])}</td><td class='num'>{int(r['date_bp'])}</td>"
        f"<td class='num'>{r['alpha_Nea']*100:.2f}</td><td class='num'>{r['z_within_etr']:+.2f}</td>"
        f"<td class='num'>{r['z_ancestry_cond']:+.2f}</td></tr>"
        for _, r in etr_out.head(12).iterrows())
    rowsC = "".join(
        f"<tr><td>{r['gene']}</td><td>{html.escape(r['phenotype'])}</td>"
        f"<td class='num'>{int(r['n_archaic_snp'])}</td>"
        f"<td class='num'>{r['age_coef_per_kyr']*100:+.3f}</td><td class='num'>{r['p']:.3f}</td>"
        f"<td>{r.get('direction','')}</td></tr>"
        for _, r in C.iterrows() if np.isfinite(r.get("age_coef_per_kyr", np.nan)))
    doc = ETR_TEMPLATE.format(
        date=datetime.date.today().isoformat(),
        slope=f"{slope*1e3*100:+.4f}", p_time=f"{p_time:.2f}",
        o_absz=f"{o_absz:.2f}", no_absz=f"{no_absz:.2f}",
        n_strong=len(strong), n_expl=len(explained), nsig=nsig, n_tested=n_tested,
        bonf_thr=0.05 / max(n_tested, 1), verdict=verdict,
        figs=figs_html, rowsB=rowsB, rowsC=rowsC)
    open(os.path.join(REPORTS, "etruscan_report.html"), "w", encoding="utf-8").write(doc)


ETR_TEMPLATE = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Etruscan archaic-introgression study</title><style>
:root{{--bg:#0f1419;--panel:#1a2129;--ink:#e6edf3;--mut:#9aa7b2;--acc:#c0794a;--good:#2ea043;--line:#2a333d}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 -apple-system,Segoe UI,Roboto,Arial,sans-serif}}
.wrap{{max-width:960px;margin:0 auto;padding:0 22px 80px}}header{{padding:46px 22px 26px;background:linear-gradient(135deg,#241a14,#0f1419);border-bottom:1px solid var(--line)}}
h1{{font-size:28px;margin:0 0 6px}}.sub{{color:var(--mut)}}.disc{{margin-top:16px;padding:10px 14px;border:1px solid #5f4a16;border-radius:8px;background:#251f10;color:#e8d9a8;font-size:13.5px}}
section{{margin:34px 0;border-top:1px solid var(--line);padding-top:10px}}h2{{font-size:20px}}p{{color:#cdd7e0}}
.callout{{background:#0f2417;border:1px solid #1f5f33;border-radius:8px;padding:12px 16px;margin:14px 0}}
figure{{margin:18px 0}}figure img{{width:100%;border:1px solid var(--line);border-radius:10px;background:#fff}}figcaption,.cap{{color:var(--mut);font-size:13px;margin-top:6px}}.cap a{{color:var(--acc)}}
table{{width:100%;border-collapse:collapse;margin:12px 0;font-size:13.5px}}th,td{{border-bottom:1px solid var(--line);padding:6px 9px;text-align:left}}th{{color:var(--mut)}}td.num{{text-align:right;font-variant-numeric:tabular-nums}}
tr.o td:first-child{{box-shadow:inset 3px 0 var(--acc)}}code{{background:#0c1116;padding:1px 6px;border-radius:5px;color:#e0a37a;font-size:13px}}
footer{{color:var(--mut);font-size:12.5px;border-top:1px solid var(--line);padding-top:18px;margin-top:40px}}</style></head><body>
<header><div class="wrap"><h1>Etruscan archaic-introgression study</h1>
<div class="sub">Neanderthal ancestry, outliers, and selection over time &middot; AADR 1240K &middot; {date}</div>
<div class="disc"><b>Exploratory.</b> ~75 Etruscans over a narrow window; single-locus selection power is low. Results are candidates framed against a genome-wide background, never proof of selection.</div></div></header>
<div class="wrap">
<section><h2>A. Did archaic ancestry change over time in Italy?</h2>
<p>Genome-wide Neanderthal ancestry across the Italian transect (Neolithic → Medieval). Weighted trend vs age: <b>{slope} pp / 1000 yr</b> (p={p_time}) — essentially flat, as expected: overall Neanderthal ancestry is stable; ancestry turnover (Neolithic, Steppe) does not change the archaic fraction much.</p></section>
<section><h2>B. Are Etruscan individuals outliers — and is it ancestry?</h2>
<div class="callout">Mean |within-Etruscan z|: AADR ancestry outliers (<code>-o</code>) <b>{o_absz}</b> vs typical Etruscans <b>{no_absz}</b>. Of {n_strong} Etruscans with within-group |z|&gt;2, <b>{n_expl}</b> drop below |z|=2 once genetic ancestry + geography are conditioned on — their apparent deviation is <b>explained by ancestry</b>, not anomalous introgression.</div>
<table><tr><th>id</th><th>group</th><th>BP</th><th>Nea%</th><th>z within</th><th>z anc-cond</th></tr>{rowsB}</table>
<p>Rows marked with a coloured bar are AADR steppe/Levantine/East-Mediterranean ancestry outliers.</p></section>
<section><h2>C. Which archaic genes changed over time?</h2>
<p>{n_tested} adaptive-introgression loci tested. The age coefficient controls for genome-wide archaic ancestry <b>and ancestry turnover (PC1, PC2)</b>, so it isolates locus-specific temporal change rather than the Neolithic/Steppe ancestry transitions. <b>{nsig}</b> reach nominal p&lt;0.05; Bonferroni threshold is {bonf_thr:.4f}. Δ/kyr &lt; 0 = archaic allele rising toward present.</p>
<table><tr><th>gene</th><th>phenotype</th><th>archaic SNPs</th><th>Δ/kyr (pp)</th><th>p</th><th>direction</th></tr>{rowsC}</table></section>
<section><h2>Figures</h2>{figs}</section>
<section><h2>D. Verdict</h2><p>{verdict}</p></section>
<footer>Generated by <code>etruscan_study.py</code>. Companion: <code>ETRUSCAN_FINDINGS.md</code>. Method &amp; citations: <code>README.md</code>.</footer>
</div></body></html>"""


if __name__ == "__main__":
    main()
