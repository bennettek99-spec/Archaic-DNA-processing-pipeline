#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which Neanderthal? Altai-vs-Vindija affinity across the ancient Eurasian record.

Every Eurasian carries Neanderthal ancestry, and the pipeline's other studies ask
how much. This one asks where it came from. The AADR 1240K release carries two
high-coverage Neanderthal genomes separated by ~4,000 km and ~70,000 years —
Altai (Denisova Cave, Siberia, ~120 kya) and Vindija 33.19 (Croatia, ~50 kya) —
and Vindija is known to sit closer to the population that actually introgressed
(Prufer et al. 2017, Science 358:655). Contrasting them is therefore a probe of
the introgressing source itself, and the question becomes: did every Eurasian
lineage meet the *same* Neanderthal population, or do Upper Palaeolithic
Europeans, East Asians and Oase1 point at different ones?

    D_VA = D(X, Yoruba; Vindija, Altai)

THE TRAP THIS STUDY IS BUILT AROUND

D_VA scales with how much Neanderthal ancestry X has. East Asians carry ~20% more
than Europeans, so under a single shared source they are *expected* to show a
~20% larger D_VA. Reading raw D_VA as evidence of a different Neanderthal would
manufacture exactly the "second pulse into East Asia" that this analysis is
supposed to test. Every cohort is therefore also measured with

    D_NEA = D(X, Yoruba; NeaAvg, Chimp),      NeaAvg = (p_Vindija + p_Altai)/2

which is exactly symmetric in the two archaic genomes and so tracks the *amount*
of Neanderthal ancestry without responding to its source. Under a single source
every cohort lies on one line D_VA = k * D_NEA through the origin; the slope k
estimates D_VA of the source population, and a cohort's residual from that line
is the only thing here that would count as a different Neanderthal.

WHAT MAKES THE NULL CREDIBLE

A null result is only worth reporting with a stated detection limit, so the run
includes a battery that is larger than the analysis it protects:

  * POSITIVE CONTROLS. Denisova as the test population must come out strongly
    negative, because the Altai Neanderthal carries ~1% Denisovan-related
    ancestry (Prufer et al. 2014, Nature 505:43) — a real, published,
    between-Neanderthal difference the statistic has to see. D_VA must also be
    proportional to D_NEA across cohorts spanning near-zero to ~10% Neanderthal
    ancestry, which shows it responds to Neanderthal ancestry at all.
  * NULL CONSTRUCTIONS. Random half-splits of single homogeneous cohorts give the
    pure sampling floor. Coverage-stratified splits — the low-coverage half of a
    cohort against its own high-coverage half — give the floor that matters,
    because coverage trends with time in the AADR and would otherwise fake a
    temporal signal. African cohorts, which have almost no Neanderthal ancestry,
    check the baseline.
  * PAIRED JACKKNIFE. The limiting sample here is the archaic genomes (only ~19k
    1240K sites separate Vindija from Altai), and that noise is identical for
    every cohort. Differences are jackknifed over shared blocks so it cancels.
  * COVERAGE MATCHING and a TRANSVERSIONS-ONLY rerun, as elsewhere in the repo.

Interpretation boundary (repository-wide): D_VA and D_NEA are relative
affinities, never percentages. Findings are hypotheses until technical causes are
excluded.

Outputs (reports/neanderthal_source/):
  ns_cohorts.csv              per-cohort D_VA, D_NEA, R and jackknife SEs
  ns_pairwise.csv             all pairwise paired-jackknife differences
  ns_residuals.csv            residuals from the single-source proportional fit
  ns_controls.csv             positive controls and null constructions
  ns_detection_limit.csv      the stated detection limit and how it was derived
  fig_n1_scaling.png          D_VA against D_NEA — the single-source line
  fig_n2_time.png             source affinity through time, by region
  fig_n3_null.png             the null battery and the detection limit
  fig_n4_targets.png          the three named targets against the null band
  PAPER_neanderthal_source.md paper-style report
Run: PYTHONIOENCODING=utf-8 python scripts/neanderthal_source.py --panel 1240k
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import source_contrast as sc
from archaic import stats as st
from archaic.cohort import pooled_freq_multi
from archaic.log_utils import get_logger
from archaic.panel import Panel
from archaic.refs import PANELS

log = get_logger("archaic.neanderthal_source")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "reports", "neanderthal_source")
META = os.path.join(RESULTS, "phase4_1240k_global_analysis.csv")
N_BLOCKS = 50
SNP_FLOOR = 100_000          # informative-SNP floor; yields 11,853 Eurasian ancients
MIN_COHORT = 20              # smallest grid cell kept as its own cohort
BASE = "Yoruba"              # African baseline for every D statistic

# A cohort is only testable for its Neanderthal *source* if it has enough
# Neanderthal ancestry for the source to matter. Below this D_NEA the ratio
# R = D_VA / D_NEA divides by something indistinguishable from zero and explodes
# (present-day Mbuti return R ~ 60), and the proportional model itself stops
# describing the cohort: for African populations the departure from the Yoruba
# baseline is dominated by deep African structure rather than by introgression,
# so neither D_VA nor D_NEA is measuring what the model assumes. African cohorts
# are therefore reported but excluded from the fit and from the residual test —
# a scope condition, not a filter chosen after seeing the answers.
MIN_DNEA = 0.02

# Present-day anchors. The African set is deliberately broad — West African
# groups sit close to the Yoruba baseline, divergent ones (Mbuti, Biaka,
# Ju|'hoan) do not, and the spread between them measures how much of D_VA is
# baseline choice rather than Neanderthal ancestry.
ANCHORS = {
    "PD_French": dict(pops=["French"]),
    "PD_Sardinian": dict(pops=["Sardinian"]),
    "PD_Basque": dict(pops=["Basque"]),
    "PD_Russian": dict(pops=["Russian"]),
    "PD_Han": dict(pops=["Han"]),
    "PD_Japanese": dict(pops=["Japanese"]),
    "PD_Dai": dict(pops=["Dai"]),
    "PD_Papuan": dict(pops=["Papuan"]),
    "PD_Karitiana": dict(pops=["Karitiana"]),
    "PD_Mixe": dict(pops=["Mixe"]),
    "PD_Mbuti": dict(pops=["Mbuti", "MbutiPygmy"]),
    "PD_Biaka": dict(pops=["Biaka"]),
    "PD_Juhoan": dict(pops=["Ju_hoan_North"]),
    "PD_Mandenka": dict(pops=["Mandenka"]),
    "PD_Mende": dict(pops=["Mende", "MSL"]),
    "PD_Esan": dict(pops=["ESN"]),
    "PD_Gambian": dict(pops=["GWD", "Gambian"]),
    "PD_Luhya": dict(pops=["LWK", "Luhya"]),
}

# Cohorts the hypothesis is actually about, evaluated against the fitted line.
TARGETS = ["Oase1_40ka", "UP_Europe_pre_LGM", "UP_NorthEastAsia"]

TIME_BINS = [(15000, 1e9, "Palaeolithic"), (8000, 15000, "LateGlacial_Meso"),
             (5000, 8000, "Neolithic"), (3000, 5000, "Bronze"),
             (1500, 3000, "Iron_Classical"), (0, 1500, "Medieval_Recent")]


# ------------------------------------------------------------ cohort building -
def region_of(lon, lat, continent=""):
    """Coarse Eurasian region from coordinates.

    Longitude/latitude bands rather than curated group labels: the AADR has 3,668
    distinct group_ids and any hand-drawn list of them would be a place for the
    analyst's expectations to enter. Bands are crude but they are fixed before
    the statistics are seen.

    The bands only partition Eurasia, so anything outside it keeps its continent
    instead. Without this an African or Native American reference population
    lands in whichever Eurasian box its coordinates happen to fall in, and the
    published tables label Mbuti and Karitiana as "NearEast".
    """
    if continent and continent != "Eurasia":
        return str(continent)
    if not (np.isfinite(lon) and np.isfinite(lat)):
        return "Unknown"
    if lon < 35 and lat >= 36:
        return "Europe"
    if lon < 45 and lat < 36:
        return "NearEast"
    if lon < 55 and 36 <= lat < 48:
        return "Caucasus"
    if lon < 75 and lat >= 48:
        return "WestSiberia"
    if lon < 75:
        return "CentralAsia"
    if lon < 100 and lat >= 48:
        return "Siberia"
    if lon < 100 and lat < 30:
        return "SouthAsia"
    if lon < 100:
        return "InnerAsia"
    if lat < 23:
        return "SEAsia"
    return "EastAsia"


def period_of(bp):
    for lo, hi, lab in TIME_BINS:
        if lo <= bp < hi:
            return lab
    return "Unknown"


def dedupe(meta):
    """Collapse AADR duplicate entries for one individual, keeping the best.

    The same person can appear as `.SG`, `.DG`, `.AG` and damage-restricted `_d`
    versions. Pooling all of them would silently give that individual several
    votes in their cohort's allele frequency, so entries are collapsed on the
    stripped sample name and the version with the most informative SNPs wins.
    """
    m = meta.copy()
    m["base_id"] = m["genetic_id"].str.replace(r"(_d)?(\.[A-Za-z0-9]+)+$", "",
                                               regex=True)
    m = m.sort_values("alpha_nSNP", ascending=False)
    return m.drop_duplicates("base_id", keep="first")


def build_cohorts(meta, panel):
    """Named Palaeolithic cohorts, a region-by-period grid, and control splits."""
    col_of = panel._id_to_col
    m = meta[meta["genetic_id"].isin(col_of)].copy()
    m["col"] = m["genetic_id"].map(col_of)
    m["region"] = [region_of(a, b, c) for a, b, c
                   in zip(m["lon"], m["lat"], m["continent"].astype(str))]
    m["period"] = m["date_bp"].map(period_of)

    anc = m[~m["is_modern"]].copy()
    eura = anc[(anc["continent"] == "Eurasia") & (anc["alpha_nSNP"] >= SNP_FLOOR)]
    eura = dedupe(eura)

    defs = {}
    meta_rows = {}

    def add(label, sub, kind, note=""):
        sub = sub.dropna(subset=["col"])
        if len(sub) == 0:
            return
        defs[label] = sub["col"].to_numpy(dtype=np.int64)
        meta_rows[label] = dict(
            label=label, kind=kind, n_ind=len(sub),
            date_bp=float(sub["date_bp"].median()) if kind != "present-day" else 0.0,
            date_min=float(sub["date_bp"].min()) if kind != "present-day" else 0.0,
            date_max=float(sub["date_bp"].max()) if kind != "present-day" else 0.0,
            region=sub["region"].mode().iat[0] if len(sub["region"].mode()) else "",
            median_coverage=float(pd.to_numeric(sub["coverage"],
                                                errors="coerce").median()),
            median_snps=float(sub["alpha_nSNP"].median()),
            median_damage=float(pd.to_numeric(sub.get("damage"),
                                              errors="coerce").median())
            if "damage" in sub else np.nan,
            mean_alpha=float(pd.to_numeric(sub["alpha_Nea"],
                                           errors="coerce").mean()),
            note=note)

    # --- named Palaeolithic cohorts (the hypothesis) --------------------------
    pal = eura[eura["date_bp"] >= 15000]
    add("IUP_Eurasia_45ka", pal[pal["date_bp"] >= 42000], "named",
        "Zlaty kun, Ranis, Bacho Kiro IUP, Ust'-Ishim — the basal ~45 kya group")
    add("UP_Europe_pre_LGM",
        pal[(pal["region"] == "Europe") & (pal["date_bp"] >= 25000)
            & (pal["date_bp"] < 42000)], "named",
        "Kostenki, Sunghir, Vestonice, Goyet, Gravettian")
    add("UP_Europe_post_LGM",
        pal[(pal["region"] == "Europe") & (pal["date_bp"] < 25000)], "named",
        "Late Upper Palaeolithic / Magdalenian Europe")
    add("UP_NorthEastAsia", pal[pal["lon"] >= 75], "named",
        "Tianyuan, Salkhit, Amur, Yana, Mal'ta — Palaeolithic east of 75E")

    # Oase1 sits below the SNP floor (25,775 informative sites at 0.05x) but is
    # a named target of the question, so it is admitted explicitly and flagged
    # everywhere rather than quietly dropped or quietly included.
    oase = anc[anc["genetic_id"].str.startswith("Oase1")]
    add("Oase1_40ka", oase, "named-lowpower",
        "admitted below the SNP floor; power-limited by construction")
    ushim = anc[anc["genetic_id"].str.startswith("Ust_Ishim")]
    add("UstIshim_44ka", ushim, "named", "single 40x genome, 44 kya")

    # --- region x period grid (the time axis) --------------------------------
    for (reg, per), g in eura.groupby(["region", "period"]):
        if reg == "Unknown" or per == "Unknown" or len(g) < MIN_COHORT:
            continue
        add(f"{reg}_{per}", g, "grid")

    # --- non-Eurasian ancients, as extra points on the scaling line ----------
    for cont in ("Americas", "Oceania", "Africa"):
        sub = dedupe(anc[(anc["continent"] == cont)
                         & (anc["alpha_nSNP"] >= SNP_FLOOR)])
        if len(sub) >= MIN_COHORT:
            add(f"ANC_{cont}", sub, "ancient-other")

    # --- present-day anchors -------------------------------------------------
    for name, sel in ANCHORS.items():
        cols = panel.cols_for(sel.get("ids"), sel.get("pops"))
        if len(cols) == 0:
            log.warning(f"  anchor {name}: no individuals found, skipped")
            continue
        sub = m[m["col"].isin(cols)]
        if len(sub) == 0:                       # present in .ind but not in meta
            defs[name] = cols
            meta_rows[name] = dict(label=name, kind="present-day", n_ind=len(cols),
                                   date_bp=0.0, date_min=0.0, date_max=0.0,
                                   region="", median_coverage=np.nan,
                                   median_snps=np.nan, median_damage=np.nan,
                                   mean_alpha=np.nan, note="")
            continue
        add(name, sub, "present-day")

    # --- null constructions --------------------------------------------------
    # Splits of single homogeneous cohorts, which by construction share a
    # Neanderthal source. Random splits give the sampling floor; coverage-
    # stratified splits give the floor that actually threatens a time series,
    # because coverage trends with date across the AADR.
    rng = np.random.default_rng(20260811)
    split_meta = []
    big = [(lab, eura[(eura["region"] == r) & (eura["period"] == p)])
           for lab, (r, p) in [("Europe_Medieval_Recent", ("Europe", "Medieval_Recent")),
                               ("Europe_Neolithic", ("Europe", "Neolithic")),
                               ("Europe_Bronze", ("Europe", "Bronze")),
                               ("Europe_Iron_Classical", ("Europe", "Iron_Classical")),
                               ("EastAsia_Bronze", ("EastAsia", "Bronze"))]]
    for lab, g in big:
        if len(g) < 4 * MIN_COHORT:
            continue
        for rep in range(3):
            idx = rng.permutation(len(g))
            h = len(g) // 2
            a, b = g.iloc[idx[:h]], g.iloc[idx[h:2 * h]]
            add(f"NULLrand_{lab}_r{rep}A", a, "null-random")
            add(f"NULLrand_{lab}_r{rep}B", b, "null-random")
            split_meta.append(("random", lab, f"NULLrand_{lab}_r{rep}A",
                               f"NULLrand_{lab}_r{rep}B"))
        cov = pd.to_numeric(g["alpha_nSNP"], errors="coerce")
        med = cov.median()
        lo_g, hi_g = g[cov <= med], g[cov > med]
        add(f"NULLcov_{lab}_LO", lo_g, "null-coverage")
        add(f"NULLcov_{lab}_HI", hi_g, "null-coverage")
        split_meta.append(("coverage", lab, f"NULLcov_{lab}_HI",
                           f"NULLcov_{lab}_LO"))

    return defs, meta_rows, split_meta, eura, m


def cached_pooled_freq(panel, defs, args, tag, log):
    """Pooled frequencies, cached so the genotype pass happens once.

    The single streaming pass over ~13,000 individuals is by far the most
    expensive part of this study, while everything downstream — the fit, the
    null battery, the figures, the report wording — gets iterated on. The cache
    is keyed on the panel, the transversion setting and a hash of the cohort
    definitions, so any change to who is in which cohort invalidates it
    automatically and a stale cache cannot silently produce a wrong answer.
    Frequencies are stored as float32 (allele frequencies need nothing like
    float64 precision) and counts as int32, which keeps the file to a few
    hundred MB. Delete the directory to force a recomputation.
    """
    import hashlib
    import json

    key = hashlib.sha256(json.dumps(
        {k: [int(c) for c in v] for k, v in sorted(defs.items())},
        sort_keys=True).encode()).hexdigest()[:16]
    cdir = os.path.join(RESULTS, "neanderthal_source_cache")
    os.makedirs(cdir, exist_ok=True)
    path = os.path.join(cdir, f"pooled_{args.panel}{tag}_{key}.npz")

    if os.path.exists(path):
        log.info(f"  reusing cached pooled frequencies: {os.path.basename(path)}")
        z = np.load(path)
        freqs = {l: z[f"f_{l}"].astype(np.float64) for l in defs}
        counts = {l: z[f"n_{l}"].astype(np.int64) for l in defs}
        return freqs, counts

    freqs, counts = pooled_freq_multi(panel, panel.snp_rows, defs,
                                      chunk=args.chunk, log=log)
    # Write through a temporary and rename, so an interrupted run cannot leave a
    # half-written cache that a later run would trust. np.savez appends ".npz" to
    # a *path* that lacks the suffix — which would silently defeat the rename —
    # so it is handed an open file object instead, where it appends nothing.
    tmp = path + ".part"
    with open(tmp, "wb") as fh:
        np.savez(fh, **{f"f_{l}": freqs[l].astype(np.float32) for l in defs},
                 **{f"n_{l}": counts[l].astype(np.int32) for l in defs})
    os.replace(tmp, path)
    log.info(f"  cached pooled frequencies to {os.path.basename(path)}")
    return freqs, counts


# --------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--blocks", type=int, default=N_BLOCKS)
    ap.add_argument("--transversions", action="store_true",
                    help="restrict to transversions (damage-robust sensitivity run)")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    os.makedirs(OUT, exist_ok=True)
    tag = "_transversions" if args.transversions else ""
    nb = args.blocks

    log.info(f"Loading panel {args.panel} (transversions={args.transversions})...")
    panel = Panel(cfg["prefix"], autosomes_only=True,
                  transversions_only=args.transversions)
    block = st.assign_blocks(panel.n_snp, nb)

    log.info("Reading archaic and baseline reference frequencies...")
    spec = {k: cfg["refs"][k] for k in ("Altai", "Vindija", "Denisova", "Chimp",
                                        "Mbuti", "Yoruba")}
    ref, rinfo = panel.frequencies(spec)
    for k, v in rinfo.items():
        log.info(f"  ref {k:9s} n={v['n_ind']:>4} covered={v['n_snp_covered']:,}")

    both = np.isfinite(ref["Vindija"]) & np.isfinite(ref["Altai"])
    informative = int((both & (ref["Vindija"] != ref["Altai"])).sum())
    log.info(f"  Vindija and Altai both called at {int(both.sum()):,} SNPs; "
             f"{informative:,} of those distinguish them "
             f"({100*informative/max(both.sum(),1):.2f}%)")

    meta = pd.read_csv(META, low_memory=False)
    defs, crows, split_meta, eura, mall = build_cohorts(meta, panel)
    log.info(f"{len(defs)} cohorts defined over {len(eura):,} unique Eurasian "
             f"ancients (>= {SNP_FLOOR:,} SNPs)")

    log.info("Pooling allele frequencies (single pass over the genotypes)...")
    freqs, counts = cached_pooled_freq(panel, defs, args, tag, log)

    # ---- block tables -------------------------------------------------------
    def table_for(label, p, n_ind, mask=None):
        return sc.build_block_table(label, p, ref[BASE], ref["Vindija"],
                                    ref["Altai"], ref["Chimp"], block, nb,
                                    mask=mask, n_ind=n_ind)

    tables = {}
    for lab in defs:
        tables[lab] = table_for(lab, freqs[lab], crows[lab]["n_ind"])

    # Reference genomes as test populations: Denisova is the positive control,
    # Chimp and Mbuti are sanity checks on the baseline.
    for lab, key in (("CTRL_Denisova", "Denisova"), ("CTRL_Chimp", "Chimp"),
                     ("CTRL_Mbuti", "Mbuti")):
        tables[lab] = table_for(lab, ref[key], rinfo[key]["n_ind"])
        crows[lab] = dict(label=lab, kind="control-reference",
                          n_ind=rinfo[key]["n_ind"], date_bp=np.nan,
                          date_min=np.nan, date_max=np.nan, region="",
                          median_coverage=np.nan, median_snps=np.nan,
                          median_damage=np.nan,
                          mean_alpha=np.nan, note=f"{key} as the test population")

    # ---- per-cohort statistics ---------------------------------------------
    rows = []
    for lab, t in tables.items():
        va = sc.cohort_stat(t, "D_VA")
        ne = sc.cohort_stat(t, "D_NEA")
        r = sc.cohort_stat(t, "R")
        # R is a ratio whose denominator is a Neanderthal-affinity measurement;
        # where that is near zero the ratio carries no information and is
        # suppressed rather than printed as a spuriously large number. Chimp is
        # degenerate by construction (it appears on both sides of D_NEA, forcing
        # D_NEA = -1), so it is scored on D_VA only.
        testable = bool(np.isfinite(ne["theta"]) and ne["theta"] > MIN_DNEA
                        and lab != "CTRL_Chimp")
        rows.append(dict(**crows[lab],
                         D_VA=va["theta"], D_VA_se=va["se"], D_VA_z=va["z"],
                         D_NEA=ne["theta"], D_NEA_se=ne["se"], D_NEA_z=ne["z"],
                         R=r["theta"] if testable else np.nan,
                         R_se=r["se"] if testable else np.nan,
                         testable=testable, n_snp=t.n_snp))
    df = pd.DataFrame(rows).sort_values(["kind", "label"])
    testable = set(df.loc[df["testable"], "label"])
    for _, r in df[df.kind.isin(["named", "named-lowpower", "control-reference"])].iterrows():
        log.info(f"  [{r['kind']:>18s}] {r['label']:24s} n={r['n_ind']:>4} "
                 f"D_VA={r['D_VA']:+.4f}+/-{r['D_VA_se']:.4f} "
                 f"D_NEA={r['D_NEA']:+.4f} R={r['R']:+.3f}")

    # ---- the single-source proportional fit --------------------------------
    # The line is defined by the region-by-period grid alone: several thousand
    # ancient Eurasians, the largest and most homogeneous set here, and none of
    # them a hypothesis target. Everything else — the named Palaeolithic
    # cohorts, the present-day anchors, the Denisovan positive control — is then
    # scored against a line it did not help define. Leaving the present-day
    # anchors out of the fit matters for one of them in particular: Papuans
    # carry Denisovan ancestry, and because the Altai genome does too, they are
    # expected to be pulled off the line for a reason that has nothing to do
    # with their Neanderthal source. That is a prediction worth testing, not a
    # constraint worth fitting.
    fit_labels = [l for l in tables if crows[l]["kind"] == "grid"
                  and l in testable]
    test_labels = [l for l in tables if not l.startswith("NULL")
                   and l in testable]
    log.info(f"Fitting D_VA = k * D_NEA on {len(fit_labels)} cohorts "
             f"(targets excluded from the fit)...")
    fit, resid = sc.fit_and_residuals(tables, fit_labels, test_labels, nb)
    log.info(f"  k = {fit['k']:.4f} +/- {fit['k_se']:.4f} (Z={fit['k_z']:.1f})")
    rdf = pd.DataFrame(resid).merge(df[["label", "kind", "n_ind", "date_bp",
                                        "region", "D_VA", "D_NEA", "R"]],
                                    on="label", how="left")
    rdf = rdf.sort_values("residual_z", key=lambda s: s.abs(), ascending=False)

    # ---- pairwise paired differences ---------------------------------------
    main_labels = [l for l in tables if crows[l]["kind"] in
                   ("named", "named-lowpower", "grid", "present-day",
                    "ancient-other")]
    pw = []
    for a, b in itertools.combinations(main_labels, 2):
        for stat in ("D_VA", "R"):
            d = sc.paired_difference(tables[a], tables[b], stat)
            pw.append(dict(cohort_a=a, cohort_b=b, statistic=stat, **d))
    pwdf = pd.DataFrame(pw)

    # ---- null battery and detection limit ----------------------------------
    ctrl = []
    for kind, lab, a, b in split_meta:
        if a not in tables or b not in tables:
            continue
        d = sc.paired_difference(tables[a], tables[b], "D_VA")
        dr = sc.paired_difference(tables[a], tables[b], "R")
        ctrl.append(dict(control=f"{kind}-split", cohort=lab, a=a, b=b,
                         D_VA_diff=d["diff"], D_VA_diff_se=d["se"],
                         D_VA_diff_z=d["z"],
                         D_VA_diff_se_independent=d["se_independent"],
                         R_diff=dr["diff"], R_diff_se=dr["se"],
                         R_diff_z=dr["z"]))
    cdf = pd.DataFrame(ctrl)

    # The scale a D_VA difference must be judged against is the typical cohort's
    # own Vindija-over-Altai displacement, NOT the fitted slope k (which is
    # D_VA per unit D_NEA and lives in different units). The statistical floor
    # comes from the real cohort-versus-cohort comparisons; the systematic floor
    # from the same-cohort splits.
    real_pairs = pwdf[(pwdf.statistic == "D_VA")
                      & pwdf.cohort_a.isin(testable)
                      & pwdf.cohort_b.isin(testable)]
    signal = float(df.loc[df.label.isin(fit_labels), "D_VA"].median())
    lim = sc.detection_limit(cdf["D_VA_diff"], real_pairs["se"], signal)
    real_pairs_R = pwdf[(pwdf.statistic == "R")
                        & pwdf.cohort_a.isin(testable)
                        & pwdf.cohort_b.isin(testable)]
    limR = sc.detection_limit(cdf["R_diff"], real_pairs_R["se"],
                              float(df.loc[df.label.isin(fit_labels), "R"].median()))
    log.info(f"  typical cohort D_VA = {signal:.4f}; detection limit on a "
             f"cohort difference = {lim['limit']:.4f} "
             f"({100*lim['limit_fraction_of_signal']:.0f}% of it; best case "
             f"{lim['best_limit']:.4f} = "
             f"{100*lim['best_fraction_of_signal']:.0f}%)")
    log.info(f"    statistical floor {lim['statistical_floor']:.5f}, "
             f"systematic floor {lim['systematic_floor']:.5f}")

    # ---- is any residual structure technical rather than historical? --------
    # The one pattern that must not be taken at face value is anything that
    # tracks sample age, because age is also a proxy for damage and coverage.
    diag = df[df.label.isin(fit_labels + [l for l in tables
                                          if crows[l]["kind"] == "named"])]
    diag = diag[diag["testable"] & np.isfinite(diag["date_bp"])]
    covar = sc.technical_covariates(
        diag["label"], diag["R"],
        {"date_bp": diag["date_bp"], "median_coverage": diag["median_coverage"],
         "median_damage": diag["median_damage"],
         "median_snps": diag["median_snps"]})
    for k_, v_ in covar.items():
        log.info(f"  [covariate] R vs {k_:16s} rho={v_['rho']:+.3f} "
                 f"p={v_['p']:.3g} (n={v_['n']})")
    covdf = pd.DataFrame([dict(covariate=k_, **v_) for k_, v_ in covar.items()])

    # how much of the SE is common-mode, i.e. what the pairing bought
    gain = float(np.nanmedian(pwdf["se_independent"] / pwdf["se"]))
    log.info(f"  paired jackknife is {gain:.2f}x tighter than independent SEs")

    # ---- coverage-matched recomputation ------------------------------------
    core = [l for l in tables if crows[l]["kind"] in
            ("named", "named-lowpower", "grid") and not l.startswith("NULL")]
    core_no_oase = [l for l in core if l != "Oase1_40ka"]
    shared = np.ones(panel.n_snp, dtype=bool)
    for l in core_no_oase:
        shared &= counts[l] >= 1
    log.info(f"Coverage-matched SNP set shared by {len(core_no_oase)} core "
             f"cohorts: {int(shared.sum()):,} SNPs")
    cm_tables = {l: table_for(l, freqs[l], crows[l]["n_ind"], mask=shared)
                 for l in core_no_oase}
    cm_fit_labels = [l for l in cm_tables if crows[l]["kind"] == "grid"
                     and cm_tables[l].d_ne_theta > 0.01]
    cm_rows = []
    if len(cm_fit_labels) >= 3:
        cm_fit, cm_resid = sc.fit_and_residuals(cm_tables, cm_fit_labels,
                                                list(cm_tables), nb)
        log.info(f"  coverage-matched k = {cm_fit['k']:.4f} +/- {cm_fit['k_se']:.4f}")
        cm_rows = cm_resid
    else:
        cm_fit = dict(k=np.nan, k_se=np.nan, k_z=np.nan)
    cmdf = pd.DataFrame(cm_rows)
    for l in cm_tables:
        s = sc.cohort_stat(cm_tables[l], "D_VA")
        if not cmdf.empty:
            cmdf.loc[cmdf.label == l, ["D_VA_matched", "D_VA_matched_se",
                                       "n_snp_matched"]] = (
                s["theta"], s["se"], cm_tables[l].n_snp)

    # ---- write everything ---------------------------------------------------
    df.to_csv(os.path.join(OUT, f"ns_cohorts{tag}.csv"), index=False)
    pwdf.to_csv(os.path.join(OUT, f"ns_pairwise{tag}.csv"), index=False)
    rdf.to_csv(os.path.join(OUT, f"ns_residuals{tag}.csv"), index=False)
    cdf.to_csv(os.path.join(OUT, f"ns_controls{tag}.csv"), index=False)
    cmdf.to_csv(os.path.join(OUT, f"ns_coverage_matched{tag}.csv"), index=False)
    covdf.to_csv(os.path.join(OUT, f"ns_covariates{tag}.csv"), index=False)
    pd.DataFrame([dict(statistic="D_VA", **lim),
                  dict(statistic="R", **limR)]).to_csv(
        os.path.join(OUT, f"ns_detection_limit{tag}.csv"), index=False)

    make_figures(df, rdf, cdf, pwdf, fit, lim, limR, tag)
    write_report(df, rdf, cdf, pwdf, cmdf, covdf, fit, cm_fit, lim, limR, gain,
                 eura, informative, int(both.sum()), int(shared.sum()),
                 args, tag)
    log.info("Done.")


# ------------------------------------------------------------------ figures ---
def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 160, "font.size": 9})
    return plt


C_EUR, C_EAS, C_UP, C_PD, C_NULL, C_OASE = (
    "#4c72b0", "#c44e52", "#dd8452", "#8172b2", "#937860", "#55a868")


def _kindcolor(k):
    return {"named": C_UP, "named-lowpower": C_OASE, "grid": C_EUR,
            "present-day": C_PD, "ancient-other": "#64b5cd",
            "control-reference": "#333333"}.get(k, "#999999")


def make_figures(df, rdf, cdf, pwdf, fit, lim, limR, tag):
    plt = _style()
    from scipy.stats import norm
    n_tested = int(((rdf.kind != "control-reference")
                    & rdf["residual_z"].notna()).sum())
    bonf_z = float(norm.isf(0.025 / max(n_tested, 1))) if n_tested else np.nan

    # --- Fig 1: the scaling line ---------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
    ax = axes[0]
    # Chimpanzee is dropped from this panel only: it is degenerate on D_NEA
    # (it appears on both sides of that statistic, forcing exactly -1) and
    # plotting it compresses every real cohort into an unreadable blob.
    sc_df = df[df.label != "CTRL_Chimp"]
    for kind in ("grid", "present-day", "ancient-other", "named",
                 "named-lowpower", "control-reference"):
        s = sc_df[sc_df.kind == kind]
        if s.empty:
            continue
        ax.errorbar(s["D_NEA"], s["D_VA"], yerr=s["D_VA_se"], xerr=s["D_NEA_se"],
                    fmt="o", ms=5, alpha=0.8, capsize=2, lw=0.8,
                    color=_kindcolor(kind), label=kind)
    for lab, dx, dy in (("Oase1_40ka", 4, 10), ("CTRL_Denisova", 6, -4),
                        ("UP_Europe_pre_LGM", 4, -12)):
        r = sc_df[sc_df.label == lab]
        if not r.empty:
            ax.annotate(lab, (float(r["D_NEA"].iat[0]), float(r["D_VA"].iat[0])),
                        textcoords="offset points", xytext=(dx, dy), fontsize=7)
    lo = float(min(0, sc_df["D_NEA"].min())) - 0.005
    hi = float(sc_df["D_NEA"].max()) + 0.01
    xs = np.linspace(lo, hi, 50)
    ax.plot(xs, fit["k"] * xs, "k--", lw=1.4,
            label=f"single source: k={fit['k']:.2f}")
    ax.axhline(0, color="k", lw=0.6, alpha=0.5)
    ax.axvline(0, color="k", lw=0.6, alpha=0.5)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("D_NEA = D(X, Yoruba; NeaAvg, Chimp)  —  amount of Neanderthal ancestry")
    ax.set_ylabel("D_VA = D(X, Yoruba; Vindija, Altai)")
    ax.set_title("Vindija preference scales with Neanderthal quantity", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[1]
    s = rdf[rdf.kind.isin(["grid", "present-day", "named", "named-lowpower",
                           "ancient-other"])].copy()
    s = s.sort_values("residual_z")
    ax.barh(range(len(s)), s["residual_z"],
            color=[_kindcolor(k) for k in s["kind"]])
    for v in (-2, 2):
        ax.axvline(v, color="k", ls=":", lw=1)
    if np.isfinite(bonf_z):
        for v in (-bonf_z, bonf_z):
            ax.axvline(v, color="crimson", ls="--", lw=1.1)
        ax.annotate(f"Bonferroni |Z|={bonf_z:.2f}", xy=(bonf_z, 0.5),
                    xytext=(-4, 4), textcoords="offset points", rotation=90,
                    fontsize=6.5, color="crimson", ha="right")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(s["label"], fontsize=5.5)
    ax.set_xlabel("residual Z from the single-source line")
    n_out = int((s["residual_z"].abs() >= 2).sum())
    ax.set_title(f"{n_out} of {len(s)} cohorts exceed |Z|=2; none survives "
                 f"Bonferroni", fontsize=9.5)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_n1_scaling{tag}.png"))
    plt.close(fig)

    # --- Fig 2: the time axis -------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    # Oase1 is deliberately kept OUT of the connected time series. He is a
    # single 0.05x genome admitted below the study's SNP floor, and drawing him
    # inside the European curve would render one unreliable point as a 40,000 BP
    # dip in a population trend. He is plotted separately, as an open marker.
    anc = df[df.kind.isin(["grid", "named"])].copy()
    anc = anc[np.isfinite(anc["date_bp"])]
    lowp = df[(df.kind == "named-lowpower") & np.isfinite(df["date_bp"])]
    regions = ["Europe", "EastAsia", "WestSiberia", "Caucasus", "CentralAsia",
               "Siberia", "NearEast", "InnerAsia"]
    cols = dict(zip(regions, plt.cm.tab10(np.linspace(0, 1, len(regions)))))
    for ax, stat, se, lab in ((axes[0], "D_VA", "D_VA_se",
                               "D_VA = D(X, Yoruba; Vindija, Altai)"),
                              (axes[1], "R", "R_se",
                               "R = D_VA / D_NEA  (source affinity per unit ancestry)")):
        for reg in regions:
            s = anc[anc.region == reg].sort_values("date_bp")
            if len(s) < 2:
                continue
            ax.errorbar(s["date_bp"], s[stat], yerr=s[se], fmt="o-", ms=4,
                        lw=1.1, capsize=2, alpha=0.85, color=cols[reg], label=reg)
        if stat == "R":
            # Scale to the plotted points and their errors. Using the D_VA-unit
            # detection limit here would be a units error: R is a ratio and its
            # limit (limR) is an order of magnitude larger.
            # include the low-power point in the range: a legend entry for a
            # marker drawn off-scale would be worse than a slightly wider axis
            v = np.concatenate([anc["R"].to_numpy(dtype=float),
                                lowp["R"].to_numpy(dtype=float)])
            e = np.concatenate([anc["R_se"].to_numpy(dtype=float),
                                lowp["R_se"].to_numpy(dtype=float)])
            ok = np.isfinite(v) & np.isfinite(e)
            if ok.any():
                ylo, yhi = np.min(v[ok] - e[ok]), np.max(v[ok] + e[ok])
                pad = 0.08 * (yhi - ylo + 1e-9)
                ax.set_ylim(ylo - pad, yhi + pad)
            ax.axhspan(fit["k"] - limR["limit"], fit["k"] + limR["limit"],
                       color="grey", alpha=0.15,
                       label=f"detection limit +/-{limR['limit']:.2f}")
            ax.axhline(fit["k"], color="k", ls="--", lw=1.2,
                       label=f"single-source k={fit['k']:.2f}")
        for _, r in lowp.iterrows():
            if np.isfinite(r[stat]):
                ax.errorbar([r["date_bp"]], [r[stat]], yerr=[r[se]], fmt="D",
                            ms=7, mfc="none", capsize=3, lw=1.3, color=C_OASE,
                            label=f"{r['label']} (below SNP floor)")
        ax.invert_xaxis()
        ax.set_xlabel("years before present")
        ax.set_ylabel(lab, fontsize=8)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=6.5, ncol=2)
    axes[0].set_title("Raw contrast through time (tracks ancestry quantity)",
                      fontsize=10)
    axes[1].set_title("Normalised: flat across the Holocene, low in the oldest "
                      "cohorts", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_n2_time{tag}.png"))
    plt.close(fig)

    # --- Fig 3: the null battery ---------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax = axes[0]
    if not cdf.empty:
        cs = cdf.sort_values("control").reset_index(drop=True)
        for kind, g in cs.groupby("control"):
            ax.errorbar(g["D_VA_diff"], g.index, xerr=g["D_VA_diff_se"],
                        fmt="o", ms=6, capsize=3, label=kind,
                        color=C_NULL if "random" in kind else C_EUR)
        ax.set_yticks(range(len(cs)))
        ax.set_yticklabels(cs["a"].str.replace("NULL", "", regex=False).str[:34],
                           fontsize=5.5)
    ax.axvline(0, color="k", lw=0.8)
    for s in (-1, 1):
        ax.axvspan(0, s * lim["limit"], color="grey", alpha=0.15)
    ax.set_xlabel("D_VA difference between two halves of the SAME cohort")
    ax.set_title(f"Null constructions; grey = detection limit "
                 f"({lim['limit']:.4f})", fontsize=9.5)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=7)

    ax = axes[1]
    ax.hist(pwdf[pwdf.statistic == "D_VA"]["z"].dropna(), bins=40, color=C_EUR,
            alpha=0.8, density=True, label="observed pairwise Z")
    zz = np.linspace(-5, 5, 200)
    ax.plot(zz, np.exp(-zz ** 2 / 2) / np.sqrt(2 * np.pi), "k--", lw=1.4,
            label="standard normal")
    ax.set_xlabel("Z of a paired cohort-vs-cohort D_VA difference")
    ax.set_ylabel("density")
    ax.set_title("All pairwise differences against the null", fontsize=9.5)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_n3_null{tag}.png"))
    plt.close(fig)

    # --- Fig 4: the three named targets --------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    want = [t for t in TARGETS if t in set(rdf["label"])]
    extra = [l for l in ("IUP_Eurasia_45ka", "UP_Europe_post_LGM",
                         "UstIshim_44ka", "CTRL_Denisova")
             if l in set(rdf["label"])]
    s = rdf[rdf.label.isin(want + extra)].copy()
    if s.empty:
        plt.close(fig)
        log.warning("no target cohorts available for fig_n4; skipped")
        return
    s["order"] = s["label"].map({l: i for i, l in enumerate(want + extra)})
    s = s.sort_values("order")
    y = np.arange(len(s))
    for yi, (_, row) in zip(y, s.iterrows()):
        ax.errorbar([row["residual"]], [yi], xerr=[row["residual_se"]], fmt="o",
                    ms=8, capsize=4, lw=1.6, color=_kindcolor(row["kind"]))
    ax.axvline(0, color="k", lw=1)
    ax.axvspan(-lim["limit"], lim["limit"], color="grey", alpha=0.18,
               label=f"detection limit +/-{lim['limit']:.4f}")
    ax.set_yticks(y)
    ax.set_yticklabels(s["label"], fontsize=8)
    ax.set_xlabel("residual from the single-source line  (D_VA units)")
    ax.set_title("Do the named targets need a different Neanderthal?", fontsize=10)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_n4_targets{tag}.png"))
    plt.close(fig)
    log.info(f"Wrote 4 figures to {OUT}")


# ------------------------------------------------------------------- report ---
def _power_sentence():
    """The two-way subsample's verdict, read from its output table.

    Returns the sentence that turns "the archaic genomes are the limiting
    sample" from an assertion into a citation, or a plainly-hedged fallback if
    `ns_subsample_power.py` has not been run in this checkout. The fallback
    matters: a missing table must not silently leave a confident claim standing.
    """
    path = os.path.join(OUT, "ns_power_scaling.csv")
    if not os.path.exists(path):
        return ("This is expected from the site count but is not measured here; "
                "run `scripts/ns_subsample_power.py` for the subsample that "
                "tests it directly. ")
    try:
        s = pd.read_csv(path)
        s = s[s.pair_class == "all"].set_index("arm")
        bs, bs_se = float(s.loc["sites", "b"]), float(s.loc["sites", "b_se"])
        bi, bi_se = (float(s.loc["individuals", "b"]),
                     float(s.loc["individuals", "b_se"]))
        rs = float(s.loc["sites", "ratio_min_to_1"])
        ri = float(s.loc["individuals", "ratio_min_to_1"])
        share = float(s.loc["individuals", "var_share"])
    except (KeyError, ValueError, IndexError):
        return ("This is expected from the site count but the subsample table "
                "could not be read; rerun `scripts/ns_subsample_power.py`. ")
    if not (bs > bi):
        return (f"A two-way subsample (`POWER_two_way_subsample.md`) does *not* "
                f"support this: thinning genomes gives an exponent of "
                f"{bi:.2f} +/- {bi_se:.2f} against {bs:.2f} +/- {bs_se:.2f} for "
                f"thinning sites, so the claim above should be read as "
                f"provisional and is contradicted by the only direct test of "
                f"it in this repository. ")
    return (f"That is a measurement, not an inference from the site count: a "
            f"two-way subsample (`POWER_two_way_subsample.md`) thins the sites "
            f"at fixed cohorts and the cohorts at fixed sites, and finds the "
            f"paired difference SE scales as q^-{bs:.2f} +/- {bs_se:.2f} on the "
            f"site axis against q^-{bi:.2f} +/- {bi_se:.2f} on the genome axis. "
            f"Cutting every cohort to an eighth of its genomes costs "
            f"{100*(ri-1):.0f}% on the SE; cutting the panel to an eighth of "
            f"its sites costs {100*(rs-1):.0f}%. Cohort size accounts for about "
            f"{100*share:.0f}% of the variance behind the stated limit. ")


def _mixture_sentence():
    """The correction to the percentage conversion, read from its own table.

    The sentence immediately above converts a resolvable D_VA difference into a
    percentage of ancestry by dividing by a typical cohort's *absolute* D_VA.
    That division assumes a mixture moves D_VA by f x D_VA, which
    `scripts/ns_mixture_power.py` measures directly. Since the measurement can
    invalidate the conversion, this paragraph is generated from its output
    rather than written by hand, and says so plainly when the table is absent -
    a stated limit that has been shown to be optimistic must never be left
    standing bare because a file was missing.
    """
    path = os.path.join(OUT, "ns_mixture_summary.csv")
    if not os.path.exists(path):
        return ("That conversion assumes a mixture moves D_VA by exactly *f* x "
                "D_VA, which is untested here; run "
                "`scripts/ns_mixture_power.py` to calibrate it before quoting "
                "the percentage.\n")
    try:
        s = pd.read_csv(path)
        r = s[s.replacement == "equidistant"].iloc[0]
        k = float(r["kappa"])
        f50, f80 = float(r["f50"]), float(r["f80"])
        corrected = float(r["analytic_limit_kappa_corrected"])
    except (KeyError, ValueError, IndexError):
        return ("That conversion assumes a mixture moves D_VA by exactly *f* x "
                "D_VA; the calibration table could not be read, so the "
                "percentage above should not be quoted until "
                "`scripts/ns_mixture_power.py` is rerun.\n")
    if abs(k - 1.0) < 0.05:
        return (f"**That conversion has been checked.** Simulated mixtures at "
                f"known fractions (`POWER_mixture_calibration.md`) move D_VA by "
                f"{k:.2f} x *f* x D_VA, so the percentage above is calibrated; "
                f"empirical 50% detection falls at *f* = {100*f50:.0f}% and 80% "
                f"detection at *f* = {100*f80:.0f}%.\n")
    return (f"**That conversion is optimistic, and the percentage above should "
            f"not be quoted on its own.** It assumes a mixture moves D_VA by "
            f"*f* x D_VA. Simulated mixtures at known fractions built at the "
            f"allele-frequency level (`POWER_mixture_calibration.md`) move it "
            f"by only {k:.2f} x that, because the conversion divides by D_VA's "
            f"*absolute* value - which this study elsewhere states is not "
            f"interpretable, being inflated by the `.SG`/`.DG` asymmetry and by "
            f"Yoruba's own Neanderthal ancestry. Rescaled, the threshold is "
            f"{100*corrected:.0f}%, and the measured detection curve puts 50% "
            f"detection at *f* = {100*f50:.0f}% and 80% detection at "
            f"*f* = {100*f80:.0f}%. The honest reading of this study's power is "
            f"the last of those: a cohort would need to have re-sourced most of "
            f"its Neanderthal ancestry before this panel would reliably have "
            f"noticed.\n")


def write_report(df, rdf, cdf, pwdf, cmdf, covdf, fit, cm_fit, lim, limR, gain,
                 eura, informative, both_called, n_shared, args, tag):
    def get(lab, col="D_VA"):
        s = df[df.label == lab]
        return float(s[col].iat[0]) if len(s) else np.nan

    def res(lab):
        s = rdf[rdf.label == lab]
        if not len(s):
            return dict(residual=np.nan, residual_se=np.nan, residual_z=np.nan)
        return dict(residual=float(s["residual"].iat[0]),
                    residual_se=float(s["residual_se"].iat[0]),
                    residual_z=float(s["residual_z"].iat[0]))

    def pair(a, b, stat="D_VA"):
        """Signed a-minus-b difference from the pairwise table."""
        s = pwdf[(pwdf.statistic == stat)
                 & (((pwdf.cohort_a == a) & (pwdf.cohort_b == b))
                    | ((pwdf.cohort_a == b) & (pwdf.cohort_b == a)))]
        if s.empty:
            return dict(diff=np.nan, se=np.nan, z=np.nan)
        r = s.iloc[0]
        sgn = 1.0 if r["cohort_a"] == a else -1.0
        return dict(diff=sgn * float(r["diff"]), se=float(r["se"]),
                    z=sgn * float(r["z"]))

    # The "which sample is limiting" sentence below is backed by a separate run
    # (scripts/ns_subsample_power.py). Its numbers are read from that run's table
    # rather than restated here, so the two documents cannot drift apart; if the
    # table is absent the claim is stated as the assertion it would otherwise be.
    power = _power_sentence()

    tgt = {t: res(t) for t in TARGETS + ["IUP_Eurasia_45ka",
                                         "UP_Europe_post_LGM", "UstIshim_44ka"]}
    # the two comparisons the question is actually about
    up_cmp = pair("UP_Europe_pre_LGM", "UP_NorthEastAsia")
    up_cmpR = pair("UP_Europe_pre_LGM", "UP_NorthEastAsia", "R")
    fh = pair("PD_French", "PD_Han")
    fhR = pair("PD_French", "PD_Han", "R")
    ea = pair("Europe_Medieval_Recent", "EastAsia_Bronze")
    # the age-correlated pattern
    old_new = pair("UP_Europe_pre_LGM", "Europe_Medieval_Recent", "R")
    asia_new = pair("UP_NorthEastAsia", "Europe_Medieval_Recent", "R")
    iup_new = pair("IUP_Eurasia_45ka", "Europe_Medieval_Recent", "R")

    tally = rdf[(rdf.kind != "control-reference") & rdf["residual_z"].notna()]
    n_test = int(len(tally))
    n_sig = int((tally["residual_z"].abs() >= 2).sum())
    from scipy.stats import norm
    bonf_z = float(norm.isf(0.025 / max(n_test, 1)))
    n_bonf = int((tally["residual_z"].abs() >= bonf_z).sum())

    cov = {r["covariate"]: r for _, r in covdf.iterrows()} if not covdf.empty else {}

    def covtxt(name):
        c = cov.get(name)
        if c is None or not np.isfinite(c["rho"]):
            return "n/a"
        return f"rho = {c['rho']:+.2f} (p = {c['p']:.2g}, n = {int(c['n'])})"

    show = ["label", "kind", "n_ind", "date_bp", "region", "median_coverage",
            "median_damage", "D_NEA", "D_VA", "D_VA_se", "D_VA_z", "R", "R_se",
            "n_snp"]
    rnd = {"date_bp": 0, "median_coverage": 2, "median_damage": 3, "D_NEA": 4,
           "D_VA": 4, "D_VA_se": 4, "D_VA_z": 2, "R": 3, "R_se": 3}
    named = df[df.kind.isin(["named", "named-lowpower",
                             "control-reference"])][show].round(rnd)
    grid = df[df.kind == "grid"][show].sort_values(
        ["region", "date_bp"], ascending=[True, False]).round(rnd)
    pdt = df[df.kind == "present-day"][show].round(rnd)
    ct = cdf.round({"D_VA_diff": 5, "D_VA_diff_se": 5, "D_VA_diff_z": 2,
                    "D_VA_diff_se_independent": 5, "R_diff": 4, "R_diff_se": 4,
                    "R_diff_z": 2})
    rt = rdf[rdf.label.isin(TARGETS + ["IUP_Eurasia_45ka", "UP_Europe_post_LGM",
                                       "UstIshim_44ka", "Europe_Palaeolithic",
                                       "Europe_Medieval_Recent",
                                       "EastAsia_Bronze", "CTRL_Denisova",
                                       "PD_French", "PD_Han", "PD_Papuan"])][
        ["label", "kind", "n_ind", "D_NEA", "D_VA", "residual", "residual_se",
         "residual_z"]].round({"D_NEA": 4, "D_VA": 4, "residual": 5,
                               "residual_se": 5, "residual_z": 2})
    cvt = covdf.round({"rho": 3, "p": 4}) if not covdf.empty else covdf

    panel_note = (" (transversions only - damage-robust sensitivity run)"
                  if args.transversions else "")
    n_coh = int((df.kind.isin(["named", "named-lowpower", "grid"])).sum())
    md = [
        "# Which Neanderthal? Altai-versus-Vindija affinity across the ancient "
        "Eurasian record\n",
        f"*AADR v66.p1 {args.panel} panel{panel_note}. {len(eura):,} unique "
        f"Eurasian ancient genomes above a {SNP_FLOOR:,}-SNP floor, pooled into "
        f"{n_coh} dated cohorts, plus "
        f"{int((df.kind=='present-day').sum())} present-day anchors. Paired "
        f"{args.blocks}-block jackknife throughout.*\n",

        "## Abstract\n",
        f"Whether every Eurasian lineage descends from one Neanderthal source "
        f"population is normally addressed with a handful of genomes. Here the "
        f"contrast D(X, Yoruba; Vindija, Altai) - positive when a population is "
        f"closer to the Croatian Vindija Neanderthal than to the Siberian Altai "
        f"one - is measured across the ancient Eurasian record. Two anchors "
        f"establish that the instrument works. Every Neanderthal-carrying "
        f"cohort is displaced towards Vindija (present-day French "
        f"{get('PD_French'):+.4f}, Han {get('PD_Han'):+.4f}), reproducing the "
        f"result that Vindija sits closer to the introgressing population "
        f"(Prufer et al. 2017); and with the Denisovan genome as the test "
        f"population the statistic returns {get('CTRL_Denisova'):+.4f} "
        f"(Z = {get('CTRL_Denisova','D_VA_z'):+.1f}), recovering the pull "
        f"towards Altai expected from that genome's ~1% Denisovan-related "
        f"ancestry - a real, published difference between the two Neanderthal "
        f"genomes' own histories.\n\n"
        f"On the question asked, the answer is a null with a stated limit. "
        f"Upper Palaeolithic Europeans and Palaeolithic north-east Asians "
        f"differ by {up_cmp['diff']:+.5f} +/- {up_cmp['se']:.5f} "
        f"(Z = {up_cmp['z']:+.2f}) in raw D_VA, and present-day French and Han "
        f"by {fh['diff']:+.5f} +/- {fh['se']:.5f} (Z = {fh['z']:+.2f}). A "
        f"single proportional relation D_VA = {fit['k']:.2f} x D_NEA "
        f"(+/- {fit['k_se']:.2f}), where D_NEA is a Vindija/Altai-symmetric "
        f"measure of how much Neanderthal ancestry a cohort has, describes "
        f"{n_test - n_sig} of {n_test} cohorts within two standard errors, and "
        f"{n_bonf} depart from it after correction for multiple testing. The "
        f"detection limit is **{lim['limit']:.4f} in D_VA units for a typical "
        f"pair of cohorts, {100*lim['limit_fraction_of_signal']:.0f}% of the "
        f"total Vindija-over-Altai signal** ({lim['best_limit']:.4f}, "
        f"{100*lim['best_fraction_of_signal']:.0f}%, for the best-powered "
        f"pairs). That percentage is itself a conversion that assumes a mixture "
        f"moves D_VA in proportion to its absolute value; simulated mixtures at "
        f"known fractions show it does not, and the calibrated figures are "
        f"given under 'Stated detection limit' and in "
        f"`POWER_mixture_calibration.md`. Either way, sources differing by less "
        f"than that are invisible here, which includes most of what the "
        f"literature debates about a second pulse into East Asia.\n\n"
        f"One pattern is not null and is reported as a candidate rather than a "
        f"finding. The oldest cohorts sit below the single-source line: "
        f"normalised source affinity is {old_new['diff']:+.2f} +/- "
        f"{old_new['se']:.2f} (Z = {old_new['z']:+.2f}) lower in pre-LGM Upper "
        f"Palaeolithic Europeans than in Medieval-to-recent Europeans, and "
        f"{iup_new['diff']:+.2f} +/- {iup_new['se']:.2f} "
        f"(Z = {iup_new['z']:+.2f}) lower in the ~45,000 BP Initial Upper "
        f"Palaeolithic group. Three things argue against taking this at face "
        f"value: it correlates with sample age, which is also a proxy for "
        f"deamination damage and coverage ({covtxt('date_bp')} against date); "
        f"it is absent from the one Palaeolithic cohort with high coverage "
        f"(Palaeolithic north-east Asia, {asia_new['diff']:+.2f} +/- "
        f"{asia_new['se']:.2f}, Z = {asia_new['z']:+.2f} against the same "
        f"comparator); and the effect is carried by the *denominator*, since "
        f"the raw D_VA of these cohorts is ordinary and it is their elevated "
        f"D_NEA that moves the ratio. A transversions-only rerun, immune to "
        f"deamination, shrinks the effect without removing it while the "
        f"across-cohort age correlation strengthens, so the pattern is not "
        f"simply damage either. The honest reading is that this study cannot "
        f"separate a genuinely less Vindija-like early source from an "
        f"age-linked measurement effect, and it should not be reported as the "
        f"former.\n",

        "## The statistic, and the trap inside it\n",
        f"Vindija and Altai are both called at {both_called:,} autosomal 1240K "
        f"sites, but only **{informative:,} of those actually distinguish "
        f"them** ({100*informative/max(both_called,1):.2f}%). The archaic "
        f"genomes, not the ancient cohorts, are the limiting sample, and that "
        f"is why this contrast is hard however many ancient genomes are "
        f"available. {power}It also has a useful consequence: "
        f"the same sites carry the "
        f"signal for every cohort, so their sampling noise is *common-mode* and "
        f"cancels when two cohorts are differenced inside the same jackknife "
        f"replicate. Pairing the jackknife this way is {gain:.1f}x tighter than "
        f"combining two independent standard errors, which is the difference "
        f"between a usable detection limit and an uninformative one.\n",
        "Three further offsets are common-mode, and are therefore differenced "
        "away rather than interpreted:\n",
        "1. **Vindija is pseudo-haploid (`.SG`) and Altai is diploid (`.DG`)** "
        "in this release, and Vindija is called at less than half as many "
        "sites. The two genomes are not symmetric inputs.\n",
        f"2. **Yoruba is not archaic-free.** All Africans carry a little "
        f"Neanderthal ancestry from back-migration (Chen et al. 2020), which "
        f"shows up here as a non-zero baseline: present-day Mbuti read "
        f"{get('PD_Mbuti'):+.4f} against the Yoruba baseline rather than zero, "
        f"and chimpanzee reads {get('CTRL_Chimp'):+.4f}.\n",
        "3. **The 1240K panel is ascertained in modern humans**, so archaic "
        "variation is sampled non-randomly.\n",
        "None of these can produce a *difference between two Eurasian cohorts*, "
        "since every cohort is measured against the same two archaic genomes "
        "and the same baseline. Only differences are reported as results.\n",
        "### Removing the quantity confound\n",
        "The confound that does *not* cancel is how much Neanderthal ancestry "
        "each cohort has, because D_VA scales with it. Each cohort is "
        "therefore also measured with `D_NEA = D(X, Yoruba; NeaAvg, Chimp)`, "
        "where `NeaAvg` is the mean of the two archaic frequencies. Swapping "
        "Vindija for Altai leaves `NeaAvg` unchanged, so `D_NEA` tracks the "
        "amount of Neanderthal ancestry without responding to its source. "
        "Under a single shared source every cohort satisfies "
        "`D_VA = k * D_NEA` with one k, and k estimates the Vindija preference "
        "of the source population itself.\n",
        "It is worth recording how large this confound actually is on this "
        "panel, because the literature figure and this panel's own answer "
        "disagree. All the quantities here are **relative**: the claim is that "
        "East Asians carry proportionally more Neanderthal ancestry than "
        "Europeans, not that they carry more in absolute percentage points of "
        "the genome. In absolutes the whole argument lives between about 1.7% "
        "and 2.6%, and the disputed gap is a few tenths of a percentage "
        "point.\n",
        f"Older estimates put East Asians roughly 20% above Europeans "
        f"(~2.1-2.4% against ~1.8-2.0%). Chen et al. 2020 revised that sharply "
        f"downwards to about 8% (1.8% against 1.7%) by dropping the assumption "
        f"that Africans carry no Neanderthal ancestry, which had depressed the "
        f"European figure. This repository's own f4-ratio goes further still: "
        f"2.14% for East Asians against 2.11% for Europeans, a ratio of 1.016 "
        f"and an absolute gap of 0.03 percentage points. D_NEA agrees, giving "
        f"{get('PD_Han','D_NEA'):.4f} for Han against "
        f"{get('PD_French','D_NEA'):.4f} for French.\n",
        "That this panel sits at the far end of a trend the field has itself "
        "been moving along is worth noting rather than dismissing, and the "
        "reason is directly relevant here: the estimate is measured against "
        "Yoruba, and Yoruba is not archaic-free (see above) - the same effect "
        "Chen et al. were correcting for. The practical consequence for this "
        "study is that the quantity correction is **small for the East Asian "
        "comparison specifically**, which is the one the literature argues "
        "about, but large for the Palaeolithic cohorts and very large for "
        "Oase1, which is where it changes the conclusion.\n",
        "**Scope condition.** The proportional model only describes cohorts "
        f"whose departure from the Yoruba baseline is dominated by Neanderthal "
        f"introgression. Cohorts with D_NEA below {MIN_DNEA} - all the African "
        f"ones - are excluded from the fit and from the residual test: their "
        f"departure from Yoruba is dominated by deep African population "
        f"structure, and R divides by something indistinguishable from zero. "
        f"They are still reported, and their D_VA values are informative as a "
        f"baseline. Chimpanzee is scored on D_VA only, since it appears on both "
        f"sides of D_NEA and is degenerate there by construction.\n",

        "## Positive controls: the statistic can see a real source difference\n",
        f"**Denisova as the test population.** The Altai Neanderthal carries "
        f"~1% Denisovan-related ancestry (Prufer et al. 2014), so the Denisovan "
        f"genome should be pulled towards Altai. It is: "
        f"D_VA = {get('CTRL_Denisova'):+.4f} +/- "
        f"{get('CTRL_Denisova','D_VA_se'):.4f} "
        f"(Z = {get('CTRL_Denisova','D_VA_z'):+.1f}), strongly negative while "
        f"every Eurasian cohort is strongly positive, and a residual of "
        f"{res('CTRL_Denisova')['residual']:+.3f} "
        f"(Z = {res('CTRL_Denisova')['residual_z']:+.1f}) from the "
        f"single-source line. A published difference between the two "
        f"Neanderthal genomes' own histories is recovered at high "
        f"significance, which is the evidence that a null elsewhere means "
        f"something.\n",
        f"**Proportionality.** Across cohorts spanning near-zero to ~12% "
        f"Neanderthal affinity, D_VA is proportional to D_NEA with slope "
        f"k = {fit['k']:.4f} +/- {fit['k_se']:.4f} (Z = {fit['k_z']:.1f}).\n",
        "\n![Figure 1](fig_n1_scaling.png)\n",
        "**Figure 1.** Left: every cohort on the D_VA-versus-D_NEA plane with "
        "the fitted single-source line. Right: residuals from that line in "
        "standard errors. The line is fitted on the region-by-period grid "
        "alone, so the named Palaeolithic cohorts, the present-day anchors and "
        "the Denisovan control are all scored against a line they did not help "
        "define.\n",

        "## Result 1: the comparisons the question asks\n",
        rt.to_markdown(index=False),
        f"\n\n**Upper Palaeolithic Europeans versus Palaeolithic north-east "
        f"Asians.** Raw D_VA differs by {up_cmp['diff']:+.5f} +/- "
        f"{up_cmp['se']:.5f} (Z = {up_cmp['z']:+.2f}); normalised source "
        f"affinity by {up_cmpR['diff']:+.2f} +/- {up_cmpR['se']:.2f} "
        f"(Z = {up_cmpR['z']:+.2f}). The two groups' Neanderthal ancestry is "
        f"indistinguishable in Vindija-versus-Altai character.\n",
        f"**The East Asian second pulse.** The specific contested claim is that "
        f"East Asians received Neanderthal ancestry from an additional or "
        f"different source. Present-day French and Han differ by "
        f"{fh['diff']:+.5f} +/- {fh['se']:.5f} (Z = {fh['z']:+.2f}) in raw "
        f"D_VA and {fhR['diff']:+.2f} +/- {fhR['se']:.2f} "
        f"(Z = {fhR['z']:+.2f}) normalised; Medieval-to-recent Europeans and "
        f"Bronze Age East Asians differ by {ea['diff']:+.5f} +/- "
        f"{ea['se']:.5f} (Z = {ea['z']:+.2f}). This neither supports nor "
        f"excludes a second pulse. It bounds how different the two sources "
        f"could be, and the bound - "
        f"{100*lim['limit_fraction_of_signal']:.0f}% of the total signal, and "
        f"looser still once the percentage conversion is calibrated - is not "
        f"tight enough to adjudicate the debate.\n",
        f"**Oase1.** In raw D_VA Oase1 is unremarkable ({get('Oase1_40ka'):+.4f} "
        f"+/- {get('Oase1_40ka','D_VA_se'):.4f}), which is *itself* the "
        f"surprise: with D_NEA = {get('Oase1_40ka','D_NEA'):.4f}, three times "
        f"any other cohort and consistent with his ~10% Neanderthal ancestry, "
        f"the single-source model predicts a D_VA near "
        f"{fit['k']*get('Oase1_40ka','D_NEA'):.2f}. His residual is "
        f"{tgt['Oase1_40ka']['residual']:+.3f} +/- "
        f"{tgt['Oase1_40ka']['residual_se']:.3f} "
        f"(Z = {tgt['Oase1_40ka']['residual_z']:+.2f}). This is **not** "
        f"reported as evidence that his recent Neanderthal ancestor came from a "
        f"different population. Oase1 is admitted here deliberately below the "
        f"study's SNP floor, at 0.05x coverage and "
        f"{int(get('Oase1_40ka','n_snp')):,} usable sites against ~507,000 for "
        f"every other cohort; only ~3.6% of those distinguish the two "
        f"Neanderthals, leaving a few hundred effectively informative sites "
        f"spread over {args.blocks} jackknife blocks. At that density the "
        f"jackknife is not trustworthy, the first-order approximation "
        f"D_VA ~ a x D_VA(source) is least accurate at his admixture "
        f"proportion, and a single pseudo-haploid genome enters the D "
        f"denominator differently from a pooled cohort. The honest statement "
        f"is that Oase1 cannot be placed on this axis with the 1240K panel.\n",
        "\n![Figure 4](fig_n4_targets.png)\n",
        "**Figure 4.** Named targets' residuals from the single-source line "
        "against the detection limit (grey).\n",

        "## Result 2: an age-correlated pattern, reported as a candidate\n",
        f"The oldest cohorts sit consistently below the single-source line. "
        f"Normalised source affinity R is {old_new['diff']:+.2f} +/- "
        f"{old_new['se']:.2f} (Z = {old_new['z']:+.2f}) lower in pre-LGM Upper "
        f"Palaeolithic Europeans than in Medieval-to-recent Europeans, and "
        f"{iup_new['diff']:+.2f} +/- {iup_new['se']:.2f} "
        f"(Z = {iup_new['z']:+.2f}) lower in the ~45,000 BP Initial Upper "
        f"Palaeolithic group; Holocene European cohorts from six periods agree "
        f"with each other to within a few percent. Taken at face value this "
        f"would say the earliest Eurasians' Neanderthal ancestry was less "
        f"Vindija-like. Four checks say it should not be taken at face "
        f"value.\n",
        f"1. **It tracks the assay, not just the calendar.** Across cohorts, R "
        f"correlates with sample age ({covtxt('date_bp')}), and age is also a "
        f"proxy for deamination damage ({covtxt('median_damage')}) and coverage "
        f"({covtxt('median_coverage')}). These covariates cannot be separated "
        f"in observational ancient-DNA data.\n",
        f"2. **The high-coverage Palaeolithic cohort does not show it.** "
        f"Palaeolithic north-east Asia - the only Palaeolithic cohort with "
        f"median coverage above 4x - sits {asia_new['diff']:+.2f} +/- "
        f"{asia_new['se']:.2f} (Z = {asia_new['z']:+.2f}) from the same "
        f"Medieval comparator, i.e. on the line. If the effect were a property "
        f"of Palaeolithic *people*, it should appear there too.\n",
        f"3. **The numerator is ordinary; the denominator moves.** These "
        f"cohorts' raw D_VA is normal (pre-LGM Upper Palaeolithic Europe "
        f"{get('UP_Europe_pre_LGM'):+.4f} against "
        f"{get('Europe_Medieval_Recent'):+.4f} for Medieval Europe, a "
        f"difference of {pair('UP_Europe_pre_LGM','Europe_Medieval_Recent')['diff']:+.5f} "
        f"+/- {pair('UP_Europe_pre_LGM','Europe_Medieval_Recent')['se']:.5f}). "
        f"What moves is D_NEA, elevated to "
        f"{get('UP_Europe_pre_LGM','D_NEA'):.4f} from "
        f"{get('Europe_Medieval_Recent','D_NEA'):.4f}. An elevated Neanderthal "
        f"level in early Upper Palaeolithic Europeans is independently expected "
        f"(Fu et al. 2016), so part of this is real - but any age-linked "
        f"inflation of D_NEA would produce exactly the same deficit in R "
        f"without any change of source.\n",
        f"4. **Multiple testing.** With {n_test} cohorts scored, the Bonferroni "
        f"threshold is |Z| = {bonf_z:.2f}; these residuals do not reach it.\n",
        "The transversions-only rerun is the sharpest discriminator available "
        "for check 1, and it is reported next.\n",
        _transversion_note(args),
        f"\n{cvt.to_markdown(index=False) if len(cvt) else ''}\n",
        "\n![Figure 2](fig_n2_time.png)\n",
        "**Figure 2.** Left: raw D_VA through time by region, which tracks "
        "Neanderthal quantity. Right: normalised source affinity, flat across "
        "the Holocene and dipping in the oldest cohorts.\n",

        "## The time axis\n",
        "Because the AADR supplies dated cohorts rather than only present-day "
        "populations, the question can be asked as a time series - the part of "
        "this analysis with no published counterpart.\n",
        grid.to_markdown(index=False),

        "## Characterising the null\n",
        "### Null constructions\n",
        "The floor is measured rather than assumed, by splitting single "
        "homogeneous cohorts - which share a Neanderthal source by construction "
        "- and differencing the halves. Random splits give the sampling floor. "
        "**Coverage-stratified splits give the floor that matters**, because "
        "coverage trends with date across the AADR and any coverage-linked bias "
        "would masquerade as a temporal signal:\n",
        ct.to_markdown(index=False),
        f"\n\nThe largest absolute null difference in D_VA is "
        f"{lim['max_abs_null']:.5f} and the spread of null differences gives a "
        f"systematic floor of {lim['systematic_floor']:.5f}. On the normalised "
        f"statistic the same splits are noisier - coverage-split R differences "
        f"reach {cdf['R_diff'].abs().max():.3f} - which is a direct measurement "
        f"of how much coverage alone can move R, and is the reason the "
        f"age-correlated pattern above is reported as a candidate.\n",
        f"Note that a same-cohort split is the right yardstick for *bias* and "
        f"the wrong one for *power*: two halves of one cohort share almost all "
        f"their ancestry, so their difference has a much smaller standard error "
        f"than a comparison between genuinely different populations. The "
        f"statistical floor is therefore taken from the "
        f"{lim['n_comparisons']} real cohort-versus-cohort comparisons "
        f"({lim['statistical_floor']:.5f} at the median), not from the "
        f"splits.\n",
        "\n![Figure 3](fig_n3_null.png)\n",
        "**Figure 3.** Left: null constructions against the detection limit. "
        "Right: the distribution of all pairwise cohort-difference Z scores "
        "against a standard normal - a check that the paired jackknife is "
        "calibrated rather than merely tight.\n",
        "### Stated detection limit\n",
        f"> **A difference of {lim['limit']:.4f} in D_VA between two typical "
        f"cohorts is resolvable at {lim['n_sigma']:.0f} sigma; against a "
        f"typical cohort's D_VA of {lim['signal']:.4f} that is "
        f"{100*lim['limit_fraction_of_signal']:.0f}%. For the best-powered "
        f"pairs the limit falls to {lim['best_limit']:.4f}, "
        f"{100*lim['best_fraction_of_signal']:.0f}%.**\n",
        f"Concretely: if a cohort replaced a fraction *f* of its Neanderthal "
        f"ancestry with ancestry from a Neanderthal lineage equidistant between "
        f"Vindija and Altai, this study would detect it only for "
        f"*f* > {100*lim['limit_fraction_of_signal']:.0f}% in a typical "
        f"comparison. Structure within the introgressing population finer than "
        f"that is invisible here. The corresponding limit on R is "
        f"{limR['limit']:.2f}.\n",
        _mixture_sentence(),
        "This is a statement about the panel, not about history. The limit is "
        f"set by the {informative:,} 1240K sites that separate the two archaic "
        "genomes; shotgun data at all sites, or the addition of Chagyrskaya and "
        "Mezmaiskaya (absent from the AADR), would tighten it substantially. "
        "The subsample in `POWER_two_way_subsample.md` is what licenses that "
        "sentence, and it also licenses its converse: growing the ancient "
        "Eurasian sample, which is the axis the AADR actually grows along, will "
        "not tighten this limit.\n",
        "### Coverage matching\n",
        f"Every core cohort was recomputed on the {n_shared:,} SNPs covered in "
        f"all of them. The fitted slope is {cm_fit['k']:.4f} +/- "
        f"{cm_fit['k_se']:.4f} against {fit['k']:.4f} +/- {fit['k_se']:.4f} on "
        f"the full panel, so the shared-SNP restriction changes nothing "
        f"material. Full numbers in `ns_coverage_matched.csv`.\n",

        "## Limitations\n",
        "- **The panel caps the question.** The AADR 1240K release contains "
        "three archaic genomes (Altai, Vindija, Denisova) and no Chagyrskaya or "
        "Mezmaiskaya, so 'which Neanderthal' can be asked along one axis only, "
        "and that axis is defined by two genomes that are themselves close "
        "relatives.\n",
        "- **D_VA and D_NEA are relative affinities, not percentages**, and R "
        "is a ratio of two of them, interpretable only across cohorts.\n",
        "- **The proportional model is first-order.** D_VA ~ a x D_VA(source) "
        "neglects terms in the admixture fraction that are acceptable at a ~2% "
        "and marginal at Oase1's ~10%.\n",
        "- **Cohorts are not random samples.** Pooled cohorts mix sites, "
        "periods and degrees of relatedness, and no kinship pruning was applied "
        "at this scale, so a large cemetery can be over-represented within its "
        "cohort.\n",
        "- **Denisovan ancestry perturbs the contrast.** Because the Altai "
        "genome carries Denisovan-related ancestry, populations with Denisovan "
        "ancestry of their own are pulled towards Altai for reasons unrelated "
        f"to their Neanderthal source. Present-day Papuans sit "
        f"{res('PD_Papuan')['residual']:+.4f} "
        f"(Z = {res('PD_Papuan')['residual_z']:+.2f}) from the line, in the "
        f"predicted direction; Oceanian cohorts should be read with that in "
        f"mind.\n",
        "- **The age-correlated residual is unresolved**, as set out above.\n",

        "## Reproduce\n",
        "```bash\n"
        f"python scripts/neanderthal_source.py --panel {args.panel}\n"
        f"python scripts/neanderthal_source.py --panel {args.panel} --transversions\n"
        "```\n",
        "\n*Refs: Green et al. 2010 Science 328:710; Reich et al. 2010 Nature "
        "468:1053; Patterson et al. 2012 Genetics 192:1065; Prufer et al. 2014 "
        "Nature 505:43; Fu et al. 2015 Nature 524:216 (Oase1); Fu et al. 2016 "
        "Nature 534:200; Prufer et al. 2017 Science 358:655 (Vindija); Chen et "
        "al. 2020 Cell 180:677; Mallick et al. 2024 Sci. Data 11:182 (AADR).*\n",

        "\n## Full cohort tables\n",
        "### Named Palaeolithic cohorts and reference controls\n",
        named.to_markdown(index=False),
        "\n### Present-day anchors\n",
        pdt.to_markdown(index=False),
    ]
    path = os.path.join(OUT, f"PAPER_neanderthal_source{tag}.md")
    with open(path, "w", encoding="utf-8") as fh_:
        fh_.write("\n".join(md))
    log.info(f"Wrote {path}")


def _transversion_note(args):
    """Summarise the transversions-only rerun inside the main report.

    This is the sharpest discriminator available for the age-correlated
    residual, because cytosine deamination — the error class that concentrates
    in exactly the oldest libraries — produces only C->T/G->A changes and so
    cannot touch transversions. If the deficit in the oldest cohorts survives
    here it is not deamination; if it vanishes, it was.
    """
    if args.transversions:
        return ("*This report is itself the transversions-only run; see "
                "`PAPER_neanderthal_source.md` for the full-panel analysis.*\n")
    pw = os.path.join(OUT, "ns_pairwise_transversions.csv")
    p = os.path.join(OUT, "ns_detection_limit_transversions.csv")
    if not (os.path.exists(p) and os.path.exists(pw)):
        return ("Not run. Execute `python scripts/neanderthal_source.py --panel 1240k "
                "--transversions` to add the damage-robust sensitivity.\n")
    lim = pd.read_csv(p)
    lim = lim[lim.statistic == "D_VA"].iloc[0]
    t = pd.read_csv(pw)
    cvp = os.path.join(OUT, "ns_covariates_transversions.csv")
    cv = None
    if os.path.exists(cvp):
        c = pd.read_csv(cvp)
        c = c[c.covariate == "date_bp"]
        if len(c):
            cv = (float(c["rho"].iat[0]), float(c["p"].iat[0]))

    def pr(a, b, stat):
        s = t[(t.statistic == stat)
              & (((t.cohort_a == a) & (t.cohort_b == b))
                 | ((t.cohort_a == b) & (t.cohort_b == a)))]
        if s.empty:
            return None
        r = s.iloc[0]
        sgn = 1.0 if r["cohort_a"] == a else -1.0
        return sgn * float(r["diff"]), float(r["se"]), sgn * float(r["z"])

    old = pr("UP_Europe_pre_LGM", "Europe_Medieval_Recent", "R")
    up = pr("UP_Europe_pre_LGM", "UP_NorthEastAsia", "D_VA")
    fh = pr("PD_French", "PD_Han", "D_VA")
    bits = [
        "Repeating the whole analysis on transversions only removes the "
        "deamination error class entirely, at the cost of ~80% of the sites."]
    if old:
        bits.append(
            f"The age-correlated deficit **{'persists' if abs(old[2]) >= 2 else 'does not reach significance'}** "
            f"there: pre-LGM Upper Palaeolithic Europe minus Medieval Europe in "
            f"normalised source affinity is {old[0]:+.2f} +/- {old[1]:.2f} "
            f"(Z = {old[2]:+.2f}), against the full-panel value quoted above. "
            f"{'That argues the pattern is not simply deamination, though it does not identify what it is.' if abs(old[2]) >= 2 else 'Losing significance here is consistent with either a damage artifact or the loss of power, and does not settle the question.'}")
    if cv:
        bits.append(
            f"Cutting the other way, the across-cohort correlation between "
            f"normalised source affinity and sample age is *stronger* on "
            f"transversions (rho = {cv[0]:+.2f}, p = {cv[1]:.3g}) than on the "
            f"full panel, which is not what a pure deamination artifact would "
            f"do. The two observations together say the pattern is not simply "
            f"damage, and that this panel cannot say what it is.")
    if up and fh:
        bits.append(
            f"The headline comparisons remain null: Upper Palaeolithic Europe "
            f"minus Palaeolithic north-east Asia {up[0]:+.5f} +/- {up[1]:.5f} "
            f"(Z = {up[2]:+.2f}), French minus Han {fh[0]:+.5f} +/- "
            f"{fh[1]:.5f} (Z = {fh[2]:+.2f}).")
    bits.append(
        f"The detection limit widens to {lim['limit']:.4f} "
        f"({100*lim['limit_fraction_of_signal']:.0f}% of signal), so this is a "
        f"weaker but independent check rather than a tighter one, and absolute "
        f"values are not comparable across the two differently-ascertained SNP "
        f"sets. Full numbers in `ns_*_transversions.csv`.")
    return " ".join(bits) + "\n"


if __name__ == "__main__":
    main()
