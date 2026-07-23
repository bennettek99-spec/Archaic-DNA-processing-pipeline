#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_build_profiles.py — read AADR 1240K, build cohort + reference profiles.

Optimised: caps cohort reads at 200 individuals (subsamples larger pools),
processes cohort profiles first (the essential output), then individual S1/S2.
All prints are flushed so progress is visible.

Outputs:
  results/profiles/cohorts.npz         cohort+ref frequency profiles, masks, chrom
  results/tables/individual_s1_s2.tsv   per-individual S1/S2 (if --with-individuals)
  results/tables/cohort_info.tsv        cohort metadata
"""
from __future__ import annotations
import os, sys, argparse, warnings, random
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _c in (os.path.join(os.path.dirname(_REPO), "archaic-introgression"), os.path.dirname(_REPO)):
    if os.path.isdir(os.path.join(_c, "archaic")):
        sys.path.insert(0, _c); break

from archaic.panel import Panel
from archaic import stats as st, snp_filters
from archaic.config import panel_prefix
from archaic.anno import load_anno

N_BLOCKS = 50
SNP_FLOOR = 30000
MIN_COHORT = 3
MAX_COHORT = 200    # cap pooled-frequency reads at this many individuals
BATCH = 150

REF_SPEC = {
    "Denisova": {"ids": ["Denisova.SG"]},
    "Altai": {"ids": ["AltaiNeanderthal.DG"]},
    "Vindija": {"ids": ["VindijaG1_final.SG"]},
    "Chimp": {"ids": ["Chimp.REF"]},
    "Mbuti": {"pops": ["Mbuti"]},
    "Yoruba": {"pops": ["Yoruba", "YRI", "YRI-Discovery"]},
}
PRESENT_POOLS = {
    "Papuan": ["Papuan"], "Han": ["Han"], "Japanese": ["Japanese"],
    "Dai": ["Dai"], "Karitiana": ["Karitiana"], "French": ["French"],
}


def p(msg):
    print(msg, flush=True)


def categorize_ancient(row):
    gid = str(row.get("genetic_id", "")).strip()
    grp = str(row.get("group_id", "")).strip().lower()
    country = str(row.get("country", "")).strip().lower()
    locality = str(row.get("locality", "")).strip().lower()
    date = row.get("date_bp")
    if not np.isfinite(date) or date < 500:
        return None
    sib = any(k in (grp + " " + locality) for k in
              ("siberia","sakha","malta","afontova","yana","kolyma","baikal","irkutsk","altai"))
    if sib or country == "russia":
        if date > 25000: return "Sib_UP"
        if date > 15000: return "Sib_LGM"
        if date > 10000: return "Sib_Hol_early"
        return "Sib_Hol_late"
    am = country in ("usa","canada","alaska","greenland","mexico","guatemala","belize",
                     "honduras","el salvador","nicaragua","costa rica","panama","colombia",
                     "venezuela","ecuador","peru","bolivia","brazil","paraguay","chile",
                     "argentina","uruguay","cuba","dominican republic","haiti","bahamas","puerto rico")
    if am:
        if country in ("alaska","greenland") or "aleut" in grp: return "AmBeringian"
        if country in ("usa","canada"): return "AmN_early" if date > 8000 else "AmN_late"
        if country in ("mexico","guatemala","belize","honduras","el salvador","nicaragua","costa rica","panama"): return "AmMeso"
        if country in ("cuba","dominican republic","haiti","bahamas","puerto rico"): return "AmCarib"
        return "AmS_early" if date > 8000 else "AmS_late"
    ea = country in ("china","japan","south korea","korea","north korea","mongolia","taiwan","vietnam","thailand","cambodia","laos","myanmar")
    if ea or "jomon" in grp:
        if "jomon" in grp: return "Jomon"
        return "EAsia_early" if date > 8000 else "EAsia_late"
    return None


ARCHAIC_IDS = {
    "Denisova.SG", "Denisova3.DG", "Denisova3_snpAD.DG", "Denisova11.SG",
    "Denisova25.SG", "AltaiNeanderthal.DG", "VindijaG1_final.SG",
    "Chagyrskaya8.DG", "Chimp.REF", "Chimp_HO.HO",
}


def define_cohorts(anno):
    cohorts = {}
    for name, pops in PRESENT_POOLS.items():
        ids = anno[anno["group_id"].isin(pops) & anno["snps_1240k"].notna()
                   & ~anno["genetic_id"].isin(ARCHAIC_IDS)]["genetic_id"].tolist()
        if len(ids) >= MIN_COHORT: cohorts[name] = ids
    anc = anno[anno["snps_1240k"].notna() & (anno["snps_1240k"] >= SNP_FLOOR)
               & ~anno["genetic_id"].isin(ARCHAIC_IDS)].copy()
    for _, row in anc.iterrows():
        cat = categorize_ancient(row)
        if cat: cohorts.setdefault(cat, []).append(row["genetic_id"])
    for _, row in anc.iterrows():
        date = row.get("date_bp")
        country = str(row.get("country", "")).strip().lower()
        if not np.isfinite(date) or date < 500: continue
        relevant = any(country == c for c in (
            "usa","canada","alaska","greenland","mexico","guatemala","belize","honduras",
            "el salvador","nicaragua","costa rica","panama","colombia","venezuela","ecuador",
            "peru","bolivia","brazil","paraguay","chile","argentina","uruguay","cuba",
            "dominican republic","haiti","bahamas","puerto rico","russia","china","japan",
            "south korea","korea","north korea","mongolia","taiwan","vietnam","thailand",
            "cambodia","laos","myanmar"))
        if not relevant: continue
        if date > 40000: tt = "TT_pre40k"
        elif date > 25000: tt = "TT_40_25k"
        elif date > 15000: tt = "TT_25_15k"
        elif date > 10000: tt = "TT_15_10k"
        elif date > 5000: tt = "TT_10_5k"
        else: tt = "TT_post5k"
        cohorts.setdefault(tt, []).append(row["genetic_id"])
    # also exclude archaic IDs from time-transect bins (they are >40k years old)
    for tt in list(cohorts.keys()):
        cohorts[tt] = [sid for sid in cohorts[tt] if sid not in ARCHAIC_IDS]
    return {k: sorted(set(v)) for k, v in cohorts.items() if len(v) >= MIN_COHORT}


def pooled_freq_batch(panel, cols, n_snp):
    """Compute pooled freq reading in sub-batches of MAX_COHORT to bound memory."""
    if len(cols) > MAX_COHORT:
        random.seed(42)
        cols = np.array(sorted(random.sample(list(cols), MAX_COHORT)), dtype=np.int64)
    G = panel.pg.read(panel.snp_rows, cols).astype(np.float32)
    G[G < 0] = np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        p = np.nanmean(G, axis=1) / 2.0 if len(cols) > 1 else (G[:, 0] / 2.0)
    return p.astype(np.float64), int(len(cols)), int(np.isfinite(p).sum())


def build_set_a(refs):
    p_den = refs["Denisova"]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        p_nea = np.nanmean(np.vstack([refs["Altai"], refs["Vindija"]]), axis=0)
        p_afr = np.nanmean(np.vstack([refs["Mbuti"], refs["Yoruba"]]), axis=0)
    den_is_a1 = p_den > 0.5
    den_ext = np.where(den_is_a1, p_den, 1.0 - p_den)
    afr = np.where(den_is_a1, p_afr, 1.0 - p_afr)
    nea = np.where(den_is_a1, p_nea, 1.0 - p_nea)
    fin = np.isfinite(p_den) & np.isfinite(p_afr) & np.isfinite(p_nea)
    set_a = fin & (den_ext >= 0.90) & (afr <= 0.10) & (nea <= 0.10)
    set_b = fin & (den_ext >= 0.90) & (afr <= 0.10) & (nea <= 0.50)
    return set_a, set_b, den_is_a1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", default="1240k")
    ap.add_argument("--config", default=None)
    ap.add_argument("--with-individuals", action="store_true",
                    help="also compute per-individual S1/S2 (slower)")
    args = ap.parse_args(argv)

    outp = os.path.join(_REPO, "results", "profiles")
    outt = os.path.join(_REPO, "results", "tables")
    os.makedirs(outp, exist_ok=True); os.makedirs(outt, exist_ok=True)

    prefix = panel_prefix(args.panel, args.config)
    p(f"Loading panel {args.panel}...")
    panel = Panel(prefix, autosomes_only=True)
    n_snp = panel.n_snp
    p(f"  autosomal SNPs: {n_snp:,}")

    anno = load_anno(prefix + ".anno")
    cohorts = define_cohorts(anno)
    p(f"\nCohorts ({len(cohorts)}):")
    for name, ids in sorted(cohorts.items(), key=lambda x: -len(x[1])):
        p(f"  {name:<20} n={len(ids):>5}")

    # --- reference frequencies ---
    p("\nReference frequencies...")
    ref_freq, ref_info = panel.frequencies(REF_SPEC)
    for k in REF_SPEC:
        p(f"  {k}: n_ind={ref_info[k]['n_ind']}, snp={ref_info[k]['n_snp_covered']:,}")

    # --- diagnostic sets ---
    p("\nBuilding diagnostic SNP sets...")
    set_a, set_b, den_is_a1 = build_set_a(ref_freq)
    p(f"  Set A: {int(set_a.sum()):,}   Set B: {int(set_b.sum()):,}")

    # --- cohort profiles ---
    p("\nComputing cohort profiles...")
    all_freq = dict(ref_freq)
    info_rows = []
    for name in REF_SPEC:
        info_rows.append({"cohort": name, "kind": "reference",
                          "n_ind": ref_info[name]["n_ind"],
                          "n_snp": ref_info[name]["n_snp_covered"]})
    for name, ids in sorted(cohorts.items()):
        cols = panel.cols_for(ids=ids)
        if len(cols) < MIN_COHORT: continue
        f, n_ind, n_snp_cov = pooled_freq_batch(panel, cols, n_snp)
        all_freq[name] = f
        kind = "present_day" if name in PRESENT_POOLS else ("transect" if name.startswith("TT_") else "ancient")
        dates = anno[anno["genetic_id"].isin(ids)]["date_bp"].dropna()
        info_rows.append({"cohort": name, "kind": kind, "n_ind": n_ind,
                          "n_snp": n_snp_cov,
                          "date_min": dates.min() if len(dates) else np.nan,
                          "date_max": dates.max() if len(dates) else np.nan})
        p(f"  {name:<20} n={n_ind:>5}  snp={n_snp_cov:>8,}")

    # --- masks ---
    chrom = panel.snp.loc[panel.snp_rows, "chrom"].to_numpy()
    tv_mask = np.asarray(snp_filters.transversion_mask(
        panel.snp.loc[panel.snp_rows]), dtype=bool)

    # --- save ---
    names = list(all_freq.keys())
    freq_matrix = np.vstack([all_freq[n] for n in names]).astype(np.float32)
    p(f"\nSaving {len(names)} profiles ({freq_matrix.shape})...")
    np.savez_compressed(
        os.path.join(outp, "cohorts.npz"),
        names=np.array(names, dtype=object),
        freqs=freq_matrix,
        snp_rows=panel.snp_rows,
        chrom=chrom, tv_mask=tv_mask,
        set_a=set_a, set_b=set_b, den_is_a1=den_is_a1,
    )
    pd.DataFrame(info_rows).to_csv(
        os.path.join(outt, "cohort_info.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"Saved cohort profiles + masks to {outp}/cohorts.npz")
    p(f"Saved cohort info to {outt}/cohort_info.tsv")

    # --- individual S1/S2 (optional) ---
    if args.with_individuals:
        p("\nComputing individual-level S1/S2...")
        anc_ids = sorted(set(sid for name, ids in cohorts.items()
                             if name not in PRESENT_POOLS and not name.startswith("TT_")
                             for sid in ids))
        all_cols = panel.cols_for(ids=anc_ids)
        p(f"  {len(anc_ids)} ancient individuals to process")
        id_by_col = {panel._id_to_col.get(sid): sid for sid in anc_ids
                     if panel._id_to_col.get(sid) is not None}

        pMbuti = ref_freq["Mbuti"]
        pDen = ref_freq["Denisova"]
        pChimp = ref_freq["Chimp"]
        pFrench = all_freq.get("French", ref_freq.get("Mbuti"))
        ok1 = np.isfinite(pMbuti) & np.isfinite(pDen) & np.isfinite(pChimp)
        ok2 = ok1 & np.isfinite(pFrench)
        idx1 = np.where(ok1)[0]
        idx2 = np.where(ok2)[0]
        starts1 = st.block_starts(len(idx1), N_BLOCKS)
        starts2 = st.block_starts(len(idx2), N_BLOCKS)
        ydiff1 = (pDen[idx1] - pChimp[idx1])
        yden1 = (pDen[idx1] + pChimp[idx1] - 2.0 * pDen[idx1] * pChimp[idx1])
        ydiff2 = (pDen[idx2] - pChimp[idx2])
        yden2 = (pDen[idx2] + pChimp[idx2] - 2.0 * pDen[idx2] * pChimp[idx2])

        rows = []
        for c_start in range(0, len(all_cols), BATCH):
            c_end = min(c_start + BATCH, len(all_cols))
            cols_chunk = all_cols[c_start:c_end]
            for stat, idx, starts, ydiff, yden, pBase in [
                ("S1", idx1, starts1, ydiff1, yden1, pMbuti),
                ("S2", idx2, starts2, ydiff2, yden2, pFrench),
            ]:
                G = panel.pg.read(panel.snp_rows[idx], cols_chunk).astype(np.float32)
                G[G < 0] = np.nan
                pX = G / 2.0  # per-individual: (n_snp, n_ind) float32
                pB = pBase[idx, None]
                num = (pX - pB) * ydiff[:, None]
                den = (pX + pB - 2 * pX * pB) * yden[:, None]
                num = np.where(np.isfinite(num), num, np.nan)
                den = np.where(np.isfinite(den), den, np.nan)
                theta, se, z, nsnp = st.batch_jackknife_ratio(num, den, starts)
                if stat == "S1":
                    for i, col in enumerate(cols_chunk):
                        sid = id_by_col.get(col, f"col_{col}")
                        r = anno[anno["genetic_id"] == sid]
                        grp = r.iloc[0]["group_id"] if len(r) else ""
                        date = r.iloc[0]["date_bp"] if len(r) else np.nan
                        rows.append({"genetic_id": sid, "group_id": grp, "date_bp": date,
                                     "S1_D": theta[i], "S1_SE": se[i], "S1_Z": z[i], "S1_nSNP": nsnp[i],
                                     "S2_D": np.nan, "S2_SE": np.nan, "S2_Z": np.nan, "S2_nSNP": np.nan})
                else:
                    for i, col in enumerate(cols_chunk):
                        sid = id_by_col.get(col, f"col_{col}")
                        # find matching row
                        for row in rows:
                            if row["genetic_id"] == id_by_col.get(col):
                                row["S2_D"] = theta[i]; row["S2_SE"] = se[i]
                                row["S2_Z"] = z[i]; row["S2_nSNP"] = nsnp[i]
                                break
            p(f"  batch {c_start}-{c_end}/{len(all_cols)} done")

        idf = pd.DataFrame(rows)
        idf.to_csv(os.path.join(outt, "individual_s1_s2.tsv"), sep="\t", index=False, na_rep="NA")
        p(f"  {len(idf)} individuals -> {outt}/individual_s1_s2.tsv")

    p("\nDone.")


if __name__ == "__main__":
    raise SystemExit(main())
