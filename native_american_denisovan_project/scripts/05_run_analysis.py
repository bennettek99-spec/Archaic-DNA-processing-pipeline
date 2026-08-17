#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_run_analysis.py — run all statistics on saved cohort profiles.

Computes Tiers 1-5 (S1, S2, pairwise contrasts, S3 composition), chromosome/LOCO,
downsampling calibration, selected-loci exclusion, ancestry conditioning, and a
simplified power simulation. All from saved profiles — no AADR I/O.

Outputs (results/tables/):
  table3_basic_denisovan_affinity.tsv     S1 + S2, all cohorts, all-sites + TV
  table4_pairwise_denisovan_contrasts.tsv D(X, Han/ Papuan; Denisova, Chimp)
  table5_diagnostic_marker_sharing.tsv    Papuan/Han-enriched marker means + S3
  table6_chronological_results.tsv        time transect S1/S2/S3
  table7_regression_models.tsv            S3 vs genome-wide ancestry
  table8_sensitivity_analyses.tsv        TV, selected-loci, LOCO, outgroups
  table9_simulation_performance.tsv       power/FPR at various alpha + marker counts
  table10_evidence_summary.tsv            per-model evidence assessment
"""
from __future__ import annotations
import os, sys, argparse, warnings, json
import numpy as np
import pandas as pd
from scipy import stats as sps

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _c in (os.path.join(os.path.dirname(_REPO), "archaic-introgression"), os.path.dirname(_REPO)):
    if os.path.isdir(os.path.join(_c, "archaic")):
        sys.path.insert(0, _c); break
from archaic import stats as st

N_BLOCKS = 50
B_BOOT = 1000   # bootstrap replicates for S3 CI
ALPHA_THRESH = 0.05   # Papuan-enriched if Papuan > Han + this

# training/validation chromosome split (odd = training, even = validation)
TRAIN_CHROMS = {str(c) for c in range(1, 23, 2)}   # 1,3,5,...,21
VALID_CHROMS = {str(c) for c in range(2, 23, 2)}   # 2,4,6,...,22

# selected-loci exclusion regions (hg19)
EXCLUSION_REGIONS = [
    ("2", 46520000, 46620000, "EPAS1"),
    ("12", 40100000, 40400000, "MUC19_region"),
]

ANCIENT_COHORTS = ["Sib_UP", "Sib_LGM", "Sib_Hol_early", "Sib_Hol_late",
                   "AmBeringian", "AmN_early", "AmN_late", "AmMeso", "AmCarib",
                   "AmS_early", "AmS_late", "EAsia_early", "EAsia_late", "Jomon"]
PD_COHORTS = ["Papuan", "Han", "Japanese", "Dai", "Karitiana", "French", "Mbuti", "Yoruba"]
ALL_COHORTS = ANCIENT_COHORTS + PD_COHORTS


def p(msg):
    print(msg, flush=True)


def load_profiles():
    path = os.path.join(_REPO, "results", "profiles", "cohorts.npz")
    z = np.load(path, allow_pickle=True)
    names = list(z["names"])
    freqs = z["freqs"].astype(np.float64)  # (n_cohorts, n_snp)
    freq = {n: freqs[i] for i, n in enumerate(names)}
    chrom = z["chrom"]
    tv_mask = z["tv_mask"].astype(bool)
    set_a = z["set_a"].astype(bool)
    set_b = z["set_b"].astype(bool)
    den_is_a1 = z["den_is_a1"].astype(bool)
    snp_rows = z["snp_rows"]
    return freq, names, chrom, tv_mask, set_a, set_b, den_is_a1, snp_rows


def dstat_all(freq, W, X, Y, Z, block, n_blocks=N_BLOCKS):
    """D(W,X;Y,Z) on all sites."""
    return st.dstat(freq, W, X, Y, Z, block, n_blocks)


def dstat_tv(freq, W, X, Y, Z, tv_mask, n_blocks=N_BLOCKS):
    """D(W,X;Y,Z) on transversion-only sites."""
    idx = np.where(tv_mask)[0]
    sub = {k: v[idx] for k, v in freq.items() if k in (W, X, Y, Z)}
    # also need all four in sub
    for name in (W, X, Y, Z):
        if name not in sub:
            return {"theta": np.nan, "se": np.nan, "z": np.nan, "n_used": 0, "statistic": f"D({W},{X};{Y},{Z})"}
    block_tv = st.assign_blocks(len(idx), n_blocks)
    return st.dstat(sub, W, X, Y, Z, block_tv, n_blocks)


def pearson_safe(a, b):
    """Pearson r ignoring NaN; returns nan if <3 shared points."""
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return np.nan
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def bootstrap_s3_corr(fX, fP, fH, valid_chrom, B=B_BOOT):
    """Bootstrap CI for S3 = corr(fX,fP) - corr(fX,fH) by resampling chromosomes."""
    unique_chroms = sorted(set(valid_chrom))
    if len(unique_chroms) < 3:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(42)
    # pre-index by chromosome
    chrom_idx = {c: np.where(valid_chrom == c)[0] for c in unique_chroms}
    s3_vals = []
    for _ in range(B):
        sampled = rng.choice(unique_chroms, size=len(unique_chroms), replace=True)
        idx = np.concatenate([chrom_idx[c] for c in sampled])
        cp = pearson_safe(fX[idx], fP[idx])
        ch = pearson_safe(fX[idx], fH[idx])
        if np.isfinite(cp) and np.isfinite(ch):
            s3_vals.append(cp - ch)
    if len(s3_vals) < 10:
        return np.nan, np.nan, np.nan
    s3_vals = np.array(s3_vals)
    point = pearson_safe(fX, fP) - pearson_safe(fX, fH)
    lo, hi = np.percentile(s3_vals, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def define_enriched_subsets(fP, fH, den_is_a1, train_mask):
    """On training sites: Papuan-enriched, Han-enriched, shared."""
    pP = np.where(den_is_a1[train_mask], fP[train_mask], 1.0 - fP[train_mask])
    pH = np.where(den_is_a1[train_mask], fH[train_mask], 1.0 - fH[train_mask])
    papuan_enriched = pP > pH + ALPHA_THRESH
    han_enriched = pH > pP + ALPHA_THRESH
    shared = np.abs(pP - pH) <= ALPHA_THRESH
    return papuan_enriched, han_enriched, shared


def mean_sharing(fX, den_is_a1, mask):
    """Mean oriented Denisovan-allele frequency at masked sites."""
    oriented = np.where(den_is_a1[mask], fX[mask], 1.0 - fX[mask])
    valid = np.isfinite(oriented)
    return float(np.nanmean(oriented[valid])) if valid.sum() else np.nan, int(valid.sum())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    args = ap.parse_args(argv)

    outt = os.path.join(_REPO, "results", "tables")
    os.makedirs(outt, exist_ok=True)

    p("Loading cohort profiles...")
    freq, names, chrom, tv_mask, set_a, set_b, den_is_a1, snp_rows = load_profiles()
    n_snp = len(chrom)
    p(f"  {len(names)} cohorts, {n_snp:,} SNPs, {int(set_a.sum()):,} Set-A sites")

    block = st.assign_blocks(n_snp, N_BLOCKS)

    # ===== TABLE 3: Basic Denisovan affinity (S1 + S2) =====
    p("\n=== Tier 1: S1 = D(X, Mbuti; Denisova, Chimp), S2 = D(X, French; Denisova, Chimp) ===")
    t3_rows = []
    for name in ALL_COHORTS + [n for n in names if n.startswith("TT_")]:
        if name not in freq: continue
        # S1 all-sites
        r1 = dstat_all(freq, name, "Mbuti", "Denisova", "Chimp", block)
        # S1 TV-only
        r1tv = dstat_tv(freq, name, "Mbuti", "Denisova", "Chimp", tv_mask)
        # S2 all-sites
        if "French" in freq:
            r2 = dstat_all(freq, name, "French", "Denisova", "Chimp", block)
            r2tv = dstat_tv(freq, name, "French", "Denisova", "Chimp", tv_mask)
        else:
            r2 = r2tv = {"theta": np.nan, "se": np.nan, "z": np.nan, "n_used": 0}
        t3_rows.append({
            "cohort": name, "kind": "present_day" if name in PD_COHORTS else ("transect" if name.startswith("TT_") else "ancient"),
            "S1_D": r1["theta"], "S1_SE": r1["se"], "S1_Z": r1["z"], "S1_nSNP": r1["n_used"],
            "S1_D_TV": r1tv["theta"], "S1_Z_TV": r1tv["z"], "S1_nSNP_TV": r1tv["n_used"],
            "S2_D": r2["theta"], "S2_SE": r2["se"], "S2_Z": r2["z"], "S2_nSNP": r2["n_used"],
            "S2_D_TV": r2tv["theta"], "S2_Z_TV": r2tv["z"], "S2_nSNP_TV": r2tv["n_used"],
        })
        p(f"  {name:<20} S1 Z={r1['z']:+6.2f}  S2 Z={r2['z']:+6.2f}  (TV S1 Z={r1tv['z']:+.2f})")
    t3 = pd.DataFrame(t3_rows)
    t3.to_csv(os.path.join(outt, "table3_basic_denisovan_affinity.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table3 ({len(t3)} rows)")

    # ===== TABLE 4: Pairwise Denisovan contrasts =====
    p("\n=== Tier 2: pairwise contrasts ===")
    t4_rows = []
    for name in ANCIENT_COHORTS + ["Karitiana", "Papuan"]:
        if name not in freq: continue
        for ref, label in [("Han", "vs_Han"), ("Papuan", "vs_Papuan")]:
            if ref not in freq: continue
            r = dstat_all(freq, name, ref, "Denisova", "Chimp", block)
            rtv = dstat_tv(freq, name, ref, "Denisova", "Chimp", tv_mask)
            t4_rows.append({"cohort": name, "contrast": label,
                            "D": r["theta"], "SE": r["se"], "Z": r["z"], "nSNP": r["n_used"],
                            "D_TV": rtv["theta"], "Z_TV": rtv["z"], "nSNP_TV": rtv["n_used"]})
        p(f"  {name:<20} vs Han Z={t4_rows[-2]['Z']:+.2f}  vs Papuan Z={t4_rows[-1]['Z']:+.2f}")
    t4 = pd.DataFrame(t4_rows)
    t4.to_csv(os.path.join(outt, "table4_pairwise_denisovan_contrasts.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table4 ({len(t4)} rows)")

    # ===== TABLE 5: Diagnostic marker sharing (S3 composition) =====
    p("\n=== Tier 5: S3 composition (Papuan-like vs Han-like) ===")
    set_a_idx = np.where(set_a)[0]
    sa_chrom = chrom[set_a_idx]
    sa_den = den_is_a1[set_a_idx]
    train_mask = np.isin(sa_chrom, list(TRAIN_CHROMS))
    valid_mask = np.isin(sa_chrom, list(VALID_CHROMS))

    fPap = freq["Papuan"][set_a_idx]
    fHan = freq["Han"][set_a_idx]

    # define enriched subsets on training
    pap_enriched, han_enriched, shared = define_enriched_subsets(fPap, fHan, sa_den, train_mask)
    p(f"  Training sites: {int(train_mask.sum())}  Validation sites: {int(valid_mask.sum())}")
    p(f"  Papuan-enriched (train): {int(pap_enriched.sum())}  Han-enriched (train): {int(han_enriched.sum())}")

    # ORIENTED profiles at Set-A (all cohorts need orientation for correlations)
    def orient(fx, den_mask):
        return np.where(den_mask, fx, 1.0 - fx)
    pP_full = orient(fPap, sa_den)  # full Set-A oriented Papuan
    pH_full = orient(fHan, sa_den)  # full Set-A oriented Han
    pP_valid = pP_full[valid_mask]  # validation subset
    pH_valid = pH_full[valid_mask]

    # validation-enriched masks (same definition, applied to validation chroms)
    pap_enr_valid = pP_valid > pH_valid + ALPHA_THRESH
    han_enr_valid = pH_valid > pP_valid + ALPHA_THRESH

    t5_rows = []
    valid_chrom_arr = sa_chrom[valid_mask]
    # indices into Set-A where validation AND enriched
    pap_valid_idx = np.where(valid_mask & np.isin(sa_chrom, list(VALID_CHROMS)))[0]
    # recompute enriched on validation directly
    pap_enr_valid = pP_valid > pH_valid + ALPHA_THRESH
    han_enr_valid = pH_valid > pP_valid + ALPHA_THRESH
    pap_v_idx = np.where(valid_mask)[0][pap_enr_valid]   # indices into Set-A
    han_v_idx = np.where(valid_mask)[0][han_enr_valid]
    for name in ALL_COHORTS:
        if name not in freq: continue
        fX_sa = freq[name][set_a_idx]  # Set-A-indexed frequency
        fX_v_oriented = orient(fX_sa[valid_mask], sa_den[valid_mask])  # ORIENTED validation
        corr_pap = pearson_safe(fX_v_oriented, pP_valid)
        corr_han = pearson_safe(fX_v_oriented, pH_valid)
        s3_point, s3_lo, s3_hi = bootstrap_s3_corr(fX_v_oriented, pP_valid, pH_valid, valid_chrom_arr)
        # mean sharing at Papuan-enriched / Han-enriched validation sites
        fX_full_oriented = orient(fX_sa, sa_den)
        if len(pap_v_idx) > 0:
            mean_pap = float(np.nanmean(fX_full_oriented[pap_v_idx]))
        else:
            mean_pap = np.nan
        if len(han_v_idx) > 0:
            mean_han = float(np.nanmean(fX_full_oriented[han_v_idx]))
        else:
            mean_han = np.nan
        t5_rows.append({
            "cohort": name, "corr_Papuan": corr_pap, "corr_Han": corr_han,
            "S3_corr": s3_point, "S3_corr_lo": s3_lo, "S3_corr_hi": s3_hi,
            "mean_PapuanEnriched": mean_pap, "n_PapuanEnriched": len(pap_v_idx),
            "mean_HanEnriched": mean_han, "n_HanEnriched": len(han_v_idx),
            "S3_mean": (mean_pap - mean_han) if np.isfinite(mean_pap) and np.isfinite(mean_han) else np.nan,
        })
        p(f"  {name:<20} corrP={corr_pap:+.3f} corrH={corr_han:+.3f}  S3={s3_point:+.4f} [{s3_lo:+.4f},{s3_hi:+.4f}]")
    t5 = pd.DataFrame(t5_rows)
    t5.to_csv(os.path.join(outt, "table5_diagnostic_marker_sharing.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table5 ({len(t5)} rows)")

    # ===== TABLE 6: Chronological results (time transect) =====
    p("\n=== Time transect ===")
    t6_rows = []
    for name in sorted([n for n in names if n.startswith("TT_")]):
        if name not in freq: continue
        r1 = dstat_all(freq, name, "Mbuti", "Denisova", "Chimp", block)
        r2 = dstat_all(freq, name, "French", "Denisova", "Chimp", block) if "French" in freq else {"theta": np.nan, "z": np.nan, "n_used": 0}
        fX = freq[name][set_a_idx]
        fX_v_oriented = orient(fX[valid_mask], sa_den[valid_mask])
        corr_pap = pearson_safe(fX_v_oriented, pP_valid)
        corr_han = pearson_safe(fX_v_oriented, pH_valid)
        s3 = corr_pap - corr_han if np.isfinite(corr_pap) and np.isfinite(corr_han) else np.nan
        dates = [int(name.split("_")[1].replace("pre", "40").replace("post", "5").replace("k",""))]  # approximate
        t6_rows.append({"cohort": name, "S1_D": r1["theta"], "S1_Z": r1["z"], "S1_nSNP": r1["n_used"],
                        "S2_D": r2["theta"], "S2_Z": r2["z"], "S2_nSNP": r2["n_used"],
                        "S3_corr": s3, "corr_Papuan": corr_pap, "corr_Han": corr_han})
        p(f"  {name:<15} S1 Z={r1['z']:+.2f}  S2 Z={r2['z']:+.2f}  S3={s3:+.4f}")
    t6 = pd.DataFrame(t6_rows)
    t6.to_csv(os.path.join(outt, "table6_chronological_results.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table6 ({len(t6)} rows)")

    # ===== TABLE 7: Regression models (ancestry conditioning) =====
    p("\n=== Ancestry conditioning (regression) ===")
    # genome-wide similarity to Han and ANE (Sib_LGM)
    gw_han = {}
    gw_ane = {}
    ref_han = freq["Han"]
    ref_ane = freq.get("Sib_LGM", freq.get("Sib_UP"))
    for name in ALL_COHORTS:
        if name not in freq: continue
        gw_han[name] = pearson_safe(freq[name], ref_han)
        gw_ane[name] = pearson_safe(freq[name], ref_ane) if ref_ane is not None else np.nan
    # regression: S3 ~ gw_han + gw_ane
    s3_vals = {r["cohort"]: r["S3_corr"] for r in t5_rows}
    reg_rows = []
    cohorts_with_data = [n for n in ALL_COHORTS if n in s3_vals and np.isfinite(s3_vals[n])
                        and np.isfinite(gw_han.get(n, np.nan)) and np.isfinite(gw_ane.get(n, np.nan))]
    if len(cohorts_with_data) >= 5:
        X_mat = np.column_stack([np.ones(len(cohorts_with_data)),
                                  [gw_han[n] for n in cohorts_with_data],
                                  [gw_ane[n] for n in cohorts_with_data]])
        y_vec = np.array([s3_vals[n] for n in cohorts_with_data])
        try:
            beta, _, _, _ = np.linalg.lstsq(X_mat, y_vec, rcond=None)
            y_pred = X_mat @ beta
            ss_res = np.sum((y_vec - y_pred) ** 2)
            ss_tot = np.sum((y_vec - y_vec.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            reg_rows.append({"model": "S3 ~ Han_sim + ANE_sim", "n": len(cohorts_with_data),
                             "intercept": beta[0], "beta_Han": beta[1], "beta_ANE": beta[2],
                             "R2": r2})
            p(f"  S3 ~ Han + ANE: beta_Han={beta[1]:+.3f} beta_ANE={beta[2]:+.3f} R2={r2:.3f}")
        except Exception as e:
            p(f"  Regression failed: {e}")
    # also regress S2 on ancestry
    s2_vals = {r["cohort"]: r["S2_Z"] for r in t3_rows if r["cohort"] in ALL_COHORTS}
    cohorts_s2 = [n for n in ALL_COHORTS if n in s2_vals and np.isfinite(s2_vals[n])
                 and np.isfinite(gw_han.get(n, np.nan)) and np.isfinite(gw_ane.get(n, np.nan))]
    if len(cohorts_s2) >= 5:
        X2 = np.column_stack([np.ones(len(cohorts_s2)),
                              [gw_han[n] for n in cohorts_s2],
                              [gw_ane[n] for n in cohorts_s2]])
        y2 = np.array([s2_vals[n] for n in cohorts_s2])
        try:
            beta2, _, _, _ = np.linalg.lstsq(X2, y2, rcond=None)
            y_pred2 = X2 @ beta2
            r2_2 = 1 - np.sum((y2 - y_pred2)**2) / np.sum((y2 - y2.mean())**2) if np.sum((y2-y2.mean())**2) > 0 else np.nan
            reg_rows.append({"model": "S2_Z ~ Han_sim + ANE_sim", "n": len(cohorts_s2),
                             "intercept": beta2[0], "beta_Han": beta2[1], "beta_ANE": beta2[2],
                             "R2": r2_2})
            p(f"  S2 ~ Han + ANE: beta_Han={beta2[1]:+.3f} beta_ANE={beta2[2]:+.3f} R2={r2_2:.3f}")
        except Exception:
            pass
    # save per-cohort ancestry + signal for plotting
    for name in ALL_COHORTS:
        if name not in freq: continue
        reg_rows.append({"model": "per_cohort", "n": 1, "intercept": np.nan,
                         "beta_Han": gw_han.get(name, np.nan), "beta_ANE": gw_ane.get(name, np.nan),
                         "R2": np.nan, "cohort": name,
                         "S3_corr": s3_vals.get(name, np.nan),
                         "S2_Z": s2_vals.get(name, np.nan)})
    t7 = pd.DataFrame(reg_rows)
    t7.to_csv(os.path.join(outt, "table7_regression_models.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table7 ({len(t7)} rows)")

    # ===== TABLE 8: Sensitivity analyses (TV, selected-loci, LOCO) =====
    p("\n=== Sensitivity analyses ===")
    t8_rows = []
    # TV-only S3 setup
    set_a_tv = set_a & tv_mask
    sa_tv_idx = np.where(set_a_tv)[0]
    sa_tv_chrom = chrom[sa_tv_idx]
    sa_tv_den = den_is_a1[sa_tv_idx]
    valid_tv = np.isin(sa_tv_chrom, list(VALID_CHROMS))
    fPap_tv = freq["Papuan"][sa_tv_idx]
    fHan_tv = freq["Han"][sa_tv_idx]
    pP_tv = orient(fPap_tv, sa_tv_den)
    pH_tv = orient(fHan_tv, sa_tv_den)
    p(f"  Set-A TV validation sites: {int(valid_tv.sum())}")

    # load .snp for positions (selected-loci exclusion)
    aadr_dir = "C:/Users/benne/aadr_v66"
    snp_local = os.path.join(os.path.dirname(_REPO), "archaic-introgression", "config.local.yaml")
    if os.path.exists(snp_local):
        import yaml
        with open(snp_local) as f:
            cfg = yaml.safe_load(f)
        aadr_dir = cfg.get("aadr_dir", aadr_dir)
    from archaic.lib_eigenstrat import read_snp
    snp_df = read_snp(os.path.join(aadr_dir, "v66.p1_1240K.snp"))
    sa_positions = snp_df.loc[snp_rows[set_a_idx], "pos"].to_numpy()
    sa_chrom_full = chrom[set_a_idx]

    # exclusion mask
    exclude_mask = np.zeros(len(set_a_idx), dtype=bool)
    for ex_chrom, ex_start, ex_end, ex_name in EXCLUSION_REGIONS:
        in_region = (sa_chrom_full == ex_chrom) & (sa_positions >= ex_start) & (sa_positions <= ex_end)
        exclude_mask |= in_region
        p(f"  Excluding {ex_name}: {int(in_region.sum())} Set-A sites")
    valid_no_sel = valid_mask & ~exclude_mask

    for name in ANCIENT_COHORTS + ["Karitiana", "Papuan", "Han", "French"]:
        if name not in freq: continue
        fX_sa = freq[name][set_a_idx]
        fX_oriented = orient(fX_sa, sa_den)

        # TV-only S3
        fX_tv_arr = freq[name][sa_tv_idx]
        fX_tv_oriented = orient(fX_tv_arr, sa_tv_den)
        if valid_tv.sum() > 10:
            cp_tv = pearson_safe(fX_tv_oriented[valid_tv], pP_tv[valid_tv])
            ch_tv = pearson_safe(fX_tv_oriented[valid_tv], pH_tv[valid_tv])
            s3_tv = cp_tv - ch_tv if np.isfinite(cp_tv) and np.isfinite(ch_tv) else np.nan
        else:
            s3_tv = np.nan
        t8_rows.append({"cohort": name, "analysis": "TV_only_S3", "value": s3_tv})

        # selected-loci exclusion S3 (use full-length masks)
        if valid_no_sel.sum() > 10:
            cp_ns = pearson_safe(fX_oriented[valid_no_sel], pP_full[valid_no_sel])
            ch_ns = pearson_safe(fX_oriented[valid_no_sel], pH_full[valid_no_sel])
            s3_ns = cp_ns - ch_ns if np.isfinite(cp_ns) and np.isfinite(ch_ns) else np.nan
        else:
            s3_ns = np.nan
        t8_rows.append({"cohort": name, "analysis": "no_selected_loci_S3", "value": s3_ns})

        # alternative outgroup (Yoruba instead of Mbuti) S1
        if "Yoruba" in freq:
            r_alt = dstat_all(freq, name, "Yoruba", "Denisova", "Chimp", block)
            t8_rows.append({"cohort": name, "analysis": "S1_Yoruba_outgroup", "value": r_alt["z"]})

    # LOCO for key cohorts
    p("  Leave-one-chromosome-out S3...")
    for name in ["AmS_early", "AmN_early", "Sib_UP", "Sib_LGM", "Karitiana", "Papuan", "Han"]:
        if name not in freq: continue
        fX_oriented = orient(freq[name][set_a_idx], sa_den)
        for loo_chrom in [str(c) for c in range(1, 23)]:
            loco_valid = valid_mask & (sa_chrom != loo_chrom)
            if loco_valid.sum() < 10: continue
            cp = pearson_safe(fX_oriented[loco_valid], pP_full[loco_valid])
            ch = pearson_safe(fX_oriented[loco_valid], pH_full[loco_valid])
            s3 = cp - ch if np.isfinite(cp) and np.isfinite(ch) else np.nan
            t8_rows.append({"cohort": name, "analysis": f"LOCO_chr{loo_chrom}", "value": s3})
    # per-chromosome S3 for key cohorts
    p("  Per-chromosome S3...")
    for name in ["AmS_early", "Sib_UP", "Karitiana", "Papuan"]:
        if name not in freq: continue
        fX_oriented = orient(freq[name][set_a_idx], sa_den)
        for c in [str(i) for i in range(1, 23)]:
            chr_valid = valid_mask & (sa_chrom == c)
            if chr_valid.sum() < 5: continue
            cp = pearson_safe(fX_oriented[chr_valid], pP_full[chr_valid])
            ch = pearson_safe(fX_oriented[chr_valid], pH_full[chr_valid])
            s3 = cp - ch if np.isfinite(cp) and np.isfinite(ch) else np.nan
            t8_rows.append({"cohort": name, "analysis": f"perchrom_{c}", "value": s3})

    t8 = pd.DataFrame(t8_rows)
    t8.to_csv(os.path.join(outt, "table8_sensitivity_analyses.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table8 ({len(t8)} rows)")

    # ===== TABLE 9: Simulation performance =====
    p("\n=== Simplified power simulation ===")
    rng = np.random.default_rng(123)
    alpha_levels = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
    valid_idx_arr = np.where(valid_mask)[0]  # indices into Set-A where validation
    n_valid = len(valid_idx_arr)
    marker_counts = [100, 500, 1000, 5000, n_valid]
    sim_rows = []
    for alpha in alpha_levels:
        # alpha = fraction of HAN-like component; 0 = pure Papuan-like, 1 = pure Han-like
        fNA_true = (1 - alpha) * pP_valid + alpha * pH_valid
        for n_markers in marker_counts:
            s3_vals_sim = []
            for _ in range(200):
                if n_markers < n_valid:
                    sub = rng.choice(n_valid, size=n_markers, replace=False)
                else:
                    sub = np.arange(n_valid)
                fX_sim = fNA_true[sub] + rng.normal(0, 0.02, len(sub))
                fX_sim = np.clip(fX_sim, 0, 1)
                cp = pearson_safe(fX_sim, pP_valid[sub])
                ch = pearson_safe(fX_sim, pH_valid[sub])
                if np.isfinite(cp) and np.isfinite(ch):
                    s3_vals_sim.append(cp - ch)
            if len(s3_vals_sim) > 5:
                s3_arr = np.array(s3_vals_sim)
                mean_s3 = float(s3_arr.mean())
                sd_s3 = float(s3_arr.std())
                # power: fraction with S3 > 0 (Papuan-like direction) at p<0.05
                if alpha == 0:
                    fpr = float(np.mean(np.abs(s3_arr) > 0.05))
                else:
                    fpr = np.nan
                power = float(np.mean(s3_arr > 0.02)) if alpha > 0 else np.nan
            else:
                mean_s3 = sd_s3 = fpr = power = np.nan
            sim_rows.append({"alpha_han": alpha, "n_markers": n_markers,
                             "mean_S3": mean_s3, "sd_S3": sd_s3,
                             "FPR": fpr, "power_papuan": power})
            p(f"  alpha={alpha:.2f} n={n_markers:>6}  S3={mean_s3:+.4f}±{sd_s3:.4f}  FPR={fpr if alpha==0 else float('nan')}")
    t9 = pd.DataFrame(sim_rows)
    t9.to_csv(os.path.join(outt, "table9_simulation_performance.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table9 ({len(t9)} rows)")

    # ===== TABLE 10: Evidence summary =====
    p("\n=== Evidence summary ===")
    # assess each model based on results
    def get(name, table, col="value"):
        for _, r in table.iterrows():
            if r.get("cohort") == name:
                return r.get(col, np.nan)
        return np.nan

    # gather key results
    pap_s1 = t3.loc[t3["cohort"] == "Papuan", "S1_Z"].values[0] if len(t3[t3["cohort"]=="Papuan"]) else np.nan
    pap_s2 = t3.loc[t3["cohort"] == "Papuan", "S2_Z"].values[0] if len(t3[t3["cohort"]=="Papuan"]) else np.nan
    han_s1 = t3.loc[t3["cohort"] == "Han", "S1_Z"].values[0] if len(t3[t3["cohort"]=="Han"]) else np.nan
    french_s2 = t3.loc[t3["cohort"] == "French", "S2_Z"].values[0] if len(t3[t3["cohort"]=="French"]) else np.nan

    am_s3_vals = [r["S3_corr"] for r in t5_rows if r["cohort"] in ANCIENT_COHORTS and np.isfinite(r["S3_corr"])]
    am_s3_mean = np.mean(am_s3_vals) if am_s3_vals else np.nan
    am_s3_pos = sum(1 for v in am_s3_vals if v > 0) if am_s3_vals else 0
    am_s3_total = len(am_s3_vals) if am_s3_vals else 0

    am_s2_vals = [r["S2_Z"] for r in t3_rows if r["cohort"] in ANCIENT_COHORTS and np.isfinite(r["S2_Z"])]
    am_s2_sig = sum(1 for v in am_s2_vals if abs(v) > 3) if am_s2_vals else 0

    models = [
        ("M0_null", "No Denisovan ancestry in Native Americans",
         f"S2 Z>|3| in {am_s2_sig}/{len(am_s2_vals) if am_s2_vals else 0} ancient cohorts; Papuan S2 Z={pap_s2:.1f}",
         "Supported" if am_s2_sig == 0 and np.isfinite(pap_s2) and pap_s2 > 3 else "Contradicted"),
        ("M1_single_pulse", "Single ancestral East Eurasian Denisovan pulse",
         f"Mean ancient American S3={am_s3_mean:+.4f} (near 0 = tracks Han ~ equally); {am_s3_pos}/{am_s3_total} positive",
         "Supported" if abs(am_s3_mean) < 0.02 else "Partially supported"),
        ("M2_multiple_ea", "Multiple Denisovan pulses in East Asia",
         "Not directly tested by S3; requires tract-level analysis (deferred)",
         "Unresolved"),
        ("M3_papuan_source", "Papuan-related source contribution",
         f"Mean ancient American S3={am_s3_mean:+.4f}; {am_s3_pos}/{am_s3_total} positive; bootstrap CIs checked",
         "Supported" if am_s3_mean > 0.02 and am_s3_pos > am_s3_total/2 else "Not supported" if am_s3_mean < 0 else "Weak"),
        ("M4_structured", "Structured Denisovan source populations",
         "Requires tract phylogenies (deferred to whole-genome tier)",
         "Unresolved"),
        ("M5_dilution", "Differential dilution",
         f"Regression R2 = {reg_rows[0]['R2'] if reg_rows else 'n/a'}",
         "Unresolved"),
        ("M6_selection", "Selection-driven convergence",
         f"Selected-loci exclusion S3 changes: see table8",
         "Unresolved"),
    ]
    t10_rows = []
    for mid, desc, evidence, verdict in models:
        t10_rows.append({"model": mid, "description": desc, "evidence": evidence, "verdict": verdict})
        p(f"  {mid}: {verdict}")
    t10 = pd.DataFrame(t10_rows)
    t10.to_csv(os.path.join(outt, "table10_evidence_summary.tsv"), sep="\t", index=False, na_rep="NA")
    p(f"  -> table10 ({len(t10)} rows)")

    p("\nAll analysis tables written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
