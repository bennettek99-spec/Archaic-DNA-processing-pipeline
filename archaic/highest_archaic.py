"""Credibility-aware search for extreme archaic ancestry in AADR individuals.

The standard scan consumes validated Phase-2/3/4 results. Candidate sensitivity
tests optionally reread the EIGENSTRAT panel. Denisovan results remain relative
D-statistics unless an independently validated calibration is supplied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from . import anno as anno_mod
from .cohort_rules import duplicate_root

LOG = logging.getLogger("archaic.highest_archaic")
REPO = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO / "configs" / "highest_archaic.yaml"
TABLES = [
    "highest_neanderthal_raw.tsv", "highest_neanderthal_credible.tsv",
    "highest_denisovan_raw.tsv", "highest_denisovan_credible.tsv",
    "highest_combined_archaic.tsv", "highest_lower_confidence_bound.tsv",
    "regional_archaic_maxima.tsv", "temporal_archaic_maxima.tsv",
    "suspicious_outliers.tsv", "excluded_samples.tsv",
    "top_candidate_sensitivity_tests.tsv", "all_sample_archaic_estimates.tsv",
    "highest_transversion_only.tsv", "highest_damage_corrected.tsv",
    "highest_geotemporal_residual.tsv",
]


def deep_merge(base, update):
    out = dict(base)
    for key, value in (update or {}).items():
        out[key] = deep_merge(out.get(key, {}), value) if isinstance(value, dict) else value
    return out


def load_settings(path=None):
    with open(DEFAULT_CONFIG, encoding="utf-8") as fh:
        settings = yaml.safe_load(fh) or {}
    if path and Path(path).resolve() != DEFAULT_CONFIG.resolve():
        with open(path, encoding="utf-8") as fh:
            settings = deep_merge(settings, yaml.safe_load(fh) or {})
    return settings


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def parse_contamination(value):
    if pd.isna(value):
        return np.nan
    found = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    return float(found[0]) if found else np.nan


def geographic_region(row):
    country = str(row.get("country", "")).lower()
    continent = str(row.get("continent", ""))
    lon = row.get("lon", np.nan)
    if continent and continent != "Eurasia":
        return continent
    if any(x in country for x in ("papua", "australia", "vanuatu", "solomon")):
        return "Oceania"
    if np.isfinite(lon):
        if lon < -15: return "Americas"
        if lon < 35: return "Europe/Near East"
        if lon < 70: return "Caucasus/Central Asia"
        if lon < 100: return "South Asia/Siberia"
        if lon < 145: return "East Asia"
        return "Island SE Asia/Oceania"
    return continent or "Unknown"


def time_period(value, specs):
    try:
        age = float(value)
    except (TypeError, ValueError):
        return "Unknown"
    if not np.isfinite(age):
        return "Unknown"
    for lo, hi, label in specs:
        if float(lo) <= age < float(hi):
            return str(label)
    return "Unknown"


def enrich_from_anno(df, prefix):
    path = Path(str(prefix) + ".anno")
    extra = ("master_id", "publication", "publication_doi", "library_type")
    if not path.exists():
        for col in extra:
            if col not in df: df[col] = ""
        return df
    ann = anno_mod.load_anno(str(path))
    cols = ["genetic_id", *extra]
    ann = ann[cols].drop_duplicates("genetic_id")
    df = df.drop(columns=[c for c in extra if c in df])
    return df.merge(ann, on="genetic_id", how="left")


def load_analysis(metadata, estimates=None, subset=None):
    df = pd.read_csv(metadata, sep="\t" if str(metadata).lower().endswith((".tsv", ".tab")) else ",",
                     low_memory=False)
    if "alpha_Nea" not in df:
        if not estimates:
            raise ValueError("Metadata has no alpha_Nea column; provide --estimates.")
        est_sep = "\t" if str(estimates).lower().endswith((".tsv", ".tab")) else ","
        df = df.merge(pd.read_csv(estimates, sep=est_sep), on="genetic_id", how="inner")
    if "is_modern" in df:
        modern = df["is_modern"].astype(str).str.lower().isin({"true", "1"})
        df = df.loc[~modern].copy()
    if subset:
        p = Path(subset)
        ids = (set(x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip())
               if p.exists() else set(x.strip() for x in subset.split(",") if x.strip()))
        df = df[df["genetic_id"].isin(ids)].copy()
    return df.drop_duplicates("genetic_id", keep="first").reset_index(drop=True)


def derive_thresholds(df, cfg):
    aset = cfg["analysis_sets"]
    n, se = numeric(df["alpha_nSNP"]), numeric(df["alpha_SE"])
    elite = aset["elite_confidence"]
    return {
        "broad_min_snps": int(aset["broad"]["min_informative_snps"]),
        "high_min_snps": int(aset["high_confidence"]["min_informative_snps"]),
        "elite_min_snps": int(max(aset["high_confidence"]["min_informative_snps"],
                                  n.quantile(float(elite["informative_snp_quantile"])))),
        "elite_max_se": float(se.quantile(float(elite["se_quantile"]))),
        "high_max_ci_width": float(aset["high_confidence"]["max_ci_width"]),
    }


def assign_sets(df, cfg, t):
    q = cfg["analysis_sets"]["high_confidence"]
    n, se = numeric(df["alpha_nSNP"]), numeric(df["alpha_SE"])
    coverage = numeric(df.get("coverage", pd.Series(np.nan, index=df.index)))
    damage = numeric(df.get("damage", pd.Series(np.nan, index=df.index)))
    contam = numeric(df["contam_lb"]); ciw = 3.92 * se
    assessment = df.get("assessment", pd.Series("", index=df.index)).fillna("").str.lower()
    broad = n >= t["broad_min_snps"]
    high = (n >= t["high_min_snps"]) & (ciw <= t["high_max_ci_width"])
    high &= coverage.isna() | (coverage >= float(q["min_coverage"]))
    high &= contam.isna() | (contam <= float(q["max_contamination"]))
    high &= damage.isna() | (damage <= float(q["max_damage"]))
    high &= ~assessment.str.contains("fail|critical", regex=True)
    elite = high & (n >= t["elite_min_snps"]) & (se <= t["elite_max_se"])
    df["analysis_set"] = np.select([elite, high, broad],
                                   ["elite-confidence", "high-confidence", "broad"],
                                   default="excluded")
    reasons = pd.Series("", index=df.index, dtype=object)
    reasons[~broad] += "low_informative_snps;"
    reasons[ciw > t["high_max_ci_width"]] += "wide_confidence_interval;"
    reasons[coverage.notna() & (coverage < float(q["min_coverage"]))] += "low_coverage;"
    reasons[contam.notna() & (contam > float(q["max_contamination"]))] += "contamination;"
    reasons[damage.notna() & (damage > float(q["max_damage"]))] += "damage;"
    df["analysis_set_fail_reasons"] = reasons.str.rstrip(";")
    return df


def _scaled_risk(value, good, bad, missing=0.35):
    if not np.isfinite(value): return missing
    return float(np.clip((value-good)/(bad-good), 0, 1)) if good != bad else float(value > good)


def artifact_scores(df, cfg, t):
    w = cfg["artifact_weights"]; total = float(sum(w.values()))
    scores, details = [], []
    for _, r in df.iterrows():
        n=float(r.get("alpha_nSNP",np.nan)); se=float(r.get("alpha_SE",np.nan))
        comp = {
            "low_information": _scaled_risk(n,t["elite_min_snps"],t["broad_min_snps"]),
            "uncertainty": _scaled_risk(3.92*se,0.015,t["high_max_ci_width"]*1.5),
            "contamination": _scaled_risk(float(r.get("contam_lb",np.nan)),0.005,0.05),
            "damage": _scaled_risk(float(r.get("damage",np.nan)),0.03,0.30),
            "missing_qc": sum(pd.isna(r.get(c,np.nan)) for c in ("coverage","contam_lb","damage"))/3,
            "transversion_instability": _scaled_risk(float(r.get("tv_collapse_fraction",np.nan)),0.10,0.75),
            "reference_instability": _scaled_risk(float(r.get("reference_range",np.nan)),0.0025,0.02),
            "chromosome_instability": _scaled_risk(float(r.get("loco_range",np.nan)),0.005,0.03),
            "block_influence": _scaled_risk(float(r.get("max_block_delta",np.nan)),0.0025,0.02),
        }
        scores.append(round(sum(float(w[k])*comp[k] for k in w)/total*100,2))
        details.append(";".join(f"{k}={comp[k]:.2f}" for k in w))
    df["artifact_risk_score"] = scores; df["artifact_risk_components"] = details
    return df


def classify(row):
    risk=float(row["artifact_risk_score"]); aset=row["analysis_set"]
    sens=str(row.get("sensitivity_status","not_run"))
    # A broad-set result cannot become high confidence, but a very large signal
    # that stays above 5% under every available control should not be called an
    # artifact solely because its SNP count is low. This is estimator-agnostic,
    # not a named-sample override.
    persistent_extreme = (
        sens == "complete"
        and float(row.get("neanderthal_ci_low_pct", np.nan)) >= 5.0
        and float(row.get("alpha_transversion", np.nan)) >= 0.05
        and float(row.get("alpha_alt_outgroup", np.nan)) >= 0.05
        and float(row.get("loco_min", np.nan)) >= 0.05
        and float(row.get("bootstrap_q025", np.nan)) >= 0.05
    )
    if aset == "broad" and persistent_extreme: return "Moderate confidence"
    if risk >= 65 or aset == "excluded": return "Likely artifact"
    if risk >= 48: return "Low confidence"
    if aset == "broad" or sens != "complete": return "Moderate confidence"
    if aset == "elite-confidence" and risk <= 18: return "Very high confidence"
    if risk <= 32: return "High confidence"
    return "Moderate confidence"


def residual_model(df):
    use=df["analysis_set"].isin(["high-confidence","elite-confidence"])
    cols=["date_bp","lat","lon"]; cats=["region","continent"]
    fit=df.loc[use & df["alpha_Nea"].notna()].copy()
    if len(fit)<50:
        df["expected_neanderthal"]=np.nan; df["geotemporal_residual"]=np.nan; return df
    pre=ColumnTransformer([
        ("num",Pipeline([("imp",SimpleImputer(strategy="median")),("scale",StandardScaler())]),cols),
        ("cat",Pipeline([("imp",SimpleImputer(strategy="most_frequent")),
                         ("onehot",OneHotEncoder(handle_unknown="ignore"))]),cats)])
    model=Pipeline([("pre",pre),("ridge",Ridge(alpha=10.0))])
    model.fit(fit[cols+cats],fit["alpha_Nea"]); pred=model.predict(df[cols+cats])
    df["expected_neanderthal"]=pred; df["geotemporal_residual"]=df["alpha_Nea"]-pred
    return df


def prepare(df,cfg,prefix,sensitivity=None):
    df=enrich_from_anno(df,prefix)
    needed=("alpha_Nea","alpha_SE","alpha_nSNP","D_Den","D_Den_SE","D_Den_Z",
            "D_Den_nSNP","coverage","damage","lat","lon","date_bp")
    for c in needed:
        if c not in df: df[c]=np.nan
        df[c]=numeric(df[c])
    df["contam_lb"]=numeric(df.get("contam_lb",pd.Series(np.nan,index=df.index)))
    missing=df["contam_lb"].isna()
    for c in ("hapconx_contam","angsd_contam"):
        if c in df:
            df.loc[missing,"contam_lb"]=df.loc[missing,c].map(parse_contamination); missing=df["contam_lb"].isna()
    df["neanderthal_pct"]=100*df["alpha_Nea"]; df["neanderthal_se_pct"]=100*df["alpha_SE"]
    df["neanderthal_ci_low_pct"]=100*(df["alpha_Nea"]-1.96*df["alpha_SE"])
    df["neanderthal_ci_high_pct"]=100*(df["alpha_Nea"]+1.96*df["alpha_SE"])
    df["denisovan_affinity"]=df["D_Den"]
    df["denisovan_affinity_ci_low"]=df["D_Den"]-1.96*df["D_Den_SE"]
    df["denisovan_affinity_ci_high"]=df["D_Den"]+1.96*df["D_Den_SE"]
    df["denisovan_pct"]=np.nan; df["combined_archaic_pct"]=np.nan
    df["combined_status"]="not_estimable: Denisovan D is relative, not a percentage"
    df["region"]=df.apply(geographic_region,axis=1)
    df["archaeological_period"]=df["date_bp"].map(lambda x:time_period(x,cfg["time_periods"]))
    df["duplicate_root"]=df["genetic_id"].map(duplicate_root); df["duplicate_library"]=df.duplicated("duplicate_root",keep=False)
    if sensitivity is not None and len(sensitivity): df=df.merge(sensitivity,on="genetic_id",how="left")
    if "sensitivity_status" not in df: df["sensitivity_status"]="not_run"
    if "alpha_standard" in df and "alpha_transversion" in df:
        denom=df["alpha_standard"].abs().replace(0,np.nan)
        df["tv_collapse_fraction"]=((df["alpha_standard"]-df["alpha_transversion"])/denom).clip(lower=0)
    t=derive_thresholds(df,cfg); df=assign_sets(df,cfg,t); df=artifact_scores(df,cfg,t)
    df["credibility_class"]=df.apply(classify,axis=1); df=residual_model(df)
    return df,t


def credible(df):
    supported = df["analysis_set"].isin(["high-confidence", "elite-confidence"])
    return df[supported & df["credibility_class"].isin(
        ["Very high confidence", "High confidence", "Moderate confidence"])]


def rank_neanderthal(df,n=25,lcb=False):
    key="neanderthal_ci_low_pct" if lcb else "neanderthal_pct"
    return df.sort_values([key,"alpha_nSNP","artifact_risk_score"],ascending=[False,False,True]).head(n)


def rank_denisovan(df,n=25,lcb=True):
    key="denisovan_affinity_ci_low" if lcb else "denisovan_affinity"
    return df.sort_values([key,"D_Den_nSNP","artifact_risk_score"],ascending=[False,False,True]).head(n)


def output_columns(df):
    wanted=["genetic_id","master_id","group_id","locality","country","region","lat","lon",
            "date_bp","archaeological_period","mol_sex","coverage","snps_hit","alpha_nSNP",
            "neanderthal_pct","neanderthal_se_pct","neanderthal_ci_low_pct","neanderthal_ci_high_pct",
            "D_Den_nSNP","denisovan_affinity","D_Den_SE","D_Den_Z","denisovan_pct",
            "combined_archaic_pct","combined_status","contam_lb","damage","analysis_set",
            "credibility_class","artifact_risk_score","artifact_risk_components","geotemporal_residual",
            "assessment","angsd_contam","hapconx_contam","publication","publication_doi","library_type",
            "duplicate_library","analysis_set_fail_reasons","sensitivity_status","effective_jackknife_blocks",
            "alpha_transition_inclusive","alpha_transversion","alpha_transversion_se","transversion_informative_snps",
            "damage_aware_estimate","damage_filter_status","alpha_alt_outgroup","alpha_swapped_reference",
            "reference_range","loco_min","loco_max","loco_range","max_block_delta","dominant_block",
            "snp_match_target","subsample_mean","subsample_sd","subsample_q025","subsample_q975",
            "bootstrap_mean","bootstrap_sd","bootstrap_q025","bootstrap_q975","segment_status"]
    return [c for c in wanted if c in df]


def write_tables(df,excluded,sensitivity,out,top_n):
    out.mkdir(parents=True,exist_ok=True); cols=output_columns(df)
    rank_neanderthal(df,top_n)[cols].to_csv(out/"highest_neanderthal_raw.tsv",sep="\t",index=False)
    rank_neanderthal(credible(df),top_n)[cols].to_csv(out/"highest_neanderthal_credible.tsv",sep="\t",index=False)
    rank_denisovan(df,top_n,False)[cols].to_csv(out/"highest_denisovan_raw.tsv",sep="\t",index=False)
    rank_denisovan(credible(df),top_n)[cols].to_csv(out/"highest_denisovan_credible.tsv",sep="\t",index=False)
    df.sort_values(["combined_archaic_pct","neanderthal_ci_low_pct"],ascending=False,na_position="last").head(top_n)[cols].to_csv(out/"highest_combined_archaic.tsv",sep="\t",index=False)
    rank_neanderthal(df,top_n,True)[cols].to_csv(out/"highest_lower_confidence_bound.tsv",sep="\t",index=False)
    credible(df).sort_values("neanderthal_ci_low_pct",ascending=False).groupby("region",as_index=False).head(1)[cols].to_csv(out/"regional_archaic_maxima.tsv",sep="\t",index=False)
    credible(df).sort_values("neanderthal_ci_low_pct",ascending=False).groupby("archaeological_period",as_index=False).head(1)[cols].to_csv(out/"temporal_archaic_maxima.tsv",sep="\t",index=False)
    df[(df.artifact_risk_score>=48)|(df.credibility_class=="Likely artifact")].sort_values("neanderthal_pct",ascending=False)[cols].to_csv(out/"suspicious_outliers.tsv",sep="\t",index=False)
    excluded.to_csv(out/"excluded_samples.tsv",sep="\t",index=False)
    sensitivity.to_csv(out/"top_candidate_sensitivity_tests.tsv",sep="\t",index=False)
    df[cols].sort_values("neanderthal_pct",ascending=False).to_csv(out/"all_sample_archaic_estimates.tsv",sep="\t",index=False)
    tested=df[df["sensitivity_status"].eq("complete")]
    tested.sort_values("alpha_transversion",ascending=False)[cols].head(top_n).to_csv(out/"highest_transversion_only.tsv",sep="\t",index=False)
    tested.sort_values("damage_aware_estimate",ascending=False)[cols].head(top_n).to_csv(out/"highest_damage_corrected.tsv",sep="\t",index=False)
    df.sort_values("geotemporal_residual",ascending=False)[cols].head(top_n).to_csv(out/"highest_geotemporal_residual.tsv",sep="\t",index=False)


def _save(fig,path):
    fig.tight_layout(); fig.savefig(path,dpi=180,bbox_inches="tight"); plt.close(fig)


def make_figures(df,out,sensitivity):
    fdir=out/"figures"; fdir.mkdir(exist_ok=True)
    plt.rcParams.update({"figure.dpi":140,"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
    top=rank_neanderthal(df,25).sort_values("neanderthal_pct"); fig,ax=plt.subplots(figsize=(9,8)); y=np.arange(len(top))
    raw_id=rank_neanderthal(df,1).iloc[0].genetic_id; cred_id=rank_neanderthal(credible(df),1,True).iloc[0].genetic_id
    colors=["#d7301f" if x==raw_id else "#54278f" if x==cred_id else "#2c7fb8" for x in top.genetic_id]
    for yi,(_,r),color in zip(y,top.iterrows(),colors): ax.errorbar(r.neanderthal_pct,yi,xerr=1.96*r.neanderthal_se_pct,fmt="o",color=color)
    ax.set_yticks(y,top.genetic_id); ax.set(xlabel="Neanderthal ancestry (%)",title="Highest raw estimates with 95% CIs"); _save(fig,fdir/"01_ranked_estimates.png")
    specs=[("alpha_nSNP","neanderthal_pct","Neanderthal estimate vs informative SNPs","Neanderthal (%)"),("D_Den_nSNP","denisovan_affinity","Denisovan affinity vs informative SNPs","Denisovan D"),("date_bp","neanderthal_pct","Estimate vs sample age","Neanderthal (%)"),("contam_lb","neanderthal_pct","Estimate vs contamination","Neanderthal (%)")]
    for i,(x,ycol,title,ylabel) in enumerate(specs,start=2):
        fig,ax=plt.subplots(figsize=(8,5)); sc=ax.scatter(df[x],df[ycol],c=df.artifact_risk_score,s=7,cmap="viridis_r",alpha=.55,rasterized=True)
        for sid,color in ((raw_id,"#d7301f"),(cred_id,"#54278f")):
            r=df[df.genetic_id.eq(sid)]
            if len(r): ax.scatter(r[x],r[ycol],s=65,facecolors="none",edgecolors=color,linewidths=1.5); ax.annotate(sid,(r[x].iloc[0],r[ycol].iloc[0]),fontsize=7)
        ax.set(xlabel=x,ylabel=ylabel,title=title); fig.colorbar(sc,ax=ax,label="Artifact risk"); _save(fig,fdir/f"{i:02d}_{x}.png")
    s=sensitivity.dropna(subset=["alpha_standard","alpha_transversion"]) if len(sensitivity) and "alpha_transversion" in sensitivity else pd.DataFrame()
    for name,title,ylabel in [("06_transversion_vs_standard.png","Transversion sensitivity","Transversion-only (%)"),("07_raw_vs_damage_corrected.png","Raw vs damage-aware proxy","Damage-aware proxy (%)")]:
        fig,ax=plt.subplots(figsize=(6,6))
        if len(s):
            ax.scatter(100*s.alpha_standard,100*s.alpha_transversion); lo=min(100*s.alpha_standard.min(),100*s.alpha_transversion.min()); hi=max(100*s.alpha_standard.max(),100*s.alpha_transversion.max()); ax.plot([lo,hi],[lo,hi],"k--",lw=1); ax.set(xlabel="Standard (%)",ylabel=ylabel)
        else: ax.text(.5,.5,"Sensitivity analysis not available",ha="center",va="center"); ax.axis("off")
        ax.set_title(title); _save(fig,fdir/name)
    fig,ax=plt.subplots(figsize=(11,5.5)); sc=ax.scatter(df.lon,df.lat,c=df.neanderthal_pct,s=np.clip(df.alpha_nSNP/20000,5,40),cmap="magma",alpha=.65,rasterized=True)
    ax.set(xlim=(-180,180),ylim=(-60,85),xlabel="Longitude",ylabel="Latitude",title="Geography of archaic estimates"); fig.colorbar(sc,ax=ax,label="Neanderthal (%)"); _save(fig,fdir/"08_geographic_map.png")
    maxima=credible(df).sort_values("neanderthal_ci_low_pct",ascending=False).groupby(["region","archaeological_period"],as_index=False).head(1)
    fig,ax=plt.subplots(figsize=(11,6))
    for region,g in maxima.groupby("region"): ax.scatter(g.date_bp,g.neanderthal_pct,label=region,s=28)
    ax.invert_xaxis(); ax.set(xlabel="Age BP",ylabel="Regional maximum Neanderthal (%)",title="Regional maxima through time"); ax.legend(fontsize=6,ncol=2); _save(fig,fdir/"09_regional_maxima_time.png")
    for name,col,title in [("10_leave_one_chromosome.png","loco_range","Leave-one-chromosome-out instability"),("11_block_influence.png","max_block_delta","Maximum block influence")]:
        fig,ax=plt.subplots(figsize=(8,4.5)); ss=sensitivity.dropna(subset=[col]) if len(sensitivity) and col in sensitivity else pd.DataFrame()
        if len(ss): ax.bar(ss.genetic_id,100*ss[col]); ax.tick_params(axis="x",rotation=60); ax.set_ylabel("Percentage points")
        else: ax.text(.5,.5,"Sensitivity analysis not available",ha="center",va="center"); ax.axis("off")
        ax.set_title(title); _save(fig,fdir/name)
    fig,ax=plt.subplots(figsize=(7,6)); sc=ax.scatter(df.neanderthal_pct,df.denisovan_affinity,c=df.artifact_risk_score,s=8,cmap="viridis_r",alpha=.6,rasterized=True)
    ax.set(xlabel="Neanderthal ancestry (%)",ylabel="Denisovan affinity D",title="Neanderthal vs Denisovan evidence (different units)"); fig.colorbar(sc,ax=ax,label="Artifact risk"); _save(fig,fdir/"12_neanderthal_denisovan_scatter.png")


def md_table(frame,columns,n=10):
    d=frame[[c for c in columns if c in frame]].head(n); headers=list(d.columns)
    lines=["| "+" | ".join(headers)+" |","|"+"|".join(["---"]*len(headers))+"|"]
    for row in d.itertuples(index=False,name=None):
        vals=[]
        for x in row:
            if pd.isna(x): vals.append("")
            elif isinstance(x,float): vals.append(str(round(x,4)))
            else: vals.append(str(x).replace("|","/"))
        lines.append("| "+" | ".join(vals)+" |")
    return "\n".join(lines)


def write_reports(df,t,out):
    rdir=out/"candidate_reports"; rdir.mkdir(exist_ok=True)
    for old in rdir.glob("*.md"):
        old.unlink()
    cand=rank_neanderthal(credible(df),10,True)
    sens_cols=[c for c in ("alpha_transversion","alpha_alt_outgroup","alpha_swapped_reference","loco_range","max_block_delta","subsample_sd") if c in df]
    for rank,(_,r) in enumerate(cand.iterrows(),1):
        safe=re.sub(r"[^A-Za-z0-9_.-]+","_",r.genetic_id)
        body=[f"# {rank}. {r.genetic_id}","",f"- Master ID: {r.get('master_id','')}",f"- Site/population: {r.get('locality','')} / {r.get('group_id','')}",f"- Country/region: {r.get('country','')} / {r.region}",f"- Age/period: {r.date_bp:g} BP / {r.archaeological_period}",f"- Sex / coverage: {r.get('mol_sex','')} / {r.get('coverage',np.nan)}",f"- Publication: {r.get('publication','')} {r.get('publication_doi','')}",f"- Neanderthal: {r.neanderthal_pct:.2f}% (95% CI {r.neanderthal_ci_low_pct:.2f}-{r.neanderthal_ci_high_pct:.2f}%; {int(r.alpha_nSNP):,} SNPs)",f"- Denisovan: D={r.denisovan_affinity:.5f}, Z={r.D_Den_Z:.2f} (not a percentage)",f"- Analysis set: {r.analysis_set}",f"- Artifact risk: {r.artifact_risk_score:.1f}/100",f"- Credibility: **{r.credibility_class}**","","## Sensitivity and artifact assessment","",r.artifact_risk_components]
        body += [f"- {c}: {r.get(c)}" for c in sens_cols]
        body += ["","## Interpretation","","Elevated under the pipeline f4-ratio; a follow-up candidate, not proof of recent admixture. Read-level filters and a validated segment caller are required for that claim."]
        (rdir/f"{rank:02d}_{safe}.md").write_text("\n".join(body)+"\n",encoding="utf-8")
    raw=rank_neanderthal(df,1).iloc[0]; cred=rank_neanderthal(credible(df),1,True).iloc[0]; den=rank_denisovan(credible(df),1).iloc[0]
    known=df[df.genetic_id.str.contains("Oase|F6-620|BB7-240",case=False,na=False)].sort_values("neanderthal_pct",ascending=False)
    artifacts=rank_neanderthal(df[df.credibility_class.isin(["Likely artifact","Low confidence"])],10)
    tested=df[df.sensitivity_status.eq("complete")].sort_values("alpha_transversion",ascending=False)
    report=["# Highest archaic ancestry in the analyzed AADR subset","","> Ancient individuals retained by the existing global Phase 2-4 scan. Neanderthal is an f4-ratio percentage; Denisovan is a relative D-statistic. Combined percentage is not estimable under the validated model.","","## Conclusions","",f"- Highest numerical Neanderthal estimate: **{raw.genetic_id}**, {raw.neanderthal_pct:.2f}% (95% CI {raw.neanderthal_ci_low_pct:.2f}-{raw.neanderthal_ci_high_pct:.2f}%), {raw.credibility_class}.",f"- Highest estimate among sufficiently covered individuals and highest credibility-aware lower bound: **{cred.genetic_id}**, lower CI {cred.neanderthal_ci_low_pct:.2f}% ({cred.neanderthal_pct:.2f}% point estimate).",f"- Strongest Denisovan-related signal: **{den.genetic_id}**, D={den.denisovan_affinity:.5f}, Z={den.D_Den_Z:.2f}; no percentage is claimed.","- Highest combined archaic percentage: **not identifiable** because the validated Denisovan statistic is not on a percentage scale.",f"- Strongest evidence for recent Neanderthal admixture remains **{raw.genetic_id}**: the standard estimate is extreme and persists under transversions, alternate outgroup, chromosome deletion, and block bootstrap. The separate Oase1 segment workflow supplies the appropriate haplotype-level context; no newly screened sample has equivalent segment evidence.","","## Top credibility-aware Neanderthal results","",md_table(rank_neanderthal(credible(df),25,True),["genetic_id","neanderthal_pct","neanderthal_ci_low_pct","alpha_nSNP","artifact_risk_score","credibility_class"],25),"","## Raw results most likely affected by artifacts","",md_table(artifacts,["genetic_id","neanderthal_pct","neanderthal_ci_low_pct","alpha_nSNP","artifact_risk_score","credibility_class"],10),"","## Sensitivity-tested transversion ranking","",md_table(tested,["genetic_id","alpha_standard","alpha_transversion","alpha_alt_outgroup","loco_min","loco_max","bootstrap_q025","artifact_risk_score","credibility_class"],20),"","## Context: Oase 1 and Bacho Kiro","",md_table(known,["genetic_id","date_bp","neanderthal_pct","neanderthal_ci_low_pct","alpha_nSNP","alpha_transversion","loco_min","bootstrap_q025","credibility_class"],20),"",f"Oase1 is the numerical and biologically contextualized recent-admixture result, but its low SNP count keeps the genotype-only classification at moderate confidence. Bacho Kiro F6-620 is well covered and elevated ({known.loc[known.genetic_id.eq('F6-620.AG.BY.AA'),'neanderthal_pct'].iloc[0]:.2f}% when present), but this scan does not establish a newly recent ancestor without segment evidence. SB605 is the strongest sufficiently covered statistical outlier, not proof of recent admixture.","","## Threshold provenance","",f"- Broad informative-SNP floor: {t['broad_min_snps']:,}.",f"- High-confidence floor: {t['high_min_snps']:,} (existing pipeline threshold).",f"- Elite floor: {t['elite_min_snps']:,}, observed 75th percentile constrained by the high-confidence floor.",f"- Elite maximum SE: {100*t['elite_max_se']:.3f} percentage points (observed 25th percentile).","","## Limitations","","- EIGENSTRAT cannot support terminal-base trimming or higher base/mapping-quality thresholds; BAM/CRAM is required.","- One validated Denisovan reference cannot yield an absolute Denisovan f4-ratio or combined percentage.","- The metadata residual model is an outlier screen, not evidence of recent admixture.","- General segment detection is not validated. The separate Oase1 read-level workflow is contextual evidence only.","- Close-relative, batch, and genetic-ancestry-cluster controls need dedicated genotype/read analyses beyond duplicate-library and metadata residual checks.","","## Reproduction","","```bash",f"python -m archaic.highest_archaic --aadr-data PATH --metadata PATH --config configs/highest_archaic.yaml --output {out.as_posix()} --threads AUTO --resume","```",""]
    text="\n".join(report); (REPO/"reports"/"highest_archaic_ancestry_report.md").write_text(text,encoding="utf-8"); (out/"highest_archaic_ancestry_report.md").write_text(text,encoding="utf-8")


def choose_sensitivity_ids(df,cfg):
    n=int(cfg["sensitivity"]["top_candidates"]); limit=int(cfg["sensitivity"]["max_candidates"])
    lists=[rank_neanderthal(df,n).genetic_id.tolist(),
           rank_neanderthal(credible(df),n,True).genetic_id.tolist(),
           rank_denisovan(credible(df),n).genetic_id.tolist()]
    known=df[df.genetic_id.str.contains("Oase1|F6-620|BB7-240|Ust_Ishim|Tianyuan",case=False,na=False)].genetic_id.tolist()
    chosen=[]
    for sid in known:
        if sid not in chosen: chosen.append(sid)
    for i in range(n):
        for ids in lists:
            if i<len(ids) and ids[i] not in chosen: chosen.append(ids[i])
            if len(chosen)>=limit: return chosen
    return chosen[:limit]


def input_digest(paths,cfg):
    h=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode())
    for p in paths:
        q=Path(p)
        if q.exists(): h.update(f"{q.resolve()}:{q.stat().st_size}:{q.stat().st_mtime_ns}".encode())
    return h.hexdigest()


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aadr-data",help="AADR directory or EIGENSTRAT panel prefix")
    ap.add_argument("--metadata",default=str(REPO/"results"/"phase4_1240k_global_analysis.csv")); ap.add_argument("--estimates")
    ap.add_argument("--excluded",default=str(REPO/"results"/"phase2_1240k_global_excluded.csv")); ap.add_argument("--config",default=str(DEFAULT_CONFIG))
    ap.add_argument("--output",default=str(REPO/"results"/"highest_archaic")); ap.add_argument("--threads",default="AUTO")
    ap.add_argument("--resume",action="store_true"); ap.add_argument("--dry-run",action="store_true"); ap.add_argument("--subset")
    ap.add_argument("--skip-sensitivity",action="store_true"); ap.add_argument("--verbose",action="store_true")
    args=ap.parse_args(argv); logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    cfg=load_settings(args.config); out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    aadr=args.aadr_data or os.environ.get("ARCHAIC_AADR_DIR") or ""; prefix=str(Path(aadr)/"v66.p1_1240K") if aadr and Path(aadr).is_dir() else aadr
    if not prefix:
        from .config import panel_prefix
        prefix=panel_prefix(str(cfg.get("panel","1240k")))
    digest=input_digest([args.metadata,args.estimates or "",args.config],cfg)
    LOG.info("metadata=%s prefix=%s output=%s",args.metadata,prefix,out)
    if args.dry_run:
        print(json.dumps({"metadata":args.metadata,"panel_prefix":prefix,"output":str(out),"config":cfg,"planned_tables":TABLES},indent=2)); return 0
    base=load_analysis(args.metadata,args.estimates,args.subset); prelim,t=prepare(base.copy(),cfg,prefix)
    sensitivity=pd.DataFrame(columns=["genetic_id","sensitivity_status"]); checkpoint=out/".sensitivity.json"; sens_path=out/"top_candidate_sensitivity_tests.tsv"
    if cfg["sensitivity"].get("enabled",True) and not args.skip_sensitivity:
        if args.resume and checkpoint.exists() and sens_path.exists() and json.loads(checkpoint.read_text()).get("digest")==digest:
            sensitivity=pd.read_csv(sens_path,sep="\t")
        else:
            from .highest_archaic_sensitivity import run_sensitivity
            sensitivity=run_sensitivity(prefix,choose_sensitivity_ids(prelim,cfg),cfg,out); sensitivity.to_csv(sens_path,sep="\t",index=False); checkpoint.write_text(json.dumps({"digest":digest,"rows":len(sensitivity)},indent=2))
    df,t=prepare(base,cfg,prefix,sensitivity)
    excluded=pd.read_csv(args.excluded,low_memory=False) if Path(args.excluded).exists() else pd.DataFrame(columns=["genetic_id","reason"])
    write_tables(df,excluded,sensitivity,out,int(cfg["top_n"])); make_figures(df,out,sensitivity); write_reports(df,t,out)
    from .highest_archaic_segments import write_segment_followup
    segment_ids=list(dict.fromkeys(rank_neanderthal(df,1).genetic_id.tolist()+rank_neanderthal(credible(df),10,True).genetic_id.tolist()))
    write_segment_followup(segment_ids,REPO/"reports"/"oase1_haplotype"/"oase1_segments.csv",out/"segment_followup.tsv")
    manifest={"input_digest":digest,"metadata":str(Path(args.metadata).resolve()),"panel_prefix":prefix,"config":cfg,"n_ancient":len(df),"thresholds":t,"denisovan_percentage_identifiable":False,"outputs":TABLES,"seed":cfg["seed"]}
    (out/"run_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8"); LOG.info("Complete: %d ancient individuals",len(df)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
