"""Genome-wide sensitivity controls for selected highest-archaic candidates."""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .panel import Panel
from .refs import REF_SAMPLES
from . import snp_filters
from . import stats as st

LOG = logging.getLogger("archaic.highest_archaic.sensitivity")


def _ratio(num, den, block, n_blocks, mask=None):
    if mask is not None:
        num, den, block = num[mask], den[mask], block[mask]
    return st.jackknife_ratio(num, den, block, n_blocks)


def _leave_one_group(num, den, groups):
    ok = np.isfinite(num) & np.isfinite(den)
    tn, td = np.nansum(np.where(ok, num, 0.0)), np.nansum(np.where(ok, den, 0.0))
    full = tn / td if td else np.nan
    rows = []
    for group in pd.unique(groups):
        take = ok & (groups == group)
        gn, gd = np.nansum(num[take]), np.nansum(den[take])
        estimate = (tn-gn)/(td-gd) if td != gd else np.nan
        rows.append((group, estimate, estimate-full))
    return full, rows


def _subsample(num, den, target, replicates, rng):
    ok = np.flatnonzero(np.isfinite(num) & np.isfinite(den))
    if not len(ok): return np.full(replicates, np.nan)
    target = min(int(target), len(ok)); out = np.empty(replicates)
    for i in range(replicates):
        take = rng.choice(ok, target, replace=False)
        d = den[take].sum(); out[i] = num[take].sum()/d if d else np.nan
    return out


def _block_bootstrap(num, den, block, n_blocks, replicates, rng):
    ok = np.isfinite(num) & np.isfinite(den)
    bn = np.bincount(block, weights=np.where(ok,num,0.0), minlength=n_blocks)
    bd = np.bincount(block, weights=np.where(ok,den,0.0), minlength=n_blocks)
    present = np.flatnonzero(np.bincount(block,weights=ok.astype(float),minlength=n_blocks)>0)
    if not len(present): return np.full(replicates,np.nan)
    out=np.empty(replicates)
    for i in range(replicates):
        take=rng.choice(present,len(present),replace=True); d=bd[take].sum(); out[i]=bn[take].sum()/d if d else np.nan
    return out


def _summary(values, prefix):
    values=np.asarray(values,dtype=float); values=values[np.isfinite(values)]
    if not len(values):
        return {f"{prefix}_mean":np.nan,f"{prefix}_sd":np.nan,f"{prefix}_q025":np.nan,f"{prefix}_q975":np.nan}
    return {f"{prefix}_mean":float(np.mean(values)),f"{prefix}_sd":float(np.std(values,ddof=1)) if len(values)>1 else 0.0,
            f"{prefix}_q025":float(np.quantile(values,.025)),f"{prefix}_q975":float(np.quantile(values,.975))}


def run_sensitivity(prefix, ids, cfg, output):
    """Recompute top candidates directly from genotypes under robust controls."""
    panel_name=str(cfg.get("panel","1240k")); n_blocks=int(cfg["n_jackknife_blocks"])
    scfg=cfg["sensitivity"]; seed=int(cfg["seed"]); outdir=Path(output)
    LOG.info("Loading %s panel for %d candidate sensitivity tests",panel_name,len(ids))
    panel=Panel(prefix,autosomes_only=True)
    refs=REF_SAMPLES.get(panel_name,REF_SAMPLES["1240k"])
    needed={k:refs[k] for k in ("Altai","Vindija","Denisova","Chimp","Mbuti","Yoruba")}
    freq,info=panel.frequencies(needed)
    present=[sid for sid in ids if sid in panel._id_to_col]; missing=[sid for sid in ids if sid not in panel._id_to_col]
    cols=np.array([panel._id_to_col[sid] for sid in present],dtype=np.int64)
    G=panel.pg.read(panel.snp_rows,cols); X=G.astype(np.float64)/2.0; X[G<0]=np.nan
    snp=panel.snp.iloc[panel.snp_rows]
    chrom=snp["chrom"].astype(str).to_numpy(); tv=np.asarray(snp_filters.transversion_mask(snp),dtype=bool)
    block=st.assign_blocks(panel.n_snp,n_blocks)
    palt,pvin,pchi,pmb,pyo=freq["Altai"],freq["Vindija"],freq["Chimp"],freq["Mbuti"],freq["Yoruba"]
    axis=palt-pchi; den_standard=axis*(pvin-pmb); axis_swap=pvin-pchi; den_swap=axis_swap*(palt-pmb)
    # Match all selected samples to the lowest callable count, with the configured
    # floor acting only when every candidate exceeds it.
    counts=[]
    for j in range(X.shape[1]):
        counts.append(int((np.isfinite(axis*(X[:,j]-pmb)) & np.isfinite(den_standard)).sum()))
    target=min(counts) if counts else int(scfg["common_snp_floor"])
    if target>=int(scfg["common_snp_floor"]): target=int(scfg["common_snp_floor"])
    rows=[]; chrom_rows=[]; block_rows=[]
    for j,sid in enumerate(present):
        rng=np.random.default_rng(seed+j); px=X[:,j]
        num=axis*(px-pmb); standard=_ratio(num,den_standard,block,n_blocks)
        trans=_ratio(num,den_standard,block,n_blocks,tv)
        num_alt=axis*(px-pyo); den_alt=axis*(pvin-pyo); alt=_ratio(num_alt,den_alt,block,n_blocks)
        swapped=_ratio(axis_swap*(px-pmb),den_swap,block,n_blocks)
        full,loco=_leave_one_group(num,den_standard,chrom)
        for c,estimate,delta in loco:
            chrom_rows.append({"genetic_id":sid,"chromosome":c,"alpha_leave_out":estimate,"delta_from_full":delta})
        _,lob=_leave_one_group(num,den_standard,block)
        for b,estimate,delta in lob:
            block_rows.append({"genetic_id":sid,"block":int(b),"alpha_leave_out":estimate,"delta_from_full":delta})
        loco_vals=np.array([x[1] for x in loco],float); block_delta=np.array([x[2] for x in lob],float)
        subs=_subsample(num,den_standard,target,int(scfg["subsampling_replicates"]),rng)
        boots=_block_bootstrap(num,den_standard,block,n_blocks,int(scfg["bootstrap_replicates"]),rng)
        ref_values=np.array([standard["theta"],alt["theta"],swapped["theta"]],float)
        row={
            "genetic_id":sid,"sensitivity_status":"complete","alpha_standard":standard["theta"],"alpha_standard_se":standard["se"],
            "standard_informative_snps":standard["n_used"],"effective_jackknife_blocks":standard["n_blocks_used"],
            "alpha_transition_inclusive":standard["theta"],"alpha_transversion":trans["theta"],"alpha_transversion_se":trans["se"],
            "transversion_informative_snps":trans["n_used"],"tv_delta":trans["theta"]-standard["theta"],
            "tv_collapse_fraction":max((standard["theta"]-trans["theta"])/abs(standard["theta"]),0.0) if standard["theta"] else np.nan,
            "damage_aware_estimate":trans["theta"],"damage_filter_status":"transversion_proxy_only; terminal trimming requires reads",
            "alpha_alt_outgroup":alt["theta"],"alt_outgroup":"Yoruba","alpha_swapped_reference":swapped["theta"],
            "reference_test":"swap Altai/Vindija statistic and scale roles","reference_range":float(np.nanmax(ref_values)-np.nanmin(ref_values)),
            "loco_min":float(np.nanmin(loco_vals)),"loco_max":float(np.nanmax(loco_vals)),"loco_range":float(np.nanmax(loco_vals)-np.nanmin(loco_vals)),
            "max_block_delta":float(np.nanmax(np.abs(block_delta))),"dominant_block":int(lob[int(np.nanargmax(np.abs(block_delta)))][0]),
            "snp_match_target":target,"segment_status":"not_run: no validated general segment caller for EIGENSTRAT",
        }
        row.update(_summary(subs,"subsample")); row.update(_summary(boots,"bootstrap")); rows.append(row)
        LOG.info("%s standard=%.4f tv=%.4f Yoruba=%.4f n=%d",sid,standard["theta"],trans["theta"],alt["theta"],standard["n_used"])
    rows.extend({"genetic_id":sid,"sensitivity_status":"missing_from_panel"} for sid in missing)
    pd.DataFrame(chrom_rows).to_csv(outdir/"leave_one_chromosome_out.tsv",sep="\t",index=False)
    pd.DataFrame(block_rows).to_csv(outdir/"block_influence.tsv",sep="\t",index=False)
    return pd.DataFrame(rows)
