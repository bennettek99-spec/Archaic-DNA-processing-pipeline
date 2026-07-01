#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — dataset preparation.

Load every ancient EURASIAN genome that actually has genotypes in the panel,
attach its metadata from the .anno (site, culture, coordinates, country, age,
coverage, SNP overlap, contamination, quality assessment, ancestry label), apply
transparent QC thresholds, and emit:

  results/phase2_<panel>_metadata.csv   retained samples + all metadata + flags
  results/phase2_<panel>_excluded.csv   excluded samples, one reason each
  results/phase2_<panel>_summary.txt    counts by continent / era / exclusion

Panel-agnostic: `--panel ho` or `--panel 1240k` (see archaic/refs.py).

Inclusion logic, in order (first failing gate wins as the exclusion reason):
  1. has genotypes in the .ind                         (else: not_in_geno)
  2. not an archaic / outgroup / non-human reference    (else: reference_or_archaic)
  3. ancient, not present-day                            (else: present_day)
  4. continent == Eurasia                                (else: non_eurasian:<continent>)
  5. quality assessment not CRITICAL/FAIL                (else: qc_assessment:<value>)
  6. SNPs hit >= panel floor                             (else: low_snps:<n>)
Retained samples additionally carry non-fatal flags: low_power, questionable,
no_date, no_coords.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic import lib_eigenstrat as le, anno as anno_mod
from archaic.refs import PANELS, NONHUMAN_OR_REF

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)

# ---- continent classification (country-first; coordinates as fallback) -------
_OCEANIA = ("australia", "new guinea", "papua", "vanuatu", "solomon", "fiji",
            "tonga", "samoa", "new zealand", "guam", "micronesia", "polynesia",
            "palau", "new caledonia", "new britain", "new ireland", "bismarck",
            "tuvalu", "kiribati", "marquesas")
_AMERICAS = ("united states", "usa", "canada", "mexico", "guatemala", "belize",
             "honduras", "salvador", "nicaragua", "costa rica", "panama",
             "colombia", "venezuela", "ecuador", "peru", "brazil", "bolivia",
             "chile", "argentina", "uruguay", "paraguay", "cuba", "dominican",
             "puerto rico", "bahamas", "greenland", "haiti", "jamaica")
_AFRICA = ("morocco", "algeria", "tunisia", "libya", "egypt", "sudan",
           "ethiopia", "eritrea", "somalia", "djibouti", "kenya", "tanzania",
           "uganda", "rwanda", "burundi", "malawi", "mozambique", "zambia",
           "zimbabwe", "botswana", "namibia", "south africa", "lesotho",
           "angola", "congo", "cameroon", "nigeria", "niger", "chad", "mali",
           "mauritania", "senegal", "gambia", "guinea-bissau", "ghana",
           "ivory", "cote d", "burkina", "benin", "togo", "sierra leone",
           "liberia", "gabon", "central african", "swahili")


def continent(country, lat, lon):
    c = (country or "").lower()
    for kw in _OCEANIA:
        if kw in c:
            return "Oceania"
    for kw in _AMERICAS:
        if kw in c:
            return "Americas"
    for kw in _AFRICA:
        if kw in c:
            return "Africa"
    if np.isfinite(lat) and np.isfinite(lon):
        if lon <= -28:
            return "Americas"
        if lat <= -10 and lon >= 100:
            return "Oceania"
        return "Eurasia"
    # no country match and no coords: leave to caller (Unknown)
    if c and c not in ("nan", "..", ""):
        return "Eurasia"   # country given but unrecognised -> assume Eurasia, flag
    return "Unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--scope", choices=["eurasia", "global"], default="eurasia",
                    help="'eurasia' = original QC (ancient Eurasians only); "
                         "'global' = keep present-day + all continents, "
                         "record continent/is_modern (for whole-AADR survey)")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    prefix = cfg["prefix"]
    snps_col = cfg["snps_col"]
    is_global = args.scope == "global"
    tag = f"{args.panel}_global" if is_global else args.panel

    print(f"Phase 2 — panel={args.panel}  scope={args.scope}  prefix={prefix}")
    ind = le.read_ind(prefix + ".ind")
    have_geno = set(ind["id"].values)
    print(f"  individuals with genotypes: {len(have_geno):,}")

    ann = anno_mod.load_anno(prefix + ".anno")
    print(f"  anno rows: {len(ann):,}")
    print(f"  snps column in use: {snps_col} "
          f"(non-null: {ann[snps_col].notna().sum():,})")

    rows, excluded = [], []
    for _, r in ann.iterrows():
        gid = r["genetic_id"]
        grp = (r["group_id"] or "")
        lat, lon = r["lat"], r["lon"]
        date_bp = r["date_bp"]
        nsnp = r[snps_col]

        def drop(reason):
            excluded.append(dict(genetic_id=gid, group_id=grp, country=r["country"],
                                 date_bp=date_bp, snps=nsnp, reason=reason))

        # 1. genotypes present
        if gid not in have_geno:
            drop("not_in_geno"); continue
        # 2. references / archaics / non-human
        gl = grp.lower()
        if any(k in gl for k in NONHUMAN_OR_REF):
            drop("reference_or_archaic"); continue
        # 3. ancient vs present-day
        is_modern = (np.isfinite(date_bp) and date_bp < 50) or \
                    ((not np.isfinite(date_bp)) and gid.endswith(".HO"))
        if is_modern and not is_global:
            drop("present_day"); continue
        # 4. continent
        cont = continent(r["country"], lat, lon)
        if cont != "Eurasia" and not is_global:
            drop(f"non_eurasian:{cont}"); continue
        # 5. quality assessment
        assess = (r["assessment"] or "").strip()
        al = assess.lower()
        if ("critical" in al) or ("fail" in al):
            drop(f"qc_assessment:{assess[:20]}"); continue
        # 6. SNP floor
        if not np.isfinite(nsnp) or nsnp < cfg["snp_floor"]:
            drop(f"low_snps:{'' if not np.isfinite(nsnp) else int(nsnp)}"); continue

        # ---- retained: collect metadata + non-fatal flags
        flags = []
        if np.isfinite(nsnp) and nsnp < cfg["snp_lowpower"]:
            flags.append("low_power")
        if "questionable" in al:
            flags.append("questionable")
        if not np.isfinite(date_bp):
            flags.append("no_date")
        if not (np.isfinite(lat) and np.isfinite(lon)):
            flags.append("no_coords")
        rows.append(dict(
            genetic_id=gid, group_id=grp, locality=r["locality"],
            country=r["country"], lat=lat, lon=lon,
            continent=cont, is_modern=bool(is_modern),
            date_bp=date_bp, date_sd=r["date_sd"], full_date=r["full_date"],
            coverage=r["coverage"], snps_hit=nsnp, mol_sex=r["mol_sex"],
            assessment=assess, angsd_contam=r["angsd"], hapconx_contam=r["hapconx"],
            damage=r["damage"], y_hap=r["y_hap"], mt_hap=r["mt_hap"],
            flags=";".join(flags),
        ))

    meta = pd.DataFrame(rows)
    exc = pd.DataFrame(excluded)

    # ----------------------------------------------------------------- outputs
    mpath = os.path.join(RESULTS, f"phase2_{tag}_metadata.csv")
    epath = os.path.join(RESULTS, f"phase2_{tag}_excluded.csv")
    meta.to_csv(mpath, index=False)
    exc.to_csv(epath, index=False)

    # summary
    lines = []
    lines.append(f"Phase 2 dataset preparation — panel {args.panel} scope {args.scope}")
    if is_global and len(meta):
        lines.append("Retained — by continent:")
        for ct, n in meta["continent"].value_counts().items():
            nmod = int(meta[(meta.continent == ct)]["is_modern"].sum())
            lines.append(f"  {str(ct):20s} {n:6,d}  (present-day: {nmod:,})")
        lines.append("")
    lines.append(f"anno rows: {len(ann):,}   genotyped: {len(have_geno):,}")
    lines.append(f"RETAINED (Eurasian ancient, QC-pass): {len(meta):,}")
    lines.append(f"EXCLUDED: {len(exc):,}")
    lines.append("")
    lines.append("Exclusion reasons (grouped):")
    er = exc["reason"].str.replace(r":.*$", "", regex=True)
    for reason, n in er.value_counts().items():
        lines.append(f"  {reason:24s} {n:6,d}")
    lines.append("")
    if len(meta):
        lines.append("Retained — non-fatal flags:")
        allflags = meta["flags"].str.split(";").explode()
        allflags = allflags[allflags != ""]
        for fl, n in allflags.value_counts().items():
            lines.append(f"  {fl:24s} {n:6,d}")
        lines.append("")
        lines.append("Retained — by era (BP):")
        bins = [0, 3000, 5000, 8000, 12000, 1e9]
        labels = ["0-3k (hist/IA-MBA)", "3-5k (Neo-EBA)", "5-8k (Neolithic)",
                  "8-12k (Meso/EpiPal)", "12k+ (UP)"]
        era = pd.cut(meta["date_bp"], bins=bins, labels=labels, right=False)
        for e, n in era.value_counts().sort_index().items():
            lines.append(f"  {str(e):24s} {n:6,d}")
        lines.append(f"  {'no_date':24s} {meta['date_bp'].isna().sum():6,d}")
        lines.append("")
        lines.append("Retained — top 15 countries:")
        for ctry, n in meta["country"].value_counts().head(15).items():
            lines.append(f"  {str(ctry):24s} {n:6,d}")
        lines.append("")
        lines.append("Retained — SNP coverage quartiles:")
        q = meta["snps_hit"].quantile([0, .25, .5, .75, 1.0])
        for k, v in q.items():
            lines.append(f"  {int(k*100):3d}%  {int(v):,}")

    summary = "\n".join(lines)
    spath = os.path.join(RESULTS, f"phase2_{tag}_summary.txt")
    with open(spath, "w", encoding="utf-8") as fh:
        fh.write(summary + "\n")
    print("\n" + summary)
    print(f"\nWrote:\n  {mpath}\n  {epath}\n  {spath}")


if __name__ == "__main__":
    main()
