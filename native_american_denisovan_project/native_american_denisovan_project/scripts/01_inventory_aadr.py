#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_inventory_aadr.py
====================
Scan the local AADR v66.p1 ``.anno`` metadata and categorize every individual
into the population groups required by the Native-American Denisovan study:
archaic references, ancient Native Americans, Siberian/ANE, ancient & present-day
East Asians (incl. Jomon), Papuan/Australasian comparisons, African outgroups,
and western-Eurasian controls.

This script is READ-ONLY with respect to the AADR: it touches only the ``.anno``
text file (no genotype I/O). It reuses the validated ``archaic.anno`` parser by
adding the main pipeline directory to ``sys.path``; it does NOT import any
genotype data or modify the production package.

Outputs (under ``native_american_denisovan_project/``):
  data/manifests/aadr_inventory_full.tsv     every AADR individual + category
  data/manifests/aadr_inventory_relevant.tsv the subset relevant to this study
  results/tables/table1_sample_inventory.tsv per-individual manifest (Table 1)
  results/tables/table2_population_groups.tsv group counts + SNP/coverage stats (Table 2)
  results/logs/inventory_summary.txt         human-readable summary

Usage:
  python scripts/01_inventory_aadr.py
  python scripts/01_inventory_aadr.py --anno "C:/Users/benne/aadr_v66/v66.p1_1240K.anno"
"""
from __future__ import annotations

import argparse
import os
import sys
import statistics as stats_mod

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _c in (os.path.join(os.path.dirname(_REPO), "archaic-introgression"), os.path.dirname(_REPO)):
    if os.path.isdir(os.path.join(_c, "archaic")):
        sys.path.insert(0, _c)
        break

from archaic.anno import load_anno


ARCHAIC_IDS = {
    "Denisova.SG", "Denisova3.DG", "Denisova3_snpAD.DG", "Denisova11.SG",
    "Denisova25.SG", "AltaiNeanderthal.DG", "VindijaG1_final.SG",
    "Chagyrskaya8.DG", "Chimp.REF", "Chimp_HO.HO",
}

ANCIENT_BP = 500.0

AMERICAN_COUNTRIES = {
    "united states", "usa", "u.s.a.", "alaska", "canada", "mexico", "guatemala",
    "belize", "honduras", "el salvador", "nicaragua", "costa rica", "panama",
    "colombia", "venezuela", "ecuador", "peru", "bolivia", "brazil",
    "paraguay", "chile", "argentina", "uruguay", "greenland", "cuba",
    "dominican republic", "haiti", "bahamas", "puerto rico",
}

EAST_ASIAN_COUNTRIES = {
    "china", "japan", "south korea", "korea", "north korea", "mongolia",
    "taiwan", "vietnam", "thailand", "cambodia", "laos", "myanmar",
    "philippines", "malaysia", "indonesia",
}

SIBERIA_KEYWORDS = (
    "siberia", "altai", "baikal", "yana", "kolyma", "mal'ta", "malta",
    "afontova", "ust'", "ust-ishim", "irkutsk", "yakutia", "kamchatka",
    "far east", "amur",
)

JOMON_KEYWORDS = ("jomon",)
NEGRITO_KEYWORDS = ("aeta", "ati", "mamanwa", "agta", "batak", "onota", "iranun")

PAPUAN_AUSTRAL_COUNTRIES = {
    "papua new guinea", "australia", "solomon islands", "vanuatu", "fiji",
    "new caledonia", "new zealand", "samoa", "tonga",
}
ANDAMAN_KEYWORDS = ("andaman", "ongea", "onge", "jarawa", "sentinel", "great andaman")

AFRICAN_POPS = {
    "mbuti", "yoruba", "mende", "dinka", "jola", "esán", "brong",
    "mandenka", "bulala", "danane", "kanembu", "ju_hoan_north", "ju_hoan_south",
    "ju/'hoan", "san", "khomani", "gumuz", "chimp", "biaka", "mozabite",
}

WESTEUR_POPS = {"french", "sardinian", "basque", "spanish", "italian", "english",
                "russian", "norsk", "norwegian", "czech", "greek", "tuscan",
                "bergamo", "albanian", "arbereshe", "cimbrianoi", "cypriot",
                "german", "icelandic", "finnish", "estonian", "lithuanian",
                "latvian", "scottish", "irish", "welsh", "hungarian",
                "bulgarian", "romanian", "serbian", "croatian", "polish",
                "ukrainian", "belarusian", "turkish", "armenian", "georgian",
                "lebanese", "syrian", "jordanian", "palestinian", "samaritan",
                "bedouin", "mozabite", "saharawi", "tunisian", "algerian",
                "moroccan", "egyptian", "ethiopian"}

EASTASIAN_POPS = {"han", "dai", "she", "tujia", "japanese", "korean", "koreans",
                  "naxi", "yaoxishuangbanna", "yizu", "buyi", "miao",
                  "mongola", "xibo", "yugur", "hezhen", "daur", "orqen",
                  "ukrainianoutlier", "tu", "laihui", "chinae", "chinesesouth",
                  "chinesenorth", "amihof", "atayal", "ami", "igorot",
                  "philippines_kankanaey", "kinh_vietnamese", "thai",
                  "cambodian", "malay", "bidayuh", "jehai", "kensiu",
                  "temuan", "classicnegrito", "aeta", "agta", "atin", "ati",
                  "mamanwa", "irigweeta", "chinese", "southchinesehan",
                  "northchinesehan", "hccnhan", "hccnhan_bj", "tibetan",
                  "tibetansherpa", "sherpa", "tulacomm", "nasueh", "dzongkha",
                  "lai", "hlabisa", "bunun", "hakka"}

AMERICAN_POPS = {"karitiana", "surui", "mayan", "maya", "zapotec", "mixe",
                 "mixtec", "pima", "nahua", "kaingang", "guarani", "wayku",
                 "quetzaltepec", "quechua", "aymara", "chipewa", "chippewa",
                 "ojibwa", "cree", "inuit", "aleut", "tlingit", "haida",
                 "kaqchikel", "k'iche", "tzotzil", "totonac", "tepehuan",
                 "huichol", "cora", "tarahumara", "purepecha", "mixtec",
                 "triqui", "amuzgo", "chatino", "mazatec", "chinantec",
                 "chontal", "popoloca", "mazahua", "matlatzinca", "ocuiltec",
                 "chichimec", "cucapa", "kiliwa", "paipai", "cochimi",
                 "seri", "tequistlatec", "georgianmayor", "northamerican",
                 "southamerican", "peruvian", "bolivian", "colombian",
                 "muisca", "chibcha", "embera", "waunana", "kogi", "arsario",
                 "yukpa", "bari", "warao", "macushi", "waiwai", "wapishana",
                 "yamana", "kawesqar", "selknam", "tehuelche", "mapuche",
                 "huelche", "pehuenche", "puelche", "teushente", "charrua",
                 "guarani", "chane", "charrua", "diaguita", "omaguaca",
                 "atacama", "colla", "quechua", "amazonian"}

PAPUAN_POPS = {"papuan", "nasioi", "bougainville", "baining", "aulua",
               "maewo", "barki", "sori", "manus", "karkar", "madak",
               "tongan", "samoa", "fijian", "rotuman", "niue", "wallis",
               "maori", "australian", "aboriginal", "wongathi", "ngaanyatjarra",
               "pijantjatjara", "western_desert", "arnhem", "tiwi", "gagadju"}


def _norm(s):
    if s is None:
        return ""
    return str(s).strip().lower()


def categorize(row):
    gid = _norm(row.get("genetic_id"))
    country = _norm(row.get("country"))
    locality = _norm(row.get("locality"))
    group = _norm(row.get("group_id"))
    date_bp = row.get("date_bp")
    is_ancient = np.isfinite(date_bp) and date_bp >= ANCIENT_BP

    if gid in {x.lower() for x in ARCHAIC_IDS}:
        return "Archaic_reference"

    if any(k in locality or k in group or k in gid for k in SIBERIA_KEYWORDS) and is_ancient:
        if "yan" in locality or "yan" in gid or "yana" in locality:
            return "Ancient_Paleo_Siberian"
        if any(k in (locality + " " + gid) for k in ("mal'ta", "malta", "afontova", "ag1", "ag2", "ag3")):
            return "Ancient_North_Eurasian"
        return "Ancient_Siberian"

    if country in AMERICAN_COUNTRIES or group in AMERICAN_POPS:
        if is_ancient:
            if country in ("alaska", "greenland") or "aleut" in (locality + group):
                return "Ancient_Arctic_Beringian"
            if country in ("united states", "usa", "canada"):
                return "Ancient_North_America"
            if country in ("mexico", "guatemala", "belize", "honduras",
                           "el salvador", "nicaragua", "costa rica", "panama"):
                return "Ancient_Mesoamerica"
            if country in ("cuba", "dominican republic", "haiti", "bahamas",
                           "puerto rico"):
                return "Ancient_Caribbean"
            return "Ancient_South_America"
        return "Modern_American_admixed"

    if is_ancient and (country in EAST_ASIAN_COUNTRIES or
                       any(k in (locality + " " + group) for k in JOMON_KEYWORDS)):
        if any(k in (locality + " " + group + " " + gid) for k in JOMON_KEYWORDS):
            return "Ancient_Jomon_Japan"
        return "Ancient_East_Asian"

    if any(k in (group + " " + locality) for k in ANDAMAN_KEYWORDS) or "andaman" in country:
        return "Andamanese"
    if country in PAPUAN_AUSTRAL_COUNTRIES or group in PAPUAN_POPS:
        if group in PAPUAN_POPS or "papua" in country or "australia" in country:
            return "Papuan_Australasian"
        return "Oceanian"
    if any(k in group for k in NEGRITO_KEYWORDS):
        return "Philippine_Negrito"

    if group in EASTASIAN_POPS or country in EAST_ASIAN_COUNTRIES:
        return "Present_day_East_Asian"

    if group in AFRICAN_POPS or "africa" in country or country in {
            "nigeria", "kenya", "tanzania", "ethiopia", "senegal", "gambia",
            "sierra leone", "cameroon", "democratic republic of congo",
            "congo", "ghana", "ivory coast", "south africa", "botswana",
            "namibia", "angola", "zambia", "zimbabwe", "mozambique",
            "madagascar", "sudan", "south sudan", "eritrea", "djibouti",
            "somalia", "uganda", "rwanda", "burundi", "central african republic",
            "chad", "niger", "mali", "burkina faso", "guinea", "guinea-bissau",
            "liberia", "togo", "benin", "gabon", "equatorial guinea",
            "sao tome and principe", "mauritania", "western sahara"}:
        return "African_outgroup"

    if group in WESTEUR_POPS or country in {
            "france", "italy", "spain", "portugal", "germany", "united kingdom",
            "uk", "ireland", "netherlands", "belgium", "switzerland", "austria",
            "sweden", "norway", "denmark", "finland", "estonia", "latvia",
            "lithuania", "poland", "czech republic", "slovakia", "hungary",
            "slovenia", "croatia", "bosnia and herzegovina", "serbia",
            "romania", "bulgaria", "greece", "albania", "macedonia",
            "montenegro", "kosovo", "russia", "ukraine", "belarus",
            "moldova", "iceland", "luxembourg", "liechtenstein", "andorra",
            "malta", "cyprus", "turkey", "georgia", "armenia", "azerbaijan",
            "israel", "jordan", "lebanon", "syria", "iraq", "iran",
            "saudi arabia", "yemen", "oman", "uae", "qatar", "kuwait",
            "bahrain", "egypt", "libya", "tunisia", "algeria", "morocco",
            "western sahara", "sudan"}:
        return "West_Eurasian_control"

    if is_ancient:
        return "Ancient_other"
    return "Other_present_day"


RELEVANT_CATEGORIES = {
    "Archaic_reference", "Ancient_Paleo_Siberian", "Ancient_North_Eurasian",
    "Ancient_Siberian", "Ancient_Arctic_Beringian", "Ancient_North_America",
    "Ancient_Mesoamerica", "Ancient_Caribbean", "Ancient_South_America",
    "Modern_American_admixed", "Ancient_Jomon_Japan", "Ancient_East_Asian",
    "Present_day_East_Asian", "Andamanese", "Papuan_Australasian", "Oceanian",
    "Philippine_Negrito", "African_outgroup", "West_Eurasian_control",
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anno", default=r"C:/Users/benne/aadr_v66/v66.p1_1240K.anno")
    ap.add_argument("--min-snps", type=int, default=30000,
                    help="1240K SNP floor for 'usable' ancient samples (matches pipeline)")
    args = ap.parse_args(argv)

    out = os.path.join(_REPO, "results")
    for sub in ("tables", "logs", "reports"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    manif = os.path.join(_REPO, "data", "manifests")
    os.makedirs(manif, exist_ok=True)

    print(f"Loading AADR annotation: {args.anno}")
    anno = load_anno(args.anno)
    print(f"  {len(anno):,} individuals in AADR 1240K release")

    cols = ["genetic_id", "group_id", "locality", "country", "lat", "lon",
            "date_bp", "date_sd", "coverage", "snps_1240k", "mol_sex",
            "assessment", "assess_warn", "damage", "publication",
            "publication_doi", "library_type"]
    cols = [c for c in cols if c in anno.columns]

    df = anno[cols].copy()
    df["category"] = [categorize(r) for _, r in anno.iterrows()]
    df["is_ancient"] = np.isfinite(df["date_bp"].to_numpy()) & (df["date_bp"].to_numpy() >= ANCIENT_BP)
    df["usable_ancient"] = (df["is_ancient"] & np.isfinite(df["snps_1240k"].to_numpy())
                            & (df["snps_1240k"].to_numpy() >= args.min_snps))

    df.to_csv(os.path.join(manif, "aadr_inventory_full.tsv"), sep="\t", index=False)
    rel = df[df["category"].isin(RELEVANT_CATEGORIES)].copy()
    rel.to_csv(os.path.join(manif, "aadr_inventory_relevant.tsv"), sep="\t", index=False)

    rel = rel.rename(columns={"genetic_id": "sample_id", "group_id": "population_label"})
    keep_cols = [c for c in ["sample_id", "population_label", "category", "locality",
                            "country", "lat", "lon", "date_bp", "date_sd",
                            "coverage", "snps_1240k", "mol_sex", "assessment",
                            "publication", "publication_doi", "usable_ancient"]
                if c in rel.columns]
    rel[keep_cols].to_csv(os.path.join(out, "tables", "table1_sample_inventory.tsv"),
                          sep="\t", index=False)

    g_rows = []
    for cat, grp in df.groupby("category"):
        snp = grp["snps_1240k"].to_numpy()
        cov = grp["coverage"].to_numpy()
        n_anc = int(grp["usable_ancient"].sum()) if "usable_ancient" in grp else 0
        g_rows.append({
            "category": cat,
            "n_individuals": len(grp),
            "n_ancient_usable": n_anc,
            "median_snps_1240k": (float(np.nanmedian(snp)) if np.isfinite(snp).any() else np.nan),
            "median_coverage": (float(np.nanmedian(cov)) if np.isfinite(cov).any() else np.nan),
        })
    gdf = (pd.DataFrame(g_rows).sort_values("n_individuals", ascending=False))
    gdf.to_csv(os.path.join(out, "tables", "table2_population_groups.tsv"),
               sep="\t", index=False)

    log_lines = []
    def w(s=""):
        print(s); log_lines.append(s)

    w("=" * 78)
    w("AADR v66.p1 1240K INVENTORY — Native-American Denisovan study")
    w("=" * 78)
    w(f"Total individuals in release: {len(anno):,}")
    w(f"SNP floor for usable ancient: {args.min_snps:,} 1240K SNPs")
    w("")
    w("Relevant categories (Table 2):")
    w(f"{'category':<32} {'n_ind':>7} {'n_anc_use':>10} {'med_SNP':>10} {'med_cov':>8}")
    w("-" * 78)
    for _, r in gdf[gdf["category"].isin(RELEVANT_CATEGORIES)].iterrows():
        w(f"{r['category']:<32} {r['n_individuals']:>7} {r['n_ancient_usable']:>10} "
          f"{r['median_snps_1240k']:>10.0f} {r['median_coverage']:>8.2f}")
    w("")
    w("All categories (for completeness):")
    for _, r in gdf.iterrows():
        w(f"  {r['category']:<32} n={r['n_individuals']:>6}  anc_usable={r['n_ancient_usable']:>6}  "
          f"med_SNP={r['median_snps_1240k']:.0f}  med_cov={r['median_coverage']:.2f}")

    sub_americas = ["Ancient_Arctic_Beringian", "Ancient_North_America",
                    "Ancient_Mesoamerica", "Ancient_Caribbean",
                    "Ancient_South_America"]
    am = df[df["category"].isin(sub_americas) & df["usable_ancient"]]
    w("")
    w("Ancient AMERICAS usable individuals (SNP >= floor):")
    w(f"  total usable ancient American = {len(am):,}")
    for cat in sub_americas:
        sub = am[am["category"] == cat]
        if len(sub):
            w(f"  {cat:<28} n={len(sub):>4}  date range "
              f"{sub['date_bp'].min():.0f}-{sub['date_bp'].max():.0f} BP  "
              f"med_SNP={sub['snps_1240k'].median():.0f}  med_cov={sub['coverage'].median():.2f}x")

    for name, countries in [("Papuan/Australasian present-day", None),
                            ("African outgroup present-day", None)]:
        pass

    w("")
    w("Key named ancient individuals (presence check):")
    named = ["MA-1", "MA1", "Mal'ta1", "AfontovaGora2", "AfontovaGora3",
            "Yana1", "YanaRHS", "Kolyma1", "USR1", "Anzick-1", "Anzick1",
            "SpiritCave", "LagoaSanta", "Sumidouro", "UpwardSunRiver11",
            "USRS11", "Beringian", "PimaSGDP"]
    ann_id = anno.set_index("genetic_id")
    for nm in named:
        hits = [i for i in ann_id.index if nm.lower() in i.lower()]
        for h in hits[:3]:
            r = ann_id.loc[h]
            w(f"  {h:<30} group={str(r.get('group_id','')):<22} "
              f"country={str(r.get('country','')):<22} "
              f"date={r.get('date_bp')} snps={r.get('snps_1240k')}")

    with open(os.path.join(out, "logs", "inventory_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    print(f"\nWrote: data/manifests/aadr_inventory_*.tsv, results/tables/table[1,2]*.tsv, "
          f"results/logs/inventory_summary.txt")


if __name__ == "__main__":
    raise SystemExit(main())
