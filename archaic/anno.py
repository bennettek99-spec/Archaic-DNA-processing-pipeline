"""
anno.py — parse the AADR .anno metadata file into a tidy DataFrame.

The .anno is a tab-separated file with a single, very long header row whose column
names contain commas, quotes and parentheses. We therefore split on TAB only
(csv.reader, quotechar='"') and locate the columns we need by matching header
text with predicates (robust to the exact wording differing slightly between the
HO and 1240K releases). Numeric fields ('..' = not computed) are coerced to NaN.
"""
from __future__ import annotations
import csv
import numpy as np
import pandas as pd


def _find(headers, pred):
    """Index of the first header satisfying pred(lowercased_header), else None."""
    for i, h in enumerate(headers):
        if pred(h.strip().lower()):
            return i
    return None


# canonical field -> predicate on the lowercased header
_FIELDS = {
    "genetic_id":   lambda h: h.startswith("genetic id"),
    "master_id":    lambda h: h.startswith("master id") or h.startswith("version id"),
    "group_id":     lambda h: h.startswith("group id") or h == "group id",
    "locality":     lambda h: h.startswith("locality"),
    "country":      lambda h: h.startswith("political entity"),
    "lat":          lambda h: h.startswith("lat"),
    "lon":          lambda h: h.startswith("long"),
    "date_bp":      lambda h: h.startswith("date mean in bp"),
    "date_sd":      lambda h: h.startswith("date standard deviation"),
    "full_date":    lambda h: h.startswith("full date"),
    "coverage":     lambda h: h.startswith("mean coverage on"),
    "snps_2m":      lambda h: ("snps hit" in h) and ("2m capture" in h or "enhance 2m" in h),
    "snps_1240k":   lambda h: ("snps hit" in h) and ("1240k snpset" in h),
    "snps_ho":      lambda h: ("snps hit" in h) and ("ho snpset" in h) and "compatibility" not in h,
    "mol_sex":      lambda h: h.startswith("molecular sex"),
    "assessment":   lambda h: h == "assessment",
    "assess_warn":  lambda h: h.startswith("assessment warnings"),
    "angsd":        lambda h: h.startswith("angsd mom"),
    "hapconx":      lambda h: h.startswith("hapconx"),
    "damage":       lambda h: h.startswith("damage rate in first"),
    "y_hap":        lambda h: h.startswith("y haplogroup in terminal"),
    "mt_hap":       lambda h: h.startswith("mtdna haplogroup"),
}

_NUM = ["lat", "lon", "date_bp", "date_sd", "coverage",
        "snps_2m", "snps_1240k", "snps_ho", "damage"]


def _num(x):
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in ("", "..", "...", "n/a", "na", "NA"):
        return np.nan
    try:
        return float(s.replace(",", ""))
    except ValueError:
        return np.nan


def load_anno(path: str) -> pd.DataFrame:
    """Return a DataFrame with one row per sample and the canonical columns above."""
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quotechar='"')
        headers = next(reader)
        idx = {f: _find(headers, pred) for f, pred in _FIELDS.items()}
        rows = []
        for r in reader:
            if not r:
                continue
            rec = {}
            for f, i in idx.items():
                rec[f] = r[i] if (i is not None and i < len(r)) else None
            rows.append(rec)

    df = pd.DataFrame(rows)
    for c in _NUM:
        if c in df:
            df[c] = df[c].map(_num)
    # tidy strings
    for c in ("genetic_id", "group_id", "locality", "country", "assessment"):
        if c in df:
            df[c] = df[c].astype(str).str.strip()
    df["_col_index"] = idx  # provenance (which header each field came from)
    df.attrs["col_index"] = idx
    df.attrs["headers"] = headers
    return df
