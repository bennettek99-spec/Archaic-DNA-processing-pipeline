"""
manifest.py - frozen AADR sample-manifest helpers.

The manifest is the audit trail for every sample considered by Phase 2.  It
records AADR release/panel, archaeological and technical metadata, the exact
inclusion/exclusion decision, and a content hash next to the CSV so downstream
analyses can say which frozen sample universe they used.
"""
from __future__ import annotations

import hashlib
import os
import re
from typing import Mapping

import numpy as np
import pandas as pd


MANIFEST_COLUMNS = [
    "aadr_release",
    "panel",
    "panel_prefix",
    "genetic_id",
    "master_id",
    "published_group_id",
    "archaeological_context",
    "locality",
    "country",
    "lat",
    "lon",
    "date_bp",
    "date_sd",
    "full_date",
    "snp_count",
    "coverage",
    "damage",
    "contamination_angsd",
    "contamination_hapconx",
    "sex",
    "assessment",
    "status",
    "exclusion_reason",
    "notes",
]


def infer_aadr_release(prefix: str) -> str:
    """Infer a stable release label from an AADR prefix such as v66.p1_1240K."""
    stem = os.path.basename(str(prefix))
    m = re.search(r"(v\d+(?:\.p\d+)?)", stem, flags=re.I)
    return m.group(1) if m else "unknown"


def _clean(value):
    if value is None:
        return ""
    try:
        if isinstance(value, float) and not np.isfinite(value):
            return ""
    except TypeError:
        pass
    return value


def manifest_record(
    row: Mapping,
    *,
    panel: str,
    panel_prefix: str,
    snps_col: str,
    status: str,
    exclusion_reason: str = "",
    notes: str = "",
) -> dict:
    """Create one canonical manifest row from an AADR .anno row."""
    group_id = _clean(row.get("group_id"))
    locality = _clean(row.get("locality"))
    country = _clean(row.get("country"))
    return {
        "aadr_release": infer_aadr_release(panel_prefix),
        "panel": panel,
        "panel_prefix": os.path.basename(str(panel_prefix)),
        "genetic_id": _clean(row.get("genetic_id")),
        "master_id": _clean(row.get("master_id")),
        "published_group_id": group_id,
        "archaeological_context": group_id,
        "locality": locality,
        "country": country,
        "lat": _clean(row.get("lat")),
        "lon": _clean(row.get("lon")),
        "date_bp": _clean(row.get("date_bp")),
        "date_sd": _clean(row.get("date_sd")),
        "full_date": _clean(row.get("full_date")),
        "snp_count": _clean(row.get(snps_col)),
        "coverage": _clean(row.get("coverage")),
        "damage": _clean(row.get("damage")),
        "contamination_angsd": _clean(row.get("angsd")),
        "contamination_hapconx": _clean(row.get("hapconx")),
        "sex": _clean(row.get("mol_sex")),
        "assessment": _clean(row.get("assessment")),
        "status": status,
        "exclusion_reason": exclusion_reason,
        "notes": notes,
    }


def write_frozen_manifest(records, path: str) -> str:
    """Write a deterministic manifest CSV and a sibling .sha256 file.

    Returns the hex digest of the CSV bytes.
    """
    df = pd.DataFrame(records)
    for col in MANIFEST_COLUMNS:
        if col not in df:
            df[col] = ""
    df = df[MANIFEST_COLUMNS].sort_values(
        ["status", "country", "published_group_id", "genetic_id"],
        na_position="last",
    )
    df.to_csv(path, index=False, lineterminator="\n")
    with open(path, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    with open(path + ".sha256", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{digest}  {os.path.basename(path)}\n")
    return digest
