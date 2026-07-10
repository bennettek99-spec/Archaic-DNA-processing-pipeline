"""
cohort_rules.py - predefined archaeological cohorts and population-test filters.

Published AADR labels are useful metadata, but they are not treated as results.
These rules assign broad cohorts before inspecting genetic outcomes and provide a
simple duplicate-library pruning mask for population-level tests.  Scripts that
need stricter READ-style relative pruning can apply archaic.kinship on top.
"""
from __future__ import annotations

import re

import pandas as pd


def _gid(row) -> str:
    return str(row.get("group_id", row.get("published_group_id", ""))).lower()


def _country(row) -> str:
    return str(row.get("country", "")).lower()


def _bp(row):
    try:
        return float(row.get("date_bp"))
    except Exception:
        return float("nan")


COHORT_RULES = [
    {
        "cohort": "Etruscan_context",
        "rule": "group_id contains Etruscan; 2300-3200 BP preferred for Iron Age Etruria",
        "predicate": lambda r: "etruscan" in _gid(r) and 1700 <= _bp(r) <= 3500,
    },
    {
        "cohort": "Latin_context",
        "rule": "Latini/Latin/Lazio Iron Age labels excluding Etruscan; 2300-3200 BP",
        "predicate": lambda r: (
            ("latini" in _gid(r) or "latin" in _gid(r) or "lazio_ia" in _gid(r))
            and "etruscan" not in _gid(r)
            and 2300 <= _bp(r) <= 3200
        ),
    },
    {
        "cohort": "Preceding_Bronze_Age_Italy",
        "rule": "Italy Bronze Age labels (_BA/_EBA/_MBA/_LBA) or Italy 3200-4500 BP",
        "predicate": lambda r: (
            "italy" in _gid(r)
            and (
                any(tag in _gid(r) for tag in ("_ba", "_eba", "_mba", "_lba"))
                or 3200 <= _bp(r) <= 4500
            )
        ),
    },
    {
        "cohort": "Imperial_Roman_context",
        "rule": "Imperial Roman / Roman Italy labels; 1700-2300 BP",
        "predicate": lambda r: (
            ("imperialroman" in _gid(r) or "roman" in _gid(r))
            and "italy" in (_gid(r) + " " + _country(r))
            and 1700 <= _bp(r) <= 2300
        ),
    },
    {
        "cohort": "Early_Medieval_Italy",
        "rule": "Italy Late Antique / Medieval labels or Italy 800-1700 BP",
        "predicate": lambda r: (
            ("italy" in _gid(r) or "italy" in _country(r))
            and (
                any(tag in _gid(r) for tag in ("medieval", "lateantique", "langobard"))
                or 800 <= _bp(r) < 1700
            )
        ),
    },
]


def assign_archaeological_cohort(row) -> tuple[str, str]:
    """Return (cohort_name, inclusion_rule); empty strings mean no study cohort."""
    for spec in COHORT_RULES:
        try:
            if spec["predicate"](row):
                return spec["cohort"], spec["rule"]
        except Exception:
            continue
    return "", ""


def apply_cohort_rules(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    assigned = out.apply(assign_archaeological_cohort, axis=1, result_type="expand")
    out["archaeological_cohort"] = assigned[0]
    out["cohort_inclusion_rule"] = assigned[1]
    out["published_label_role"] = "metadata_not_truth"
    return out


def duplicate_root(genetic_id: str) -> str:
    """Collapse likely duplicate libraries by stripping the terminal data suffix."""
    s = str(genetic_id)
    s = re.sub(r"\.BY\.AA$", "", s)
    return s.rsplit(".", 1)[0] if "." in s else s


def add_population_test_keep(df: pd.DataFrame) -> pd.DataFrame:
    """Flag one highest-SNP library per duplicate root inside each cohort.

    This is the cheap deterministic pre-filter.  It removes duplicate libraries
    before population tests and leaves close-relative READ pruning to kinship.py
    in analyses that load genotype data.
    """
    out = df.copy()
    if "archaeological_cohort" not in out:
        out = apply_cohort_rules(out)
    out["duplicate_root"] = out["genetic_id"].map(duplicate_root)
    out["population_test_keep"] = True
    out["population_test_exclusion"] = ""
    snp_col = "alpha_nSNP" if "alpha_nSNP" in out else "snps_hit"
    cohorts = out["archaeological_cohort"].fillna("")
    for cohort, sub in out[cohorts != ""].groupby("archaeological_cohort"):
        ranked = sub.sort_values(snp_col, ascending=False)
        keep_idx = ranked.groupby("duplicate_root", sort=False).head(1).index
        drop_idx = sub.index.difference(keep_idx)
        out.loc[drop_idx, "population_test_keep"] = False
        out.loc[drop_idx, "population_test_exclusion"] = (
            "duplicate_library_lower_snp_count"
        )
    return out


def rules_table() -> pd.DataFrame:
    return pd.DataFrame(
        {"cohort": r["cohort"], "inclusion_rule": r["rule"]} for r in COHORT_RULES
    )
