"""
snp_filters.py - small SNP-set filters for robustness analyses.
"""
from __future__ import annotations

import pandas as pd


TRANSVERSIONS = {frozenset(x) for x in (("A", "C"), ("A", "T"), ("C", "G"), ("G", "T"))}


def transversion_mask(snp: pd.DataFrame):
    """Return True for A/C, A/T, C/G, and G/T SNPs."""
    a1 = snp["a1"].astype(str).str.upper()
    a2 = snp["a2"].astype(str).str.upper()
    return [frozenset((x, y)) in TRANSVERSIONS for x, y in zip(a1, a2)]
