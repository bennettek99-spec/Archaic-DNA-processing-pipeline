"""Transparent tract filters and interval-mask diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _append_flag(existing: pd.Series, flag: str) -> pd.Series:
    values = existing.fillna("").astype(str)
    return values.where(values.eq(""), values + ";") + flag


def load_bed(path: str | Path) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 1, 2])
    bed.columns = ["chromosome", "start_bp", "end_bp"]
    bed["chromosome"] = bed["chromosome"].astype(str).str.replace(r"^chr", "", regex=True).str.upper()
    return bed


def overlap_mask(tracts: pd.DataFrame, intervals: pd.DataFrame) -> np.ndarray:
    result = np.zeros(len(tracts), dtype=bool)
    for chromosome, index in tracts.groupby("chromosome").groups.items():
        regions = intervals.loc[intervals["chromosome"] == chromosome]
        if regions.empty:
            continue
        starts = tracts.loc[index, "start_bp"].to_numpy()
        ends = tracts.loc[index, "end_bp"].to_numpy()
        hit = np.zeros(len(index), dtype=bool)
        for region in regions.itertuples(index=False):
            hit |= (starts < region.end_bp) & (ends > region.start_bp)
        result[np.asarray(index)] = hit
    return result


def filter_tracts(
    tracts: pd.DataFrame,
    *,
    minimum_length_cm: float,
    minimum_confidence: float,
    minimum_callable_fraction: float,
    masks: dict[str, pd.DataFrame] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = tracts.copy()
    reasons: list[list[str]] = [[] for _ in range(len(work))]

    def flag(mask: np.ndarray | pd.Series, reason: str) -> None:
        for idx in np.flatnonzero(np.asarray(mask, dtype=bool)):
            reasons[int(idx)].append(reason)

    flag(work["length_cm"].lt(minimum_length_cm), "below_minimum_length")
    confidence_known = work["posterior_denisovan"].notna()
    flag(confidence_known & work["posterior_denisovan"].lt(minimum_confidence), "below_confidence")
    callable_known = work["callable_fraction"].notna()
    flag(callable_known & work["callable_fraction"].lt(minimum_callable_fraction), "low_callable_fraction")
    for label, intervals in (masks or {}).items():
        flag(overlap_mask(work, intervals), f"overlaps_{label}")

    text = [";".join(values) for values in reasons]
    excluded_mask = np.asarray([bool(value) for value in text])
    work["_exclusion_reason"] = text
    excluded = work.loc[excluded_mask].copy()
    retained = work.loc[~excluded_mask].drop(columns=["_exclusion_reason"]).copy()
    return retained.reset_index(drop=True), excluded.reset_index(drop=True)


def overlapping_pairs(tracts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    ordered = tracts.sort_values(["sample_id", "chromosome", "start_bp", "end_bp"])
    for (sample, chromosome), group in ordered.groupby(["sample_id", "chromosome"], sort=False):
        previous_index = None
        previous_end = -1
        previous_caller = None
        for index, row in group.iterrows():
            if row["start_bp"] < previous_end:
                rows.append(
                    {
                        "sample_id": sample,
                        "chromosome": chromosome,
                        "left_index": previous_index,
                        "right_index": index,
                        "same_caller": previous_caller == row["caller"],
                        "exact_duplicate": (
                            row["start_bp"] == group.loc[previous_index, "start_bp"]
                            and row["end_bp"] == group.loc[previous_index, "end_bp"]
                            and previous_caller == row["caller"]
                        ),
                    }
                )
            if row["end_bp"] > previous_end:
                previous_index = index
                previous_end = row["end_bp"]
                previous_caller = row["caller"]
    return pd.DataFrame(rows)
