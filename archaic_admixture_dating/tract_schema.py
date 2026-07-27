"""Standard tract schema, normalization, and fail-closed validation."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

STANDARD_COLUMNS = [
    "sample_id",
    "population",
    "chromosome",
    "start_bp",
    "end_bp",
    "start_cm",
    "end_cm",
    "length_bp",
    "length_cm",
    "posterior_archaic",
    "posterior_denisovan",
    "caller",
    "source_class",
    "source_class_probability",
    "callable_fraction",
    "qc_flags",
]
VALID_CHROMS = {str(value) for value in range(1, 23)} | {"X", "Y", "MT"}


class TractValidationError(ValueError):
    pass


def normalize_chromosome(value: object) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.upper()
    if text == "M":
        text = "MT"
    return text


def _reason_strings(reason_sets: list[list[str]]) -> list[str]:
    return [";".join(dict.fromkeys(values)) for values in reason_sets]


def validate_tracts(frame: pd.DataFrame, *, reject_invalid: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("tracts must be a pandas DataFrame")
    work = frame.copy()
    missing_core = [column for column in ("sample_id", "population", "chromosome", "start_bp", "end_bp") if column not in work]
    if missing_core:
        raise TractValidationError(f"Missing required tract columns: {', '.join(missing_core)}")
    defaults: dict[str, object] = {
        "start_cm": np.nan,
        "end_cm": np.nan,
        "length_bp": np.nan,
        "length_cm": np.nan,
        "posterior_archaic": np.nan,
        "posterior_denisovan": np.nan,
        "caller": "unknown",
        "source_class": "unresolved",
        "source_class_probability": np.nan,
        "callable_fraction": np.nan,
        "qc_flags": "",
    }
    for column, value in defaults.items():
        if column not in work:
            work[column] = value
    work["chromosome"] = work["chromosome"].map(normalize_chromosome)
    numeric = [
        "start_bp",
        "end_bp",
        "start_cm",
        "end_cm",
        "length_bp",
        "length_cm",
        "posterior_archaic",
        "posterior_denisovan",
        "source_class_probability",
        "callable_fraction",
    ]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    missing_bp_length = work["length_bp"].isna()
    work.loc[missing_bp_length, "length_bp"] = (
        work.loc[missing_bp_length, "end_bp"] - work.loc[missing_bp_length, "start_bp"]
    )
    missing_cm_length = work["length_cm"].isna() & work["start_cm"].notna() & work["end_cm"].notna()
    work.loc[missing_cm_length, "length_cm"] = (
        work.loc[missing_cm_length, "end_cm"] - work.loc[missing_cm_length, "start_cm"]
    )

    reasons: list[list[str]] = [[] for _ in range(len(work))]

    def flag(mask: Iterable[bool], message: str) -> None:
        for idx in np.flatnonzero(np.asarray(mask, dtype=bool)):
            reasons[int(idx)].append(message)

    flag(~work["chromosome"].isin(VALID_CHROMS), "invalid_chromosome")
    flag(work["sample_id"].astype(str).str.strip().eq(""), "missing_sample_id")
    flag(work["population"].astype(str).str.strip().eq(""), "missing_population")
    flag(work["start_bp"].isna() | work["end_bp"].isna(), "missing_bp_coordinate")
    flag(work["start_bp"].lt(0) | work["end_bp"].lt(0), "negative_bp_coordinate")
    flag(work["end_bp"].le(work["start_bp"]), "end_not_greater_than_start")
    flag(work["length_bp"].lt(0), "negative_length_bp")
    flag(work["length_cm"].isna(), "missing_length_cm")
    flag(work["length_cm"].le(0), "nonpositive_length_cm")
    for column in (
        "posterior_archaic",
        "posterior_denisovan",
        "source_class_probability",
        "callable_fraction",
    ):
        finite = work[column].notna()
        flag(finite & ~work[column].between(0, 1), f"{column}_outside_0_1")

    reason_text = _reason_strings(reasons)
    invalid = np.asarray([bool(value) for value in reason_text])
    work["_exclusion_reason"] = reason_text
    excluded = work.loc[invalid].copy()
    valid = work.loc[~invalid].drop(columns=["_exclusion_reason"]).copy()
    valid = valid[STANDARD_COLUMNS + [column for column in valid.columns if column not in STANDARD_COLUMNS]]
    if reject_invalid and len(excluded):
        counts = excluded["_exclusion_reason"].value_counts().to_dict()
        raise TractValidationError(f"Rejected {len(excluded)} malformed tracts: {counts}")
    return valid.reset_index(drop=True), excluded.reset_index(drop=True)


def read_tracts(path: str | Path, **kwargs) -> pd.DataFrame:
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix in {".tsv", ".txt", ".bed"}:
        return pd.read_csv(source, sep="\t", **kwargs)
    if suffix == ".csv":
        return pd.read_csv(source, **kwargs)
    if suffix in {".parquet", ".pq"}:
        try:
            return pd.read_parquet(source, **kwargs)
        except ImportError as error:
            raise RuntimeError("Parquet support requires the optional pyarrow dependency") from error
    raise ValueError(f"Unsupported tract format: {source.suffix}")


def write_tracts(frame: pd.DataFrame, path: str | Path) -> None:
    from .checkpointing import atomic_write_text

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.suffix.lower() in {".parquet", ".pq"}:
        temporary = target.with_name(f".{target.name}.tmp.parquet")
        try:
            frame.to_parquet(temporary, index=False)
            temporary.replace(target)
        except ImportError as error:
            raise RuntimeError("Parquet support requires the optional pyarrow dependency") from error
        finally:
            if temporary.exists():
                temporary.unlink()
    else:
        separator = "," if target.suffix.lower() == ".csv" else "\t"
        atomic_write_text(target, frame.to_csv(index=False, sep=separator))
