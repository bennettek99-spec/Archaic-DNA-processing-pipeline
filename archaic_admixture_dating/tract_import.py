"""Adapters from established caller or generic tables to the V1 tract schema."""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from .tract_schema import read_tracts, validate_tracts, write_tracts

IBDMIX_ALIASES = {
    "ID": "sample_id",
    "sample": "sample_id",
    "pop": "population",
    "chr": "chromosome",
    "chrom": "chromosome",
    "start": "start_bp",
    "end": "end_bp",
    "Start": "start_bp",
    "End": "end_bp",
    "LOD": "caller_score",
    "posterior": "posterior_denisovan",
}


def import_tracts(
    source: str | Path,
    *,
    caller: str = "generic",
    column_map: dict[str, str] | None = None,
    population: str | None = None,
    output: str | Path | None = None,
    excluded_output: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = read_tracts(source)
    aliases: dict[str, str] = {}
    if caller.lower() == "ibdmix":
        aliases.update({key: value for key, value in IBDMIX_ALIASES.items() if key in frame.columns})
    aliases.update(column_map or {})
    frame = frame.rename(columns=aliases)
    if population is not None and "population" not in frame:
        frame["population"] = population
    if "caller" not in frame:
        frame["caller"] = caller
    if "caller_metadata" not in frame:
        known = {
            "sample_id", "population", "chromosome", "start_bp", "end_bp",
            "start_cm", "end_cm", "length_bp", "length_cm",
            "posterior_denisovan", "caller", "source_class",
            "source_class_probability", "callable_fraction", "qc_flags",
        }
        extras = [column for column in frame.columns if column not in known]
        if extras:
            frame["caller_metadata"] = frame[extras].apply(
                lambda row: row.dropna().to_json(), axis=1
            )
    valid, excluded = validate_tracts(frame)
    if output:
        write_tracts(valid, output)
    if excluded_output:
        write_tracts(excluded, excluded_output)
    return valid, excluded


def parse_column_map(values: list[str] | None) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for item in values or []:
        if "=" not in item:
            raise ValueError(f"Column mapping must be SOURCE=STANDARD, got {item!r}")
        source, target = item.split("=", 1)
        mapping[source.strip()] = target.strip()
    return mapping
