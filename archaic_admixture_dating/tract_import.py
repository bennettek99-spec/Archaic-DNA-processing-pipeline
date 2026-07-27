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

SKOV_2018_ALIASES = {
    "name": "sample_id",
    "pop": "population",
    "chrom": "chromosome",
    "start": "start_bp",
    "end": "end_bp",
    "length": "length_bp",
    "MeanProb": "posterior_archaic",
}
SKOV_2018_RECOMBINATION_RATE = 1.2e-8


def _adapt_skov_hmm(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "name",
        "pop",
        "chrom",
        "start",
        "end",
        "length",
        "MeanProb",
        "Shared_with_Altai",
        "Shared_with_Denisova",
        "Shared_with_Vindija",
        "method",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Skov HMM input is missing columns: {', '.join(missing)}")
    papuan = frame.loc[
        frame["pop"].astype(str).str.casefold().eq("papuans")
        & frame["method"].astype(str).str.casefold().eq("hmm")
    ].copy()
    papuan = papuan.rename(columns=SKOV_2018_ALIASES)
    for column in (
        "start_bp",
        "end_bp",
        "length_bp",
        "posterior_archaic",
        "Shared_with_Altai",
        "Shared_with_Denisova",
        "Shared_with_Vindija",
    ):
        papuan[column] = pd.to_numeric(papuan[column], errors="coerce")
    scale_to_cm = SKOV_2018_RECOMBINATION_RATE * 100.0
    papuan["start_cm"] = papuan["start_bp"] * scale_to_cm
    papuan["end_cm"] = papuan["end_bp"] * scale_to_cm
    papuan["length_cm"] = papuan["length_bp"] * scale_to_cm
    papuan["posterior_denisovan"] = float("nan")
    denisova = papuan["Shared_with_Denisova"]
    vindija = papuan["Shared_with_Vindija"]
    altai = papuan["Shared_with_Altai"]
    papuan["source_class"] = "unresolved"
    papuan.loc[vindija.gt(denisova), "source_class"] = "neanderthal_affinity"
    papuan.loc[denisova.gt(vindija), "source_class"] = "denisovan_affinity"
    papuan.loc[
        denisova.gt(vindija) & denisova.gt(altai),
        "source_class",
    ] = "denisovan_affinity_strict"
    papuan["caller"] = "skov_hmm_2018"
    papuan["qc_flags"] = "published_repeatmask_and_1000g_strictmask_applied_upstream"
    papuan["genetic_length_method"] = "constant_1.2e-8_recombination_per_bp_per_generation"
    papuan["source_class_rule"] = (
        "published_broad=Denisova>Vindija; strict=Denisova>Vindija_and_Altai"
    )
    return papuan


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
    if caller.lower() == "skov_hmm":
        frame = _adapt_skov_hmm(frame)
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
            "posterior_archaic", "posterior_denisovan", "caller", "source_class",
            "source_class_probability", "callable_fraction", "qc_flags",
            "snps", "country", "region", "Shared_with_Altai",
            "Shared_with_Denisova", "Shared_with_Vindija", "outgroup", "method",
            "genetic_length_method", "source_class_rule",
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
