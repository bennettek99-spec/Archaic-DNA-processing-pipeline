"""Required tract-date sensitivity axes with explicit unavailable states."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .dating_single_pulse import fit_single_pulse


def _fit_row(
    subset: pd.DataFrame,
    *,
    axis: str,
    level: str,
    minimum_length_cm: float,
    generation_time_years: float,
    bounds: tuple[float, float],
    length_scale: float = 1.0,
) -> dict[str, object]:
    if len(subset) < 10:
        return {
            "axis": axis,
            "level": level,
            "n_tracts": len(subset),
            "status": "insufficient_tracts",
        }
    try:
        fit = fit_single_pulse(
            subset["length_cm"].to_numpy(dtype=float) * length_scale,
            minimum_length_cm=minimum_length_cm * length_scale,
            generation_time_years=generation_time_years,
            bounds_generations=bounds,
        )
    except ValueError as error:
        return {
            "axis": axis,
            "level": level,
            "n_tracts": len(subset),
            "status": "not_estimable",
            "warning_flags": str(error),
        }
    return {
        "axis": axis,
        "level": level,
        "n_tracts": len(subset),
        "status": "ok",
        "generations": fit["generations"],
        "kya": fit["kya"],
        "ci_low_kya": fit["ci_low_kya"],
        "ci_high_kya": fit["ci_high_kya"],
        "warning_flags": ";".join(fit["warning_flags"]),
    }


def run_sensitivity(tracts: pd.DataFrame, config: dict) -> pd.DataFrame:
    dating = config["dating"]
    tract_config = config["tracts"]
    bounds = tuple(dating["single_pulse"]["bounds_generations"])
    default_min = float(tract_config["minimum_length_cm"])
    default_gen = float(config["project"]["generation_time_years"])
    rows: list[dict[str, object]] = []

    for generation_time in dating["generation_time_sensitivity"]:
        rows.append(
            _fit_row(
                tracts,
                axis="generation_time",
                level=str(generation_time),
                minimum_length_cm=default_min,
                generation_time_years=float(generation_time),
                bounds=bounds,
            )
        )
    for threshold in dating["minimum_tract_length_sensitivity_cm"]:
        subset = tracts.loc[tracts["length_cm"] >= float(threshold)]
        rows.append(
            _fit_row(
                subset,
                axis="minimum_length_cm",
                level=str(threshold),
                minimum_length_cm=float(threshold),
                generation_time_years=default_gen,
                bounds=bounds,
            )
        )
    for confidence in dating["confidence_sensitivity"]:
        known = tracts["posterior_denisovan"].notna()
        subset = tracts.loc[~known | (tracts["posterior_denisovan"] >= float(confidence))]
        rows.append(
            _fit_row(
                subset,
                axis="minimum_confidence",
                level=str(confidence),
                minimum_length_cm=default_min,
                generation_time_years=default_gen,
                bounds=bounds,
            )
        )
    for scale in dating.get("recombination_map_scale_sensitivity", [0.95, 1.0, 1.05]):
        rows.append(
            _fit_row(
                tracts,
                axis="recombination_map_scale",
                level=str(scale),
                minimum_length_cm=default_min,
                generation_time_years=default_gen,
                bounds=bounds,
                length_scale=float(scale),
            )
        )
    for chromosome in sorted(tracts["chromosome"].astype(str).unique(), key=lambda x: (len(x), x)):
        subset = tracts.loc[tracts["chromosome"].astype(str) != chromosome]
        rows.append(
            _fit_row(
                subset,
                axis="leave_one_chromosome_out",
                level=chromosome,
                minimum_length_cm=default_min,
                generation_time_years=default_gen,
                bounds=bounds,
            )
        )
    for fraction in (0.01, 0.05):
        keep = max(1, int(np.floor(len(tracts) * (1 - fraction))))
        subset = tracts.nsmallest(keep, "length_cm")
        rows.append(
            _fit_row(
                subset,
                axis="exclude_longest_fraction",
                level=str(fraction),
                minimum_length_cm=default_min,
                generation_time_years=default_gen,
                bounds=bounds,
            )
        )
    selected = tracts["qc_flags"].fillna("").astype(str).str.contains("selected", case=False)
    if selected.any():
        rows.append(
            _fit_row(
                tracts.loc[~selected],
                axis="exclude_candidate_selected_loci",
                level="qc_flags_selected",
                minimum_length_cm=default_min,
                generation_time_years=default_gen,
                bounds=bounds,
            )
        )
    else:
        rows.append(
            {
                "axis": "exclude_candidate_selected_loci",
                "level": "unavailable",
                "n_tracts": len(tracts),
                "status": "not_available_no_selected_locus_annotation",
            }
        )
    return pd.DataFrame(rows)
