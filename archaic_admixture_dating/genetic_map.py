"""Genetic-map loading and tract-coordinate interpolation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .tract_schema import normalize_chromosome


DEFAULT_MAP_PATTERN = "genetic_map_GRCh37_chr{chromosome}.txt.gz"


def load_genetic_map(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a SHAPEIT/IMPUTE-style map and validate monotonic coordinates."""
    source = Path(path)
    frame = pd.read_csv(source, sep=r"\s+")
    normalized = {
        str(column).strip().casefold().replace(" ", ""): column
        for column in frame.columns
    }
    position_column = normalized.get("position(bp)") or normalized.get("position")
    map_column = normalized.get("map(cm)") or normalized.get("map")
    if position_column is None or map_column is None:
        raise ValueError(
            f"Genetic map {source.name} must contain Position(bp) and Map(cM) columns"
        )
    positions = pd.to_numeric(frame[position_column], errors="coerce").to_numpy(float)
    map_cm = pd.to_numeric(frame[map_column], errors="coerce").to_numpy(float)
    finite = np.isfinite(positions) & np.isfinite(map_cm)
    positions = positions[finite]
    map_cm = map_cm[finite]
    if len(positions) < 2:
        raise ValueError(f"Genetic map {source.name} has fewer than two usable positions")
    if np.any(np.diff(positions) <= 0):
        raise ValueError(f"Genetic map {source.name} positions are not strictly increasing")
    if np.any(np.diff(map_cm) < -1e-9):
        raise ValueError(f"Genetic map {source.name} cM coordinates are not monotonic")
    return positions, map_cm


def apply_genetic_map(
    frame: pd.DataFrame,
    map_directory: str | Path,
    *,
    pattern: str = DEFAULT_MAP_PATTERN,
    build: str = "GRCh37",
) -> pd.DataFrame:
    """Interpolate tract endpoints onto a chromosome-specific genetic map.

    Tracts extending beyond a map's covered interval receive missing genetic
    coordinates and therefore fail closed during standard tract validation.
    """
    work = frame.copy()
    directory = Path(map_directory)
    if not directory.is_dir():
        raise ValueError(f"Genetic-map directory does not exist: {directory}")

    chromosomes = work["chromosome"].map(normalize_chromosome)
    starts = pd.to_numeric(work["start_bp"], errors="coerce")
    ends = pd.to_numeric(work["end_bp"], errors="coerce")
    work["chromosome"] = chromosomes
    work["start_cm"] = np.nan
    work["end_cm"] = np.nan
    work["length_cm"] = np.nan
    work["genetic_map_status"] = "not_interpolated"

    for chromosome in sorted(chromosomes.dropna().unique()):
        map_path = directory / pattern.format(chromosome=chromosome)
        if not map_path.is_file():
            raise ValueError(
                f"Missing genetic map for chromosome {chromosome}: {map_path.name}"
            )
        positions, map_cm = load_genetic_map(map_path)
        selected = chromosomes.eq(chromosome)
        covered = (
            selected
            & starts.ge(positions[0])
            & ends.le(positions[-1])
            & starts.notna()
            & ends.notna()
        )
        work.loc[selected & ~covered, "genetic_map_status"] = "outside_map_range"
        if covered.any():
            start_cm = np.interp(starts.loc[covered].to_numpy(float), positions, map_cm)
            end_cm = np.interp(ends.loc[covered].to_numpy(float), positions, map_cm)
            work.loc[covered, "start_cm"] = start_cm
            work.loc[covered, "end_cm"] = end_cm
            work.loc[covered, "length_cm"] = end_cm - start_cm
            work.loc[covered, "genetic_map_status"] = "interpolated"

    work["genetic_length_method"] = (
        f"linear_interpolation_{build}_chromosome_specific_genetic_map"
    )
    work["genetic_map_build"] = build
    work["genetic_map_pattern"] = pattern
    return work
