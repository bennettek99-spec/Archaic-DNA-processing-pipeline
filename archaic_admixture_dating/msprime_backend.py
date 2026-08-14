"""Coalescent pilot backend that extracts introgressed migration tracts.

The backend deliberately models a compact chromosome-sized sequence and records
backward-time migrations into one or two Denisovan-related populations.
Intervals carried by those migration records are mapped to descendant Papuan
sample haplotypes and merged per diploid individual. Selection and caller error
remain explicit post-coalescent sensitivity approximations.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np
import pandas as pd


def _demography(model: dict[str, Any], generation_time_years: float, introgression_fraction: float):
    try:
        import msprime
    except ImportError as error:
        raise RuntimeError("The msprime engine requires the repository's 'sim' extra") from error

    demo = msprime.Demography()
    for name, size in (
        ("PAP", 10_000),
        ("EAS", 10_000),
        ("AFR", 15_000),
        ("HUM", 10_000),
        ("DEN1", 3_000),
        ("DEN2", 3_000),
        ("ARCH", 4_000),
        ("ROOT", 10_000),
    ):
        demo.add_population(name=name, initial_size=size)

    dates = [float(value) * 1000.0 / generation_time_years for value in model["dates_kya"]]
    kind = model["kind"]
    if kind in {"single", "bottleneck", "selection"}:
        demo.add_mass_migration(
            time=dates[0], source="PAP", dest="DEN1", proportion=introgression_fraction
        )
    elif kind in {"two", "divergent_sources"}:
        weights = model.get("weights", [0.5, 0.5])
        destinations = ["DEN1", "DEN2"] if kind == "divergent_sources" else ["DEN1", "DEN1"]
        for date, weight, destination in zip(dates, weights, destinations):
            demo.add_mass_migration(
                time=date,
                source="PAP",
                dest=destination,
                proportion=introgression_fraction * float(weight),
            )
    elif kind == "continuous":
        older, younger = max(dates), min(dates)
        for date in np.linspace(younger, older, 7):
            demo.add_mass_migration(
                time=float(date),
                source="PAP",
                dest="DEN1",
                proportion=introgression_fraction / 7.0,
            )
    elif kind == "modern_mixing":
        older, recent = dates
        demo.add_mass_migration(
            time=older, source="PAP", dest="DEN1", proportion=introgression_fraction
        )
        demo.add_mass_migration(
            time=older, source="EAS", dest="DEN1", proportion=introgression_fraction * 0.25
        )
        demo.add_mass_migration(
            time=recent, source="PAP", dest="EAS", proportion=float(model.get("weights", [0.75, 0.25])[1])
        )
    else:
        raise ValueError(f"Unsupported msprime model kind {kind!r}")

    if kind == "bottleneck":
        demo.add_population_parameters_change(time=700, initial_size=800, population="PAP")
        demo.add_population_parameters_change(time=1400, initial_size=10_000, population="PAP")

    demo.add_population_split(
        time=70_000 / generation_time_years,
        derived=["PAP", "EAS", "AFR"],
        ancestral="HUM",
    )
    demo.add_population_split(
        time=350_000 / generation_time_years,
        derived=["DEN1", "DEN2"],
        ancestral="ARCH",
    )
    demo.add_population_split(
        time=600_000 / generation_time_years,
        derived=["HUM", "ARCH"],
        ancestral="ROOT",
    )
    demo.sort_events()
    return demo


def _merge(intervals: list[tuple[float, float]]) -> list[tuple[float, float]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for left, right in ordered[1:]:
        old_left, old_right = merged[-1]
        if left <= old_right + 1:
            merged[-1] = (old_left, max(old_right, right))
        else:
            merged.append((left, right))
    return merged


def _extract(ts, archaic_population_ids: set[int]) -> dict[int, list[tuple[float, float]]]:
    """Introgressed intervals per individual, from migration records.

    .. warning::

       This over-attributes. A migration record's interval is the span of the
       ancestral lineage when it moved, which is wider than the segment any one
       modern sample inherits from it, so the intervals come out too long and
       the measured archaic fraction depends on how much sequence was
       simulated. Measured directly: fraction 0.158 at 10 Mb against 0.093 at
       30 Mb, and a 1400-generation pulse recovered as 921 and 1180
       generations respectively.

       It is left in place because the tract-level M1-M10 workflow was
       validated against it and its outputs are already published as such, but
       it must not be used as truth for anything quantitative.
       :func:`archaic_admixture_dating.genotype_simulation._census_intervals`
       is the corrected implementation: it recovers the simulated proportion to
       0.0400 +/- 0.0089 against a target of 0.0400, with no dependence on
       sequence length.
    """
    intervals: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for migration in ts.migrations():
        if migration.dest not in archaic_population_ids:
            continue
        midpoint = migration.left + (migration.right - migration.left) / 2.0
        tree = ts.at(midpoint)
        try:
            descendants = tree.samples(migration.node)
        except ValueError:
            continue
        for sample_node in descendants:
            individual = ts.node(sample_node).individual
            if individual >= 0:
                intervals[individual].append((migration.left, migration.right))
    return {individual: _merge(values) for individual, values in intervals.items()}


def simulate_msprime_tracts(
    model_id: str,
    model: dict[str, Any],
    *,
    seed: int,
    generation_time_years: float,
    sequence_length_bp: int,
    recombination_rate: float,
    sample_individuals: int,
    introgression_fraction: float,
    chromosome: str = "21",
    minimum_length_cm: float = 0.02,
    false_negative_rate: float = 0.0,
    length_noise_sd_fraction: float = 0.0,
) -> pd.DataFrame:
    try:
        import msprime
    except ImportError as error:
        raise RuntimeError("The msprime engine requires the repository's 'sim' extra") from error

    demo = _demography(model, generation_time_years, introgression_fraction)
    ts = msprime.sim_ancestry(
        samples={"PAP": int(sample_individuals)},
        ploidy=2,
        demography=demo,
        sequence_length=int(sequence_length_bp),
        recombination_rate=float(recombination_rate),
        random_seed=int(seed),
        record_migrations=True,
    )
    populations = {population.metadata["name"]: population.id for population in ts.populations()}
    intervals = _extract(ts, {populations["DEN1"], populations["DEN2"]})
    rng = np.random.default_rng(seed + 1)
    rows: list[dict[str, Any]] = []
    for individual, values in intervals.items():
        for left, right in values:
            length_bp = max(1, int(round(right - left)))
            length_cm = length_bp * recombination_rate * 100.0
            if model["kind"] == "selection" and rng.random() < float(model.get("selected_fraction", 0.08)):
                length_cm *= rng.uniform(2.0, 4.0)
                length_bp = int(round(length_cm / (recombination_rate * 100.0)))
            if length_noise_sd_fraction:
                scale = rng.lognormal(0, length_noise_sd_fraction)
                length_cm *= scale
                length_bp = int(round(length_bp * scale))
            if length_cm < minimum_length_cm or rng.random() < false_negative_rate:
                continue
            start_bp = int(round(left))
            rows.append(
                {
                    "sample_id": f"PAPUAN_MSPRIME_{individual + 1:02d}",
                    "population": "Papuan_synthetic",
                    "chromosome": str(chromosome),
                    "start_bp": start_bp,
                    "end_bp": min(int(sequence_length_bp), start_bp + length_bp),
                    "start_cm": start_bp * recombination_rate * 100.0,
                    "end_cm": start_bp * recombination_rate * 100.0 + length_cm,
                    "length_bp": length_bp,
                    "length_cm": length_cm,
                    "posterior_denisovan": 1.0,
                    "caller": "msprime_migration_truth",
                    "source_class": "unresolved",
                    "source_class_probability": np.nan,
                    "callable_fraction": 1.0,
                    "qc_flags": "",
                    "simulation_model": model_id,
                    "simulation_seed": seed,
                    "simulator_version": msprime.__version__,
                }
            )
    return pd.DataFrame(rows)
