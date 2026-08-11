"""Deterministic tract-level simulation calibration for models M1-M10."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .checkpointing import CheckpointStore, Deadline
from .tract_schema import write_tracts
from .tract_summary import summarize_tracts


def derived_seed(master_seed: int, model_id: str, replicate: int) -> int:
    payload = f"{int(master_seed)}:{model_id}:{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little") or 1


def _rates_for_model(model: dict[str, Any], generation_time_years: float, count: int, rng) -> np.ndarray:
    dates = np.asarray(model["dates_kya"], dtype=float)
    rates = dates * 1000.0 / generation_time_years
    kind = model["kind"]
    if kind == "single":
        return np.repeat(rates[0], count)
    if kind in {"two", "modern_mixing", "divergent_sources"}:
        weights = np.asarray(model.get("weights", [0.5, 0.5]), dtype=float)
        weights /= weights.sum()
        return rng.choice(rates, size=count, p=weights)
    if kind == "continuous":
        older, younger = max(rates), min(rates)
        return rng.uniform(younger, older, size=count)
    if kind in {"bottleneck", "selection"}:
        return np.repeat(rates[0], count)
    raise ValueError(f"Unsupported simulation model kind {kind!r}")


def simulate_tracts(
    model_id: str,
    model: dict[str, Any],
    *,
    n_tracts: int,
    generation_time_years: float,
    minimum_length_cm: float,
    seed: int,
    false_negative_rate: float = 0.0,
    length_noise_sd_fraction: float = 0.0,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rates = _rates_for_model(model, generation_time_years, n_tracts, rng)
    excess_cm = rng.exponential(1.0 / rates) * 100.0
    kind = model["kind"]
    if kind == "bottleneck":
        long_tail = rng.random(n_tracts) < 0.18
        excess_cm[long_tail] *= rng.uniform(1.3, 2.0, long_tail.sum())
    elif kind == "selection":
        selected = rng.random(n_tracts) < float(model.get("selected_fraction", 0.08))
        excess_cm[selected] *= rng.uniform(2.0, 4.0, selected.sum())
    if length_noise_sd_fraction > 0:
        excess_cm *= rng.lognormal(mean=0, sigma=length_noise_sd_fraction, size=n_tracts)
    lengths_cm = minimum_length_cm + excess_cm
    retained = rng.random(n_tracts) >= false_negative_rate
    lengths_cm = lengths_cm[retained]
    count = len(lengths_cm)
    chromosome = rng.integers(1, 23, size=count).astype(str)
    sample_number = rng.integers(1, 21, size=count)
    start_bp = rng.integers(1_000_000, 200_000_000, size=count, dtype=np.int64)
    length_bp = np.maximum(1000, np.rint(lengths_cm * 1_000_000).astype(np.int64))
    start_cm = start_bp / 1_000_000.0
    frame = pd.DataFrame(
        {
            "sample_id": [f"PAPUAN_SYNTH_{value:02d}" for value in sample_number],
            "population": "Papuan_synthetic",
            "chromosome": chromosome,
            "start_bp": start_bp,
            "end_bp": start_bp + length_bp,
            "start_cm": start_cm,
            "end_cm": start_cm + lengths_cm,
            "length_bp": length_bp,
            "length_cm": lengths_cm,
            "posterior_denisovan": rng.uniform(0.82, 1.0, size=count),
            "caller": "synthetic_truth",
            "source_class": "unresolved",
            "source_class_probability": np.nan,
            "callable_fraction": rng.uniform(0.90, 1.0, size=count),
            "qc_flags": "",
            "simulation_model": model_id,
            "simulation_seed": seed,
        }
    )
    return frame


def simulate_model_set(
    config: dict[str, Any],
    output_dir: str | Path,
    *,
    replicates: int,
    resume: bool,
    deadline: Deadline,
    config_digest: str,
) -> tuple[pd.DataFrame, bool]:
    output = Path(output_dir)
    simulation_dir = output / "simulations"
    simulation_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = CheckpointStore(
        output / "checkpoints" / "simulations.json",
        "simulations",
        config_digest,
        int(config["project"]["random_seed"]),
    )
    state = checkpoint.load(resume=resume)
    rows: list[dict[str, Any]] = []
    paused = False
    settings = config["simulation"]
    error = settings.get("caller_error", {})
    for model_id, model in config["models"].items():
        for replicate in range(int(replicates)):
            unit = f"{model_id}:{replicate}"
            target = simulation_dir / f"{model_id}_replicate_{replicate:04d}.tsv"
            if resume and checkpoint.unit_valid(state, unit):
                tracts = pd.read_csv(target, sep="\t")
            else:
                if deadline.should_stop(estimated_next_seconds=2):
                    paused = True
                    break
                seed = derived_seed(config["project"]["random_seed"], model_id, replicate)
                try:
                    if settings.get("engine") == "msprime":
                        from .msprime_backend import simulate_msprime_tracts

                        backend = settings.get("msprime", {})
                        chromosomes = settings.get("pilot_chromosomes", [21, 22])
                        chromosome = str(chromosomes[replicate % len(chromosomes)])
                        tracts = simulate_msprime_tracts(
                            model_id,
                            model,
                            seed=seed,
                            generation_time_years=float(config["project"]["generation_time_years"]),
                            sequence_length_bp=int(backend.get("sequence_length_bp", 20_000_000)),
                            recombination_rate=float(backend.get("recombination_rate", 1e-8)),
                            sample_individuals=int(backend.get("sample_individuals", 20)),
                            introgression_fraction=float(backend.get("introgression_fraction", 0.04)),
                            chromosome=chromosome,
                            minimum_length_cm=float(config["tracts"]["minimum_length_cm"]),
                            false_negative_rate=float(error.get("false_negative_rate", 0)),
                            length_noise_sd_fraction=float(error.get("length_noise_sd_fraction", 0)),
                        )
                        if tracts.empty:
                            raise RuntimeError(
                                f"{model_id} replicate {replicate} produced no detectable "
                                "msprime tracts; increase sequence length or sample count"
                            )
                    else:
                        tracts = simulate_tracts(
                            model_id,
                            model,
                            n_tracts=int(settings["observed_tract_count"]),
                            generation_time_years=float(config["project"]["generation_time_years"]),
                            minimum_length_cm=float(config["tracts"]["minimum_length_cm"]),
                            seed=seed,
                            false_negative_rate=float(error.get("false_negative_rate", 0)),
                            length_noise_sd_fraction=float(error.get("length_noise_sd_fraction", 0)),
                        )
                    write_tracts(tracts, target)
                    checkpoint.mark_completed(state, unit, [target])
                except Exception as exc:
                    checkpoint.mark_failed(state, unit, exc)
                    raise
            summary = summarize_tracts(tracts)["overall"].iloc[0].to_dict()
            rows.append({"model_id": model_id, "replicate": replicate, **summary})
        if paused:
            break
    summary_frame = pd.DataFrame(rows)
    summary_path = simulation_dir / "simulation_summaries.tsv"
    if not summary_frame.empty:
        write_tracts(summary_frame, summary_path)
    state["state"] = "paused" if paused else "complete"
    checkpoint.save(state)
    return summary_frame, paused
