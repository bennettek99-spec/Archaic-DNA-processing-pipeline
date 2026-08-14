"""Caller-aware calibration: what decay does the caller return for a known date?

The dating problem this addresses is not that the estimator is wrong but that
the *observation process* distorts what it sees. Posterior decoding merges
adjacent archaic runs and drops short ones, so the decoded tract-length
distribution decays more slowly than the true one and implies a date more
recent than the truth. The distortion cannot be divided out of the real data,
because measuring it there requires knowing the answer.

Simulation breaks the circle. For a known pulse time, this module runs the
whole chain -- coalescent, genotypes, outgroup-private window counts, HMM fit,
posterior decoding, tract extraction, decay estimate -- and records what comes
out. Sweeping the pulse time gives a calibration curve from true date to
observed decay, which is then inverted on the real measurement.

Three numbers are tracked per replicate and they are different things:

``true_decay``
    Exponential decay of the true introgressed intervals, recovered from a
    census placed above every archaic pulse. Should recover the input pulse
    time. Only populated when ``record_truth`` is set; the sweep does not need
    it, because the true pulse time is the input.
``fitted_generations``
    The HMM's own admixture-time parameter, from the fitted transition rate.
``decoded_decay``
    Decay of the runs posterior decoding actually returns. This is the
    quantity the real analysis measured, and the one the curve inverts.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .genotype_simulation import (
    PulseConfig,
    interval_lengths_morgans,
    simulate_replicate,
)
from .skov_hmm import call_individual, decay_generations

# Anchors measured on the real 89 Papuan individuals, for validation and for
# the inversion target. See SEGMENT_STRUCTURE_RESULTS.md.
REAL_DECODED_DECAY = 655.3
REAL_FITTED_PARAMETER = 1019.8
REAL_ARCHAIC_FRACTION = 0.0682
REAL_RATE_CONTRAST = 8.8
REAL_MODERN_RATE = 0.0256
REAL_ARCHAIC_RATE = 0.2245

DEFAULT_MIN_LENGTH_MORGANS = 5e-4      # 0.05 cM, matching the real analysis

# Variant density is scaled because the real callset is filtered for
# callability and the simulation is not. Fitted by matching the two Poisson
# rates and the decoded archaic fraction to the real analysis; see
# CALLER_CALIBRATION_RESULTS.md. Both entry points share it so a run cannot
# accidentally sit at a different operating point from the curve it inverts on.
CALIBRATED_MUTATION_SCALE = 0.40
BASE_MUTATION_RATE = 1.4e-8            # the Jacobs model's published value


@dataclass
class ReplicateResult:
    pulse_generations: float
    seed: int
    true_decay: float
    true_denisovan_decay: float
    fitted_generations: float
    decoded_decay: float
    decoded_over_fitted: float
    decoded_over_true: float
    archaic_fraction_true: float
    archaic_fraction_decoded: float
    modern_rate: float
    archaic_rate: float
    rate_contrast: float
    n_true_tracts: int
    n_decoded_tracts: int
    n_individuals: int
    sequence_length: int

    def to_row(self) -> dict[str, Any]:
        return self.__dict__.copy()


def _true_window_mask(intervals: dict[int, list[tuple[float, float]]],
                      shape: tuple[int, int], window_bp: int) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    for i, values in intervals.items():
        for left, right in values:
            mask[i, int(left) // window_bp:int(right) // window_bp] = True
    return mask


def run_scenario(
    pulse: PulseConfig,
    *,
    seed: int,
    sequence_length: int = 10_000_000,
    n_papuan: int = 20,
    n_outgroup: int = 100,
    recombination_rate: float = 1.2e-8,
    window_bp: int = 1000,
    min_length_morgans: float = DEFAULT_MIN_LENGTH_MORGANS,
    mutation_scale: float = CALIBRATED_MUTATION_SCALE,
    max_iter: int = 30,
) -> dict[str, Any]:
    """Run any pulse configuration through the caller and summarise it.

    Unlike :func:`run_replicate`, which exists to build the single-pulse
    calibration curve, this accepts mixtures and prolonged flow so that
    histories which are not a single pulse can be put on the same axis.
    """
    sim = simulate_replicate(
        pulse,
        seed=seed,
        sequence_length=sequence_length,
        n_papuan=n_papuan,
        n_outgroup=n_outgroup,
        recombination_rate=recombination_rate,
        window_bp=window_bp,
        mutation_rate=BASE_MUTATION_RATE * mutation_scale,
        record_truth=False,
    )
    counts = sim["counts"]

    fitted: list[float] = []
    decoded_lengths: list[np.ndarray] = []
    modern_rates: list[float] = []
    archaic_rates: list[float] = []
    decoded_windows = 0

    for i in range(counts.shape[0]):
        result = call_individual(
            counts[i],
            window_bp=window_bp,
            recombination_rate=recombination_rate,
            min_length_morgans=min_length_morgans,
            archaic_fraction=0.06,
            admixture_generations=1000.0,
            max_iter=max_iter,
        )
        fitted.append(result["fitted_generations"])
        decoded_lengths.append(result["lengths_morgans"])
        decoded_windows += int((result["ends"] - result["starts"]).sum())
        modern_rates.append(float(result["fit"].rates[0]))
        archaic_rates.append(float(result["fit"].rates[1]))

    pooled = np.concatenate(decoded_lengths) if decoded_lengths else np.empty(0)
    kept = pooled[pooled >= min_length_morgans]
    decoded_decay = decay_generations(pooled, min_length_morgans)
    fitted_median = float(np.median(fitted)) if fitted else float("nan")
    modern_rate = float(np.median(modern_rates))
    archaic_rate = float(np.median(archaic_rates))

    row: dict[str, Any] = dict(pulse.describe())
    row.update(
        seed=int(seed),
        fitted_generations=fitted_median,
        decoded_decay=decoded_decay,
        decoded_over_fitted=decoded_decay / fitted_median if fitted_median else float("nan"),
        archaic_fraction_decoded=decoded_windows / (counts.shape[0] * counts.shape[1]),
        modern_rate=modern_rate,
        archaic_rate=archaic_rate,
        rate_contrast=archaic_rate / modern_rate if modern_rate else float("nan"),
        n_decoded_tracts=int(kept.size),
        n_individuals=int(counts.shape[0]),
        sequence_length=int(sequence_length),
    )
    return row


def scenario_sweep(
    scenarios: Sequence[PulseConfig],
    *,
    replicates: int = 4,
    base_seed: int = 20260813,
    progress: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    rows = []
    total = len(scenarios) * replicates
    done = 0
    for index, pulse in enumerate(scenarios):
        for r in range(replicates):
            row = run_scenario(pulse, seed=base_seed + index * 977 + r, **kwargs)
            rows.append(row)
            done += 1
            if progress:
                print(
                    f"[{done:3d}/{total}] {row['pulse_label']:<28s} "
                    f"fitted={row['fitted_generations']:7.1f}  "
                    f"decoded={row['decoded_decay']:7.1f}  "
                    f"frac={row['archaic_fraction_decoded']:.4f}  "
                    f"n={row['n_decoded_tracts']:5d}",
                    flush=True,
                )
    return pd.DataFrame(rows)


def run_replicate(
    pulse_generations: float,
    *,
    seed: int,
    sequence_length: int = 10_000_000,
    n_papuan: int = 20,
    n_outgroup: int = 200,
    recombination_rate: float = 1.2e-8,
    window_bp: int = 1000,
    min_length_morgans: float = DEFAULT_MIN_LENGTH_MORGANS,
    mutation_scale: float = CALIBRATED_MUTATION_SCALE,
    max_iter: int = 30,
    record_truth: bool = False,
) -> ReplicateResult:
    """Simulate one replicate and push it through the caller."""
    pulse = PulseConfig(mode="single", generations=float(pulse_generations))
    sim = simulate_replicate(
        pulse,
        seed=seed,
        sequence_length=sequence_length,
        n_papuan=n_papuan,
        n_outgroup=n_outgroup,
        recombination_rate=recombination_rate,
        window_bp=window_bp,
        mutation_rate=BASE_MUTATION_RATE * mutation_scale,
        record_truth=record_truth,
    )
    counts = sim["counts"]

    if record_truth:
        true_lengths = interval_lengths_morgans(sim["true_archaic"], recombination_rate)
        true_den = interval_lengths_morgans(sim["true_denisovan"], recombination_rate)
        true_decay = decay_generations(true_lengths, min_length_morgans)
        true_den_decay = decay_generations(true_den, min_length_morgans)
        true_fraction = float(
            _true_window_mask(sim["true_archaic"], counts.shape, window_bp).mean()
        )
        n_true = int((true_lengths >= min_length_morgans).sum())
    else:
        true_decay = float("nan")
        true_den_decay = float("nan")
        true_fraction = float("nan")
        n_true = 0

    fitted: list[float] = []
    decoded_lengths: list[np.ndarray] = []
    modern_rates: list[float] = []
    archaic_rates: list[float] = []
    decoded_windows = 0

    for i in range(counts.shape[0]):
        result = call_individual(
            counts[i],
            window_bp=window_bp,
            recombination_rate=recombination_rate,
            min_length_morgans=min_length_morgans,
            archaic_fraction=0.06,
            admixture_generations=1000.0,
            max_iter=max_iter,
        )
        fitted.append(result["fitted_generations"])
        decoded_lengths.append(result["lengths_morgans"])
        decoded_windows += int((result["ends"] - result["starts"]).sum())
        modern_rates.append(float(result["fit"].rates[0]))
        archaic_rates.append(float(result["fit"].rates[1]))

    # The real anchors (0.0256 / 0.2245 per kb, 8.8x) were themselves estimated
    # by fitting a two-component Poisson mixture, so comparing against the
    # fitted rates rather than the simulation's truth is the like-for-like check.
    modern_rate = float(np.median(modern_rates))
    archaic_rate = float(np.median(archaic_rates))

    pooled = np.concatenate(decoded_lengths) if decoded_lengths else np.empty(0)
    kept = pooled[pooled >= min_length_morgans]
    decoded_decay = decay_generations(pooled, min_length_morgans)
    fitted_median = float(np.median(fitted)) if fitted else float("nan")

    total_windows = counts.shape[0] * counts.shape[1]
    return ReplicateResult(
        pulse_generations=float(pulse_generations),
        seed=int(seed),
        true_decay=true_decay,
        true_denisovan_decay=true_den_decay,
        fitted_generations=fitted_median,
        decoded_decay=decoded_decay,
        decoded_over_fitted=decoded_decay / fitted_median if fitted_median else float("nan"),
        decoded_over_true=decoded_decay / true_decay if true_decay else float("nan"),
        archaic_fraction_true=true_fraction,
        archaic_fraction_decoded=decoded_windows / total_windows,
        modern_rate=modern_rate,
        archaic_rate=archaic_rate,
        rate_contrast=archaic_rate / modern_rate if modern_rate else float("nan"),
        n_true_tracts=n_true,
        n_decoded_tracts=int(kept.size),
        n_individuals=int(counts.shape[0]),
        sequence_length=int(sequence_length),
    )


def sweep(
    pulse_times: Sequence[float],
    *,
    replicates: int = 3,
    base_seed: int = 20260813,
    progress: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run the calibration sweep over true pulse times."""
    rows = []
    total = len(pulse_times) * replicates
    done = 0
    for t in pulse_times:
        for r in range(replicates):
            seed = base_seed + int(t) * 1000 + r
            result = run_replicate(t, seed=seed, **kwargs)
            rows.append(result.to_row())
            done += 1
            if progress:
                print(
                    f"[{done:3d}/{total}] t={t:7.1f}  true={result.true_decay:7.1f}  "
                    f"fitted={result.fitted_generations:7.1f}  "
                    f"decoded={result.decoded_decay:7.1f}  "
                    f"contrast={result.rate_contrast:4.1f}x  "
                    f"n={result.n_decoded_tracts:5d}",
                    flush=True,
                )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# inversion
# --------------------------------------------------------------------------


def fit_curve(table: pd.DataFrame, x: str = "pulse_generations",
              y: str = "decoded_decay") -> np.ndarray:
    """Least-squares line through the calibration points.

    A straight line is used deliberately. The relationship between true date
    and decoded decay is close to proportional over the tested range, and a
    flexible fit on a handful of noisy points would invent structure the
    replicates cannot support.
    """
    grouped = table.groupby(x)[y].mean()
    return np.polyfit(grouped.index.to_numpy(float), grouped.to_numpy(float), 1)


def invert(
    table: pd.DataFrame,
    observed_decay: float,
    *,
    n_boot: int = 2000,
    rng_seed: int = 7,
) -> dict[str, Any]:
    """Invert the calibration curve on an observed decoded decay.

    Uncertainty comes from resampling replicates within each pulse time, which
    propagates the simulation's own scatter into the returned interval. It does
    not include uncertainty in the real measurement or any error in the
    demographic model.
    """
    coeffs = fit_curve(table)
    slope, intercept = float(coeffs[0]), float(coeffs[1])

    # A slope near zero is not merely numerically awkward: it means the caller's
    # output does not respond to the pulse time, so nothing can be inverted.
    # Test the span the curve covers rather than the slope alone, since the
    # slope's scale depends on the units of both axes.
    times = table["pulse_generations"].to_numpy(float)
    span = abs(slope) * (times.max() - times.min())
    scale = float(table["decoded_decay"].mean())
    if not np.isfinite(span) or span < 0.01 * abs(scale):
        raise ValueError(
            "degenerate calibration curve: decoded decay does not respond to "
            f"pulse time (span {span:.4g} over decoded scale {scale:.4g})"
        )
    point = (observed_decay - intercept) / slope

    rng = np.random.default_rng(rng_seed)
    groups = {t: g for t, g in table.groupby("pulse_generations")}
    draws = np.empty(n_boot)
    for b in range(n_boot):
        rows = []
        for t, g in groups.items():
            idx = rng.integers(0, len(g), len(g))
            rows.append(g.iloc[idx])
        resampled = pd.concat(rows)
        c = fit_curve(resampled)
        draws[b] = (observed_decay - c[1]) / c[0] if c[0] != 0 else np.nan

    draws = draws[np.isfinite(draws)]
    return {
        "observed_decay": float(observed_decay),
        "slope": slope,
        "intercept": intercept,
        "point_estimate_generations": float(point),
        "ci95_generations": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "n_bootstrap": int(draws.size),
    }


def required_inflation(
    table: pd.DataFrame,
    observed_decay: float,
    candidate_generations: float,
) -> dict[str, Any]:
    """How much more distortion would the caller need to explain the data?

    If the real measurement lies outside the simulated range, the useful
    question is not "what date does the curve return" but "what would have to
    be true of the observation process for a given candidate date to produce
    this measurement". That converts a failed inversion into a quantified
    requirement on the forward model, which is testable.

    Returns the extra compression factor needed on top of whatever the
    simulation already produces, and the total decoder inflation it implies.
    """
    coeffs = fit_curve(table)
    predicted = float(np.polyval(coeffs, candidate_generations))
    simulated_ratio = float(table["decoded_over_fitted"].mean())
    extra = observed_decay / predicted if predicted else float("nan")
    return {
        "candidate_generations": float(candidate_generations),
        "candidate_kya": float(candidate_generations) * 29.0 / 1000.0,
        "predicted_decoded_decay": predicted,
        "observed_decoded_decay": float(observed_decay),
        "extra_compression_needed": extra,
        "simulated_inflation": 1.0 / simulated_ratio if simulated_ratio else float("nan"),
        "required_total_inflation": (1.0 / simulated_ratio) / extra
        if simulated_ratio and extra
        else float("nan"),
        "real_measured_inflation": REAL_FITTED_PARAMETER / REAL_DECODED_DECAY,
    }


def save(table: pd.DataFrame, inversion: dict[str, Any], directory: Path) -> None:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    table.to_csv(directory / "calibration_sweep.tsv", sep="\t", index=False)
    (directory / "inversion.json").write_text(
        json.dumps(inversion, indent=2), encoding="utf-8"
    )
