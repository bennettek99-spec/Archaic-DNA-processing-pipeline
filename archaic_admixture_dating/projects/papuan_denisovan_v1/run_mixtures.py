"""Can a mixture or prolonged gene flow reproduce the observed Papuan decay?

The single-pulse calibration
([CALLER_CALIBRATION_RESULTS.md](CALLER_CALIBRATION_RESULTS.md)) put the real
decoded decay of 655.3 generations below every single pulse from 14.5 to
51.6 kya. Mixtures are the obvious next candidate, and there is a specific
reason to doubt them: for a mixture of exponentials the mean tract length is
the weighted mean of the component means, so the implied decay is bounded
between the components. A mixture of pulses inside the swept range should not
be able to fall below the youngest of them.

That argument is about *true* tract lengths, and the quantity being measured is
a *decoded* one. A mixture lays down more archaic segments in close proximity,
which gives posterior decoding more opportunity to bridge and merge them, and
merging is length-creating. So the bound can in principle be beaten by the
observation process even though it holds for the underlying tracts. That is
what this run tests.

Each scenario is also inverted through the single-pulse calibration curve, so
the output says what apparent single-pulse date each mixture would be mistaken
for -- the identifiability question, on the same axis.

Usage::

    python archaic_admixture_dating/projects/papuan_denisovan_v1/run_mixtures.py \
        --output archaic_admixture_dating/outputs/papuan_mixtures_v1
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from archaic_admixture_dating.caller_calibration import (
    REAL_DECODED_DECAY,
    fit_curve,
    scenario_sweep,
)
from archaic_admixture_dating.genotype_simulation import PulseConfig

GENERATION_TIME_YEARS = 29.0

OPERATING_POINT = dict(
    sequence_length=10_000_000,
    n_papuan=20,
    n_outgroup=100,
    mutation_scale=0.40,
    recombination_rate=1.2e-8,
    window_bp=1000,
)


def _kya(generations: float) -> float:
    return generations * GENERATION_TIME_YEARS / 1000.0


SCENARIOS = [
    PulseConfig(mode="single", generations=1550.0, label="single 45 kya"),
    PulseConfig(mode="single", generations=500.0, label="single 14.5 kya (floor)"),
    PulseConfig(mode="published", label="published Jacobs two-pulse"),
    PulseConfig(mode="two", generations=(1550.0, 1000.0), weights=(0.5, 0.5),
                label="two 45+29 kya 50/50"),
    PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5),
                label="two 45+17 kya 50/50"),
    PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.25, 0.75),
                label="two 45+17 kya 25/75"),
    PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.75, 0.25),
                label="two 45+17 kya 75/25"),
    PulseConfig(mode="two", generations=(1724.0, 1000.0), weights=(0.5, 0.5),
                label="two 50+29 kya 50/50"),
    PulseConfig(mode="continuous", generations=(1550.0, 600.0), n_bins=7,
                label="continuous 45->17 kya"),
    PulseConfig(mode="continuous", generations=(1724.0, 1035.0), n_bins=7,
                label="continuous 50->30 kya"),
]


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--base-seed", type=int, default=20260814)
    parser.add_argument(
        "--calibration",
        type=Path,
        default=Path("archaic_admixture_dating/outputs/papuan_caller_calibration_v1")
        / "calibration_sweep.tsv",
        help="single-pulse sweep used to convert decoded decay to apparent date",
    )
    args = parser.parse_args(argv)

    table = scenario_sweep(
        SCENARIOS,
        replicates=args.replicates,
        base_seed=args.base_seed,
        **OPERATING_POINT,
    )

    grouped = (
        table.groupby("pulse_label", sort=False)
        .agg(
            decoded_decay=("decoded_decay", "mean"),
            decoded_sd=("decoded_decay", "std"),
            fitted=("fitted_generations", "mean"),
            fraction=("archaic_fraction_decoded", "mean"),
            n=("n_decoded_tracts", "mean"),
        )
        .reset_index()
    )

    apparent = None
    if args.calibration.exists():
        curve = fit_curve(pd.read_csv(args.calibration, sep="\t"))
        slope, intercept = float(curve[0]), float(curve[1])
        apparent = (grouped["decoded_decay"] - intercept) / slope
        grouped["apparent_single_pulse_gen"] = apparent
        grouped["apparent_single_pulse_kya"] = apparent * GENERATION_TIME_YEARS / 1000.0

    grouped["distance_from_real"] = grouped["decoded_decay"] - REAL_DECODED_DECAY

    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "mixture_sweep.tsv", sep="\t", index=False)
    grouped.to_csv(args.output / "mixture_summary.tsv", sep="\t", index=False,
                   float_format="%.4f")

    reaches = grouped[grouped["decoded_decay"] <= REAL_DECODED_DECAY]
    verdict = {
        "real_decoded_decay": REAL_DECODED_DECAY,
        "lowest_scenario": grouped.loc[grouped["decoded_decay"].idxmin(), "pulse_label"],
        "lowest_decoded_decay": float(grouped["decoded_decay"].min()),
        "any_scenario_reaches_real": bool(len(reaches) > 0),
        "scenarios_reaching_real": reaches["pulse_label"].tolist(),
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "git_commit": _git_commit(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "operating_point": OPERATING_POINT,
            "replicates": args.replicates,
            "base_seed": args.base_seed,
        },
    }
    (args.output / "verdict.json").write_text(json.dumps(verdict, indent=2),
                                              encoding="utf-8")

    summary_path = Path(__file__).parent / "MIXTURE_SUMMARY.tsv"
    grouped.to_csv(summary_path, sep="\t", index=False, float_format="%.4f")

    print()
    print("=" * 88)
    print("MIXTURE AND PROLONGED-FLOW SCENARIOS")
    print("=" * 88)
    header = f"{'scenario':<30} {'decoded':>9} {'+/-':>7} {'fitted':>9} {'frac':>7}"
    if apparent is not None:
        header += f" {'apparent kya':>13}"
    print(header)
    for _, row in grouped.iterrows():
        line = (
            f"{row['pulse_label']:<30} {row['decoded_decay']:9.1f} "
            f"{row['decoded_sd']:7.1f} {row['fitted']:9.1f} {row['fraction']:7.4f}"
        )
        if apparent is not None:
            line += f" {row['apparent_single_pulse_kya']:13.1f}"
        print(line)

    print()
    print(f"real decoded decay: {REAL_DECODED_DECAY}")
    print(f"lowest scenario:    {verdict['lowest_scenario']} "
          f"at {verdict['lowest_decoded_decay']:.1f}")
    if verdict["any_scenario_reaches_real"]:
        print("VERDICT: at least one mixture reaches the real decoded decay:")
        for name in verdict["scenarios_reaching_real"]:
            print(f"  - {name}")
    else:
        print("VERDICT: no tested mixture or prolonged-flow history reaches the "
              "real decoded decay.")
    print()
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
