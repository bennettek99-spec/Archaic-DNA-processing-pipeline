"""Run the caller-aware calibration sweep and invert it on the real data.

Usage::

    python -m archaic_admixture_dating.projects.papuan_denisovan_v1.run_calibration \
        --output archaic_admixture_dating/outputs/papuan_caller_calibration_v1

The operating point below was chosen by matching the simulation's *measurable*
observation-process anchors to the real 89-individual analysis: the two fitted
Poisson rates and the decoded archaic genome fraction. Variant density is
scaled because the real callset is filtered for callability and the simulation
is not. The Denisovan pulse time is never tuned; it is the swept parameter.

Pulse times stay at or below the Papuan/Ghost merge (1784 generations) so that
every point in the sweep delivers gene flow to the same recipient population.
Crossing that boundary changes the recipient and would put a discontinuity in
the middle of the curve.
"""

from __future__ import annotations

import argparse
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from archaic_admixture_dating.caller_calibration import (
    REAL_ARCHAIC_FRACTION,
    REAL_DECODED_DECAY,
    REAL_FITTED_PARAMETER,
    REAL_MODERN_RATE,
    REAL_ARCHAIC_RATE,
    REAL_RATE_CONTRAST,
    invert,
    save,
    sweep,
)

# The lower points are calibration anchors, not candidate histories. The first
# sweep put the real decoded decay (655.3) below every simulated point, so the
# inversion extrapolated off the bottom of the curve and returned an interval
# spanning negative time. Extending downwards turns that into interpolation and
# makes the resulting number interpretable, whether or not it is credible.
PULSE_TIMES = [
    500.0, 600.0, 700.0, 800.0,
    900.0, 1000.0, 1100.0, 1250.0, 1400.0, 1550.0, 1700.0, 1780.0,
]

OPERATING_POINT = dict(
    sequence_length=10_000_000,
    n_papuan=20,
    n_outgroup=100,
    mutation_scale=0.40,
    recombination_rate=1.2e-8,
    window_bp=1000,
)

GENERATION_TIME_YEARS = 29.0


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
    parser.add_argument("--base-seed", type=int, default=20260813)
    args = parser.parse_args(argv)

    table = sweep(
        PULSE_TIMES,
        replicates=args.replicates,
        base_seed=args.base_seed,
        **OPERATING_POINT,
    )

    inversion = invert(table, REAL_DECODED_DECAY)
    point = inversion["point_estimate_generations"]
    low, high = inversion["ci95_generations"]
    inversion["point_estimate_kya"] = point * GENERATION_TIME_YEARS / 1000.0
    inversion["ci95_kya"] = [
        low * GENERATION_TIME_YEARS / 1000.0,
        high * GENERATION_TIME_YEARS / 1000.0,
    ]
    inversion["generation_time_years"] = GENERATION_TIME_YEARS

    grouped = table.groupby("pulse_generations").mean(numeric_only=True)
    fidelity = {
        "modern_rate": [float(grouped["modern_rate"].mean()), REAL_MODERN_RATE],
        "archaic_rate": [float(grouped["archaic_rate"].mean()), REAL_ARCHAIC_RATE],
        "rate_contrast": [float(grouped["rate_contrast"].mean()), REAL_RATE_CONTRAST],
        "archaic_fraction": [
            float(grouped["archaic_fraction_decoded"].mean()),
            REAL_ARCHAIC_FRACTION,
        ],
        "decoded_over_fitted": [
            float(grouped["decoded_over_fitted"].mean()),
            REAL_DECODED_DECAY / REAL_FITTED_PARAMETER,
        ],
    }
    inversion["observation_process_fidelity"] = {
        key: {"simulated": value[0], "real": value[1], "ratio": value[0] / value[1]}
        for key, value in fidelity.items()
    }
    inversion["provenance"] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "operating_point": OPERATING_POINT,
        "pulse_times_generations": PULSE_TIMES,
        "replicates": args.replicates,
        "base_seed": args.base_seed,
    }

    save(table, inversion, args.output)

    # Full run directories are gitignored generated products; the compact
    # summary lives beside the other project result tables so the numbers in
    # the write-up have a tracked source.
    summary = grouped.reset_index()[
        [
            "pulse_generations",
            "fitted_generations",
            "decoded_decay",
            "decoded_over_fitted",
            "archaic_fraction_decoded",
            "modern_rate",
            "archaic_rate",
            "rate_contrast",
            "n_decoded_tracts",
        ]
    ].copy()
    summary.insert(1, "pulse_kya", summary["pulse_generations"] * GENERATION_TIME_YEARS / 1000.0)
    summary.to_csv(
        Path(__file__).parent / "CALLER_CALIBRATION_SUMMARY.tsv",
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    print()
    print("=" * 68)
    print("CALIBRATION CURVE  (decoded decay as a function of true pulse time)")
    print("=" * 68)
    print(f"{'true gen':>9} {'true kya':>9} {'fitted':>9} {'decoded':>9} "
          f"{'d/f':>6} {'n':>6}")
    for t, row in grouped.iterrows():
        print(f"{t:9.0f} {t*GENERATION_TIME_YEARS/1000:9.1f} "
              f"{row['fitted_generations']:9.1f} {row['decoded_decay']:9.1f} "
              f"{row['decoded_over_fitted']:6.3f} {row['n_decoded_tracts']:6.0f}")

    print()
    print("Observation-process fidelity (simulated vs real):")
    for key, value in inversion["observation_process_fidelity"].items():
        print(f"  {key:22s} {value['simulated']:9.4f} vs {value['real']:9.4f}"
              f"   ratio {value['ratio']:.2f}")

    print()
    print(f"slope {inversion['slope']:.4f}  intercept {inversion['intercept']:.1f}")
    print(f"real decoded decay: {REAL_DECODED_DECAY} generations")
    print(f"INVERTED ESTIMATE:  {point:.0f} generations "
          f"({inversion['point_estimate_kya']:.1f} kya)")
    print(f"  95% CI (simulation scatter only): {low:.0f}-{high:.0f} generations "
          f"({inversion['ci95_kya'][0]:.1f}-{inversion['ci95_kya'][1]:.1f} kya)")
    print()
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
