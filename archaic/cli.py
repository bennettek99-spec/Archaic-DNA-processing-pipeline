"""
cli.py — `archaic-pipeline` console-script entry point.

Thin wrapper around the existing top-level phase scripts (phase1_validate.py,
phase2_prepare.py, ... run_pipeline.py) so the pipeline can be invoked as an
installed command instead of `python phase3_estimate.py`:

    archaic-pipeline validate --panel 1240k
    archaic-pipeline all --panel 1240k
    archaic-pipeline smoke-test

Each subcommand just execs the corresponding script with the remaining
arguments passed through unchanged, so behaviour is identical to calling the
script directly — this file adds no new logic beyond dispatch. This assumes an
editable/in-place install (`pip install -e .`): the phase scripts live at the
repo root next to the `archaic/` package, which is where this module resolves
them from.
"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# subcommand -> script at the repo root
SCRIPTS = {
    "validate": "phase1_validate.py",
    "prepare": "phase2_prepare.py",
    "estimate": "phase3_estimate.py",
    "normalize": "phase4_normalize.py",
    "pca": "phase5_pca.py",
    "outliers": "phase6_outliers.py",
    "reports": "phase7_reports.py",
    "figures": "phase8_figures.py",
    "robustness": "phase9_robustness.py",
    "report": "generate_report.py",
    "all": "run_pipeline.py",
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="archaic-pipeline",
        description="Dispatch to the archaic-introgression pipeline's phase scripts.")
    sub = ap.add_subparsers(dest="command", required=True)
    for name in SCRIPTS:
        sub.add_parser(name, help=f"run {SCRIPTS[name]}",
                       add_help=False)
    sub.add_parser("smoke-test", help="synthetic-data plumbing check (no AADR data needed)")

    args, rest = ap.parse_known_args(argv)

    if args.command == "smoke-test":
        from .smoke import run_smoke_test
        sys.exit(0 if run_smoke_test() else 1)

    script = os.path.join(REPO_ROOT, SCRIPTS[args.command])
    if not os.path.exists(script):
        print(f"archaic-pipeline: {script} not found (expected an editable "
              f"install alongside the repo's phase scripts).", file=sys.stderr)
        sys.exit(1)
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    r = subprocess.run([sys.executable, script] + rest, env=env)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
