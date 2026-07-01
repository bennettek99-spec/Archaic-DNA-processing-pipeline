#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_pipeline.py — orchestrate the full pipeline for a panel.

Runs Phase 2 -> Phase 9 + the HTML report in order (Phase 3 is the long step).
Each phase is resumable/idempotent, so re-running skips finished work.

    python run_pipeline.py --panel 1240k
    python run_pipeline.py --panel 1240k --from 4   # resume at Phase 4
"""
import os, sys, subprocess, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.log_utils import get_logger, StepTimer

log = get_logger("archaic.run_pipeline")

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = [
    (2, "phase2_prepare.py",  ["--panel"]),
    (3, "phase3_estimate.py", ["--panel"]),
    (4, "phase4_normalize.py",["--panel"]),
    (5, "phase5_pca.py",      ["--panel"]),
    (6, "phase6_outliers.py", ["--panel"]),
    (7, "phase7_reports.py",  ["POS"]),          # positional panel arg
    (8, "phase8_figures.py",  ["POS"]),
    (9, "phase9_robustness.py",["--panel"]),
    (99, "generate_report.py",["--panel"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="1240k")
    ap.add_argument("--from", dest="start", type=int, default=2)
    args = ap.parse_args()
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    for ph, script, flags in STEPS:
        if ph < args.start and ph != 99:
            continue
        cmd = [sys.executable, os.path.join(HERE, script)]
        cmd += ([args.panel] if flags == ["POS"] else ["--panel", args.panel])
        with StepTimer(log, f"Phase {ph}: {script}"):
            r = subprocess.run(cmd, env=env)
            if r.returncode != 0:
                log.error(f"Phase {ph} failed (exit {r.returncode}); stopping.")
                sys.exit(r.returncode)
    log.info("Pipeline complete. See reports/ and results/.")


if __name__ == "__main__":
    main()
