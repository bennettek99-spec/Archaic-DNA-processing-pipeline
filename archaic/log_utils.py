"""
log_utils.py — shared logging setup for the orchestrator (run_pipeline.py) and
any long-running phase that wants timestamped progress output.

Phase/report scripts (phase1_validate.py, phase7_reports.py, etc.) intentionally
keep using plain print() for their tabular results and markdown/report bodies —
that IS their output, not a diagnostic log, and should stay undecorated. This
module is for progress/status messages: "starting phase N", "chunk done, ETA",
that kind of thing.

Level is controlled by the ARCHAIC_LOG_LEVEL env var (default INFO), so a run
can be made quieter/louder without touching code:
    ARCHAIC_LOG_LEVEL=DEBUG python run_pipeline.py --panel 1240k
"""
from __future__ import annotations
import logging
import os
import time


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s", "%H:%M:%S"))
        logger.addHandler(handler)
        level = os.environ.get("ARCHAIC_LOG_LEVEL", "INFO").upper()
        logger.setLevel(getattr(logging, level, logging.INFO))
        logger.propagate = False
    return logger


class StepTimer:
    """Context manager: logs start and elapsed time for a named step.

        with StepTimer(log, "Phase 3: phase3_estimate.py"):
            ...
    """
    def __init__(self, logger: logging.Logger, label: str):
        self.logger = logger
        self.label = label

    def __enter__(self):
        self.t0 = time.time()
        self.logger.info(f"start: {self.label}")
        return self

    def __exit__(self, exc_type, exc, tb):
        elapsed = time.time() - self.t0
        if exc is None:
            self.logger.info(f"done:  {self.label} ({elapsed/60:.1f} min)")
        else:
            self.logger.error(f"failed: {self.label} ({elapsed/60:.1f} min) — {exc}")
        return False
