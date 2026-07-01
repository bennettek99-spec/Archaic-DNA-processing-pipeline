"""
Runs the synthetic-data plumbing smoke test (archaic/smoke.py) as part of the
normal test suite, so a broken reader/estimator is caught in CI without needing
the real AADR data. This is secondary to — and does not replace — the
AADR-based validation (phase1_validate.py, VALIDATION.md, SIMULATION_VALIDATION.md).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.smoke import run_smoke_test


def test_synthetic_smoke():
    assert run_smoke_test(seed=1, verbose=False)
