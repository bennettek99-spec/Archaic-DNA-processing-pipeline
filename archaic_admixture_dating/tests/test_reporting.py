from __future__ import annotations

import pandas as pd

from archaic_admixture_dating.reporting import _dating_estimability


def test_dating_fails_closed_when_fits_hit_lower_bounds():
    model_table = pd.DataFrame(
        [
            {
                "model_id": "two_pulse",
                "parameter_recovery_quality": "poor",
                "simulation_classification_accuracy": 0.0,
            }
        ]
    )
    fits = {
        "single_pulse": {
            "warning_flags": [
                "estimate_at_parameter_bound",
                "poor_single_exponential_fit",
            ]
        },
        "two_pulse": {"younger_generations": 50.0000001},
        "continuous_flow": {"younger_generations": 50.0000002},
    }

    estimable, reasons = _dating_estimability(model_table, fits)

    assert estimable is False
    assert any("single-pulse estimate hit" in reason for reason in reasons)
    assert any("two_pulse younger component hit" in reason for reason in reasons)
    assert any("poor simulation parameter recovery" in reason for reason in reasons)


def test_dating_can_remain_estimable_when_guardrails_pass():
    model_table = pd.DataFrame(
        [
            {
                "model_id": "single_pulse",
                "parameter_recovery_quality": "good",
                "simulation_classification_accuracy": 0.9,
            }
        ]
    )
    fits = {
        "single_pulse": {"warning_flags": []},
        "two_pulse": {"younger_generations": 800},
        "continuous_flow": {"younger_generations": 700},
    }

    estimable, reasons = _dating_estimability(model_table, fits)

    assert estimable is True
    assert reasons == []
