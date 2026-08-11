from __future__ import annotations

import numpy as np
import pandas as pd

from archaic_admixture_dating.observation_calibration import (
    FEATURES,
    leave_one_out_classification,
    posterior_predictive_compatibility,
    summarize_lengths,
)


def test_summarize_lengths_recovers_exponential_rate():
    rng = np.random.default_rng(7)
    generations = 1200.0
    lengths = 0.02 + rng.exponential(1.0 / generations, 100_000) * 100.0
    result = summarize_lengths(lengths)

    assert abs(result["effective_generations_unbounded"] - generations) / generations < 0.02


def _synthetic_summaries() -> pd.DataFrame:
    rows = []
    for model_index, model_id in enumerate(["M1", "M2"]):
        for replicate in range(5):
            row = {
                "scenario": "configured_error",
                "model_id": model_id,
                "model_name": model_id,
                "replicate": replicate,
            }
            for feature_index, feature in enumerate(FEATURES):
                row[feature] = model_index * 100 + feature_index + replicate * 0.01
            rows.append(row)
    return pd.DataFrame(rows)


def test_leave_one_out_classification_separates_distinct_models():
    _, accuracy = leave_one_out_classification(_synthetic_summaries())

    assert accuracy["classification_accuracy"].eq(1.0).all()


def test_posterior_predictive_check_rejects_distant_observation():
    table = _synthetic_summaries()
    observed = {feature: 10_000.0 for feature in FEATURES}
    summary, features = posterior_predictive_compatibility(table, observed)

    assert summary["compatibility"].eq("rejected").all()
    assert not features["within_95pct_envelope"].any()
