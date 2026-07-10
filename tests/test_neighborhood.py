import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.neighborhood import local_residual_stats


def test_effective_neighbor_count_controls_expected_se():
    df = pd.DataFrame({
        "genetic_id": ["target", "near_hi", "near_lo", "far"],
        "alpha_adj": [0.03, 0.02, 0.04, 0.05],
        "alpha_SE": [0.01, 0.001, 0.05, 0.02],
    })
    X = np.array([[0.0], [0.1], [0.2], [4.0]])
    ref = np.array([False, True, True, True])
    out = local_residual_stats(df, X, ref, df["genetic_id"].to_numpy(dtype=object), K=3)

    w = np.array([1 / 0.001 ** 2, 1 / 0.05 ** 2, 1 / 0.02 ** 2])
    a = np.array([0.02, 0.04, 0.05])
    expected = np.average(a, weights=w)
    var_obs = np.average((a - expected) ** 2, weights=w)
    n_eff = w.sum() ** 2 / (w ** 2).sum()

    assert np.isclose(out["expected"][0], expected)
    assert np.isclose(out["n_neighbors_eff"][0], n_eff)
    assert np.isclose(out["se_expected"][0], np.sqrt(var_obs / n_eff))
    assert out["n_neighbors_eff"][0] < 2.0
