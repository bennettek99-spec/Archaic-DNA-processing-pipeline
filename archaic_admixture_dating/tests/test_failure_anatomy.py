from __future__ import annotations

import numpy as np
import pandas as pd

from archaic_admixture_dating.failure_anatomy import (
    _analysis_sets,
    _effective_decay,
    _selected_locus_sensitivity,
    _trim_sensitivity,
    build_s4_bridge,
)


def test_effective_decay_changes_with_nonexponential_tail():
    lengths = np.concatenate([np.linspace(0.02, 0.20, 200), np.linspace(1.0, 4.0, 40)])
    low = _effective_decay(lengths, 0.02)
    high = _effective_decay(lengths, 0.5)

    assert low["n_tracts"] == 240
    assert high["n_tracts"] == 40
    assert low["effective_generations_unbounded"] != high["effective_generations_unbounded"]


def test_analysis_sets_keep_affinity_semantics_separate():
    frame = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C", "D"],
            "source_class": [
                "denisovan_affinity",
                "denisovan_affinity_strict",
                "neanderthal_affinity",
                "unresolved",
            ],
            "length_cm": [0.1, 0.2, 0.3, 0.4],
        }
    )
    sets = _analysis_sets(frame)

    assert len(sets["denisovan_broad"]) == 2
    assert len(sets["denisovan_strict"]) == 1
    assert len(sets["neanderthal_affinity"]) == 1


def test_s4_bridge_matches_all_individuals_and_is_reproducible(tmp_path, monkeypatch):
    sample_ids = [f"P{i:02d}" for i in range(89)]
    s4 = pd.DataFrame(
        {
            "name": sample_ids,
            "Dataset": ["A"] * 30 + ["B"] * 30 + ["C"] * 29,
            "Outgroup": "Whole~world",
            "Admixture_time": np.linspace(900, 1200, 89),
            "Admixture_proportion": np.linspace(2, 5, 89),
            "Coalescence_with_Human": 1500,
            "Coalescence_with_Archaic": 30000,
        }
    )
    workbook = tmp_path / "s4.xlsx"
    monkeypatch.setattr(pd, "read_excel", lambda *args, **kwargs: s4.copy())
    features = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "denisovan_broad_median_length_cm": np.linspace(0.5, 0.1, 89),
            "denisovan_broad_tract_count": np.arange(89) + 10,
        }
    )

    bridge, correlations, cohorts = build_s4_bridge(workbook, features, seed=42)

    assert len(bridge) == 89
    assert len(cohorts) == 3
    median = correlations.loc[
        correlations["feature"].eq("denisovan_broad_median_length_cm")
    ].iloc[0]
    assert median["spearman_rho"] < -0.99
    assert median["spearman_bootstrap_high"] < -0.99


def test_longest_tract_trim_reports_requested_levels():
    frame = pd.DataFrame({"length_cm": np.linspace(0.02, 5.0, 1000)})
    result = _trim_sensitivity(frame)

    assert result["longest_fraction_removed"].tolist() == [0.0, 0.001, 0.005, 0.01, 0.05]
    assert result["n_tracts"].is_monotonic_decreasing


def test_selected_locus_sensitivity_accepts_noncontiguous_index(tmp_path):
    frame = pd.DataFrame(
        {
            "chromosome": ["1"] * 12,
            "start_bp": np.arange(12) * 100,
            "end_bp": np.arange(12) * 100 + 50,
            "length_cm": np.linspace(0.02, 0.2, 12),
        },
        index=np.arange(12) * 5,
    )
    bed = tmp_path / "candidates.bed"
    bed.write_text("1\t0\t60\n", encoding="utf-8")

    result = _selected_locus_sensitivity(frame, bed).iloc[0]

    assert result["status"] == "complete"
    assert result["tracts_overlapping_selected_loci"] == 1
    assert result["tracts_retained"] == 11
