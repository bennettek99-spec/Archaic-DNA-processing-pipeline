from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from archaic_admixture_dating.segment_structure import (
    audit_tiling,
    ascertainment_table,
    classify_source_affinity,
    effective_decay,
    fit_state_mixture,
    load_decoded_segments,
    stability_summary,
    subsampled_gof,
    threshold_table,
    validate_state_calls,
)


def _tiled_genome(individuals: int = 4, blocks: int = 60, seed: int = 7) -> pd.DataFrame:
    """Alternating modern-human/archaic segments that abut on a 1-kb grid."""
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for index in range(individuals):
        cursor = 1_000_000
        for block in range(blocks):
            archaic = block % 2 == 1
            length = int(rng.integers(40, 120) if archaic else rng.integers(500, 1500)) * 1000
            rate = 0.25 if archaic else 0.025
            rows.append(
                {
                    "name": f"P{index:02d}",
                    "chrom": "1",
                    "start": cursor,
                    "end": cursor + length + 1000,
                    "length": length,
                    "snps": int(rng.poisson(rate * length / 1000)) + int(archaic),
                    "pop": "Papuans",
                    "MeanProb": float(rng.uniform(0.5, 0.999)),
                    "Shared_with_Altai": int(rng.integers(0, 5)) if archaic else 0,
                    "Shared_with_Denisova": int(rng.integers(0, 9)) if archaic else 0,
                    "Shared_with_Vindija": int(rng.integers(0, 5)) if archaic else 0,
                    "method": "HMM",
                }
            )
            cursor += length
    return pd.DataFrame(rows)


def test_load_decoded_segments_reproduces_the_published_length_column(tmp_path):
    frame = _tiled_genome(individuals=2, blocks=10)
    source = tmp_path / "s5.tsv"
    frame.to_csv(source, sep="\t", index=False)

    loaded = load_decoded_segments(source, chromosomes=["1"])

    # The published `end` sits one 1-kb window past the segment, so the
    # half-open interval must be rebuilt from `start` and `length`.
    assert (loaded["end_bp"] - loaded["start_bp"] == loaded["length_bp"]).all()
    assert (loaded["published_end"] - loaded["start_bp"] == loaded["length_bp"] + 1000).all()


def test_load_decoded_segments_applies_no_posterior_filter(tmp_path):
    frame = _tiled_genome(individuals=2, blocks=10)
    frame.loc[0, "MeanProb"] = 0.5001
    source = tmp_path / "s5.tsv"
    frame.to_csv(source, sep="\t", index=False)

    loaded = load_decoded_segments(source, chromosomes=["1"])

    assert len(loaded) == len(frame)
    assert loaded["MeanProb"].min() == pytest.approx(0.5001)


def test_audit_tiling_detects_abutting_whole_genome_coverage(tmp_path):
    frame = _tiled_genome()
    source = tmp_path / "s5.tsv"
    frame.to_csv(source, sep="\t", index=False)
    loaded = load_decoded_segments(source, chromosomes=["1"])

    per_individual, summary = audit_tiling(loaded)

    assert len(per_individual) == 4
    assert summary["fraction_exactly_abutting"] == pytest.approx(1.0)
    assert summary["median_neighbour_gap_bp"] == 0.0
    assert summary["start_off_1kb_grid"] == 0
    assert summary["mean_prob_below_0_5"] == 0


def test_state_mixture_separates_the_two_emission_rates(tmp_path):
    frame = _tiled_genome(individuals=6, blocks=80)
    source = tmp_path / "s5.tsv"
    frame.to_csv(source, sep="\t", index=False)
    loaded = load_decoded_segments(source, chromosomes=["1"])

    labelled, parameters = fit_state_mixture(loaded)

    assert parameters["archaic_private_snps_per_kb"] > parameters["modern_human_private_snps_per_kb"]
    assert parameters["rate_ratio"] > 3.0
    # Archaic segments are the short, variant-dense ones.
    assert parameters["median_archaic_length_bp"] < parameters["median_modern_human_length_bp"]
    assert labelled["hidden_state"].isin(["archaic", "modern_human"]).all()


def test_state_calls_respect_hmm_path_alternation(tmp_path):
    frame = _tiled_genome(individuals=6, blocks=80)
    source = tmp_path / "s5.tsv"
    frame.to_csv(source, sep="\t", index=False)
    labelled, _ = fit_state_mixture(load_decoded_segments(source, chromosomes=["1"]))

    validation = validate_state_calls(labelled)

    assert validation["abutting_pairs_tested"] > 0
    assert validation["fraction_alternating"] > 0.9


def test_subsampled_gof_separates_exponential_from_contaminated_mixture():
    rng = np.random.default_rng(11)
    clean = rng.exponential(0.2, 20_000) + 0.02
    contaminated = np.concatenate([clean, rng.exponential(4.0, 6_000) + 0.02])

    clean_result = subsampled_gof(clean, 0.02, subsample=400, replicates=60, seed=3)
    dirty_result = subsampled_gof(contaminated, 0.02, subsample=400, replicates=60, seed=3)

    assert clean_result["rejection_rate_alpha_0_05"] < 0.2
    assert dirty_result["rejection_rate_alpha_0_05"] > 0.8


def test_subsampled_gof_reports_insufficient_tracts():
    assert subsampled_gof(np.linspace(0.02, 0.5, 50), 0.02, subsample=500)["status"] == (
        "insufficient_tracts"
    )


def test_stability_summary_flags_threshold_dependent_decay():
    rng = np.random.default_rng(5)
    clean = pd.DataFrame({"length_cm": rng.exponential(0.2, 40_000) + 0.02})
    contaminated = pd.DataFrame(
        {"length_cm": np.concatenate([clean["length_cm"], rng.exponential(5.0, 12_000) + 0.02])}
    )

    table = threshold_table({"clean": clean, "contaminated": contaminated})
    stability = stability_summary(table)

    spread = stability.set_index("analysis_set")["spread_ratio"]
    assert spread["clean"] < 1.1
    assert spread["contaminated"] > spread["clean"]


def test_effective_decay_recovers_a_known_exponential_rate():
    rng = np.random.default_rng(19)
    lengths = rng.exponential(100.0 / 600.0, 50_000) + 0.05

    result = effective_decay(lengths, 0.05)

    assert result["n_tracts"] == 50_000
    assert result["effective_generations"] == pytest.approx(600.0, rel=0.02)


def test_effective_decay_fails_closed_on_tiny_samples():
    result = effective_decay([0.1, 0.2, 0.3], 0.02)

    assert result["n_tracts"] == 3
    assert np.isnan(result["effective_generations"])


def test_source_affinity_rule_is_relative_sharing_only():
    frame = pd.DataFrame(
        {
            "Shared_with_Denisova": [5, 1, 5, 0],
            "Shared_with_Vindija": [1, 5, 1, 0],
            "Shared_with_Altai": [1, 1, 9, 0],
        }
    )

    result = classify_source_affinity(frame)

    assert result.tolist() == [
        "denisovan_affinity_strict",
        "neanderthal_affinity",
        "denisovan_affinity",
        "unresolved",
    ]


def test_ascertainment_table_exposes_length_dependent_labelling():
    rng = np.random.default_rng(23)
    lengths = rng.exponential(0.2, 4_000) + 0.01
    # Sharing counts scale with length, so short segments stay unlabelled.
    denisova = rng.poisson(lengths * 20)
    frame = pd.DataFrame(
        {
            "length_cm": lengths,
            "snps": rng.poisson(lengths * 60) + 1,
            "Shared_with_Denisova": denisova,
            "Shared_with_Vindija": rng.poisson(lengths * 20),
            "Shared_with_Altai": rng.poisson(lengths * 20),
        }
    )
    frame["source_class"] = classify_source_affinity(frame)

    table = ascertainment_table(frame)

    assert len(table) == 10
    assert table["fraction_classifiable"].iloc[0] < table["fraction_classifiable"].iloc[-1]
    assert table["median_length_cm"].is_monotonic_increasing
