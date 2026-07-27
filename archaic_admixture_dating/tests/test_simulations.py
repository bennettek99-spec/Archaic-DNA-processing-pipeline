from __future__ import annotations

import pandas as pd

from archaic_admixture_dating.config import apply_profile, load_config
from archaic_admixture_dating.msprime_backend import simulate_msprime_tracts
from archaic_admixture_dating.simulations import derived_seed, simulate_tracts


def test_deterministic_seed_and_simulation():
    config = apply_profile(load_config(), "smoke")
    seed = derived_seed(config["project"]["random_seed"], "M5", 2)
    first = simulate_tracts(
        "M5",
        config["models"]["M5"],
        n_tracts=100,
        generation_time_years=29,
        minimum_length_cm=0.02,
        seed=seed,
    )
    second = simulate_tracts(
        "M5",
        config["models"]["M5"],
        n_tracts=100,
        generation_time_years=29,
        minimum_length_cm=0.02,
        seed=seed,
    )
    pd.testing.assert_frame_equal(first, second)


def test_all_ten_competing_models_generate_valid_tracts():
    config = load_config()
    assert set(config["models"]) == {f"M{number}" for number in range(1, 11)}
    for index, (model_id, model) in enumerate(config["models"].items()):
        frame = simulate_tracts(
            model_id,
            model,
            n_tracts=40,
            generation_time_years=29,
            minimum_length_cm=0.02,
            seed=100 + index,
        )
        assert len(frame) == 40
        assert (frame["length_cm"] >= 0.02).all()


def test_msprime_backend_extracts_recorded_migration_tracts():
    config = load_config()
    frame = simulate_msprime_tracts(
        "M1",
        config["models"]["M1"],
        seed=41,
        generation_time_years=29,
        sequence_length_bp=5_000_000,
        recombination_rate=1e-8,
        sample_individuals=10,
        introgression_fraction=0.5,
        minimum_length_cm=0.001,
    )
    assert not frame.empty
    assert frame["caller"].eq("msprime_migration_truth").all()
    assert (frame["length_cm"] > 0).all()
