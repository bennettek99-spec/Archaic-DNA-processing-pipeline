from __future__ import annotations

import numpy as np
import pytest

from archaic_admixture_dating.genotype_simulation import (
    GHOST_MERGE_TIME,
    PAPUAN_MERGE_TIME,
    POP_DEN1,
    POP_GHOST,
    POP_PAPUAN,
    PUBLISHED_DEN1_TIME,
    PUBLISHED_DEN2_TIME,
    PulseConfig,
    build_demography,
    interval_lengths_morgans,
    simulate_replicate,
)

stdpopsim = pytest.importorskip("stdpopsim")
msprime = pytest.importorskip("msprime")


def _pulses(demography, dest):
    return [
        e
        for e in demography.events
        if type(e).__name__ == "MassMigration" and int(e.dest) == dest
    ]


def test_published_mode_keeps_both_denisovan_pulses():
    demography, _ = build_demography(PulseConfig(mode="published"))
    times = sorted(
        e.time
        for e in demography.events
        if type(e).__name__ == "MassMigration"
        and int(e.source) == POP_PAPUAN
        and int(e.dest) in (6, 7)
    )
    assert np.allclose(times, [PUBLISHED_DEN1_TIME, PUBLISHED_DEN2_TIME])


def test_single_mode_replaces_both_pulses_with_one():
    demography, _ = build_demography(
        PulseConfig(mode="single", generations=1400.0)
    )
    pulses = _pulses(demography, POP_DEN1)
    assert len(pulses) == 1
    assert pulses[0].time == 1400.0
    assert pulses[0].source == POP_PAPUAN
    assert not _pulses(demography, 7)


def test_single_mode_preserves_total_denisovan_proportion():
    pulse = PulseConfig(mode="single", generations=1400.0)
    demography, _ = build_demography(pulse)
    assert _pulses(demography, POP_DEN1)[0].proportion == pytest.approx(0.04)


def test_old_pulses_retarget_onto_the_papuan_ancestor():
    """Papuan does not exist beyond its merge into Ghost."""
    young, _ = build_demography(
        PulseConfig(mode="single", generations=PAPUAN_MERGE_TIME - 10)
    )
    old, _ = build_demography(
        PulseConfig(mode="single", generations=PAPUAN_MERGE_TIME + 10)
    )
    assert _pulses(young, POP_DEN1)[0].source == POP_PAPUAN
    assert _pulses(old, POP_DEN1)[0].source == POP_GHOST


def test_events_stay_time_ordered():
    demography, _ = build_demography(
        PulseConfig(mode="single", generations=1400.0)
    )
    times = [e.time for e in demography.events]
    assert times == sorted(times)


def test_pulse_time_outside_the_demography_is_rejected():
    with pytest.raises(ValueError):
        build_demography(
            PulseConfig(mode="single", generations=GHOST_MERGE_TIME + 1)
        )
    with pytest.raises(ValueError):
        build_demography(PulseConfig(mode="single", generations=-5))


def test_single_mode_requires_a_time():
    with pytest.raises(ValueError):
        build_demography(PulseConfig(mode="single", generations=None))


def test_window_counts_have_expected_shape_and_are_nonnegative():
    sequence_length = 500_000
    result = simulate_replicate(
        PulseConfig(mode="single", generations=1400.0),
        seed=4,
        sequence_length=sequence_length,
        n_papuan=3,
        n_outgroup=6,
        record_truth=False,
    )
    counts = result["counts"]
    assert counts.shape == (3, sequence_length // 1000)
    assert (counts >= 0).all()
    assert counts.sum() > 0


def test_truth_intervals_lie_inside_the_sequence():
    sequence_length = 500_000
    result = simulate_replicate(
        PulseConfig(mode="single", generations=1400.0),
        seed=5,
        sequence_length=sequence_length,
        n_papuan=3,
        n_outgroup=6,
        record_truth=True,
    )
    for intervals in result["true_archaic"].values():
        for left, right in intervals:
            assert 0 <= left < right <= sequence_length
    lengths = interval_lengths_morgans(result["true_archaic"], 1.2e-8)
    assert (lengths >= 0).all()


def test_simulation_is_deterministic_for_a_fixed_seed():
    kwargs = dict(
        seed=6, sequence_length=300_000, n_papuan=3, n_outgroup=5, record_truth=False
    )
    first = simulate_replicate(PulseConfig(mode="single", generations=1200.0), **kwargs)
    second = simulate_replicate(PulseConfig(mode="single", generations=1200.0), **kwargs)
    assert np.array_equal(first["counts"], second["counts"])
