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


def _haplotype_fraction(result, sequence_length, key="true_denisovan_haplotype"):
    covered = sum(right - left for values in result[key].values() for left, right in values)
    n_haplotypes = 2 * len(result["individuals"])
    return covered / (n_haplotypes * sequence_length)


def test_census_sits_above_every_archaic_pulse():
    """A pulse older than the census is invisible to it."""
    from archaic_admixture_dating.genotype_simulation import POP_DEN1, POP_DEN2, POP_NEA1, _census_time

    demography, _ = build_demography(
        PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5))
    )
    census = _census_time(demography)
    archaic = {POP_DEN1, POP_DEN2, POP_NEA1}
    pulse_times = [
        e.time
        for e in demography.events
        if type(e).__name__ == "MassMigration" and int(e.dest) in archaic
    ]
    assert pulse_times
    assert census > max(pulse_times)


@pytest.mark.parametrize(
    "pulse",
    [
        PulseConfig(mode="single", generations=1400.0),
        PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5)),
        PulseConfig(mode="continuous", generations=(1550.0, 600.0), n_bins=5),
    ],
    ids=["single", "two", "continuous"],
)
def test_census_truth_recovers_the_simulated_proportion(pulse):
    """Every mode should deliver the 0.04 it was configured with.

    Averaged over seeds on purpose. At this sequence length the per-replicate
    scatter is comparable to the mean -- single seeds have come out at 0.009 and
    0.064 for the same configuration -- so a one-seed assertion would be flaky
    rather than strict. The window below still excludes the old migration-record
    method, which returned 0.158 for the single-pulse case.
    """
    sequence_length = 5_000_000
    fractions = [
        _haplotype_fraction(
            simulate_replicate(
                pulse,
                seed=seed,
                sequence_length=sequence_length,
                n_papuan=10,
                n_outgroup=15,
                record_truth=True,
            ),
            sequence_length,
        )
        for seed in (23, 24, 25)
    ]
    assert 0.018 < float(np.mean(fractions)) < 0.075


def test_diploid_coverage_exceeds_haplotype_coverage():
    """Either haplotype being archaic is enough for a diploid caller."""
    sequence_length = 5_000_000
    result = simulate_replicate(
        PulseConfig(mode="single", generations=1400.0),
        seed=22,
        sequence_length=sequence_length,
        n_papuan=10,
        n_outgroup=15,
        record_truth=True,
    )
    haplotype = _haplotype_fraction(result, sequence_length)
    covered = sum(
        right - left
        for values in result["true_denisovan"].values()
        for left, right in values
    )
    individual = covered / (len(result["individuals"]) * sequence_length)
    assert individual > haplotype
    # Two independent haplotypes give 1 - (1 - p)^2, so never more than double.
    assert individual <= 2 * haplotype + 1e-9


def test_archaic_truth_includes_neanderthal_as_well_as_denisovan():
    result = simulate_replicate(
        PulseConfig(mode="published"),
        seed=24,
        sequence_length=5_000_000,
        n_papuan=10,
        n_outgroup=15,
        record_truth=True,
    )
    denisovan = sum(
        right - left for v in result["true_denisovan"].values() for left, right in v
    )
    archaic = sum(
        right - left for v in result["true_archaic"].values() for left, right in v
    )
    assert archaic >= denisovan


def test_simulation_is_deterministic_for_a_fixed_seed():
    kwargs = dict(
        seed=6, sequence_length=300_000, n_papuan=3, n_outgroup=5, record_truth=False
    )
    first = simulate_replicate(PulseConfig(mode="single", generations=1200.0), **kwargs)
    second = simulate_replicate(PulseConfig(mode="single", generations=1200.0), **kwargs)
    assert np.array_equal(first["counts"], second["counts"])
