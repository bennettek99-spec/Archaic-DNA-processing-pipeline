from __future__ import annotations

import numpy as np
import pytest

from archaic_admixture_dating.genotype_simulation import (
    POP_DEN1,
    POP_GHOST,
    POP_PAPUAN,
    PulseConfig,
    _staged_proportions,
    build_demography,
)

stdpopsim = pytest.importorskip("stdpopsim")
msprime = pytest.importorskip("msprime")

TOTAL = 0.04


def _denisovan_pulses(demography):
    return sorted(
        (
            e
            for e in demography.events
            if type(e).__name__ == "MassMigration" and int(e.dest) == POP_DEN1
        ),
        key=lambda e: e.time,
    )


def _realised_total(demography) -> float:
    """Mass migrations compose backwards in time, so proportions do not add."""
    remaining = 1.0
    for event in _denisovan_pulses(demography):
        remaining *= 1.0 - event.proportion
    return 1.0 - remaining


@pytest.mark.parametrize("weights", [(0.5, 0.5), (0.25, 0.75), (0.75, 0.25), (0.9, 0.1)])
def test_staged_proportions_compose_to_the_target_total(weights):
    p_old, p_young = _staged_proportions(TOTAL, weights)
    realised = p_young + (1.0 - p_young) * p_old
    assert realised == pytest.approx(TOTAL, rel=1e-9)


def test_staged_proportions_respect_the_requested_split():
    p_old, p_young = _staged_proportions(TOTAL, (0.25, 0.75))
    assert p_young == pytest.approx(TOTAL * 0.75)


def test_two_pulse_delivers_the_same_total_as_a_single_pulse():
    single, _ = build_demography(PulseConfig(mode="single", generations=1550.0))
    two, _ = build_demography(
        PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5))
    )
    assert _realised_total(two) == pytest.approx(_realised_total(single), rel=1e-9)


def test_continuous_flow_delivers_the_same_total():
    single, _ = build_demography(PulseConfig(mode="single", generations=1550.0))
    flow, _ = build_demography(
        PulseConfig(mode="continuous", generations=(1550.0, 600.0), n_bins=7)
    )
    assert len(_denisovan_pulses(flow)) == 7
    assert _realised_total(flow) == pytest.approx(_realised_total(single), rel=1e-9)


def test_continuous_flow_spans_the_requested_interval():
    flow, _ = build_demography(
        PulseConfig(mode="continuous", generations=(1550.0, 600.0), n_bins=5)
    )
    times = [e.time for e in _denisovan_pulses(flow)]
    assert times[0] == pytest.approx(600.0)
    assert times[-1] == pytest.approx(1550.0)
    assert np.allclose(np.diff(times), np.diff(times)[0])


def test_two_pulse_places_events_at_both_requested_times():
    two, _ = build_demography(
        PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5))
    )
    times = [e.time for e in _denisovan_pulses(two)]
    assert times == [600.0, 1550.0]


def test_components_retarget_independently_across_the_merge():
    """Only the component older than the merge moves onto the ancestor."""
    two, _ = build_demography(
        PulseConfig(mode="two", generations=(1800.0, 600.0), weights=(0.5, 0.5))
    )
    pulses = _denisovan_pulses(two)
    assert pulses[0].source == POP_PAPUAN     # 600 generations
    assert pulses[1].source == POP_GHOST      # 1800 generations


def test_reversed_interval_is_rejected():
    with pytest.raises(ValueError):
        build_demography(
            PulseConfig(mode="two", generations=(600.0, 1550.0), weights=(0.5, 0.5))
        )


def test_mixture_modes_require_a_pair_of_times():
    with pytest.raises(ValueError):
        build_demography(PulseConfig(mode="two", generations=1550.0))


def test_continuous_requires_at_least_two_bins():
    with pytest.raises(ValueError):
        build_demography(
            PulseConfig(mode="continuous", generations=(1550.0, 600.0), n_bins=1)
        )


def test_single_mode_rejects_a_pair_of_times():
    with pytest.raises(ValueError):
        build_demography(PulseConfig(mode="single", generations=(1550.0, 600.0)))


def test_unknown_mode_is_rejected():
    with pytest.raises(ValueError):
        build_demography(PulseConfig(mode="spiral", generations=1000.0))


def test_all_mixture_events_stay_time_ordered():
    for pulse in (
        PulseConfig(mode="two", generations=(1550.0, 600.0), weights=(0.5, 0.5)),
        PulseConfig(mode="continuous", generations=(1724.0, 1035.0), n_bins=7),
    ):
        demography, _ = build_demography(pulse)
        times = [e.time for e in demography.events]
        assert times == sorted(times)


def test_describe_reports_both_component_times():
    described = PulseConfig(
        mode="two", generations=(1550.0, 600.0), label="two 45+17"
    ).describe()
    assert described["pulse_generations"] == 1550.0
    assert described["pulse_generations_younger"] == 600.0
    assert described["pulse_label"] == "two 45+17"
