"""Genotype-level simulation of the Skov observation process.

The existing :mod:`msprime_backend` reads true introgressed tracts straight out
of the tree sequence's migration records. That is the right tool for asking
whether a dating estimator recovers a rate it was given, but it cannot say
anything about what a *caller* does to those tracts, because no caller is ever
run: detection error enters as post-hoc multiplicative noise, which is
length-preserving. Posterior-decoding run inflation is length-creating, so no
setting of those knobs reproduces it.

This module closes that gap. It simulates genotypes, derives the per-window
outgroup-private variant counts that the Skov HMM actually consumes, and hands
them to :mod:`skov_hmm`. Truth and observation come out of the same replicate,
so the distortion is measured rather than assumed.

Demography is the published Jacobs et al. (2019) model as distributed by
stdpopsim (``PapuansOutOfAfrica_10J19``), with one deliberate change: the
Denisovan pulse time is a free parameter. That model encodes a *young*
Denisovan pulse (1027.6 generations, 29.8 kya) as a published finding, so
running it unmodified would assume the very thing under test. Every other
parameter stays at its published value.

The Papuan population exists only back to 1784 generations, where it merges
into Ghost. Pulses older than that are retargeted onto Ghost, which is the
Papuan ancestor at that depth. This also exposes the pulse to the CHB lineage,
which is what an older shared pulse means, but is worth remembering when
reading a sweep that crosses the boundary.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

import numpy as np

STDPOPSIM_MODEL_ID = "PapuansOutOfAfrica_10J19"

# Published pulse times in the Jacobs model, in generations.
PUBLISHED_DEN1_TIME = 1027.5862068965516
PUBLISHED_DEN2_TIME = 1575.8620689655172
PUBLISHED_DEN1_PROPORTION = 0.022
PUBLISHED_DEN2_PROPORTION = 0.018

# Population indices in the stdpopsim model.
POP_YRI = 0
POP_PAPUAN = 3
POP_DEN1 = 6
POP_DEN2 = 7
POP_NEA1 = 8
POP_GHOST = 9

PAPUAN_MERGE_TIME = 1784.0
GHOST_MERGE_TIME = 2218.0

DEFAULT_WINDOW_BP = 1000
DEFAULT_RECOMBINATION_RATE = 1.2e-8


def _load_base_model():
    import stdpopsim

    species = stdpopsim.get_species("HomSap")
    return species.get_demographic_model(STDPOPSIM_MODEL_ID)


def _is_pulse(event, time: float, source: int, dest: int) -> bool:
    return (
        type(event).__name__ == "MassMigration"
        and abs(float(event.time) - time) < 1e-6
        and int(event.source) == source
        and int(event.dest) == dest
    )


@dataclass
class PulseConfig:
    """How the Denisovan gene flow is configured for one replicate.

    ``single``      one pulse at ``generations``
    ``published``   both Jacobs pulses exactly as published
    ``two``         two pulses at ``generations = (older, younger)``, split by
                    ``weights``
    ``continuous``  gene flow spread over ``generations = (older, younger)``
                    in ``n_bins`` equal steps
    """

    mode: str = "single"
    generations: float | tuple[float, float] | None = None
    proportion: float = PUBLISHED_DEN1_PROPORTION + PUBLISHED_DEN2_PROPORTION
    weights: tuple[float, float] = (0.5, 0.5)
    n_bins: int = 7
    label: str = ""

    def describe(self) -> dict[str, Any]:
        older, younger = (None, None)
        if isinstance(self.generations, (tuple, list)):
            older, younger = float(self.generations[0]), float(self.generations[1])
        elif self.generations is not None:
            older = younger = float(self.generations)
        return {
            "pulse_mode": self.mode,
            "pulse_label": self.label or self.mode,
            "pulse_generations": older,
            "pulse_generations_younger": younger,
            "pulse_proportion": self.proportion,
            "pulse_weight_older": self.weights[0],
        }


def _interval(pulse: PulseConfig) -> tuple[float, float]:
    if not isinstance(pulse.generations, (tuple, list)) or len(pulse.generations) != 2:
        raise ValueError(f"{pulse.mode!r} mode requires generations=(older, younger)")
    older, younger = float(pulse.generations[0]), float(pulse.generations[1])
    if younger >= older:
        raise ValueError(
            f"generations must be (older, younger) with older > younger, "
            f"got ({older}, {younger})"
        )
    return older, younger


def _recipient(time: float) -> int:
    return POP_PAPUAN if time <= PAPUAN_MERGE_TIME else POP_GHOST


def _check_time(time: float) -> None:
    if not 0 < time < GHOST_MERGE_TIME:
        raise ValueError(
            f"pulse time {time} must be between 0 and {GHOST_MERGE_TIME} "
            "generations for this demography"
        )


def _staged_proportions(total: float, weights: tuple[float, float]) -> tuple[float, float]:
    """Per-event proportions giving a target total archaic fraction.

    Mass migrations compose backwards in time: a lineage that already moved at
    the younger event is no longer available at the older one. Applying
    ``p_young`` then ``p_old`` leaves a realised fraction of
    ``p_young + (1 - p_young) * p_old``, so the older event must be inflated to
    hit the intended total. Passing the raw weights straight through would
    quietly under-deliver archaic ancestry and confound the mixture sweep with
    a change in total introgression.
    """
    w_old, w_young = weights
    scale = w_old + w_young
    p_young = total * (w_young / scale)
    denominator = 1.0 - p_young
    if denominator <= 0:
        raise ValueError("younger pulse proportion leaves nothing for the older pulse")
    p_old = total * (w_old / scale) / denominator
    if not 0 <= p_old < 1:
        raise ValueError(f"derived older-pulse proportion {p_old} is out of range")
    return p_old, p_young


def build_demography(pulse: PulseConfig):
    """Published Jacobs demography with the Denisovan pulse retimed.

    ``mode="published"`` leaves both Denisovan pulses exactly as published.
    ``mode="single"`` collapses them into one pulse of the same total
    proportion at ``pulse.generations``.
    """
    model = _load_base_model()
    demography = copy.deepcopy(model.model)

    if pulse.mode == "published":
        return demography, model

    events = [
        e
        for e in demography.events
        if not _is_pulse(e, PUBLISHED_DEN1_TIME, POP_PAPUAN, POP_DEN1)
        and not _is_pulse(e, PUBLISHED_DEN2_TIME, POP_PAPUAN, POP_DEN2)
    ]
    removed = len(demography.events) - len(events)
    if removed != 2:
        raise RuntimeError(
            f"expected to remove 2 Denisovan pulses, removed {removed}; "
            "the stdpopsim model definition may have changed"
        )

    import msprime

    def add(time: float, proportion: float) -> None:
        _check_time(time)
        events.append(
            msprime.MassMigration(
                time=float(time),
                source=_recipient(time),
                dest=POP_DEN1,
                proportion=float(proportion),
            )
        )

    if pulse.mode == "single":
        if pulse.generations is None:
            raise ValueError("single-pulse mode requires pulse.generations")
        if isinstance(pulse.generations, (tuple, list)):
            raise ValueError("single-pulse mode takes a scalar generations value")
        add(float(pulse.generations), pulse.proportion)

    elif pulse.mode == "two":
        older, younger = _interval(pulse)
        p_old, p_young = _staged_proportions(pulse.proportion, pulse.weights)
        add(younger, p_young)
        add(older, p_old)

    elif pulse.mode == "continuous":
        older, younger = _interval(pulse)
        if pulse.n_bins < 2:
            raise ValueError("continuous mode requires n_bins >= 2")
        # Each step removes the same fraction of what remains, so the realised
        # total is 1 - (1 - q)**n_bins. Solving for q keeps the total archaic
        # fraction identical to the single-pulse runs.
        q = 1.0 - (1.0 - pulse.proportion) ** (1.0 / pulse.n_bins)
        for time in np.linspace(younger, older, pulse.n_bins):
            add(float(time), q)

    else:
        raise ValueError(f"unsupported pulse mode {pulse.mode!r}")

    demography.events = sorted(events, key=lambda e: e.time)
    return demography, model


# --------------------------------------------------------------------------
# simulation
# --------------------------------------------------------------------------


def _census_intervals(
    ts, archaic_ids: set[int], sample_nodes: np.ndarray, census_time: float
) -> dict[int, list[tuple[float, float]]]:
    """Introgressed intervals per sample haplotype, from a census.

    Attributing migration *records* to descendants -- the obvious approach, and
    the one :mod:`msprime_backend` still uses -- is wrong for this purpose. A
    migration record's interval is the span of the ancestral lineage at the
    time it moved, which is far wider than the segment any one modern sample
    ends up inheriting from it, so the intervals over-attribute. The error is
    not subtle: measured archaic fraction moved from 0.158 at 10 Mb to 0.093 at
    30 Mb, and the recovered decay for a 1400-generation pulse came out at 921
    and 1180 respectively. A quantity that depends on how much sequence you
    simulated is not a truth.

    A census node sits on every lineage at a chosen time. Placing one older
    than every archaic pulse means introgressed lineages have already entered
    their source population, so ``link_ancestors`` returns exactly the segments
    each sample inherits from an archaic ancestor. No approximation.
    """
    from collections import defaultdict

    import msprime

    nodes = ts.tables.nodes
    is_census = (nodes.flags & msprime.NODE_IS_CEN_EVENT) != 0
    census_nodes = np.flatnonzero(is_census & (nodes.time == census_time))
    if census_nodes.size == 0:
        return {}
    archaic_census = census_nodes[
        np.isin(nodes.population[census_nodes], np.fromiter(archaic_ids, dtype=np.int32))
    ]
    if archaic_census.size == 0:
        return {}

    edges = ts.link_ancestors(
        samples=[int(n) for n in sample_nodes],
        ancestors=[int(n) for n in archaic_census],
    )
    intervals: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for left, right, child in zip(edges.left, edges.right, edges.child):
        intervals[int(child)].append((float(left), float(right)))
    return {node: _merge(values) for node, values in intervals.items()}


def _merge_to_individuals(
    per_node: dict[int, list[tuple[float, float]]],
    ingroup_nodes: np.ndarray,
) -> dict[int, list[tuple[float, float]]]:
    """Collapse haplotype intervals onto their diploid individual.

    Per-haplotype coverage is the quantity that should equal the simulated
    admixture proportion. Per-individual coverage is what a diploid caller
    sees, and is larger because either haplotype being archaic is enough.
    """
    merged: dict[int, list[tuple[float, float]]] = {}
    for index in range(ingroup_nodes.shape[0]):
        combined: list[tuple[float, float]] = []
        for node in ingroup_nodes[index]:
            combined.extend(per_node.get(int(node), []))
        if combined:
            merged[index] = _merge(combined)
    return merged


def _merge(intervals: list[tuple[float, float]]) -> list[tuple[float, float]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for left, right in ordered[1:]:
        prev_left, prev_right = merged[-1]
        if left <= prev_right:
            merged[-1] = (prev_left, max(prev_right, right))
        else:
            merged.append((left, right))
    return merged


def _window_counts(ts, ingroup_nodes: np.ndarray, outgroup_nodes: np.ndarray,
                   sequence_length: int, window_bp: int) -> np.ndarray:
    """Per-individual counts of outgroup-private derived variants per window.

    This is the Skov HMM's observation: a variant counts for an individual if
    the individual carries it on either haplotype and no outgroup sample
    carries it at all.
    """
    genotypes = ts.genotype_matrix() > 0          # sites x haplotypes
    positions = ts.tables.sites.position

    private = ~genotypes[:, outgroup_nodes].any(axis=1)
    if not private.any():
        n_ind = ingroup_nodes.shape[0]
        return np.zeros((n_ind, sequence_length // window_bp), dtype=np.int32)

    rows = np.flatnonzero(private)
    # ingroup_nodes is (n_individuals, 2)
    carried = genotypes[np.ix_(rows, ingroup_nodes.ravel())]
    carried = carried.reshape(rows.size, ingroup_nodes.shape[0], 2).any(axis=2)

    window_index = (positions[rows] // window_bp).astype(np.int64)
    n_windows = sequence_length // window_bp
    counts = np.zeros((ingroup_nodes.shape[0], n_windows), dtype=np.int32)
    for i in range(ingroup_nodes.shape[0]):
        hit = window_index[carried[:, i]]
        if hit.size:
            counts[i] = np.bincount(hit, minlength=n_windows)[:n_windows]
    return counts


def simulate_replicate(
    pulse: PulseConfig,
    *,
    seed: int,
    sequence_length: int = 20_000_000,
    n_papuan: int = 30,
    n_outgroup: int = 50,
    recombination_rate: float = DEFAULT_RECOMBINATION_RATE,
    mutation_rate: float | None = None,
    window_bp: int = DEFAULT_WINDOW_BP,
    record_truth: bool = True,
) -> dict[str, Any]:
    """One simulated chunk: window counts plus true introgressed intervals."""
    import msprime

    demography, model = build_demography(pulse)
    if mutation_rate is None:
        mutation_rate = model.mutation_rate

    census_time = None
    if record_truth:
        census_time = _census_time(demography)
        demography.events = sorted(
            list(demography.events) + [msprime.CensusEvent(time=census_time)],
            key=lambda e: e.time,
        )

    ts = msprime.sim_ancestry(
        samples={"Papuan": n_papuan, "YRI": n_outgroup},
        ploidy=2,
        demography=demography,
        sequence_length=int(sequence_length),
        recombination_rate=float(recombination_rate),
        random_seed=int(seed),
    )
    ts = msprime.sim_mutations(ts, rate=float(mutation_rate), random_seed=int(seed) + 1)

    populations = {p.metadata["name"]: p.id for p in ts.populations()}
    papuan_individuals = [
        ind.id
        for ind in ts.individuals()
        if ts.node(ind.nodes[0]).population == populations["Papuan"]
    ]
    ingroup_nodes = np.array(
        [ts.individual(i).nodes for i in papuan_individuals], dtype=np.int64
    )
    outgroup_nodes = np.array(
        ts.samples(population=populations["YRI"]), dtype=np.int64
    )

    counts = _window_counts(
        ts, ingroup_nodes, outgroup_nodes, int(sequence_length), window_bp
    )

    result: dict[str, Any] = {
        "counts": counts,
        "individuals": papuan_individuals,
        "sequence_length": int(sequence_length),
        "window_bp": window_bp,
        "recombination_rate": recombination_rate,
        "mutation_rate": mutation_rate,
        "n_sites": int(ts.num_sites),
        "seed": int(seed),
        **pulse.describe(),
    }

    if record_truth:
        den_ids = {populations["Den1"], populations["Den2"]}
        nea_ids = {populations["Nea1"]}
        flat_nodes = ingroup_nodes.ravel()

        den_by_node = _census_intervals(ts, den_ids, flat_nodes, census_time)
        nea_by_node = _census_intervals(ts, nea_ids, flat_nodes, census_time)
        archaic_by_node: dict[int, list[tuple[float, float]]] = {}
        for source in (den_by_node, nea_by_node):
            for node, values in source.items():
                archaic_by_node.setdefault(node, []).extend(values)
        archaic_by_node = {k: _merge(v) for k, v in archaic_by_node.items()}

        result["census_time"] = census_time
        # Per haplotype: should equal the simulated admixture proportion.
        result["true_denisovan_haplotype"] = den_by_node
        result["true_archaic_haplotype"] = archaic_by_node
        # Per individual: what a diploid caller can see.
        result["true_denisovan"] = _merge_to_individuals(den_by_node, ingroup_nodes)
        result["true_neanderthal"] = _merge_to_individuals(nea_by_node, ingroup_nodes)
        result["true_archaic"] = _merge_to_individuals(archaic_by_node, ingroup_nodes)
    return result


def _census_time(demography) -> float:
    """A time older than every archaic pulse in the model.

    The census has to sit above all of them, otherwise pulses older than it are
    invisible: their lineages have not yet entered the archaic population when
    the census is taken.
    """
    archaic = {POP_DEN1, POP_DEN2, POP_NEA1}
    times = [
        float(e.time)
        for e in demography.events
        if type(e).__name__ == "MassMigration" and int(e.dest) in archaic
    ]
    if not times:
        raise ValueError("no archaic pulses found; cannot place a census")
    return max(times) + 1.0


def interval_lengths_morgans(
    intervals: dict[int, list[tuple[float, float]]], recombination_rate: float
) -> np.ndarray:
    lengths = [
        (right - left) * recombination_rate
        for values in intervals.values()
        for left, right in values
    ]
    return np.asarray(lengths, dtype=np.float64)
