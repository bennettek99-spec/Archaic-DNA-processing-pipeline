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
    """How the Denisovan gene flow is configured for one replicate."""

    mode: str = "single"          # "single" | "published"
    generations: float | None = None
    proportion: float = PUBLISHED_DEN1_PROPORTION + PUBLISHED_DEN2_PROPORTION

    def describe(self) -> dict[str, Any]:
        return {
            "pulse_mode": self.mode,
            "pulse_generations": self.generations,
            "pulse_proportion": self.proportion,
        }


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

    if pulse.generations is None:
        raise ValueError("single-pulse mode requires pulse.generations")
    if not 0 < pulse.generations < GHOST_MERGE_TIME:
        raise ValueError(
            f"pulse time {pulse.generations} must be between 0 and "
            f"{GHOST_MERGE_TIME} generations for this demography"
        )

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

    recipient = POP_PAPUAN if pulse.generations <= PAPUAN_MERGE_TIME else POP_GHOST
    events.append(
        msprime.MassMigration(
            time=float(pulse.generations),
            source=recipient,
            dest=POP_DEN1,
            proportion=float(pulse.proportion),
        )
    )
    demography.events = sorted(events, key=lambda e: e.time)
    return demography, model


# --------------------------------------------------------------------------
# simulation
# --------------------------------------------------------------------------


def _true_archaic_intervals(
    ts, archaic_ids: set[int], keep_nodes: np.ndarray | None = None
) -> dict[int, list[tuple[float, float]]]:
    """Introgressed intervals per sampled individual, from migration records.

    Each migration is attributed at its own midpoint, so the tree sequence is
    walked once in order rather than rebuilt per record. Rebuilding (``ts.at``
    per migration) is quadratic in practice and makes chromosome-scale
    replicates unusable.
    """
    from collections import defaultdict

    relevant = [
        (m.left + (m.right - m.left) / 2.0, m.left, m.right, m.node)
        for m in ts.migrations()
        if int(m.dest) in archaic_ids
    ]
    if not relevant:
        return {}
    relevant.sort(key=lambda row: row[0])

    node_individual = ts.tables.nodes.individual
    keep = None if keep_nodes is None else set(int(n) for n in keep_nodes)

    intervals: dict[int, list[tuple[float, float]]] = defaultdict(list)
    cursor = 0
    n = len(relevant)
    for tree in ts.trees():
        left, right = tree.interval
        while cursor < n and relevant[cursor][0] < left:
            cursor += 1
        probe = cursor
        while probe < n and relevant[probe][0] < right:
            _, m_left, m_right, node = relevant[probe]
            try:
                descendants = tree.samples(node)
            except ValueError:
                probe += 1
                continue
            for sample_node in descendants:
                if keep is not None and int(sample_node) not in keep:
                    continue
                individual = node_individual[sample_node]
                if individual >= 0:
                    intervals[int(individual)].append((m_left, m_right))
            probe += 1
        cursor = probe
        if cursor >= n:
            break
    return {ind: _merge(vals) for ind, vals in intervals.items()}


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

    ts = msprime.sim_ancestry(
        samples={"Papuan": n_papuan, "YRI": n_outgroup},
        ploidy=2,
        demography=demography,
        sequence_length=int(sequence_length),
        recombination_rate=float(recombination_rate),
        random_seed=int(seed),
        record_migrations=record_truth,
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
        index = {ind: i for i, ind in enumerate(papuan_individuals)}
        flat_nodes = ingroup_nodes.ravel()
        denisovan = _index_intervals(
            _true_archaic_intervals(ts, den_ids, flat_nodes), index
        )
        neanderthal = _index_intervals(
            _true_archaic_intervals(ts, nea_ids, flat_nodes), index
        )
        combined: dict[int, list[tuple[float, float]]] = {}
        for source in (denisovan, neanderthal):
            for key, values in source.items():
                combined.setdefault(key, []).extend(values)
        result["true_denisovan"] = denisovan
        result["true_neanderthal"] = neanderthal
        result["true_archaic"] = {k: _merge(v) for k, v in combined.items()}
    return result


def _index_intervals(intervals: dict[int, list[tuple[float, float]]],
                     index: dict[int, int]) -> dict[int, list[tuple[float, float]]]:
    return {index[k]: v for k, v in intervals.items() if k in index}


def interval_lengths_morgans(
    intervals: dict[int, list[tuple[float, float]]], recombination_rate: float
) -> np.ndarray:
    lengths = [
        (right - left) * recombination_rate
        for values in intervals.values()
        for left, right in values
    ]
    return np.asarray(lengths, dtype=np.float64)
