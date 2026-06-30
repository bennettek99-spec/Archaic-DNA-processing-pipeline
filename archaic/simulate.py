"""
simulate.py — coalescent simulation with archaic introgression at a KNOWN rate,
for ground-truth validation of the estimator (via msprime).

Demography (a deliberately standard human/Neanderthal model; times in
generations, 29 yr/gen):
  * chimp outgroup splits ~6.5 Mya
  * Neanderthal lineage splits from modern humans ~600 kya (small Ne ~2500)
  * Altai/Vindija split within Neanderthals ~150 kya
  * African (Mbuti-like) and out-of-Africa/Test split ~70 kya
  * a single Neanderthal -> Test admixture pulse of proportion ALPHA ~55 kya

We sample diploids from Test, African, Altai, Vindija and one Chimp, simulate
unlinked-ish blocks (independent replicates = jackknife blocks), and return
per-population derived-allele frequencies + per-Test-individual dosages so the
SAME estimator used on real data can try to recover ALPHA.
"""
from __future__ import annotations
import numpy as np
import msprime

GEN = 29.0
def _g(years):
    return years / GEN


def build_demography(test_alphas, Ne_hum=10000, Ne_nea=2500, Ne_chimp=10000):
    """test_alphas: dict {test_population_name: neanderthal_admixture_proportion}.
    All test populations are sisters of AFR (split 70 kya) and each receives its
    own Neanderthal pulse at 55 kya, so they differ only in archaic dose."""
    d = msprime.Demography()
    tests = list(test_alphas)
    # INTRO = introgressing Neanderthal ghost lineage (sister to Altai/Vindija).
    base = [("CHIMP", Ne_chimp), ("ALTAI", Ne_nea), ("VINDIJA", Ne_nea),
            ("INTRO", Ne_nea), ("NEA", Ne_nea), ("AFR", Ne_hum),
            ("HUM", Ne_hum), ("HN", Ne_hum), ("ROOT", Ne_hum)]
    for name, size in base:
        d.add_population(name=name, initial_size=size)
    for t in tests:
        d.add_population(name=t, initial_size=Ne_hum)
    for t in tests:                                   # pulse per test population
        d.add_mass_migration(time=_g(55_000), source=t, dest="INTRO",
                             proportion=test_alphas[t])
    d.add_population_split(time=_g(70_000), derived=["AFR"] + tests, ancestral="HUM")
    d.add_population_split(time=_g(150_000), derived=["ALTAI", "VINDIJA", "INTRO"],
                          ancestral="NEA")
    d.add_population_split(time=_g(600_000), derived=["HUM", "NEA"], ancestral="HN")
    d.add_population_split(time=_g(6_500_000), derived=["HN", "CHIMP"], ancestral="ROOT")
    d.sort_events()
    return d


def simulate_multi(test_alphas, n_per=20, n_afr=15, n_blocks=120, seq_len=1_000_000,
                   mu=1.25e-8, rho=1e-8, seed=0):
    """Simulate one or more test populations (sharing the SAME sites/refs).

    Returns (refs, block, dosage, popfreq):
      refs[name]    : derived-allele freq per site for Altai,Vindija,Chimp,Mbuti
      block         : block id per site (replicate index) for the jackknife
      dosage[pop]   : (n_sites, n_per) diploid dosage for each test population
      popfreq[pop]  : derived-allele freq per site for each test population
    """
    tests = list(test_alphas)
    demo = build_demography(test_alphas)
    samples = {**{t: n_per for t in tests}, "AFR": n_afr,
               "ALTAI": 1, "VINDIJA": 1, "CHIMP": 1}
    refacc = {k: [] for k in ["Altai", "Vindija", "Chimp", "Mbuti"]}
    popacc = {t: [] for t in tests}
    dosacc = {t: [] for t in tests}
    blocks = []
    cols = {}; ind_cols = {}
    refmap = {"AFR": "Mbuti", "ALTAI": "Altai", "VINDIJA": "Vindija", "CHIMP": "Chimp"}

    reps = msprime.sim_ancestry(
        samples=samples, demography=demo, sequence_length=seq_len,
        recombination_rate=rho, num_replicates=n_blocks, random_seed=seed + 1, ploidy=2)
    for bi, ts in enumerate(reps):
        mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed + 1000 + bi,
                                    model=msprime.BinaryMutationModel())
        if mts.num_sites == 0:
            continue
        G = mts.genotype_matrix()
        if not cols:
            popname = {p.id: p.metadata["name"] for p in mts.populations()}
            cols = {}; byind = {}
            for col, u in enumerate(mts.samples()):
                nd = mts.node(u); pn = popname[nd.population]
                cols.setdefault(pn, []).append(col)
                byind.setdefault((pn, nd.individual), []).append(col)
            for t in tests:
                ind_cols[t] = [v for (p, _), v in sorted(byind.items()) if p == t]
        for src, dst in refmap.items():
            refacc[dst].append(G[:, cols[src]].mean(1))
        for t in tests:
            popacc[t].append(G[:, cols[t]].mean(1))
            dosacc[t].append(np.column_stack([G[:, c].sum(1) for c in ind_cols[t]]))
        blocks.append(np.full(G.shape[0], bi, dtype=np.int32))

    refs = {k: np.concatenate(v).astype(np.float64) for k, v in refacc.items()}
    popfreq = {t: np.concatenate(v).astype(np.float64) for t, v in popacc.items()}
    dosage = {t: np.concatenate(v, axis=0) for t, v in dosacc.items()}
    block = np.concatenate(blocks)
    return refs, block, dosage, popfreq


def simulate(alpha, n_test=20, n_afr=15, n_blocks=120, seq_len=1_000_000,
             mu=1.25e-8, rho=1e-8, seed=0):
    """Single test population (back-compat). Returns (freq, block, test_dosage, info);
    freq has Altai,Vindija,Chimp,Mbuti and X (the test population)."""
    refs, block, dosage, popfreq = simulate_multi(
        {"TEST": alpha}, n_per=n_test, n_afr=n_afr, n_blocks=n_blocks,
        seq_len=seq_len, mu=mu, rho=rho, seed=seed)
    freq = {**refs, "X": popfreq["TEST"]}
    info = dict(n_sites=len(block), alpha_true=alpha, n_blocks=int(block.max()) + 1)
    return freq, block, dosage["TEST"], info
