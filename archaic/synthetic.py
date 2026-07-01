"""
synthetic.py — build a small, fully synthetic AADR-shaped panel (.ind/.snp/.geno)
for smoke-testing the file-format reader and the f4-ratio/D-statistic estimator
WITHOUT the real (large, non-redistributable) AADR data.

This is deliberately NOT a scientific-validation tool — for ground-truth accuracy
testing against a real coalescent model, see archaic/simulate.py + msprime
(used by validate_simulation.py). Here we only need internally-consistent
population allele frequencies that make the estimator behave sanely, packed
into real TGENO/GENO files, so CI and new contributors can exercise
Panel + PackedGeno + stats end-to-end on-disk in a couple of seconds.

Population model (independent sites, no LD — irrelevant for a format/plumbing
smoke test): draw an ancestral frequency per SNP, then perturb it per lineage
with lineage-appropriate divergence (small for within-clade splits, larger for
Neanderthal/Denisovan/chimp), and build the ancient "Test" population's
frequency as a two-source admixture mixture (1-alpha)*modern + alpha*archaic,
mirroring the ghost-lineage model in archaic/simulate.py but at the
allele-frequency level (no coalescent needed for a plumbing test).
"""
from __future__ import annotations
import os
import numpy as np

CHROMS = [str(c) for c in (1, 2, 3)]


def _clip(p):
    return np.clip(p, 1e-3, 1 - 1e-3)


def simulate_synthetic_frequencies(n_snp=4000, alpha_true=0.02, seed=0):
    """Return per-population allele-frequency arrays (n_snp,) for a small,
    internally-consistent synthetic panel: Altai, Vindija, Denisova, Chimp,
    Mbuti (African baseline), Test (ancient target, alpha_true Neanderthal
    ancestry), and AltaiEcho (a second draw from the Altai frequency, used as
    an "archaic-as-its-own-test" sanity check — alpha should come out ~100%)."""
    rng = np.random.default_rng(seed)
    p0 = rng.uniform(0.05, 0.95, n_snp)

    p_hum_base = _clip(p0 + rng.normal(0, 0.03, n_snp))
    p_mbuti = _clip(p_hum_base + rng.normal(0, 0.02, n_snp))
    p_test_modern = _clip(p_hum_base + rng.normal(0, 0.02, n_snp))

    p_nea_anc = _clip(p0 + rng.normal(0, 0.15, n_snp))       # Neanderthal split
    p_altai = _clip(p_nea_anc + rng.normal(0, 0.04, n_snp))   # within-Nea split
    p_vindija = _clip(p_nea_anc + rng.normal(0, 0.04, n_snp))
    p_intro = _clip(p_nea_anc + rng.normal(0, 0.04, n_snp))   # introgressing ghost lineage

    p_denisova = _clip(p0 + rng.normal(0, 0.18, n_snp))
    p_chimp = rng.uniform(0.05, 0.95, n_snp)                  # deep outgroup

    p_test = _clip((1 - alpha_true) * p_test_modern + alpha_true * p_intro)

    return dict(Altai=p_altai, Vindija=p_vindija, Denisova=p_denisova,
                Chimp=p_chimp, Mbuti=p_mbuti, Test=p_test, AltaiEcho=p_altai)


def _draw_genotypes(freqs, n_per_pop, missing_rate, seed):
    """freqs: dict name -> (n_snp,) freq array. n_per_pop: dict name -> n_ind
    (1 for single reference genomes). Returns (dosage[int8, n_ind_total, n_snp],
    ind_ids[list[str]], ind_pops[list[str]])."""
    rng = np.random.default_rng(seed + 1)
    n_snp = len(next(iter(freqs.values())))
    dosage_rows, ids, pops = [], [], []
    for name, n_ind in n_per_pop.items():
        p = freqs[name]
        for i in range(n_ind):
            g = rng.binomial(2, p).astype(np.int8)
            miss = rng.random(n_snp) < missing_rate
            g[miss] = -1
            dosage_rows.append(g)
            ids.append(f"{name}_{i}" if n_ind > 1 else name)
            pops.append(name)
    return np.stack(dosage_rows), ids, pops


def _pack_codes(codes, n):
    """codes: 1-D int array, length n, values in 0..3 -> packed bytes (ceil(n/4))."""
    codes = np.asarray(codes, dtype=np.int64)
    pad = (-n) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros(pad, dtype=np.int64)])
    codes = codes.reshape(-1, 4)
    byte = (codes[:, 0] << 6) | (codes[:, 1] << 4) | (codes[:, 2] << 2) | codes[:, 3]
    return byte.astype(np.uint8).tobytes()


def write_packed_geno(path, dosage, transposed=True):
    """Write `dosage` (int8 [n_ind, n_snp], 0/1/2, -1=missing) as either the
    TGENO (individual-major) or GENO (SNP-major) packed layout that
    archaic.lib_eigenstrat.PackedGeno reads. AADR v66 HO/1240K ship TGENO only
    (transposed=True); GENO support exists purely for testing the reader's
    other branch."""
    dosage = np.asarray(dosage)
    nind, nsnp = dosage.shape
    codes = np.where(dosage < 0, 3, dosage).astype(np.int64)
    with open(path, "wb") as fh:
        if transposed:
            header = f"TGENO {nind} {nsnp} synthetic".encode("ascii")[:48].ljust(48, b"\x00")
            fh.write(header)
            for i in range(nind):
                fh.write(_pack_codes(codes[i], nsnp))
        else:
            rlen = max(48, (nind + 3) // 4)
            header = f"GENO {nind} {nsnp} synthetic".encode("ascii")[:rlen].ljust(rlen, b"\x00")
            fh.write(header)
            for s in range(nsnp):
                row = _pack_codes(codes[:, s], nind)
                fh.write(row.ljust(rlen, b"\x00"))


def write_ind(path, ind_ids, ind_pops):
    with open(path, "w", encoding="utf-8") as fh:
        for i, p in zip(ind_ids, ind_pops):
            fh.write(f"{i}\tU\t{p}\n")


def write_snp(path, n_snp):
    with open(path, "w", encoding="utf-8") as fh:
        for i in range(n_snp):
            chrom = CHROMS[i % len(CHROMS)]
            pos = 1000 + 100 * (i // len(CHROMS))
            gpos = pos / 1e8
            fh.write(f"rs{i}\t{chrom}\t{gpos:.8f}\t{pos}\tA\tG\n")


def write_synthetic_panel(out_dir, n_snp=4000, alpha_true=0.02, n_mbuti=10,
                          n_test=15, missing_rate=0.02, seed=0, transposed=True):
    """Write a complete synthetic .ind/.snp/.geno panel to out_dir/synthetic.*
    and return (prefix, info) where info has the true alpha and population
    sizes, for use by smoke tests."""
    os.makedirs(out_dir, exist_ok=True)
    freqs = simulate_synthetic_frequencies(n_snp=n_snp, alpha_true=alpha_true, seed=seed)
    n_per_pop = {"Altai": 1, "Vindija": 1, "Denisova": 1, "Chimp": 1,
                 "Mbuti": n_mbuti, "Test": n_test, "AltaiEcho": 1}
    dosage, ind_ids, ind_pops = _draw_genotypes(freqs, n_per_pop, missing_rate, seed)

    prefix = os.path.join(out_dir, "synthetic")
    write_ind(prefix + ".ind", ind_ids, ind_pops)
    write_snp(prefix + ".snp", n_snp)
    write_packed_geno(prefix + ".geno", dosage, transposed=transposed)

    info = dict(alpha_true=alpha_true, n_snp=n_snp, n_per_pop=n_per_pop,
               ind_ids=ind_ids, ind_pops=ind_pops)
    return prefix, info
