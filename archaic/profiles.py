"""
profiles.py — population / regional "mean-genome" allele-frequency profiles.

A single ancient genome is noisy (Phase 3-4: per-individual SE ~0.3-0.6% on
Neanderthal %). Averaging many genomes from one population/region/period into a
mean allele-frequency vector ("mean genome") beats that noise down, giving much
tighter GROUP-level archaic estimates and a clean way to compare populations.

This module:
  * builds mean allele-frequency profiles for named cohorts (memory-safe: one
    cohort read at a time);
  * computes group-level Neanderthal/Denisovan estimates with block-jackknife SEs
    (treating a cohort's mean frequency as the test population) — far more precise
    than any single genome;
  * gives pairwise genome-wide allele-frequency distances between cohorts (for a
    genetic-context / clustering view);
  * can save/load profiles (.npz) so they are reusable across analyses.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from . import stats as st


def cohort_frequencies(panel, name_to_cols, min_ind=1):
    """name_to_cols: cohort name -> array of .geno column indices.
    min_ind=1 by default so single-genome references (Altai, Chimp, ...) are kept.

    Returns (freq, info): freq[name] = float64 mean allele-1 frequency over the
    autosomes (NaN where the cohort has no data); info[name] = dict(n, n_snp).
    Reads one cohort at a time to bound memory.
    """
    freq, info = {}, {}
    for name, cols in name_to_cols.items():
        cols = np.asarray(cols, dtype=np.int64)
        if len(cols) < min_ind:
            freq[name] = np.full(panel.n_snp, np.nan)
            info[name] = dict(n=len(cols), n_snp=0)
            continue
        G = panel.pg.read(panel.snp_rows, cols)
        Gf = G.astype(np.float32); Gf[G < 0] = np.nan
        with np.errstate(invalid="ignore"):
            p = np.nanmean(Gf, axis=1) / 2.0
        freq[name] = p.astype(np.float64)
        info[name] = dict(n=int(len(cols)), n_snp=int(np.isfinite(p).sum()))
        del G, Gf
    return freq, info


def group_archaic(freq, cohorts, block, n_blocks=50,
                  A="Altai", O="Chimp", B="Mbuti", Ref="Vindija", Den="Denisova"):
    """Group-level archaic estimates per cohort, using the same validated stats.

    Requires reference frequencies (A,O,B,Ref,Den) to be present in `freq`.
    Returns a DataFrame: alpha_Nea, alpha_SE, D_Nea(+Z), D_Den(+Z), n_used.
    """
    rows = []
    for c in cohorts:
        a = st.f4_ratio(freq, A, O, c, B, Ref, block, n_blocks)
        dn = st.dstat(freq, c, B, A, O, block, n_blocks)
        dd = st.dstat(freq, c, B, Den, O, block, n_blocks)
        rows.append(dict(cohort=c, alpha_Nea=a["theta"], alpha_SE=a["se"],
                         alpha_nSNP=a["n_used"], D_Nea=dn["theta"], D_Nea_Z=dn["z"],
                         D_Den=dd["theta"], D_Den_Z=dd["z"]))
    return pd.DataFrame(rows)


def f4_contrast(freq, P1, P2, block, n_blocks=50, A="Altai", O="Chimp"):
    """D(P1, P2; Altai, Chimp): do two cohorts differ in Neanderthal sharing?
    Uses Yoruba-free Chimp rooting; for a powered differential test prefer an
    African outgroup via dstat(P1,P2,Altai,'Yoruba')."""
    return st.dstat(freq, P1, P2, A, O, block, n_blocks)


def distance_matrix(freq, names, min_overlap=20000):
    """Pairwise genome-wide allele-frequency distance: mean (pA - pB)^2 over SNPs
    where both cohorts have data. Returns an (n,n) DataFrame (0 on diagonal)."""
    F = {n: freq[n] for n in names}
    D = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            m = np.isfinite(F[a]) & np.isfinite(F[b])
            d = np.mean((F[a][m] - F[b][m]) ** 2) if m.sum() >= min_overlap else np.nan
            D[i, j] = D[j, i] = d
    return pd.DataFrame(D, index=names, columns=names)


def save_profiles(path, freq, info, snp_rows):
    names = list(freq)
    np.savez_compressed(path, names=np.array(names, dtype=object),
                        freqs=np.vstack([freq[n] for n in names]).astype(np.float32),
                        n=np.array([info[n]["n"] for n in names]),
                        n_snp=np.array([info[n]["n_snp"] for n in names]),
                        snp_rows=snp_rows)


def load_profiles(path):
    z = np.load(path, allow_pickle=True)
    names = list(z["names"])
    freq = {n: z["freqs"][i].astype(np.float64) for i, n in enumerate(names)}
    info = {n: dict(n=int(z["n"][i]), n_snp=int(z["n_snp"][i])) for i, n in enumerate(names)}
    return freq, info, z["snp_rows"]
