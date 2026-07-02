"""
consumer_dna.py — read a direct-to-consumer genotype file (MyHeritage / 23andMe /
AncestryDNA / build37) and align it to an AADR panel so a single personal genome
can be run through the *exact* validated archaic-affinity estimator as one more
"test population" (like a single-genome ancient sample, e.g. Ust'-Ishim).

The only thing the f4/D statistics need is each population's frequency of ONE
consistent allele per SNP (they are polarisation-invariant, see stats.py). The
panel counts the `a1` allele from column 5 of the .snp file. So for every panel
SNP that the consumer file also genotypes, we convert the diploid call into the
dosage of that same `a1` allele (0, 1 or 2 copies) -> frequency in {0, 0.5, 1},
exactly how Panel.frequencies() treats a single-genome reference.

Alignment rules (standard genotype-harmonisation practice):
  * match consumer SNP to panel SNP by (chromosome, position) — both build37 —
    which is robust to differing SNP-name schemes (rsID vs Affx-*).
  * accept a call if its two alleles are a subset of {a1, a2}; if not, try the
    reverse-complement (opposite-strand report) before giving up.
  * DROP strand-ambiguous sites (a1/a2 == A/T or C/G): their strand cannot be
    resolved by allele identity, so a wrong-strand call would silently invert the
    dosage. With >10^5 unambiguous sites remaining this costs nothing and removes
    a whole class of error.
  * DROP no-calls ("--", "00", "..") and non-autosomal sites.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
_BASES = set("ACGT")
_AMBIG = ({"A", "T"}, {"C", "G"})
_NOCALL = {"-", "0", ".", "I", "D", "N"}   # indels/no-calls we can't use for f-stats
AUTOSOMES = {str(c) for c in range(1, 23)}


def read_consumer_file(path: str) -> pd.DataFrame:
    """Parse a MyHeritage / 23andMe / AncestryDNA raw file into a tidy frame.

    Returns columns: chrom (str), pos (int), geno (2-char upper-case string).
    Handles both the comma-quoted MyHeritage layout and the tab-separated
    23andMe/Ancestry layout, skipping '#'/'##' comment and header lines.
    """
    rows_chrom, rows_pos, rows_geno = [], [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            line = line.strip()
            if not line:
                continue
            # split on comma or whitespace; strip surrounding quotes on each field
            parts = line.split(",") if "," in line else line.split()
            parts = [p.strip().strip('"').strip("'") for p in parts]
            if len(parts) < 4:
                continue
            rsid, chrom, pos, geno = parts[0], parts[1], parts[2], parts[3]
            if chrom.upper() == "CHROMOSOME" or pos.upper() == "POSITION":
                continue  # header row
            if chrom not in AUTOSOMES:
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            g = geno.upper()
            if len(g) == 1:            # haploid call on an autosome -> duplicate
                g = g + g
            if len(g) != 2:
                continue
            rows_chrom.append(chrom)
            rows_pos.append(pos_i)
            rows_geno.append(g)
    df = pd.DataFrame({"chrom": rows_chrom, "pos": rows_pos, "geno": rows_geno})
    return df


def _dosage_of_a1(b1, b2, a1, a2):
    """Copies of allele a1 in the diploid call (b1,b2), trying a strand flip.

    Returns np.nan if the call is a no-call, ambiguous-strand, or inconsistent
    with the biallelic panel site {a1, a2}.
    """
    if b1 in _NOCALL or b2 in _NOCALL or b1 not in _BASES or b2 not in _BASES:
        return np.nan
    if {a1, a2} in _AMBIG:                     # strand-ambiguous panel site -> drop
        return np.nan
    call = {b1, b2}
    site = {a1, a2}
    if call <= site:
        return float((b1 == a1) + (b2 == a1))
    # try opposite strand
    c1, c2 = _COMP.get(b1), _COMP.get(b2)
    if c1 and c2 and {c1, c2} <= site:
        return float((c1 == a1) + (c2 == a1))
    return np.nan                              # triallelic / mismatch -> drop


def align_to_panel(panel, consumer: pd.DataFrame):
    """Build a length-`panel.n_snp` frequency array for the consumer genome.

    freq[i] is the consumer's frequency of the panel's a1 allele at autosomal
    SNP `panel.snp_rows[i]` (matching the ordering used by Panel.frequencies),
    or NaN where the consumer has no usable call. Also returns a stats dict.
    """
    snp = panel.snp.iloc[panel.snp_rows].reset_index(drop=True).copy()
    snp["i"] = np.arange(len(snp), dtype=np.int64)      # row in the freq arrays
    snp["a1"] = snp["a1"].str.upper()
    snp["a2"] = snp["a2"].str.upper()

    cons = consumer.drop_duplicates(subset=["chrom", "pos"], keep="first")
    merged = snp.merge(cons, on=["chrom", "pos"], how="inner")

    freq = np.full(panel.n_snp, np.nan, dtype=np.float64)
    n_amb = n_mismatch = n_used = n_flip = 0
    b1 = merged["geno"].str[0].to_numpy()
    b2 = merged["geno"].str[1].to_numpy()
    a1 = merged["a1"].to_numpy()
    a2 = merged["a2"].to_numpy()
    ii = merged["i"].to_numpy()

    for k in range(len(merged)):
        A1, A2 = a1[k], a2[k]
        if {A1, A2} in _AMBIG:
            n_amb += 1
            continue
        d = _dosage_of_a1(b1[k], b2[k], A1, A2)
        if np.isnan(d):
            n_mismatch += 1
            continue
        freq[ii[k]] = d / 2.0
        n_used += 1
        if b1[k] not in {A1, A2} or b2[k] not in {A1, A2}:
            n_flip += 1

    stats = dict(
        consumer_snps=len(consumer),
        panel_autosomal_snps=panel.n_snp,
        matched_by_position=len(merged),
        strand_ambiguous_dropped=n_amb,
        mismatch_dropped=n_mismatch,
        strand_flipped=n_flip,
        usable_snps=n_used,
    )
    return freq, stats


# ===========================================================================
# Archaic-allele MATCH-RATE estimator — the statistic that actually works on a
# single consumer-array genome (f4/D against the archaic *genomes* need far more
# informative sites than a consumer array shares, and collapse to ~0 even for
# populations known to be ~2% Neanderthal). Here we instead define curated marker
# alleles that are (a) DERIVED (differ from Chimp), (b) carried by BOTH high-
# coverage archaics, and (c) essentially absent in Africans, then measure how
# often the target genome carries them. At such sites the archaic-allele frequency
# is a near-direct proxy for the introgression fraction; a light internal
# calibration against reference populations of known archaic ancestry converts the
# match rate to a percentage. (Same idea as 23andMe / Vernot-&-Akey tag-SNP counts,
# but built from the in-panel archaic genomes rather than a published map.)
# ===========================================================================

def archaic_markers(freq, kind="neanderthal", afr_max=0.05, arch_thresh=0.99):
    """Boolean marker mask + 'archaic allele is a1?' mask for the panel.

    kind='neanderthal': Altai AND Vindija fixed for the same DERIVED allele
                        (opposite to Chimp), that allele African-rare.
    kind='denisovan'  : Denisova fixed for a DERIVED allele that Altai & Vindija
                        do NOT carry (Neanderthal-ancestral) and that is
                        African-rare -> Denisova-specific.
    """
    A, V, D = freq["Altai"], freq["Vindija"], freq["Denisova"]
    C, Y, M = freq["Chimp"], freq["Yoruba"], freq["Mbuti"]
    lo, hi = 1.0 - arch_thresh, arch_thresh
    if kind == "neanderthal":
        fin = np.isfinite(A) & np.isfinite(V) & np.isfinite(C) & np.isfinite(Y) & np.isfinite(M)
        a1 = fin & (A > hi) & (V > hi) & (C < lo)        # archaic-derived allele = a1
        a2 = fin & (A < lo) & (V < lo) & (C > hi)        # archaic-derived allele = a2
    elif kind == "denisovan":
        fin = (np.isfinite(D) & np.isfinite(A) & np.isfinite(V) & np.isfinite(C)
               & np.isfinite(Y) & np.isfinite(M))
        a1 = fin & (D > hi) & (C < lo) & (A < lo) & (V < lo)
        a2 = fin & (D < lo) & (C > hi) & (A > hi) & (V > hi)
    else:
        raise ValueError(kind)
    arch_a1 = a1 & ~a2
    arch_a2 = a2 & ~a1
    keep = arch_a1 | arch_a2

    def arch_freq(p):
        return np.where(arch_a1, p, np.where(arch_a2, 1.0 - p, np.nan))

    afr = 0.5 * (arch_freq(Y) + arch_freq(M))
    keep = keep & (afr < afr_max)
    return keep, arch_a1, arch_a2


def archaic_allele_freq(p, arch_a1, arch_a2):
    """Frequency of the archaic-derived allele (whichever of a1/a2 it is)."""
    return np.where(arch_a1, p, np.where(arch_a2, 1.0 - p, np.nan))


def block_mean_jackknife(vals, mask, block, n_blocks=50):
    """Delete-one-block jackknife of the mean of `vals` over `mask` sites.

    Returns (theta, se, n_used, block_sum, block_count) so callers can reuse the
    per-block sums to jackknife a downstream calibration without recomputation.
    """
    m = mask & np.isfinite(vals)
    v = np.where(m, vals, 0.0)
    w = m.astype(np.float64)
    bn = np.bincount(block, weights=v, minlength=n_blocks)
    bw = np.bincount(block, weights=w, minlength=n_blocks)
    T, W = v.sum(), w.sum()
    theta = T / W if W else np.nan
    loo = np.array([(T - bn[b]) / (W - bw[b])
                    for b in range(n_blocks) if bw[b] > 0 and (W - bw[b]) > 0])
    g = len(loo)
    se = np.sqrt((g - 1) / g * np.sum((loo - loo.mean()) ** 2)) if g > 1 else np.nan
    return theta, se, int(W), bn, bw
