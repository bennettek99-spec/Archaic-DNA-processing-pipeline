"""
loci.py — locus- and gene-level archaic (Neanderthal) allele analysis.

The genome-wide pipeline answers "how much archaic ancestry"; this module answers
"WHICH archaic alleles, and did they change over time" — needed to ask whether
particular archaic-introgressed genes were under selection.

Definitions (the standard putative-introgression heuristic):
  A panel SNP is "archaic-informative" if the high-coverage archaics (Altai +
  Vindija) are ~fixed for one allele that is ~absent in Africans (Mbuti + Yoruba)
  but segregating outside Africa. The "archaic allele" is the allele the archaics
  carry. Its frequency in a cohort is our locus-level introgression signal.
  (Not every such allele is truly introgressed — incomplete lineage sorting can
  mimic it — so results are candidates, not proof.)

A curated panel of well-known adaptive-introgression loci (approximate hg19 gene
windows) is provided for targeted analysis; functions also support genome-wide
sampling for a neutral background.
"""
from __future__ import annotations
import numpy as np

# Curated adaptive-introgression / archaic-ancestry loci (approximate hg19 windows).
# Each: gene, chrom, start, end, phenotype, reference.
LOCI = [
    ("BNC2",      "9",  16_400_000, 16_880_000, "skin pigmentation",         "Vernot & Akey 2014"),
    ("OAS1-3",    "12", 113_340_000, 113_450_000, "antiviral innate immunity","Mendez et al. 2013"),
    ("TLR1-6-10", "4",  38_750_000, 38_870_000, "innate immunity (TLRs)",    "Dannemann et al. 2016"),
    ("STAT2",     "12", 56_730_000, 56_760_000, "interferon immunity",       "Mendez et al. 2012"),
    ("SLC16A11",  "17", 6_900_000,  6_980_000,  "lipid metabolism / T2D",    "SIGMA 2014"),
    ("POU2F3",    "11", 120_090_000, 120_200_000, "keratinocyte differentiation","Vernot et al. 2016"),
    ("TBX15-WARS2","1", 119_290_000, 119_630_000, "body-fat / cold response", "Racimo et al. 2017"),
    ("FADS1-2",   "11", 61_550_000, 61_660_000, "fatty-acid metabolism",     "Buckley et al. 2017"),
    ("HYAL2",     "3",  50_350_000, 50_380_000, "UV response / skin",        "Ding et al. 2014"),
    ("KRT-cluster","12",52_600_000, 53_100_000, "keratin / skin & hair",     "Vernot & Akey 2014"),
    ("OCA2-HERC2","15", 28_000_000, 28_600_000, "eye/skin pigmentation",     "archaic pigmentation"),
    ("HLA",       "6",  29_690_000, 33_500_000, "MHC / adaptive immunity",   "Abi-Rached et al. 2011"),
    ("EPAS1",     "2",  46_520_000, 46_620_000, "hypoxia (Denisovan; control)","Huerta-Sanchez 2014"),
    ("IL18RAP",   "2",  102_950_000, 103_060_000, "inflammatory immunity",   "Dannemann et al. 2016"),
]


def panel_rows_in_window(panel, chrom, start, end):
    """geno_rows of autosomal panel SNPs inside [start, end] on chrom (hg19)."""
    snp = panel.snp
    sub = snp[(snp["chrom"] == str(chrom)) & (snp["pos"] >= start) & (snp["pos"] <= end)]
    rows = sub.index.values.astype(np.int64)
    # restrict to the autosomal set the panel actually uses
    return np.intersect1d(rows, panel.snp_rows)


def archaic_informative(panel, snp_rows, refs, arch_thresh=0.9, afr_thresh=0.1):
    """For the given snp_rows, decide which are archaic-informative and which
    allele is the archaic one.

    refs: the panel's reference selector dict (refs.PANELS[..]['refs']).
    Returns dict with: rows (kept geno_rows), arch_is_a1 (bool array; True if the
    archaic allele is allele1/col5), p_arch, p_afr (allele1 freqs at kept rows).
    """
    # read Altai, Vindija, and African baseline at these SNPs
    cols_alt = panel.cols_for(**refs["Altai"])
    cols_vin = panel.cols_for(**refs["Vindija"])
    cols_afr = np.concatenate([panel.cols_for(**refs["Mbuti"]),
                               panel.cols_for(**refs["Yoruba"])])
    need = np.unique(np.concatenate([cols_alt, cols_vin, cols_afr]))
    G = panel.pg.read(snp_rows, need)
    Gf = G.astype(np.float32); Gf[G < 0] = np.nan
    pos = {c: i for i, c in enumerate(need)}

    def freq(cols):
        idx = [pos[c] for c in cols]
        with np.errstate(invalid="ignore"):
            return np.nanmean(Gf[:, idx], axis=1) / 2.0
    p_arch = np.nanmean(np.vstack([freq(cols_alt), freq(cols_vin)]), axis=0)
    p_afr = freq(cols_afr)

    arch_a1 = p_arch > 0.5                              # archaic allele = allele1?
    p_arch_allele_in_afr = np.where(arch_a1, p_afr, 1 - p_afr)
    arch_extreme = np.where(arch_a1, p_arch, 1 - p_arch)
    keep = (arch_extreme >= arch_thresh) & (p_arch_allele_in_afr <= afr_thresh) \
        & np.isfinite(p_arch) & np.isfinite(p_afr)
    return dict(rows=snp_rows[keep], arch_is_a1=arch_a1[keep],
                p_arch=p_arch[keep], p_afr=p_afr[keep])


def archaic_allele_freq(panel, snp_rows, arch_is_a1, cols):
    """Mean archaic-allele frequency per SNP for a cohort (individual cols).

    Returns (per_snp_freq, per_snp_n) where freq is the frequency of the archaic
    allele (allele1 if arch_is_a1 else allele2) among non-missing genotypes.
    """
    if len(cols) == 0 or len(snp_rows) == 0:
        return np.array([]), np.array([])
    G = panel.pg.read(snp_rows, np.asarray(cols, dtype=np.int64))
    Gf = G.astype(np.float32); Gf[G < 0] = np.nan
    with np.errstate(invalid="ignore"):
        f_a1 = np.nanmean(Gf, axis=1) / 2.0
        n = np.sum(G >= 0, axis=1)
    f_arch = np.where(arch_is_a1, f_a1, 1.0 - f_a1)
    return f_arch, n
