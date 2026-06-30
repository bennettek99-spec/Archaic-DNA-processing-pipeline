"""
refs.py — per-panel configuration: file prefix, the relevant 'SNPs hit' QC column,
and the archaic/outgroup/baseline reference samples (their IDs differ between the
HO and 1240K releases — e.g. the chimp is 'Chimp_HO.HO' in HO but 'Chimp.REF' in
1240K). One place so every phase stays panel-agnostic.
"""

PANELS = {
    "ho": dict(
        prefix=r"C:\Users\benne\aadr_v66\v66.p1_HO",
        snps_col="snps_ho",
        snp_floor=15_000,      # hard exclusion below this many autosomal SNPs hit
        snp_lowpower=100_000,  # flag (not exclude) low-power samples below this
        refs=dict(
            Altai=dict(ids=["AltaiNeanderthal.DG"]),
            Vindija=dict(ids=["VindijaG1_final.SG"]),
            Denisova=dict(ids=["Denisova.SG"]),
            Chimp=dict(ids=["Chimp_HO.HO"]),
            Mbuti=dict(pops=["Mbuti"]),
            Yoruba=dict(pops=["Yoruba"]),
        ),
    ),
    "1240k": dict(
        prefix=r"C:\Users\benne\aadr_v66\v66.p1_1240K",
        snps_col="snps_1240k",
        snp_floor=30_000,
        snp_lowpower=200_000,
        refs=dict(
            Altai=dict(ids=["AltaiNeanderthal.DG"]),
            Vindija=dict(ids=["VindijaG1_final.SG"]),
            Denisova=dict(ids=["Denisova.SG"]),
            Chimp=dict(ids=["Chimp.REF"]),
            # high-coverage shotgun Africans in 1240K (HGDP/SGDP .DG)
            Mbuti=dict(pops=["Mbuti"]),
            Yoruba=dict(pops=["Yoruba", "YRI", "YRI-Discovery"]),
        ),
    ),
}

# group_id keyword fragments that mark a NON-test reference / archaic / non-human;
# these are removed from the analysis set regardless of geography/age.
NONHUMAN_OR_REF = (
    "neanderthal", "denisova", "chimp", "gorilla", "ancestor",
    "href", "hg19ref", "_mix",
)
