"""
refs.py — per-panel configuration. File paths and QC thresholds come from
config.yaml (archaic.config), so no user paths are hard-coded here. This module
holds only the archaic/outgroup/baseline reference *samples*, whose IDs differ
between the HO and 1240K releases (e.g. the chimp is 'Chimp_HO.HO' in HO but
'Chimp.REF' in 1240K). PANELS merges the two so every phase stays panel-agnostic.
"""
from .config import load_config, panel_prefix

# archaic / outgroup / baseline reference samples per panel (structural, not paths)
REF_SAMPLES = {
    "ho": dict(
        Altai=dict(ids=["AltaiNeanderthal.DG"]),
        Vindija=dict(ids=["VindijaG1_final.SG"]),
        Denisova=dict(ids=["Denisova.SG"]),
        Chimp=dict(ids=["Chimp_HO.HO"]),
        Mbuti=dict(pops=["Mbuti"]),
        Yoruba=dict(pops=["Yoruba"]),
    ),
    "1240k": dict(
        Altai=dict(ids=["AltaiNeanderthal.DG"]),
        Vindija=dict(ids=["VindijaG1_final.SG"]),
        Denisova=dict(ids=["Denisova.SG"]),
        Chimp=dict(ids=["Chimp.REF"]),
        Mbuti=dict(pops=["Mbuti"]),
        Yoruba=dict(pops=["Yoruba", "YRI", "YRI-Discovery"]),
    ),
}


def _build_panels():
    cfg = load_config()
    panels = {}
    for name, pc in cfg["panels"].items():
        panels[name] = dict(
            prefix=panel_prefix(name),
            snps_col=pc["snps_col"],
            snp_floor=pc["snp_floor"],
            snp_lowpower=pc["snp_lowpower"],
            refs=REF_SAMPLES.get(name, REF_SAMPLES["1240k"]),
        )
    return panels


PANELS = _build_panels()

# group_id keyword fragments that mark a NON-test reference / archaic / non-human;
# these are removed from the analysis set regardless of geography/age.
NONHUMAN_OR_REF = (
    "neanderthal", "denisova", "chimp", "gorilla", "ancestor",
    "href", "hg19ref", "_mix",
)
