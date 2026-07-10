"""
refs.py - per-panel configuration.

File paths and QC thresholds come from config.yaml (archaic.config), so no user
paths are hard-coded here. This module holds only the archaic/outgroup/baseline
reference sample definitions, whose IDs differ between HO and 1240K releases.
"""
from __future__ import annotations

import os

from .config import load_config

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
    aadr_dir = str(cfg.get("aadr_dir") or "").strip()
    panels = {}
    for name, pc in cfg["panels"].items():
        panel = dict(
            prefix=os.path.join(aadr_dir, pc["prefix"]) if aadr_dir else pc["prefix"],
            snps_col=pc["snps_col"],
            snp_floor=pc["snp_floor"],
            snp_lowpower=pc["snp_lowpower"],
            refs=REF_SAMPLES.get(name, REF_SAMPLES["1240k"]),
        )
        for key in ("high_conf_min_coverage", "high_conf_max_contam_lb",
                    "high_conf_max_damage"):
            if key in pc:
                panel[key] = pc[key]
        panels[name] = panel
    return panels


PANELS = _build_panels()

NONHUMAN_OR_REF = (
    "neanderthal", "denisova", "chimp", "gorilla", "ancestor",
    "href", "hg19ref", "_mix",
)
