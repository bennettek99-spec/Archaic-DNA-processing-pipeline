"""
config.py — load pipeline configuration (paths, panel thresholds) from a YAML file
so nothing is hard-coded. Resolution order:
  1. explicit path argument
  2. $ARCHAIC_CONFIG
  3. config.yaml at the repo root
  4. built-in defaults (the original developer paths) — so existing runs still work.
"""
from __future__ import annotations
import os
import yaml

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULTS = {
    "aadr_dir": r"C:\Users\benne\aadr_v66",
    "results_dir": "results",
    "panels": {
        "ho": {"prefix": "v66.p1_HO", "snps_col": "snps_ho",
               "snp_floor": 15000, "snp_lowpower": 100000},
        "1240k": {"prefix": "v66.p1_1240K", "snps_col": "snps_1240k",
                  "snp_floor": 30000, "snp_lowpower": 200000},
    },
}
_CACHE = None


def load_config(path=None):
    global _CACHE
    if path is None and _CACHE is not None:
        return _CACHE
    path = path or os.environ.get("ARCHAIC_CONFIG") or os.path.join(_REPO, "config.yaml")
    cfg = dict(_DEFAULTS)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
        cfg.update({k: v for k, v in user.items() if k != "panels"})
        if "panels" in user:
            cfg["panels"] = {**_DEFAULTS["panels"], **user["panels"]}
    if path == os.path.join(_REPO, "config.yaml"):
        _CACHE = cfg
    return cfg


def panel_prefix(panel_name, path=None):
    cfg = load_config(path)
    return os.path.join(cfg["aadr_dir"], cfg["panels"][panel_name]["prefix"])
