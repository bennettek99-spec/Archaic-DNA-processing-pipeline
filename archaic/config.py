"""
config.py - load pipeline configuration for paths and panel thresholds.

Resolution order:
  1. explicit path argument
  2. $ARCHAIC_CONFIG
  3. config.local.yaml at the repo root
  4. config.yaml at the repo root
  5. built-in panel defaults, with aadr_dir required from config or
     $ARCHAIC_AADR_DIR
"""
from __future__ import annotations

import copy
import os

import yaml

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULTS = {
    "aadr_dir": "",
    "results_dir": "results",
    "panels": {
        "ho": {
            "prefix": "v66.p1_HO",
            "snps_col": "snps_ho",
            "snp_floor": 15000,
            "snp_lowpower": 100000,
            "high_conf_min_coverage": 0.02,
            "high_conf_max_contam_lb": 0.02,
            "high_conf_max_damage": 0.30,
        },
        "1240k": {
            "prefix": "v66.p1_1240K",
            "snps_col": "snps_1240k",
            "snp_floor": 30000,
            "snp_lowpower": 200000,
            "high_conf_min_coverage": 0.02,
            "high_conf_max_contam_lb": 0.02,
            "high_conf_max_damage": 0.30,
        },
    },
}
_CACHE = None


def load_config(path=None):
    global _CACHE
    if path is None and _CACHE is not None:
        return _CACHE

    default_path = os.path.join(_REPO, "config.yaml")
    local_path = os.path.join(_REPO, "config.local.yaml")
    path = path or os.environ.get("ARCHAIC_CONFIG") or (
        local_path if os.path.exists(local_path) else default_path
    )
    cfg = copy.deepcopy(_DEFAULTS)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
        cfg.update({k: v for k, v in user.items() if k != "panels"})
        for name, panel_cfg in (user.get("panels") or {}).items():
            cfg["panels"][name] = {**cfg["panels"].get(name, {}), **panel_cfg}

    if os.environ.get("ARCHAIC_AADR_DIR"):
        cfg["aadr_dir"] = os.environ["ARCHAIC_AADR_DIR"]

    if path in {default_path, local_path}:
        _CACHE = cfg
    return cfg


def panel_prefix(panel_name, path=None):
    cfg = load_config(path)
    aadr_dir = str(cfg.get("aadr_dir") or "").strip()
    if not aadr_dir or aadr_dir in {"/path/to/aadr", "C:/path/to/aadr"}:
        raise RuntimeError(
            "AADR data directory is not configured. Set aadr_dir in config.yaml, "
            "set ARCHAIC_CONFIG to a machine-local config, or set ARCHAIC_AADR_DIR."
        )
    return os.path.join(aadr_dir, cfg["panels"][panel_name]["prefix"])
