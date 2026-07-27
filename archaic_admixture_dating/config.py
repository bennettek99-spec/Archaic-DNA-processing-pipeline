"""Configuration loading, validation, hashing, and profile resolution."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PACKAGE_ROOT / "configs" / "papuan_denisovan_v1.yaml"


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    config_path = Path(path) if path else DEFAULT_CONFIG
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    if overrides:
        config = deep_merge(config, overrides)
    config["_config_path"] = str(config_path.resolve())
    validate_config(config)
    return config


def validate_config(config: dict[str, Any]) -> None:
    required = ("project", "runtime", "tracts", "dating", "simulation", "models")
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError(f"Missing configuration sections: {', '.join(missing)}")
    runtime = config["runtime"]
    if runtime.get("max_block_minutes", 0) <= runtime.get("stop_before_limit_minutes", 0):
        raise ValueError("max_block_minutes must exceed stop_before_limit_minutes")
    if int(runtime.get("max_threads", 0)) < 1:
        raise ValueError("runtime.max_threads must be positive")
    min_cm = float(config["tracts"].get("minimum_length_cm", 0))
    if min_cm < 0:
        raise ValueError("tracts.minimum_length_cm cannot be negative")
    gen = float(config["project"].get("generation_time_years", 0))
    if gen <= 0:
        raise ValueError("project.generation_time_years must be positive")


def apply_profile(config: dict[str, Any], profile: str | None) -> dict[str, Any]:
    name = profile or str(config["runtime"].get("profile", "laptop"))
    profiles = config.get("profiles", {})
    if name not in profiles:
        raise ValueError(f"Unknown runtime profile {name!r}; choose from {sorted(profiles)}")
    resolved = deep_merge(config, profiles[name])
    resolved["runtime"]["profile"] = name
    validate_config(resolved)
    return resolved


def config_hash(config: dict[str, Any]) -> str:
    serializable = {key: value for key, value in config.items() if not key.startswith("_")}
    payload = json.dumps(serializable, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def save_snapshot(config: dict[str, Any], path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    serializable = {key: value for key, value in config.items() if not key.startswith("_")}
    text = yaml.safe_dump(serializable, sort_keys=False, allow_unicode=True)
    from .checkpointing import atomic_write_text

    atomic_write_text(target, text)
