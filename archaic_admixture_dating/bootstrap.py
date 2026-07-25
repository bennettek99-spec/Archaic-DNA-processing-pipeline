"""Deterministic chromosome-block and sample bootstrap utilities."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd


def bootstrap_fit(
    tracts: pd.DataFrame,
    fitter: Callable[..., dict[str, Any]],
    *,
    replicates: int,
    seed: int,
    group_column: str = "chromosome",
    fitter_kwargs: dict[str, Any] | None = None,
) -> pd.DataFrame:
    if group_column not in tracts:
        raise ValueError(f"Bootstrap group column {group_column!r} is absent")
    groups = list(tracts.groupby(group_column, sort=False))
    if len(groups) < 2:
        raise ValueError("At least two bootstrap groups are required")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    kwargs = fitter_kwargs or {}
    for replicate in range(int(replicates)):
        chosen = rng.integers(0, len(groups), size=len(groups))
        sampled = pd.concat([groups[index][1] for index in chosen], ignore_index=True)
        try:
            fit = fitter(sampled["length_cm"].to_numpy(), **kwargs)
            rows.append({"replicate": replicate, "status": "ok", **fit})
        except Exception as error:
            rows.append({"replicate": replicate, "status": "failed", "error": str(error)})
    return pd.DataFrame(rows)


def interval(frame: pd.DataFrame, column: str, alpha: float = 0.05) -> dict[str, float]:
    values = pd.to_numeric(frame.loc[frame["status"] == "ok", column], errors="coerce").dropna()
    if values.empty:
        return {"low": np.nan, "median": np.nan, "high": np.nan, "successful": 0}
    return {
        "low": float(values.quantile(alpha / 2)),
        "median": float(values.median()),
        "high": float(values.quantile(1 - alpha / 2)),
        "successful": int(len(values)),
    }
