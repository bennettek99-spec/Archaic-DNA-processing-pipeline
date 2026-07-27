"""Compact, non-interactive QC and model-fit figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_qc_figure(tracts: pd.DataFrame, model_table: pd.DataFrame, output: str | Path) -> Path:
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    lengths = tracts["length_cm"].to_numpy(dtype=float)
    figure, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].hist(lengths, bins=min(40, max(10, int(np.sqrt(len(lengths))))), color="#6c5ce7", alpha=0.85)
    axes[0].set(xlabel="Tract length (cM)", ylabel="Count", title="Observed tract lengths")
    ordered = np.sort(lengths)
    survival = 1 - np.arange(len(ordered)) / len(ordered)
    axes[1].step(ordered, survival, where="post", color="#00897b")
    axes[1].set(xlabel="Tract length (cM)", ylabel="P(length ≥ x)", title="Empirical survival")
    axes[1].set_yscale("log")
    axes[2].bar(model_table["model_name"], model_table["delta_bic"], color="#ef6c00")
    axes[2].tick_params(axis="x", rotation=30)
    axes[2].set(ylabel="ΔBIC", title="Model comparison")
    figure.tight_layout()
    figure.savefig(target, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return target
