"""
validation.py - publication-style validation summaries beyond correlation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def validation_summary(df: pd.DataFrame, *, label: str = "all") -> dict:
    """MAE, bias, RMSE, CI coverage, correlation, and Bland-Altman limits."""
    d = df.dropna(subset=["my_Nea", "pub_Nea"]).copy()
    if len(d) == 0:
        return {"set": label, "n": 0}
    err = d["my_Nea"].to_numpy(float) - d["pub_Nea"].to_numpy(float)
    avg = 0.5 * (d["my_Nea"].to_numpy(float) + d["pub_Nea"].to_numpy(float))
    se = d.get("my_SE", pd.Series(np.nan, index=d.index)).to_numpy(float)
    lo = d.get("pub_lo", d["pub_Nea"]).to_numpy(float)
    hi = d.get("pub_hi", d["pub_Nea"]).to_numpy(float)
    within_pub = (d["my_Nea"].to_numpy(float) >= lo) & (d["my_Nea"].to_numpy(float) <= hi)
    within_pub_se = (d["my_Nea"].to_numpy(float) >= lo - se) & (
        d["my_Nea"].to_numpy(float) <= hi + se
    )
    corr = np.corrcoef(d["my_Nea"], d["pub_Nea"])[0, 1] if len(d) > 1 else np.nan
    bias = float(np.mean(err))
    sd = float(np.std(err, ddof=1)) if len(err) > 1 else np.nan
    return {
        "set": label,
        "n": int(len(d)),
        "mae_pct_points": float(np.mean(np.abs(err))),
        "bias_pct_points": bias,
        "rmse_pct_points": float(np.sqrt(np.mean(err ** 2))),
        "pearson_r": float(corr),
        "ci_coverage_pub_range": float(np.mean(within_pub)),
        "ci_coverage_pub_range_plus_1se": float(np.mean(within_pub_se)),
        "bland_altman_bias": bias,
        "bland_altman_lower": bias - 1.96 * sd if np.isfinite(sd) else np.nan,
        "bland_altman_upper": bias + 1.96 * sd if np.isfinite(sd) else np.nan,
        "mean_of_methods": float(np.mean(avg)),
    }


def bland_altman_table(df: pd.DataFrame) -> pd.DataFrame:
    d = df.dropna(subset=["my_Nea", "pub_Nea"]).copy()
    d["method_mean"] = 0.5 * (d["my_Nea"] + d["pub_Nea"])
    d["method_diff"] = d["my_Nea"] - d["pub_Nea"]
    return d[["name", "category", "method_mean", "method_diff", "my_SE", "source"]]


def threshold_sensitivity(
    df: pd.DataFrame,
    *,
    thresholds=(2.0, 3.0, 5.0),
    observed_col="my_Nea",
    published_col="pub_Nea",
) -> pd.DataFrame:
    """How stable are above-threshold calls under nearby thresholds?"""
    rows = []
    d = df.dropna(subset=[observed_col, published_col])
    for thr in thresholds:
        obs = d[observed_col] >= thr
        pub = d[published_col] >= thr
        rows.append(
            {
                "threshold_pct": float(thr),
                "observed_positive": int(obs.sum()),
                "published_positive": int(pub.sum()),
                "agreement": float((obs == pub).mean()) if len(d) else np.nan,
                "false_positive": int((obs & ~pub).sum()),
                "false_negative": int((~obs & pub).sum()),
            }
        )
    return pd.DataFrame(rows)
