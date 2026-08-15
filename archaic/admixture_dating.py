"""
Estimate the time of Neanderthal admixture in a single ancient human genome.

The implementation follows the single-sample ancestry-covariance statistic in
Moorjani et al. (2016, PNAS 113:5652-5657). SNPs are ascertained where the
Altai Neanderthal carries a chimp-polarized derived allele and an African panel
is fixed for the ancestral allele. For pairs of ascertained SNPs in a genetic
distance bin, the target sample's genotype covariance is calculated and fit to

    C(d) = A * exp(-lambda * d) + c

where d is in Morgans and lambda is generations from the admixture event to the
sample. Uncertainty is estimated by leave-one-autosome-out jackknife.

This module intentionally reports dates for ancient samples. The single-sample
statistic is downward biased for the much older admixture time separating the
shared event from present-day genomes; population-sample methods are required
for that use case.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import least_squares

from .panel import AUTOSOMES, Panel


METHOD_VERSION = "0.1.0"
DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "configs" / "neanderthal_dating.yaml"


@dataclass
class ExponentialFit:
    amplitude: float
    generations: float
    offset: float
    objective: float
    r_squared: float
    n_bins: int
    weighted: bool
    converged: bool


@dataclass
class DatingEstimate:
    genetic_id: str
    label: str
    role: str
    sample_age_bp: float
    sample_age_sd: float
    n_ascertained_sites: int
    n_target_sites: int
    target_missing_fraction: float
    generations: float
    generations_se: float
    generations_ci_low: float
    generations_ci_high: float
    event_date_bp: float
    event_date_ci_low: float
    event_date_ci_high: float
    amplitude: float
    offset: float
    r_squared: float
    n_fit_bins: int
    n_jackknife_replicates: int
    fit_min_cm: float
    fit_max_cm: float
    bin_cm: float
    fit_weighted: bool
    status: str
    interpretation: str


def derived_dosage(genotype: np.ndarray, chimp_genotype: np.ndarray) -> np.ndarray:
    """Return dosage of the allele opposite the chimp state.

    AADR's packed dosage orientation is handled consistently with the existing
    Oase1 workflow: when the chimp dosage is 0, the packed dosage is treated as
    derived; when the chimp dosage is 2, 2-dosage is treated as derived.
    Missing genotypes remain -1.
    """
    genotype = np.asarray(genotype)
    chimp_genotype = np.asarray(chimp_genotype)
    out = np.full(np.broadcast_shapes(genotype.shape, chimp_genotype.shape),
                  -1, dtype=np.int8)
    g = np.broadcast_to(genotype, out.shape)
    c = np.broadcast_to(chimp_genotype, out.shape)
    called = g >= 0
    chimp0 = c == 0
    chimp2 = c == 2
    out[called & chimp0] = g[called & chimp0]
    out[called & chimp2] = 2 - g[called & chimp2]
    return out


def ascertain_neanderthal_sites(
    panel: Panel,
    refs: dict,
    *,
    min_african_calls: int = 80,
    max_african_derived_frequency: float = 0.01,
    non_african_pops: Iterable[str] = ("CEU", "CHB"),
    require_vindija: bool = False,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Ascertain Altai-derived, African-ancestral autosomal SNPs.

    This is AADR's closest reproducible analogue of "ascertainment 0" from
    Moorjani et al. Africans must have at least ``min_african_calls`` observed
    genotypes and at most ``max_african_derived_frequency`` derived frequency.
    A small nonzero ceiling is useful for the pseudo-haploid AADR representation,
    where a single error would otherwise discard a site across a large African
    panel. Set it to zero for the paper's strict fixed-ancestral definition.
    At least one derived allele must also be observed in selected non-African
    present-day controls, matching the paper's sequencing-error guard.
    """
    rows = panel.snp_rows
    altai = panel.cols_for(**refs["Altai"])
    chimp = panel.cols_for(**refs["Chimp"])
    africans = np.unique(np.concatenate([
        panel.cols_for(**refs["Mbuti"]),
        panel.cols_for(**refs["Yoruba"]),
    ]))
    vindija = panel.cols_for(**refs["Vindija"])
    if len(altai) != 1 or len(chimp) != 1:
        raise ValueError("Dating requires exactly one Altai and one chimp reference.")
    if len(africans) < min_african_calls:
        raise ValueError(
            f"Only {len(africans)} African references are available; "
            f"min_african_calls={min_african_calls}.")

    needed = np.unique(np.concatenate([altai, chimp, africans, vindija]))
    G = panel.pg.read(rows, needed)
    colpos = {int(c): i for i, c in enumerate(needed)}
    g_chimp = G[:, colpos[int(chimp[0])]]
    g_altai = G[:, colpos[int(altai[0])]]
    g_afr = G[:, [colpos[int(c)] for c in africans]]

    chimp_defined = (g_chimp == 0) | (g_chimp == 2)
    d_altai = derived_dosage(g_altai, g_chimp)
    d_afr = derived_dosage(g_afr, g_chimp[:, None])
    african_called = np.sum(g_afr >= 0, axis=1)
    african_derived_frequency = np.divide(
        np.sum(np.where(d_afr >= 0, d_afr, 0), axis=1),
        2.0 * african_called,
        out=np.full(len(rows), np.nan, dtype=float),
        where=african_called > 0,
    )
    african_fixed_ancestral = (
        (african_called >= min_african_calls)
        & (african_derived_frequency <= max_african_derived_frequency)
    )
    keep = chimp_defined & (d_altai >= 1) & african_fixed_ancestral

    if require_vindija:
        if len(vindija) != 1:
            raise ValueError("Vindija-confirmed ascertainment requires one reference.")
        g_vindija = G[:, colpos[int(vindija[0])]]
        keep &= derived_dosage(g_vindija, g_chimp) >= 1

    candidate_rows = rows[keep]
    candidate_chimp = g_chimp[keep]

    non_african_cols = panel.cols_for(pops=list(non_african_pops))
    if len(non_african_cols):
        G_nonafr = panel.pg.read(candidate_rows, non_african_cols)
        d_nonafr = derived_dosage(G_nonafr, candidate_chimp[:, None])
        observed = np.any(d_nonafr > 0, axis=1)
        candidate_rows = candidate_rows[observed]
        candidate_chimp = candidate_chimp[observed]
    else:
        observed = np.ones(len(candidate_rows), dtype=bool)

    stats = {
        "autosomal_rows": int(len(rows)),
        "chimp_defined": int(np.sum(chimp_defined)),
        "altai_derived_african_fixed": int(np.sum(keep)),
        "non_african_control_individuals": int(len(non_african_cols)),
        "final_ascertained": int(len(candidate_rows)),
        "african_reference_individuals": int(len(africans)),
        "min_african_calls": int(min_african_calls),
        "max_african_derived_frequency": float(max_african_derived_frequency),
        "require_vindija": bool(require_vindija),
        "non_african_pops": list(non_african_pops),
    }
    return candidate_rows, candidate_chimp, stats


def pair_covariance_aggregates(
    chrom: np.ndarray,
    gpos_morgans: np.ndarray,
    genotype: np.ndarray,
    *,
    min_cm: float,
    max_cm: float,
    bin_cm: float,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Accumulate per-chromosome sufficient statistics for C(d).

    Each chromosome array has shape (4, n_bins), containing sum(x), sum(y),
    sum(x*y), and number of pairs. This makes leave-one-chromosome-out
    jackknifing possible without re-enumerating SNP pairs.
    """
    if not (0 <= min_cm < max_cm):
        raise ValueError("Require 0 <= min_cm < max_cm.")
    if bin_cm <= 0:
        raise ValueError("bin_cm must be positive.")
    n_bins = int(math.ceil((max_cm - min_cm) / bin_cm))
    centers = min_cm + (np.arange(n_bins, dtype=float) + 0.5) * bin_cm
    per_chrom: dict[str, np.ndarray] = {}

    chrom = np.asarray(chrom).astype(str)
    gpos_morgans = np.asarray(gpos_morgans, dtype=float)
    genotype = np.asarray(genotype, dtype=float)
    valid = (
        np.isfinite(gpos_morgans)
        & np.isfinite(genotype)
        & (genotype >= 0)
        & np.isin(chrom, sorted(AUTOSOMES, key=int))
    )

    min_m = min_cm / 100.0
    max_m = max_cm / 100.0
    for ch in sorted(AUTOSOMES, key=int):
        use = valid & (chrom == ch)
        if np.sum(use) < 2:
            continue
        pos = gpos_morgans[use]
        val = genotype[use]
        order = np.argsort(pos, kind="mergesort")
        pos = pos[order]
        val = val[order]
        agg = np.zeros((4, n_bins), dtype=np.float64)

        for i in range(len(pos) - 1):
            lo = max(i + 1, int(np.searchsorted(pos, pos[i] + min_m, side="left")))
            hi = int(np.searchsorted(pos, pos[i] + max_m, side="right"))
            if hi <= lo:
                continue
            distances_cm = (pos[lo:hi] - pos[i]) * 100.0
            bins = np.floor((distances_cm - min_cm) / bin_cm).astype(np.int64)
            inside = (bins >= 0) & (bins < n_bins)
            if not np.any(inside):
                continue
            bins = bins[inside]
            y = val[lo:hi][inside]
            x = float(val[i])
            counts = np.bincount(bins, minlength=n_bins).astype(float)
            agg[0] += x * counts
            agg[1] += np.bincount(bins, weights=y, minlength=n_bins)
            agg[2] += np.bincount(bins, weights=x * y, minlength=n_bins)
            agg[3] += counts
        per_chrom[ch] = agg
    return centers, per_chrom


def covariance_curve(
    centers_cm: np.ndarray,
    aggregates: np.ndarray,
) -> pd.DataFrame:
    """Convert sufficient statistics to the paper's binned covariance."""
    sum_x, sum_y, sum_xy, n = np.asarray(aggregates, dtype=float)
    cov = np.full(len(n), np.nan, dtype=float)
    ok = n > 1
    cov[ok] = (
        sum_xy[ok] - (sum_x[ok] * sum_y[ok] / n[ok])
    ) / (n[ok] - 1.0)
    return pd.DataFrame({
        "distance_cm": np.asarray(centers_cm, dtype=float),
        "covariance": cov,
        "n_pairs": n.astype(np.int64),
    })


def combine_aggregates(per_chrom: dict[str, np.ndarray],
                       exclude: str | None = None) -> np.ndarray:
    kept = [v for ch, v in per_chrom.items() if ch != exclude]
    if not kept:
        raise ValueError("No chromosome aggregates remain.")
    return np.sum(np.stack(kept, axis=0), axis=0)


def fit_exponential_curve(
    curve: pd.DataFrame,
    *,
    min_pairs: int = 50,
    weighted: bool = False,
) -> ExponentialFit:
    """Fit A*exp(-lambda*d)+c using bounded multi-start least squares."""
    use = (
        np.isfinite(curve["covariance"].to_numpy())
        & (curve["n_pairs"].to_numpy() >= min_pairs)
    )
    if np.sum(use) < 12:
        raise ValueError(
            f"Only {int(np.sum(use))} covariance bins have >= {min_pairs} pairs.")
    d = curve.loc[use, "distance_cm"].to_numpy(dtype=float) / 100.0
    y = curve.loc[use, "covariance"].to_numpy(dtype=float)
    n = curve.loc[use, "n_pairs"].to_numpy(dtype=float)
    weights = np.sqrt(n / np.median(n)) if weighted else np.ones_like(n)

    tail_n = max(3, len(y) // 5)
    c0 = float(np.median(y[-tail_n:]))
    a0 = float(max(y[0] - c0, np.nanpercentile(y, 90) - c0, 1e-6))
    # Genotypes are coded 0/1/2, so their covariance cannot exceed 1 in
    # magnitude. Keeping the affine term within that physical range prevents
    # near-flat fits built from a huge positive amplitude and negative offset.
    amp_upper = 1.0
    offset_bound = 1.0

    def model(par):
        return par[0] * np.exp(-par[1] * d) + par[2]

    def residual(par):
        return (model(par) - y) * weights

    best = None
    for rate0 in (20.0, 50.0, 80.0, 120.0, 200.0, 400.0, 800.0, 1500.0):
        trial = least_squares(
            residual,
            x0=np.array([a0, rate0, c0], dtype=float),
            bounds=(
                np.array([0.0, 1.0, -offset_bound]),
                np.array([amp_upper, 5000.0, offset_bound]),
            ),
            max_nfev=20000,
        )
        objective = float(np.sum(residual(trial.x) ** 2))
        if best is None or objective < best[0]:
            best = (objective, trial)
    assert best is not None
    objective, result = best
    pred = model(result.x)
    denom = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - np.sum((y - pred) ** 2) / denom) if denom > 0 else np.nan
    return ExponentialFit(
        amplitude=float(result.x[0]),
        generations=float(result.x[1]),
        offset=float(result.x[2]),
        objective=objective,
        r_squared=r2,
        n_bins=int(np.sum(use)),
        weighted=bool(weighted),
        converged=bool(result.success),
    )


def fit_with_chromosome_jackknife(
    centers_cm: np.ndarray,
    per_chrom: dict[str, np.ndarray],
    *,
    min_pairs: int = 50,
    weighted: bool = False,
) -> tuple[ExponentialFit, float, dict[str, float], pd.DataFrame]:
    total = combine_aggregates(per_chrom)
    curve = covariance_curve(centers_cm, total)
    full = fit_exponential_curve(curve, min_pairs=min_pairs, weighted=weighted)
    leave_one_out: dict[str, float] = {}
    for ch in sorted(per_chrom, key=int):
        try:
            loo_curve = covariance_curve(
                centers_cm, combine_aggregates(per_chrom, exclude=ch))
            loo_fit = fit_exponential_curve(
                loo_curve, min_pairs=min_pairs, weighted=weighted)
            if np.isfinite(loo_fit.generations):
                leave_one_out[ch] = loo_fit.generations
        except (ValueError, RuntimeError):
            continue
    vals = np.asarray(list(leave_one_out.values()), dtype=float)
    if len(vals) >= 2:
        mean = float(np.mean(vals))
        se = float(np.sqrt((len(vals) - 1.0) / len(vals)
                           * np.sum((vals - mean) ** 2)))
    else:
        se = np.nan
    return full, se, leave_one_out, curve


def calendar_interval(
    generations: float,
    generations_se: float,
    sample_age_bp: float,
    sample_age_sd: float,
    *,
    generation_years: float = 29.0,
    generation_years_low: float = 25.0,
    generation_years_high: float = 33.0,
    seed: int = 20260723,
    draws: int = 100_000,
) -> tuple[float, float, float]:
    """Propagate rate, sample-age, and generation-interval uncertainty."""
    point = float(sample_age_bp + generations * generation_years)
    if not np.isfinite(generations_se) or generations_se <= 0:
        return point, np.nan, np.nan
    rng = np.random.default_rng(seed)
    rate = rng.normal(generations, generations_se, size=draws)
    rate = rate[rate > 0]
    if len(rate) == 0:
        return point, np.nan, np.nan
    age = rng.normal(sample_age_bp, max(sample_age_sd, 0.0), size=len(rate))
    interval = rng.uniform(generation_years_low, generation_years_high,
                           size=len(rate))
    dates = age + rate * interval
    lo, hi = np.quantile(dates, [0.025, 0.975])
    return point, float(lo), float(hi)


def _git_value(repo: Path, *args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), *args],
            text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _file_fingerprint(path: Path) -> dict:
    stat = path.stat()
    # Hash the first and last MiB: enough to identify a local run input without
    # rereading the full 7.1 GB TGENO file for every analysis.
    h = hashlib.sha256()
    with path.open("rb") as fh:
        head = fh.read(1024 * 1024)
        h.update(head)
        if stat.st_size > len(head):
            fh.seek(max(0, stat.st_size - 1024 * 1024))
            h.update(fh.read(1024 * 1024))
    return {
        "filename": path.name,
        "bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256_first_last_1mib": h.hexdigest(),
    }


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _repo_relative(path: Path, repo: Path) -> str:
    try:
        return path.resolve().relative_to(repo.resolve()).as_posix()
    except ValueError:
        return path.name


def _classification(fit: ExponentialFit, se: float, n_sites: int,
                    n_jackknife: int) -> tuple[str, str]:
    rel_se = se / fit.generations if np.isfinite(se) and fit.generations > 0 else np.inf
    if n_sites < 3_000 or n_jackknife < 15:
        return (
            "inconclusive_data_limited",
            "Too few informative target sites or successful chromosome replicates.",
        )
    if not fit.converged or not np.isfinite(fit.r_squared) or fit.r_squared < 0.02:
        return (
            "inconclusive_poor_fit",
            "The exponential ancestry-covariance curve is not a stable fit.",
        )
    if rel_se > 0.75 or fit.generations - 1.96 * se <= 0:
        return (
            "exploratory_high_uncertainty",
            "The exponential signal is present but its 95% chromosome-jackknife interval reaches zero or uncertainty is high.",
        )
    return (
        "supported_exploratory_estimate",
        "The curve, site count, and chromosome jackknife support an exploratory date.",
    )


def run_from_config(
    config_path: str | os.PathLike = DEFAULT_CONFIG,
    *,
    output_dir: str | os.PathLike = "results/admixture_dating",
    transversions_only: bool = False,
    weighted_fit: bool = False,
    sample_ids: Iterable[str] | None = None,
    method_overrides: dict | None = None,
) -> list[DatingEstimate]:
    """Run the configured AADR dating analysis and write auditable outputs."""
    config_path = Path(config_path)
    with config_path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    panel_name = str(cfg.get("panel", "1240k"))

    # Import here so config.local.yaml/environment resolution has already
    # occurred in the calling process.
    from .refs import PANELS
    if panel_name not in PANELS:
        raise ValueError(f"Unknown panel {panel_name!r}.")
    panel_cfg = PANELS[panel_name]
    panel = Panel(panel_cfg["prefix"], transversions_only=transversions_only)
    method = {**(cfg.get("method", {}) or {}), **(method_overrides or {})}
    centers, _ = np.array([]), None
    rows, chimp_calls, ascertain_stats = ascertain_neanderthal_sites(
        panel,
        panel_cfg["refs"],
        min_african_calls=int(method.get("min_african_calls", 80)),
        max_african_derived_frequency=float(
            method.get("max_african_derived_frequency", 0.01)),
        non_african_pops=method.get("non_african_pops", ["CEU", "CHB"]),
        require_vindija=bool(method.get("require_vindija", False)),
    )
    if len(rows) < 100:
        raise RuntimeError(f"Only {len(rows)} sites passed ascertainment.")

    snp = panel.snp.loc[rows]
    chrom = snp["chrom"].to_numpy(dtype=str)
    gpos = snp["gpos"].to_numpy(dtype=float)
    min_cm = float(method.get("min_cm", 0.02))
    max_cm = float(method.get("max_cm", 1.0))
    bin_cm = float(method.get("bin_cm", 0.005))
    min_pairs = int(method.get("min_pairs_per_bin", 50))
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    all_curves = []
    all_jackknife = []
    estimates: list[DatingEstimate] = []
    requested = set(sample_ids or [])
    targets = [
        target for target in cfg.get("targets", [])
        if not requested or str(target["genetic_id"]) in requested
    ]
    if requested:
        found = {str(target["genetic_id"]) for target in targets}
        missing = sorted(requested - found)
        if missing:
            raise ValueError(f"Requested samples are absent from config: {missing}")
    for index, target in enumerate(targets):
        sid = str(target["genetic_id"])
        cols = panel.cols_for(ids=[sid])
        if len(cols) != 1:
            print(f"WARNING: target {sid!r} is absent; skipping.")
            continue
        g = panel.pg.read(rows, cols)[:, 0]
        dgeno = derived_dosage(g, chimp_calls)
        called = dgeno >= 0
        try:
            centers, per_chrom = pair_covariance_aggregates(
                chrom,
                gpos,
                dgeno,
                min_cm=min_cm,
                max_cm=max_cm,
                bin_cm=bin_cm,
            )
            fit, se, loo, curve = fit_with_chromosome_jackknife(
                centers,
                per_chrom,
                min_pairs=min_pairs,
                weighted=weighted_fit,
            )
        except (ValueError, RuntimeError) as exc:
            estimates.append(DatingEstimate(
                genetic_id=sid,
                label=str(target.get("label", sid)),
                role=str(target.get("role", "replication")),
                sample_age_bp=float(target["age_bp"]),
                sample_age_sd=float(target.get("age_sd", 0.0)),
                n_ascertained_sites=int(len(rows)),
                n_target_sites=int(np.sum(called)),
                target_missing_fraction=float(1.0 - np.mean(called)),
                generations=np.nan,
                generations_se=np.nan,
                generations_ci_low=np.nan,
                generations_ci_high=np.nan,
                event_date_bp=np.nan,
                event_date_ci_low=np.nan,
                event_date_ci_high=np.nan,
                amplitude=np.nan,
                offset=np.nan,
                r_squared=np.nan,
                n_fit_bins=0,
                n_jackknife_replicates=0,
                fit_min_cm=min_cm,
                fit_max_cm=max_cm,
                bin_cm=bin_cm,
                fit_weighted=weighted_fit,
                status="inconclusive_data_limited",
                interpretation=str(exc),
            ))
            print(
                f"{sid}: no estimable curve; sites={int(np.sum(called)):,}; "
                f"reason={exc}")
            continue
        ci_low = max(0.0, fit.generations - 1.96 * se) if np.isfinite(se) else np.nan
        ci_high = fit.generations + 1.96 * se if np.isfinite(se) else np.nan
        event, event_low, event_high = calendar_interval(
            fit.generations,
            se,
            float(target["age_bp"]),
            float(target.get("age_sd", 0.0)),
            generation_years=float(method.get("generation_years", 29.0)),
            generation_years_low=float(method.get("generation_years_low", 25.0)),
            generation_years_high=float(method.get("generation_years_high", 33.0)),
            seed=int(method.get("seed", 20260723)) + index,
        )
        status, interpretation = _classification(
            fit, se, int(np.sum(called)), len(loo))
        estimate = DatingEstimate(
            genetic_id=sid,
            label=str(target.get("label", sid)),
            role=str(target.get("role", "replication")),
            sample_age_bp=float(target["age_bp"]),
            sample_age_sd=float(target.get("age_sd", 0.0)),
            n_ascertained_sites=int(len(rows)),
            n_target_sites=int(np.sum(called)),
            target_missing_fraction=float(1.0 - np.mean(called)),
            generations=fit.generations,
            generations_se=float(se),
            generations_ci_low=float(ci_low),
            generations_ci_high=float(ci_high),
            event_date_bp=event,
            event_date_ci_low=event_low,
            event_date_ci_high=event_high,
            amplitude=fit.amplitude,
            offset=fit.offset,
            r_squared=fit.r_squared,
            n_fit_bins=fit.n_bins,
            n_jackknife_replicates=len(loo),
            fit_min_cm=min_cm,
            fit_max_cm=max_cm,
            bin_cm=bin_cm,
            fit_weighted=weighted_fit,
            status=status,
            interpretation=interpretation,
        )
        estimates.append(estimate)
        c = curve.copy()
        c.insert(0, "genetic_id", sid)
        c["fitted_covariance"] = (
            fit.amplitude
            * np.exp(-fit.generations * c["distance_cm"].to_numpy() / 100.0)
            + fit.offset
        )
        all_curves.append(c)
        for ch, val in loo.items():
            all_jackknife.append({
                "genetic_id": sid,
                "excluded_chromosome": ch,
                "generations": val,
            })
        print(
            f"{sid}: {fit.generations:.1f} +/- {se:.1f} generations; "
            f"event {event:,.0f} BP; sites={int(np.sum(called)):,}; "
            f"R2={fit.r_squared:.3f}; status={status}")

    estimates_df = pd.DataFrame([asdict(x) for x in estimates])
    estimates_df.to_csv(out / "estimates.tsv", sep="\t", index=False)
    if all_curves:
        pd.concat(all_curves, ignore_index=True).to_csv(
            out / "covariance_curves.tsv", sep="\t", index=False)
    pd.DataFrame(all_jackknife).to_csv(
        out / "chromosome_jackknife.tsv", sep="\t", index=False)

    repo = Path(__file__).resolve().parent.parent
    geno_path = Path(panel.prefix + ".geno")
    code_path = Path(__file__).resolve()
    manifest = {
        "method": "single-sample Neanderthal ancestry covariance decay",
        "method_version": METHOD_VERSION,
        "citation": "Moorjani et al. 2016 PNAS 113:5652-5657",
        "config": _repo_relative(config_path, repo),
        "panel": panel_name,
        "panel_prefix": Path(panel.prefix).name,
        "transversions_only": transversions_only,
        "weighted_fit": weighted_fit,
        "ascertainment": ascertain_stats,
        "fit": {
            "min_cm": min_cm,
            "max_cm": max_cm,
            "bin_cm": bin_cm,
            "min_pairs_per_bin": min_pairs,
        },
        "git_commit": _git_value(repo, "rev-parse", "HEAD"),
        "git_branch": _git_value(repo, "branch", "--show-current"),
        "git_status_short": _git_value(repo, "status", "--short"),
        "analysis_code_sha256": _sha256_file(code_path),
        "input_geno": _file_fingerprint(geno_path),
    }
    with (out / "run_manifest.json").open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
    _plot_results(estimates, all_curves, out / "admixture_date_curves.png")
    return estimates


def _plot_results(estimates: list[DatingEstimate], curves: list[pd.DataFrame],
                  path: Path) -> None:
    if not estimates or not curves:
        return
    curve_by_id = {str(c["genetic_id"].iloc[0]): c for c in curves}
    plottable = [x for x in estimates if x.genetic_id in curve_by_id]
    if not plottable:
        return
    n = len(plottable)
    ncols = 2
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.8 * nrows),
                             squeeze=False)
    for ax, estimate in zip(axes.flat, plottable):
        c = curve_by_id[estimate.genetic_id]
        usable = c["n_pairs"] >= 50
        ax.scatter(c.loc[usable, "distance_cm"],
                   c.loc[usable, "covariance"],
                   s=7, alpha=0.45, color="#365f91", label="binned covariance")
        ax.plot(c["distance_cm"], c["fitted_covariance"],
                color="#b22222", linewidth=2, label="exponential fit")
        ax.set_title(
            f"{estimate.label}\n"
            f"{estimate.generations:.0f} +/- {estimate.generations_se:.0f} generations; "
            f"{estimate.event_date_bp / 1000:.1f} ka BP")
        ax.set_xlabel("genetic distance (cM)")
        ax.set_ylabel("genotype covariance")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle(
        "Neanderthal ancestry-covariance decay in early modern humans",
        fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="archaic-admixture-date",
        description="Estimate Neanderthal admixture time in ancient AADR genomes.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--out", default="results/admixture_dating")
    parser.add_argument("--transversions-only", action="store_true")
    parser.add_argument(
        "--sample", action="append", default=[],
        help="run only this configured genetic ID (repeatable)")
    parser.add_argument(
        "--weighted-fit", action="store_true",
        help="weight residuals by sqrt(pair count); default matches unweighted paper fit")
    parser.add_argument("--min-cm", type=float)
    parser.add_argument("--max-cm", type=float)
    parser.add_argument("--bin-cm", type=float)
    parser.add_argument("--min-pairs-per-bin", type=int)
    parser.add_argument("--max-african-derived-frequency", type=float)
    parser.add_argument("--require-vindija", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    overrides = {
        key: value for key, value in {
            "min_cm": args.min_cm,
            "max_cm": args.max_cm,
            "bin_cm": args.bin_cm,
            "min_pairs_per_bin": args.min_pairs_per_bin,
            "max_african_derived_frequency": args.max_african_derived_frequency,
            "require_vindija": True if args.require_vindija else None,
        }.items()
        if value is not None
    }
    run_from_config(
        args.config,
        output_dir=args.out,
        transversions_only=args.transversions_only,
        weighted_fit=args.weighted_fit,
        sample_ids=args.sample,
        method_overrides=overrides,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
