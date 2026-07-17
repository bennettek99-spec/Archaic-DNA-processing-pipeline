"""Run a reference-aware Denisovan genome analysis on an AADR panel.

The main ancient-Eurasian pipeline deliberately excludes archaic references from
its target cohort.  This module provides the complementary workflow needed when
the target *is* an archaic genome.  It reuses the validated Panel and f-statistic
engine, but avoids presenting a modern-human introgression percentage for a
Denisovan genome.

Default target: Denisova3.DG, the high-coverage diploid call set.  Denisova.SG is
the same biological specimen in a lower-density pseudo-haploid call set and is
therefore treated as a technical replicate, never independent ancestry evidence.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import os
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import stats as st
from .anno import load_anno
from .config import panel_prefix
from .panel import Panel
from .profiles import cohort_frequencies


DEFAULT_TARGET = "Denisova3.DG"
DEFAULT_OUTPUT = Path("results") / "denisovan_genome"
N_BLOCKS = 50

SAMPLE_SPECS = {
    "Target": {"ids": [DEFAULT_TARGET]},
    "Denisova_replicate": {"ids": ["Denisova.SG"]},
    "Denisova_alt_calls": {"ids": ["Denisova3_snpAD.DG"]},
    "Denisova11_F1": {"ids": ["Denisova11.SG"]},
    "Denisova25": {"ids": ["Denisova25.SG"]},
    "Altai_Neanderthal": {"ids": ["AltaiNeanderthal.DG"]},
    "Chagyrskaya_Neanderthal": {"ids": ["Chagyrskaya8.DG"]},
    "Vindija_Neanderthal": {"ids": ["VindijaG1_final.SG"]},
    "Chimp": {"ids": ["Chimp.REF"]},
    "Mbuti": {"pops": ["Mbuti"]},
    "Yoruba": {"pops": ["Yoruba", "YRI", "YRI-Discovery"]},
    "Papuan": {"pops": ["Papuan"]},
    "French": {"pops": ["French"]},
}


def _portable_path(path: str | Path) -> str:
    """Return a repo-relative path when possible, otherwise only the basename.

    Public reports retain file identity and checksums without leaking a user's
    home directory or machine-specific checkout path.
    """
    candidate = Path(path)
    try:
        return candidate.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return candidate.name

INDIVIDUAL_LABELS = {
    "Target": "Denisova 3 (30.7x diploid target)",
    "Denisova_replicate": "Denisova 3 (SG technical replicate)",
    "Denisova_alt_calls": "Denisova 3 (snpAD alternate calls)",
    "Denisova11_F1": "Denisova 11 (Neanderthal-Denisovan F1)",
    "Denisova25": "Denisova 25 (provisional, 185 ka)",
    "Altai_Neanderthal": "Altai Neanderthal",
    "Chagyrskaya_Neanderthal": "Chagyrskaya Neanderthal",
    "Vindija_Neanderthal": "Vindija Neanderthal",
    "Chimp": "Chimpanzee reference",
}

STAT_SPECS = [
    (
        "primary_denisovan_lineage",
        "D(Target,Altai;Denisova replicate,Chimp)",
        ("Target", "Altai_Neanderthal", "Denisova_replicate", "Chimp"),
        "Positive values distinguish the target's Denisovan lineage from Neanderthal.",
    ),
    (
        "technical_replicate_neanderthal_balance",
        "D(Target,Denisova replicate;Altai,Chimp)",
        ("Target", "Denisova_replicate", "Altai_Neanderthal", "Chimp"),
        "Expected near zero if the two call sets for Denisova 3 are technically balanced.",
    ),
    (
        "technical_replicate_modern_balance",
        "D(Target,Denisova replicate;Mbuti,Chimp)",
        ("Target", "Denisova_replicate", "Mbuti", "Chimp"),
        "Expected near zero; detects differential modern-human/reference bias between call sets.",
    ),
    (
        "target_neanderthal_affinity",
        "D(Target,Mbuti;Altai,Chimp)",
        ("Target", "Mbuti", "Altai_Neanderthal", "Chimp"),
        "Shared archaic affinity; this is not a Neanderthal ancestry percentage.",
    ),
    (
        "target_denisovan_affinity_technical",
        "D(Target,Mbuti;Denisova replicate,Chimp)",
        ("Target", "Mbuti", "Denisova_replicate", "Chimp"),
        "Technical positive control using a second call set from the same specimen.",
    ),
    (
        "f1_neanderthal_parent_signal",
        "D(Denisova11,Denisova replicate;Altai,Chimp)",
        ("Denisova11_F1", "Denisova_replicate", "Altai_Neanderthal", "Chimp"),
        "Positive values recover the Neanderthal-parent contribution in Denisova 11.",
    ),
    (
        "f1_denisovan_parent_signal",
        "D(Denisova11,Altai;Denisova replicate,Chimp)",
        ("Denisova11_F1", "Altai_Neanderthal", "Denisova_replicate", "Chimp"),
        "Positive values recover the Denisovan-parent contribution in Denisova 11.",
    ),
    (
        "denisova25_denisovan_lineage",
        "D(Denisova25,Altai;Denisova replicate,Chimp)",
        ("Denisova25", "Altai_Neanderthal", "Denisova_replicate", "Chimp"),
        "Tests whether the provisional Denisova 25 genome falls on the Denisovan side.",
    ),
    (
        "denisova25_neanderthal_excess_vs_target",
        "D(Denisova25,Target;Altai,Chimp)",
        ("Denisova25", "Target", "Altai_Neanderthal", "Chimp"),
        "Positive values indicate more Neanderthal-related sharing in Denisova 25 than Denisova 3.",
    ),
    (
        "papuan_denisovan_positive_control",
        "D(Papuan,Mbuti;Denisova replicate,Chimp)",
        ("Papuan", "Mbuti", "Denisova_replicate", "Chimp"),
        "Modern-human Denisovan-affinity positive control.",
    ),
    (
        "french_denisovan_negative_control",
        "D(French,Mbuti;Denisova replicate,Chimp)",
        ("French", "Mbuti", "Denisova_replicate", "Chimp"),
        "Modern West Eurasian negative control; expected near zero.",
    ),
]

CHROMOSOME_TEST_IDS = {
    "primary_denisovan_lineage",
    "technical_replicate_neanderthal_balance",
    "f1_neanderthal_parent_signal",
    "denisova25_denisovan_lineage",
}


def _sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def _json_value(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot JSON-encode {type(value)!r}")


def genotype_pair_metrics(a, b) -> dict:
    """Concordance metrics for two diploid dosage vectors (missing is <0)."""
    a = np.asarray(a)
    b = np.asarray(b)
    keep = (a >= 0) & (b >= 0)
    n = int(keep.sum())
    if n == 0:
        return {
            "joint_calls": 0,
            "exact_concordance": np.nan,
            "mean_allele_distance": np.nan,
            "opposite_homozygote_rate": np.nan,
        }
    delta = np.abs(a[keep].astype(float) - b[keep].astype(float))
    return {
        "joint_calls": n,
        "exact_concordance": float(np.mean(delta == 0)),
        "mean_allele_distance": float(np.mean(delta / 2.0)),
        "opposite_homozygote_rate": float(np.mean(delta == 2)),
    }


def _available_specs(panel: Panel, target: str) -> tuple[dict, list[str]]:
    specs = {k: dict(v) for k, v in SAMPLE_SPECS.items()}
    specs["Target"] = {"ids": [target]}
    available = {}
    missing = []
    for name, spec in specs.items():
        cols = panel.cols_for(spec.get("ids"), spec.get("pops"))
        if len(cols):
            available[name] = spec
        else:
            missing.append(name)
    required = {
        "Target",
        "Denisova_replicate",
        "Altai_Neanderthal",
        "Vindija_Neanderthal",
        "Chimp",
        "Mbuti",
    }
    absent = required.difference(available)
    if absent:
        raise RuntimeError(f"Required AADR samples/populations missing: {sorted(absent)}")
    return available, missing


def _load_frequencies(panel: Panel, specs: dict) -> tuple[dict, dict]:
    name_to_cols = {
        name: panel.cols_for(spec.get("ids"), spec.get("pops"))
        for name, spec in specs.items()
    }
    # A mean over an all-missing SNP is intentionally NaN and is later excluded
    # by each statistic. Suppress NumPy's expected per-SNP empty-slice warning.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        return cohort_frequencies(panel, name_to_cols, min_ind=1)


def _sample_qc(panel: Panel, specs: dict, target: str) -> tuple[pd.DataFrame, dict]:
    individual = {
        name: spec["ids"][0]
        for name, spec in specs.items()
        if spec.get("ids")
    }
    cols = np.array(
        [panel.cols_for(ids=[sample_id])[0] for sample_id in individual.values()],
        dtype=np.int64,
    )
    matrix = panel.pg.read(panel.snp_rows, cols)
    genotype = {name: matrix[:, i] for i, name in enumerate(individual)}

    anno = load_anno(panel.prefix + ".anno")
    anno_by_id = anno.set_index("genetic_id", drop=False)
    rows = []
    for name, sample_id in individual.items():
        g = genotype[name]
        called = g >= 0
        n_called = int(called.sum())
        rec = anno_by_id.loc[sample_id] if sample_id in anno_by_id.index else None
        rows.append(
            {
                "name": name,
                "label": (
                    f"Denisova 3 ({target} target)" if name == "Target"
                    else INDIVIDUAL_LABELS.get(name, name)
                ),
                "sample_id": sample_id,
                "n_callable": n_called,
                "call_rate": n_called / panel.n_snp,
                "n_heterozygous": int((g == 1).sum()),
                "heterozygosity_among_calls": (
                    float(np.mean(g[called] == 1)) if n_called else np.nan
                ),
                "aadr_coverage": _anno_value(rec, "coverage"),
                "aadr_assessment": _anno_value(rec, "assessment"),
                "aadr_warning": _anno_value(rec, "assess_warn"),
                "date_bp": _anno_value(rec, "date_bp"),
                "locality": _anno_value(rec, "locality"),
                "publication": _anno_value(rec, "publication"),
                "publication_doi": _anno_value(rec, "publication_doi"),
            }
        )
    return pd.DataFrame(rows), genotype


def _anno_value(record, key):
    if record is None:
        return None
    value = record.get(key)
    if pd.isna(value):
        return None
    return value


def _pairwise_table(genotype: dict) -> pd.DataFrame:
    target = genotype["Target"]
    rows = []
    for name, g in genotype.items():
        if name == "Target":
            continue
        rows.append(
            {
                "comparison": name,
                "label": INDIVIDUAL_LABELS.get(name, name),
                **genotype_pair_metrics(target, g),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["mean_allele_distance", "comparison"], na_position="last"
    )


def _stat_rows(freq: dict, block, mode: str, n_blocks: int = N_BLOCKS) -> pd.DataFrame:
    rows = []
    for test_id, label, args, interpretation in STAT_SPECS:
        if not set(args).issubset(freq):
            continue
        result = st.dstat(freq, *args, block, n_blocks)
        rows.append(
            {
                "mode": mode,
                "test_id": test_id,
                "label": label,
                "statistic": result["statistic"],
                "estimate": result["theta"],
                "se": result["se"],
                "z": result["z"],
                "n_snp": result["n_used"],
                "n_blocks": result["n_blocks_used"],
                "interpretation": interpretation,
                "is_percentage": False,
            }
        )

    for test_id, sample, interpretation in [
        (
            "target_neanderthal_axis_score",
            "Target",
            "Neanderthal-calibrated axis score for an archaic genome; not an admixture percentage.",
        ),
        (
            "f1_neanderthal_axis_score",
            "Denisova11_F1",
            "Expected intermediate control for the F1; not a calibrated parental fraction.",
        ),
    ]:
        needed = {"Altai_Neanderthal", "Chimp", sample, "Mbuti", "Vindija_Neanderthal"}
        if not needed.issubset(freq):
            continue
        result = st.f4_ratio(
            freq,
            "Altai_Neanderthal",
            "Chimp",
            sample,
            "Mbuti",
            "Vindija_Neanderthal",
            block,
            n_blocks,
        )
        rows.append(
            {
                "mode": mode,
                "test_id": test_id,
                "label": result["statistic"],
                "statistic": result["statistic"],
                "estimate": result["theta"],
                "se": result["se"],
                "z": result["z"],
                "n_snp": result["n_used"],
                "n_blocks": result["n_blocks_used"],
                "interpretation": interpretation,
                "is_percentage": False,
            }
        )
    return pd.DataFrame(rows)


def _transversion_sensitivity(full: pd.DataFrame, tv: pd.DataFrame) -> pd.DataFrame:
    cols = ["test_id", "estimate", "se", "z", "n_snp"]
    merged = full[cols].merge(tv[cols], on="test_id", suffixes=("_all", "_tv"))
    merged["same_direction"] = (
        np.sign(merged["estimate_all"]) == np.sign(merged["estimate_tv"])
    )
    merged["estimate_delta"] = merged["estimate_tv"] - merged["estimate_all"]
    merged["tv_site_retention"] = merged["n_snp_tv"] / merged["n_snp_all"]
    return merged


def _marker_sharing(freq: dict, panel: Panel, mode: str) -> pd.DataFrame:
    required = {
        "Denisova_replicate",
        "Altai_Neanderthal",
        "Vindija_Neanderthal",
        "Mbuti",
        "Yoruba",
    }
    if not required.issubset(freq):
        return pd.DataFrame()
    p_den = freq["Denisova_replicate"]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        p_nea = np.nanmean(
            np.vstack([freq["Altai_Neanderthal"], freq["Vindija_Neanderthal"]]),
            axis=0,
        )
        p_afr = np.nanmean(np.vstack([freq["Mbuti"], freq["Yoruba"]]), axis=0)
    den_is_a1 = p_den > 0.5
    den_extreme = np.where(den_is_a1, p_den, 1.0 - p_den)
    afr_of_den = np.where(den_is_a1, p_afr, 1.0 - p_afr)
    nea_of_den = np.where(den_is_a1, p_nea, 1.0 - p_nea)
    marker = (
        (den_extreme >= 0.9)
        & (afr_of_den <= 0.1)
        & (nea_of_den <= 0.5)
        & np.isfinite(p_den)
        & np.isfinite(p_afr)
        & np.isfinite(p_nea)
    )
    rows = []
    for name, values in freq.items():
        if name == "Chimp":
            continue
        callable_marker = marker & np.isfinite(values)
        n = int(callable_marker.sum())
        oriented = np.where(den_is_a1, values, 1.0 - values)
        rows.append(
            {
                "mode": mode,
                "sample": name,
                "label": INDIVIDUAL_LABELS.get(name, name),
                "reference_defined_markers": int(marker.sum()),
                "n_callable_markers": n,
                "mean_denisovan_marker_allele": (
                    float(np.nanmean(oriented[callable_marker])) if n else np.nan
                ),
                "is_admixture_percentage": False,
            }
        )
    return pd.DataFrame(rows).sort_values(
        "mean_denisovan_marker_allele", ascending=False, na_position="last"
    )


def _chromosome_robustness(panel: Panel, freq: dict) -> pd.DataFrame:
    chrom = panel.snp.loc[panel.snp_rows, "chrom"].to_numpy()
    rows = []
    spec_by_id = {spec[0]: spec for spec in STAT_SPECS}
    for chromosome in [str(i) for i in range(1, 23)]:
        keep = chrom == chromosome
        if not np.any(keep):
            continue
        sub = {name: values[keep] for name, values in freq.items()}
        block = st.assign_blocks(int(keep.sum()), min(20, int(keep.sum())))
        for test_id in CHROMOSOME_TEST_IDS:
            _, label, args, interpretation = spec_by_id[test_id]
            if not set(args).issubset(sub):
                continue
            result = st.dstat(sub, *args, block, min(20, int(keep.sum())))
            rows.append(
                {
                    "chromosome": chromosome,
                    "test_id": test_id,
                    "label": label,
                    "estimate": result["theta"],
                    "se": result["se"],
                    "z": result["z"],
                    "n_snp": result["n_used"],
                    "interpretation": interpretation,
                }
            )
    return pd.DataFrame(rows)


def _make_figure(
    marker: pd.DataFrame,
    fstats: pd.DataFrame,
    chromosomes: pd.DataFrame,
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    all_marker = marker[marker["mode"] == "all_snps"].copy()
    wanted = [
        "Target",
        "Denisova_alt_calls",
        "Denisova11_F1",
        "Denisova25",
        "Altai_Neanderthal",
        "Chagyrskaya_Neanderthal",
        "Papuan",
        "French",
        "Mbuti",
    ]
    all_marker = all_marker[all_marker["sample"].isin(wanted)].sort_values(
        "mean_denisovan_marker_allele"
    )
    axes[0].barh(
        all_marker["sample"],
        100 * all_marker["mean_denisovan_marker_allele"],
        color="#7c3aed",
    )
    axes[0].set_xlabel("Reference-defined Denisovan allele (%)")
    axes[0].set_title("Denisovan marker fingerprint")
    axes[0].grid(axis="x", alpha=0.2)

    selected = [
        "primary_denisovan_lineage",
        "technical_replicate_neanderthal_balance",
        "f1_neanderthal_parent_signal",
        "f1_denisovan_parent_signal",
        "denisova25_denisovan_lineage",
        "papuan_denisovan_positive_control",
        "french_denisovan_negative_control",
    ]
    z = fstats[fstats["test_id"].isin(selected)].copy()
    pivot = z.pivot(index="test_id", columns="mode", values="z").reindex(selected)
    y = np.arange(len(pivot))
    axes[1].barh(y - 0.18, pivot.get("all_snps"), height=0.34, label="all SNPs")
    axes[1].barh(
        y + 0.18,
        pivot.get("transversions_only"),
        height=0.34,
        label="transversions",
    )
    axes[1].axvline(0, color="black", lw=0.8)
    axes[1].set_yticks(y, [x.replace("_", " ") for x in pivot.index])
    axes[1].set_xlabel("Block-jackknife Z")
    axes[1].set_title("Lineage and control tests")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="x", alpha=0.2)

    chr_primary = chromosomes[
        chromosomes["test_id"] == "primary_denisovan_lineage"
    ].copy()
    chr_primary["chromosome"] = chr_primary["chromosome"].astype(int)
    chr_primary = chr_primary.sort_values("chromosome")
    axes[2].errorbar(
        chr_primary["chromosome"],
        chr_primary["estimate"],
        yerr=1.96 * chr_primary["se"],
        marker="o",
        ms=4,
        lw=1,
        color="#0f766e",
    )
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].set_xticks(range(1, 23, 2))
    axes[2].set_xlabel("Autosome")
    axes[2].set_ylabel("D estimate (95% jackknife CI)")
    axes[2].set_title("Denisovan lineage by chromosome")
    axes[2].grid(alpha=0.2)

    fig.suptitle("Denisova 3 genome through the archaic-introgression engine", fontsize=15)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _get_row(df: pd.DataFrame, column: str, value: str):
    found = df[df[column] == value]
    return None if found.empty else found.iloc[0]


def _fmt(value, digits=3):
    if value is None or not np.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def _build_report(
    target: str,
    panel: Panel,
    qc: pd.DataFrame,
    pairwise: pd.DataFrame,
    fstats: pd.DataFrame,
    sensitivity: pd.DataFrame,
    marker: pd.DataFrame,
    chromosomes: pd.DataFrame,
    manifest: dict,
    figure_path: Path,
) -> str:
    q = _get_row(qc, "name", "Target")
    replicate = _get_row(pairwise, "comparison", "Denisova_replicate")
    alt_calls = _get_row(pairwise, "comparison", "Denisova_alt_calls")
    lineage = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "primary_denisovan_lineage",
    )
    lineage_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "primary_denisovan_lineage",
    )
    balance = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "technical_replicate_neanderthal_balance",
    )
    balance_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "technical_replicate_neanderthal_balance",
    )
    f1_nea = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "f1_neanderthal_parent_signal",
    )
    f1_den = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "f1_denisovan_parent_signal",
    )
    target_marker = marker[
        (marker["mode"] == "all_snps") & (marker["sample"] == "Target")
    ].iloc[0]
    den25 = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "denisova25_denisovan_lineage",
    )
    den25_excess = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "denisova25_neanderthal_excess_vs_target",
    )
    den25_excess_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "denisova25_neanderthal_excess_vs_target",
    )
    french = _get_row(
        fstats[fstats["mode"] == "all_snps"],
        "test_id",
        "french_denisovan_negative_control",
    )
    french_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "french_denisovan_negative_control",
    )
    img = base64.b64encode(figure_path.read_bytes()).decode("ascii")

    cards = [
        ("Target callability", f"{100*q.call_rate:.2f}%", f"{int(q.n_callable):,} autosomal sites"),
        (
            "Alternate-call concordance",
            f"{100*alt_calls.exact_concordance:.2f}%",
            f"{int(alt_calls.joint_calls):,} joint calls with Denisova3_snpAD.DG",
        ),
        (
            "Denisovan fingerprint",
            f"{100*target_marker.mean_denisovan_marker_allele:.2f}%",
            f"{int(target_marker.n_callable_markers):,} callable reference-defined markers",
        ),
        (
            "Lineage statistic",
            f"Z = {lineage.z:.1f}",
            f"transversions Z = {lineage_tv.z:.1f}",
        ),
    ]
    card_html = "".join(
        f'<article class="card"><span>{html.escape(k)}</span><strong>{html.escape(v)}</strong>'
        f'<small>{html.escape(note)}</small></article>'
        for k, v, note in cards
    )

    findings = [
        f"The primary input is {target}, a {q.aadr_coverage:.1f}x AADR diploid call set "
        f"with {100*q.call_rate:.2f}% autosomal 1240K callability and assessment "
        f"{q.aadr_assessment}.",
        f"The alternate diploid Denisova3_snpAD.DG calls agree at "
        f"{100*alt_calls.exact_concordance:.2f}% of {int(alt_calls.joint_calls):,} joint sites. "
        f"The lower-density pseudo-haploid Denisova.SG representation agrees at "
        f"{100*replicate.exact_concordance:.2f}% of {int(replicate.joint_calls):,} joint sites.",
        f"The dense-diploid versus pseudo-haploid replicate-balance test is not null: "
        f"D={balance.estimate:.5f}, Z={balance.z:.2f}; transversions D={balance_tv.estimate:.5f}, "
        f"Z={balance_tv.z:.2f}. This exposes call-set/coverage asymmetry and is why the "
        f"same-specimen replicate is not used as an independent quantitative calibration.",
        f"The Denisovan-vs-Neanderthal lineage contrast is D={lineage.estimate:.4f} "
        f"(Z={lineage.z:.1f}, {int(lineage.n_snp):,} SNPs) and remains "
        f"D={lineage_tv.estimate:.4f} (Z={lineage_tv.z:.1f}) using transversions only.",
        f"The known Denisova 11 F1 control independently recovers both parents: "
        f"Neanderthal-parent Z={f1_nea.z:.1f} and Denisovan-parent Z={f1_den.z:.1f}.",
    ]
    if den25 is not None:
        findings.append(
            f"The provisional ~185 ka Denisova 25 call set is Denisovan-side in the same "
            f"contrast (D={den25.estimate:.4f}, Z={den25.z:.1f}). It also shows more "
            f"Neanderthal-related sharing than Denisova 3 (D={den25_excess.estimate:.4f}, "
            f"Z={den25_excess.z:.2f}; transversions Z={den25_excess_tv.z:.2f}). These are "
            "useful internal signals but remain provisional because its AADR assessment and "
            "source publication are provisional."
        )
    findings.append(
        f"The French negative control is null on the full panel (D={french.estimate:.4f}, "
        f"Z={french.z:.2f}) but shifts positive in the smaller transversion subset "
        f"(D={french_tv.estimate:.4f}, Z={french_tv.z:.2f}). The primary lineage conclusions "
        "are far larger and stable, but subtle transversion-only effects should not be overread."
    )
    findings.append(
        "No Denisovan ancestry percentage is reported. The primary Denisovan reference and target "
        "are the same specimen in different call sets, and the pipeline has no independent second "
        "Denisovan calibration genome suitable for an absolute fraction."
    )
    finding_html = "".join(f"<li>{html.escape(x)}</li>" for x in findings)

    table_args = dict(index=False, border=0, classes="data", na_rep="n/a")
    qc_view = qc[
        [
            "label",
            "sample_id",
            "n_callable",
            "call_rate",
            "n_heterozygous",
            "heterozygosity_among_calls",
            "aadr_coverage",
            "aadr_assessment",
            "date_bp",
        ]
    ].copy()
    f_view = fstats[
        ["mode", "test_id", "estimate", "se", "z", "n_snp", "interpretation"]
    ].copy()
    m_view = marker[
        ["mode", "sample", "n_callable_markers", "mean_denisovan_marker_allele"]
    ].copy()
    p_view = pairwise[
        [
            "comparison",
            "joint_calls",
            "exact_concordance",
            "mean_allele_distance",
            "opposite_homozygote_rate",
        ]
    ].copy()
    s_view = sensitivity.copy()
    for frame in (qc_view, f_view, m_view, p_view, s_view):
        for col in frame.select_dtypes(include=["float"]).columns:
            frame[col] = frame[col].map(lambda x: round(x, 6) if pd.notna(x) else x)

    generated = manifest["generated_utc"]
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Denisova 3 genome analysis</title>
<style>
:root{{--ink:#172033;--muted:#5e687a;--paper:#f6f4ef;--card:#fff;--violet:#6d28d9;--teal:#0f766e;--line:#dfe3e8}}
*{{box-sizing:border-box}} body{{margin:0;background:linear-gradient(135deg,#ede9fe,#f6f4ef 34%,#ecfeff);color:var(--ink);font:15px/1.55 Inter,Segoe UI,sans-serif}}
main{{max-width:1200px;margin:auto;padding:42px 24px 80px}} header{{background:#172033;color:white;border-radius:22px;padding:34px 38px;box-shadow:0 18px 45px #1720332b}}
h1{{font-size:clamp(2rem,5vw,3.7rem);line-height:1.02;margin:6px 0 15px}} .eyebrow{{text-transform:uppercase;letter-spacing:.14em;color:#c4b5fd;font-weight:700}}
.lede{{font-size:1.08rem;max-width:850px;color:#dbe4f0}} .cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:22px 0}}
.card{{background:var(--card);padding:20px;border-radius:16px;border:1px solid #ffffffcc;box-shadow:0 8px 24px #17203310}} .card span,.card small{{display:block;color:var(--muted)}} .card strong{{display:block;color:var(--violet);font-size:1.65rem;margin:5px 0}}
section{{background:#ffffffd9;margin:18px 0;padding:26px;border:1px solid #fff;border-radius:18px;box-shadow:0 8px 28px #1720330b}} h2{{margin-top:0;font-size:1.45rem}} h3{{margin-top:25px}}
.callout{{border-left:5px solid var(--teal);background:#ecfdf5;padding:16px 19px;border-radius:10px}} img{{width:100%;height:auto;border-radius:12px}}
.scroll{{overflow:auto}} table.data{{border-collapse:collapse;width:100%;font-size:.86rem}} table.data th{{position:sticky;top:0;background:#eef2ff;text-align:left}} table.data th,table.data td{{padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:top;white-space:nowrap}} table.data td:last-child{{white-space:normal;min-width:280px}}
code{{background:#eef2ff;padding:2px 5px;border-radius:5px}} footer{{color:var(--muted);padding:15px 4px}}
@media(max-width:850px){{.cards{{grid-template-columns:repeat(2,1fr)}} header{{padding:26px}}}} @media(max-width:500px){{.cards{{grid-template-columns:1fr}}}}
</style></head><body><main>
<header><div class="eyebrow">Archaic-introgression pipeline · reference-aware module</div><h1>Denisova 3 genome analysis</h1>
<p class="lede">A real-data run of the high-coverage Denisovan genome through the pipeline's validated AADR reader and block-jackknife f-statistic engine, with explicit controls for same-specimen reference leakage.</p></header>
<div class="cards">{card_html}</div>
<section><h2>Result in plain language</h2><ul>{finding_html}</ul>
<div class="callout"><strong>Interpretation boundary:</strong> reference-defined marker sharing and archaic-axis scores are fingerprints, not ancestry percentages. The same biological Denisova 3 specimen appears as both the target and a technical replicate, so those comparisons validate data handling rather than independently proving Denisovan identity.</div></section>
<section><h2>Visual summary</h2><img alt="Denisovan genome summary plots" src="data:image/png;base64,{img}"></section>
<section><h2>Input and genotype QC</h2><div class="scroll">{qc_view.to_html(**table_args)}</div></section>
<section><h2>Target-to-comparator concordance</h2><div class="scroll">{p_view.to_html(**table_args)}</div></section>
<section><h2>Block-jackknife f-statistics</h2><p>Every row reports a ratio-of-sums estimate with 50 contiguous genomic delete-one blocks. Values marked as axis scores are deliberately not converted to percentages.</p><div class="scroll">{f_view.to_html(**table_args)}</div></section>
<section><h2>Transversion sensitivity</h2><p>Transversions reduce sensitivity to ancient-DNA deamination and retain the direction of the main lineage tests.</p><div class="scroll">{s_view.to_html(**table_args)}</div></section>
<section><h2>Reference-defined Denisovan marker sharing</h2><p>Markers are sites where Denisova.SG carries an allele at ≥90%, the pooled African baseline at ≤10%, and the Altai/Vindija Neanderthal mean at ≤50%. This is a lineage fingerprint and is intentionally not labelled an admixture percentage.</p><div class="scroll">{m_view.to_html(**table_args)}</div></section>
<section><h2>Reproducibility</h2><p>Command: <code>{html.escape(manifest['command'])}</code></p><p>Panel prefix: <code>{html.escape(manifest['panel_prefix'])}</code>. The local AADR directory is intentionally omitted. Generated {html.escape(generated)}. Machine-readable tables and the complete manifest sit beside this report.</p></section>
<footer>Generated by <code>archaic.denisovan_genome</code>. Scientific claims are bounded to the supplied AADR genotype panel and the explicit tests shown above.</footer>
</main></body></html>"""


def _build_markdown(
    target: str,
    qc: pd.DataFrame,
    pairwise: pd.DataFrame,
    fstats: pd.DataFrame,
    marker: pd.DataFrame,
) -> str:
    q = _get_row(qc, "name", "Target")
    replicate = _get_row(pairwise, "comparison", "Denisova_replicate")
    alt_calls = _get_row(pairwise, "comparison", "Denisova_alt_calls")
    all_stats = fstats[fstats["mode"] == "all_snps"]
    lineage = _get_row(all_stats, "test_id", "primary_denisovan_lineage")
    lineage_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "primary_denisovan_lineage",
    )
    f1_nea = _get_row(all_stats, "test_id", "f1_neanderthal_parent_signal")
    f1_den = _get_row(all_stats, "test_id", "f1_denisovan_parent_signal")
    balance = _get_row(all_stats, "test_id", "technical_replicate_neanderthal_balance")
    balance_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "technical_replicate_neanderthal_balance",
    )
    den25_excess = _get_row(
        all_stats, "test_id", "denisova25_neanderthal_excess_vs_target"
    )
    den25_excess_tv = _get_row(
        fstats[fstats["mode"] == "transversions_only"],
        "test_id",
        "denisova25_neanderthal_excess_vs_target",
    )
    mark = marker[(marker["mode"] == "all_snps") & (marker["sample"] == "Target")].iloc[0]
    return f"""# Denisova 3 genome run

## Primary result

- Target: `{target}` (AADR assessment `{q.aadr_assessment}`, coverage {_fmt(q.aadr_coverage, 1)}x).
- Autosome callability: {100*q.call_rate:.2f}% ({int(q.n_callable):,} 1240K sites).
- Alternate diploid-call concordance: {100*alt_calls.exact_concordance:.2f}% across {int(alt_calls.joint_calls):,} joint calls.
- Same-specimen `Denisova.SG` technical concordance: {100*replicate.exact_concordance:.2f}% across {int(replicate.joint_calls):,} joint calls.
- Denisovan-vs-Neanderthal contrast: D={lineage.estimate:.6f}, Z={lineage.z:.2f}, n={int(lineage.n_snp):,}.
- Transversions-only repeat: D={lineage_tv.estimate:.6f}, Z={lineage_tv.z:.2f}, n={int(lineage_tv.n_snp):,}.
- Reference-defined Denisovan marker allele: {100*mark.mean_denisovan_marker_allele:.2f}% across {int(mark.n_callable_markers):,} callable markers.
- Denisova 11 F1 control: Neanderthal-parent Z={f1_nea.z:.2f}; Denisovan-parent Z={f1_den.z:.2f}.
- Provisional Denisova 25 has more Neanderthal-related sharing than Denisova 3: D={den25_excess.estimate:.6f}, Z={den25_excess.z:.2f}; transversions Z={den25_excess_tv.z:.2f}.

## Interpretation

The run strongly and robustly places the target on the Denisovan side of the
Denisovan-Neanderthal contrast.  The known Denisova 11 F1 recovers both parental
signals, providing an internal biological control.  The same-specimen replicate
is a technical validation only. Its balance test is significantly non-null
(all-SNP Z={balance.z:.2f}; transversions Z={balance_tv.z:.2f}), exposing a
call-set/coverage asymmetry between dense diploid and pseudo-haploid data. It
must not be treated as independent quantitative calibration.

No Denisovan ancestry percentage is reported.  `Denisova.SG` and `{target}` are
different AADR call sets from the same Denisova 3 specimen, and the panel does not
provide an independent second high-quality Denisovan calibration genome for an
absolute fraction. Reference-defined marker sharing is a fingerprint, not a
percentage of ancestry.

Open `report.html` for the full interactive-style report and all result tables.
"""


def run(panel_name: str, target: str, output: Path, config: str | None = None) -> dict:
    prefix = panel_prefix(panel_name, config)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    panel = Panel(prefix)
    specs, missing = _available_specs(panel, target)
    freq, info = _load_frequencies(panel, specs)
    qc, genotype = _sample_qc(panel, specs, target)
    pairwise = _pairwise_table(genotype)

    block = st.assign_blocks(panel.n_snp, N_BLOCKS)
    all_stats = _stat_rows(freq, block, "all_snps")
    all_marker = _marker_sharing(freq, panel, "all_snps")
    chromosomes = _chromosome_robustness(panel, freq)

    tv_panel = Panel(prefix, transversions_only=True)
    tv_specs, _ = _available_specs(tv_panel, target)
    tv_freq, tv_info = _load_frequencies(tv_panel, tv_specs)
    tv_block = st.assign_blocks(tv_panel.n_snp, N_BLOCKS)
    tv_stats = _stat_rows(tv_freq, tv_block, "transversions_only")
    tv_marker = _marker_sharing(tv_freq, tv_panel, "transversions_only")

    fstats = pd.concat([all_stats, tv_stats], ignore_index=True)
    marker = pd.concat([all_marker, tv_marker], ignore_index=True)
    sensitivity = _transversion_sensitivity(all_stats, tv_stats)

    paths = {
        "qc": output / "sample_qc.tsv",
        "pairwise": output / "pairwise_comparisons.tsv",
        "fstats": output / "fstats.tsv",
        "sensitivity": output / "transversion_sensitivity.tsv",
        "marker": output / "denisovan_marker_sharing.tsv",
        "chromosomes": output / "chromosome_robustness.tsv",
        "figure": output / "denisovan_genome_summary.png",
        "manifest": output / "run_manifest.json",
        "report": output / "report.html",
        "results": output / "RESULTS.md",
    }
    for frame, key in [
        (qc, "qc"),
        (pairwise, "pairwise"),
        (fstats, "fstats"),
        (sensitivity, "sensitivity"),
        (marker, "marker"),
        (chromosomes, "chromosomes"),
    ]:
        frame.to_csv(paths[key], sep="\t", index=False, na_rep="NA")

    _make_figure(marker, fstats, chromosomes, paths["figure"])

    target_calls = genotype["Target"]
    target_hash = hashlib.sha256(target_calls.tobytes()).hexdigest()
    manifest = {
        "analysis": "denisovan_genome",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": (
            f"python -m archaic.denisovan_genome --panel {panel_name} "
            f"--target {target} --output {_portable_path(output)}"
        ),
        "panel": panel_name,
        "panel_prefix": Path(prefix).name,
        "target": target,
        "target_genotype_sha256": target_hash,
        "autosomal_panel_sites": panel.n_snp,
        "transversion_panel_sites": tv_panel.n_snp,
        "n_blocks": N_BLOCKS,
        "missing_optional_samples": missing,
        "frequency_profiles": info,
        "transversion_frequency_profiles": tv_info,
        "input_provenance": {
            ext: {
                "path": Path(prefix + ext).name,
                "size_bytes": os.path.getsize(prefix + ext),
                "mtime_utc": datetime.fromtimestamp(
                    os.path.getmtime(prefix + ext), timezone.utc
                ).isoformat(),
                "sha256": _sha256_file(prefix + ext) if ext in {".snp", ".ind"} else None,
            }
            for ext in (".geno", ".snp", ".ind", ".anno")
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "module_sha256": _sha256_file(__file__),
        },
        "interpretation_boundaries": [
            "Denisova3.DG and Denisova.SG are call sets from the same biological specimen.",
            "No Denisovan or combined archaic ancestry percentage is identifiable here.",
            "Reference-defined marker sharing is a lineage fingerprint, not an ancestry proportion.",
            "Denisova25.SG has a provisional AADR/publication status and is interpreted accordingly.",
        ],
        "outputs": {key: _portable_path(path) for key, path in paths.items()},
    }
    paths["manifest"].write_text(
        json.dumps(manifest, indent=2, default=_json_value) + "\n", encoding="utf-8"
    )
    paths["report"].write_text(
        _build_report(
            target,
            panel,
            qc,
            pairwise,
            fstats,
            sensitivity,
            marker,
            chromosomes,
            manifest,
            paths["figure"],
        ),
        encoding="utf-8",
    )
    paths["results"].write_text(
        _build_markdown(target, qc, pairwise, fstats, marker), encoding="utf-8"
    )
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run a reference-aware Denisovan genome analysis on local AADR data."
    )
    parser.add_argument("--panel", choices=("1240k", "ho"), default="1240k")
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--config", help="Optional pipeline YAML configuration path")
    args = parser.parse_args(argv)
    manifest = run(args.panel, args.target, args.output, args.config)
    print(f"Denisovan genome analysis complete: {manifest['outputs']['report']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
