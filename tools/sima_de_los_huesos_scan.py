#!/usr/bin/env python3
"""Run the AADR f4-ratio scan on published Sima de los Huesos BAMs.

This is intentionally a small adapter, not a general ancient-DNA caller.  The
Meyer et al. (2016) ENA BAMs are already duplicate-removed and filtered at
read length >=35 bp and mapping quality >=30.  We make deterministic
pseudo-haploid calls only at 1240K AADR SNPs, retain bases with base quality
>=30, and run the repository's existing Phase-3 f4/D estimators.

The results are exploratory: these ~430 ka samples have extremely sparse
nuclear coverage and are not part of the AADR release used by the main
pipeline.  Transition-only damage is checked with a transversion-only repeat.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    import bamnostic
except ImportError as exc:  # pragma: no cover - documented operational guard
    raise SystemExit(
        "This adapter needs the pure-Python 'bamnostic' package. "
        "Install it in an isolated environment, then rerun."
    ) from exc

from archaic import stats as st
from archaic.panel import Panel
from archaic.refs import PANELS

NUCLEAR_BAMS = {
    "femurXIII": ("ERR995357", "femurXIII.L35MQ30.bam", "09d66d541a63b4e93bee59ad505770dd"),
    "femur_fragment": ("ERR995361", "femur_fragment.L35MQ30.bam", "f21cd599ad206caeda0b24be5f180a8c"),
    "incisor": ("ERR995358", "incisor.L35MQ30.bam", "660f1ab7f5ced32eb32169f0b71d6492"),
    "scapula": ("ERR995359", "scapula.L35MQ30.bam", "0b70f2c1597e2d461094786e40b89c39"),
    "molar": ("ERR995360", "molar.L35MQ30.bam", "9ecb7d4c628f812d62ba15a72df03f56"),
}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pileup_at_panel_sites(path: Path, targets: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]],
                          min_mapq: int, min_baseq: int) -> dict[int, list[str]]:
    """Sequentially scan an unindexed BAM and collect AADR-compatible bases."""
    calls: dict[int, list[str]] = defaultdict(list)
    bam = bamnostic.AlignmentFile(str(path), "rb", require_index=False, check_sq=False)
    try:
        for read in bam:
            if (read.is_unmapped or read.is_secondary or read.is_supplementary or
                    read.is_duplicate or read.mapq < min_mapq):
                continue
            chrom = int(read.reference_id) + 1  # hg19 BAM ref IDs are chr1=0 ... chr22=21
            if chrom not in targets:
                continue
            positions, rows, alleles = targets[chrom]
            seq = read.query_sequence.upper()
            quals = read.query_qualities
            qpos, rpos = 0, int(read.reference_start) + 1  # BAM is 0-based; AADR is 1-based
            for op, length in read.cigartuples:
                if op in (0, 7, 8):  # M, =, X: consume both query and reference
                    left = np.searchsorted(positions, rpos, side="left")
                    right = np.searchsorted(positions, rpos + length, side="left")
                    for idx in range(left, right):
                        offset = int(positions[idx] - rpos)
                        if qpos + offset >= len(seq) or quals[qpos + offset] < min_baseq:
                            continue
                        base = seq[qpos + offset]
                        a1, a2 = alleles[idx]
                        if base == a1 or base == a2:
                            calls[int(rows[idx])].append(base)
                    qpos += length
                    rpos += length
                elif op in (1, 4):  # insertion or soft clip: query only
                    qpos += length
                elif op in (2, 3):  # deletion or skipped region: reference only
                    rpos += length
                # hard clip/pad do not advance either coordinate
    finally:
        bam.close()
    return calls


def deterministic_pseudohaploid(calls: dict[int, list[str]], alleles: np.ndarray,
                                n_snp: int, transversions_only: bool) -> tuple[np.ndarray, int]:
    """Choose one observation per locus by a stable hash, yielding allele-1 frequency."""
    out = np.full(n_snp, np.nan, dtype=np.float32)
    kept = 0
    for row, bases in calls.items():
        a1, a2 = alleles[row]
        if transversions_only and {a1, a2} not in ({"A", "C"}, {"A", "T"}, {"C", "G"}, {"G", "T"}):
            continue
        # Stable, unbiased within-site choice avoids presenting a majority-vote as a diploid call.
        base = bases[row % len(bases)]
        out[row] = 1.0 if base == a1 else 0.0
        kept += 1
    return out, kept


def estimate(panel: Panel, rf: dict[str, np.ndarray], p_x: np.ndarray) -> dict[str, float | int]:
    p_alt, p_chi, p_vin, p_den, p_mb = (rf[k].astype(np.float32) for k in
                                         ["Altai", "Chimp", "Vindija", "Denisova", "Mbuti"])
    starts = st.block_starts(panel.n_snp, 50)
    diff_xm = p_x - p_mb
    den_x = p_x + p_mb - 2 * p_x * p_mb
    ax_n = p_alt - p_chi
    ax_d = p_den - p_chi
    alpha, alpha_se, _, alpha_n = st.batch_jackknife_ratio(
        (ax_n * diff_xm)[:, None], (ax_n * (p_vin - p_mb))[:, None], starts)
    d_n, d_n_se, d_n_z, d_n_n = st.batch_jackknife_ratio(
        (diff_xm * ax_n)[:, None], (den_x * (p_alt + p_chi - 2 * p_alt * p_chi))[:, None], starts)
    d_d, d_d_se, d_d_z, d_d_n = st.batch_jackknife_ratio(
        (diff_xm * ax_d)[:, None], (den_x * (p_den + p_chi - 2 * p_den * p_chi))[:, None], starts)
    return {
        "alpha_Nea": float(alpha[0]), "alpha_SE": float(alpha_se[0]), "alpha_nSNP": int(alpha_n[0]),
        "D_Nea": float(d_n[0]), "D_Nea_SE": float(d_n_se[0]), "D_Nea_Z": float(d_n_z[0]), "D_Nea_nSNP": int(d_n_n[0]),
        "D_Den": float(d_d[0]), "D_Den_SE": float(d_d_se[0]), "D_Den_Z": float(d_d_z[0]), "D_Den_nSNP": int(d_d_n[0]),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-baseq", type=int, default=30)
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    panel = Panel(PANELS["1240k"]["prefix"], autosomes_only=True)
    print("Loading AADR reference genotypes once...", flush=True)
    reference_freqs, _ = panel.frequencies({
        k: PANELS["1240k"]["refs"][k] for k in ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]
    })
    snp = panel.snp.iloc[panel.snp_rows].copy()
    # The estimator uses compact autosomal row positions, while the BAM lookup
    # starts from source .geno row numbers.  Make that correspondence explicit.
    source_to_compact = np.full(len(panel.snp), -1, dtype=np.int64)
    source_to_compact[panel.snp_rows] = np.arange(panel.n_snp, dtype=np.int64)
    all_alleles = snp[["a1", "a2"]].to_numpy(dtype=str)
    targets = {}
    for chrom in range(1, 23):
        sub = snp[snp.chrom.astype(str).eq(str(chrom))]
        targets[chrom] = (sub.pos.to_numpy(dtype=np.int64), source_to_compact[sub.index.to_numpy(dtype=np.int64)],
                          sub[["a1", "a2"]].to_numpy(dtype=str))

    per_sample: dict[str, dict[int, list[str]]] = {}
    source = []
    for sample, (ena_run, filename, expected_md5) in NUCLEAR_BAMS.items():
        path = args.bam_dir / filename
        actual_md5 = md5(path)
        if actual_md5 != expected_md5:
            raise RuntimeError(f"{path}: MD5 {actual_md5} does not match ENA record {expected_md5}")
        print(f"Scanning {sample}: {path.name}", flush=True)
        per_sample[sample] = pileup_at_panel_sites(path, targets, args.min_mapq, args.min_baseq)
        source.append({"sample": sample, "file": filename, "bytes": path.stat().st_size,
                       "md5": actual_md5, "ena_run": ena_run})

    combined: dict[int, list[str]] = defaultdict(list)
    for calls in per_sample.values():
        for row, bases in calls.items():
            combined[row].extend(bases)
    per_sample["Sima_combined"] = combined

    rows = []
    for sample, calls in per_sample.items():
        for transversions in (False, True):
            p_x, n_called = deterministic_pseudohaploid(calls, all_alleles, panel.n_snp, transversions)
            row = {"sample": sample, "mode": "transversions_only" if transversions else "all_snps",
                   "pseudo_haploid_AADR_calls": n_called, **estimate(panel, reference_freqs, p_x)}
            rows.append(row)
            print(f"{sample} {row['mode']}: {n_called:,} calls, alpha_nSNP={row['alpha_nSNP']:,}", flush=True)

    fields = list(rows[0])
    with (args.output / "sima_aadr_f4_results.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "study": "PRJEB10597", "publication": "Meyer et al. 2016, Nature 531:504-507",
        "method": "Published L35/MQ30 BAMs; AADR 1240K-targeted deterministic pseudo-haploid calls; Phase-3 f4/D estimators",
        "limits": [
            "Not a whole-genome call set or a replacement for the AADR panel.",
            "Sparse coverage means individual estimates are exploratory.",
            "The transition/transversion comparison is a damage-sensitivity check, not full ancient-DNA authentication.",
        ],
        "filters": {"min_mapping_quality": args.min_mapq, "min_base_quality": args.min_baseq},
        "inputs": source,
    }
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
