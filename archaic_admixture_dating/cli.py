"""Command-line interface for the Papuan Denisovan admixture-dating V1 module."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd

from .bootstrap import bootstrap_fit, interval
from .calibration import calibrate_simulations
from .checkpointing import (
    RESUMABLE_EXIT_CODE,
    CheckpointStore,
    Deadline,
    atomic_write_json,
    atomic_write_text,
    file_fingerprint,
)
from .config import apply_profile, config_hash, load_config, save_snapshot
from .dating_continuous import fit_continuous_flow
from .dating_single_pulse import fit_single_pulse
from .dating_two_pulse import fit_two_pulse
from .diagnostics import collect_warnings
from .downloads import Deadline as DownloadDeadline
from .downloads import DownloadError, download, verify_record
from .logging_utils import get_logger
from .manifests import DownloadManifest, DownloadRecord
from .model_comparison import compare_models
from .plotting import make_qc_figure
from .reporting import write_report
from .simulations import simulate_model_set, simulate_tracts
from .sensitivity import run_sensitivity
from .tract_filtering import filter_tracts, load_bed, overlapping_pairs
from .tract_import import import_tracts, parse_column_map
from .tract_schema import read_tracts, validate_tracts, write_tracts
from .tract_summary import summarize_tracts

LOG = get_logger("archaic_admixture_dating")
REPOSITORY_ROOT = Path(__file__).resolve().parent.parent


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", help="YAML configuration path")
    parser.add_argument("--profile", choices=["smoke", "laptop", "full"], default=None)
    parser.add_argument("--output", help="Run output directory")
    parser.add_argument("--run-id", default=None, help="Stable run identifier for resume")
    parser.add_argument("--max-block-minutes", type=float)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--memory-gb", type=float)
    parser.add_argument("--chromosomes", nargs="+")
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--verbose", action="store_true")


def _resolved(args: argparse.Namespace) -> dict[str, Any]:
    overrides: dict[str, Any] = {"runtime": {}}
    if args.max_block_minutes is not None:
        overrides["runtime"]["max_block_minutes"] = args.max_block_minutes
    if args.threads is not None:
        overrides["runtime"]["max_threads"] = args.threads
    if args.memory_gb is not None:
        overrides["runtime"]["max_memory_gb"] = args.memory_gb
    config = load_config(args.config, overrides=overrides)
    return apply_profile(config, args.profile)


def _run_dir(args: argparse.Namespace, config: dict[str, Any]) -> Path:
    if args.output:
        return Path(args.output).resolve()
    profile = config["runtime"]["profile"]
    run_id = args.run_id or f"papuan_denisovan_v1_{profile}"
    root = Path(config["runtime"]["output_root"])
    if not root.is_absolute():
        root = REPOSITORY_ROOT / root
    return (root / run_id).resolve()


def _deadline(config: dict[str, Any]) -> Deadline:
    runtime = config["runtime"]
    return Deadline(
        float(runtime["max_block_minutes"]),
        float(runtime["stop_before_limit_minutes"]),
    )


def _tract_path(args: argparse.Namespace, run_dir: Path) -> Path:
    explicit = getattr(args, "input", None)
    return Path(explicit).resolve() if explicit else run_dir / "tracts" / "filtered_tracts.tsv"


def _select(frame: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    selected = frame
    if args.chromosomes:
        chroms = {str(value).replace("chr", "") for value in args.chromosomes}
        selected = selected.loc[selected["chromosome"].astype(str).isin(chroms)]
    if args.samples:
        selected = selected.loc[selected["sample_id"].isin(args.samples)]
    return selected.reset_index(drop=True)


def _write_frames(frames: dict[str, pd.DataFrame], directory: Path) -> list[Path]:
    paths: list[Path] = []
    directory.mkdir(parents=True, exist_ok=True)
    for name, frame in frames.items():
        path = directory / f"{name}.tsv"
        write_tracts(frame, path)
        paths.append(path)
    return paths


def _fit_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "minimum_length_cm": float(config["tracts"]["minimum_length_cm"]),
        "generation_time_years": float(config["project"]["generation_time_years"]),
    }


def _configured_masks(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    masks: dict[str, pd.DataFrame] = {}
    paths = config["tracts"].get("masks", {})
    enabled = {
        "low_mappability": config["tracts"].get("exclude_low_mappability", True),
        "centromeres": config["tracts"].get("exclude_centromeres", True),
        "telomeres": config["tracts"].get("exclude_telomere_buffers", True),
        "selected_loci": bool(config["tracts"].get("selected_loci_bed")),
    }
    if config["tracts"].get("selected_loci_bed"):
        paths = {**paths, "selected_loci": config["tracts"]["selected_loci_bed"]}
    for label, path in paths.items():
        if not path or not enabled.get(label, True):
            continue
        target = Path(path)
        if not target.exists():
            raise ValueError(f"Configured {label} mask does not exist: {target}")
        masks[label] = load_bed(target)
    return masks


def command_init(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    for name in ("checkpoints", "downloads", "tracts", "qc", "models", "simulations", "report"):
        (output / name).mkdir(parents=True, exist_ok=True)
    save_snapshot(config, output / "config.snapshot.yaml")
    atomic_write_json(
        output / "run.json",
        {
            "project": config["project"]["name"],
            "profile": config["runtime"]["profile"],
            "config_hash": config_hash(config),
            "state": "initialized",
        },
    )
    print(output)
    return 0


def command_inspect(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    frame = read_tracts(_tract_path(args, output))
    valid, excluded = validate_tracts(frame)
    valid = _select(valid, args)
    result = {
        "rows": len(frame),
        "valid_rows": len(valid),
        "excluded_rows": len(excluded),
        "samples": int(valid["sample_id"].nunique()),
        "populations": sorted(valid["population"].dropna().astype(str).unique()),
        "chromosomes": sorted(valid["chromosome"].dropna().astype(str).unique()),
        "length_cm_range": [
            float(valid["length_cm"].min()) if len(valid) else None,
            float(valid["length_cm"].max()) if len(valid) else None,
        ],
    }
    print(json.dumps(result, indent=2))
    return 0


def command_estimate_storage(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    datasets = config["downloads"].get("datasets", {})
    expected = sum(int(value.get("expected_size", 0) or 0) for value in datasets.values())
    available = shutil.disk_usage(output.parent if output.parent.exists() else REPOSITORY_ROOT).free
    estimate = {
        "configured_download_bytes": expected,
        "recommended_with_headroom_bytes": int(expected * 1.10),
        "available_bytes": available,
        "simulation_replicates": int(config["simulation"]["final_replicates_per_model"]) * len(config["models"]),
        "note": "Unknown-size or controlled-access datasets are excluded and will not auto-download.",
    }
    print(json.dumps(estimate, indent=2))
    return 0


def _dataset_record(config: dict[str, Any], dataset_id: str, output: Path) -> DownloadRecord:
    dataset = config["downloads"].get("datasets", {}).get(dataset_id)
    if not dataset:
        raise DownloadError(f"Dataset {dataset_id!r} is not configured")
    return DownloadRecord(
        dataset_id=dataset_id,
        source_url=dataset["url"],
        destination=str((output / "downloads" / dataset.get("filename", dataset_id)).resolve()),
        expected_size=dataset.get("expected_size"),
        checksum_algorithm=dataset.get("checksum_algorithm"),
        checksum=dataset.get("checksum"),
        access=dataset.get("access", "public"),
    )


def command_download(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    record = _dataset_record(config, args.dataset, output)
    if (
        record.expected_size is not None
        and record.expected_size >= 1024**3
        and not args.dry_run
        and not args.force
    ):
        raise DownloadError(
            "Multi-gigabyte download requires an explicit preflight: rerun with "
            "--dry-run, then repeat with --force to confirm the scoped transfer"
        )
    manifest = DownloadManifest(output / "downloads" / "download_manifest.json")
    settings = config["downloads"]
    result = download(
        record,
        manifest,
        resume=args.resume,
        dry_run=args.dry_run,
        chunk_size_mb=settings["chunk_size_mb"],
        retries=settings["retries"],
        retry_backoff_seconds=settings["retry_backoff_seconds"],
        bandwidth_limit_mbps=settings.get("bandwidth_limit_mbps"),
        deadline=DownloadDeadline(
            float(settings["max_block_minutes"]),
            float(config["runtime"]["stop_before_limit_minutes"]),
        ),
    )
    print(json.dumps(result.__dict__, indent=2))
    return RESUMABLE_EXIT_CODE if result.completion_state == "paused" else 0


def command_download_status(args: argparse.Namespace) -> int:
    config = _resolved(args)
    manifest = DownloadManifest(_run_dir(args, config) / "downloads" / "download_manifest.json")
    print(json.dumps(manifest.load(), indent=2))
    return 0


def command_verify_downloads(args: argparse.Namespace) -> int:
    config = _resolved(args)
    manifest = DownloadManifest(_run_dir(args, config) / "downloads" / "download_manifest.json")
    data = manifest.load()
    results = {}
    for dataset_id, value in data.get("records", {}).items():
        results[dataset_id] = verify_record(DownloadRecord(**value))
    print(json.dumps(results, indent=2))
    return 0 if all(results.values()) else 1


def command_import(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    valid, excluded = import_tracts(
        args.input,
        caller=args.caller,
        column_map=parse_column_map(args.map_column),
        population=args.population,
        output=output / "tracts" / "imported_tracts.tsv",
        excluded_output=output / "tracts" / "malformed_tracts.tsv",
    )
    atomic_write_json(
        output / "input_manifest.json",
        {
            "input": file_fingerprint(args.input),
            "caller": args.caller,
            "population_override": args.population,
            "access_authorization": "user-supplied path; authorization must be documented by operator",
        },
    )
    print(json.dumps({"imported": len(valid), "malformed": len(excluded)}, indent=2))
    return 0


def command_qc(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    source = _tract_path(args, output)
    frame, malformed = validate_tracts(read_tracts(source))
    frame = _select(frame, args)
    tracts = config["tracts"]
    retained, excluded = filter_tracts(
        frame,
        minimum_length_cm=float(tracts["minimum_length_cm"]),
        minimum_confidence=float(tracts["minimum_confidence"]),
        minimum_callable_fraction=float(tracts["minimum_callable_fraction"]),
        masks=_configured_masks(config),
    )
    write_tracts(retained, output / "tracts" / "filtered_tracts.tsv")
    write_tracts(excluded, output / "tracts" / "qc_excluded_tracts.tsv")
    if len(malformed):
        write_tracts(malformed, output / "tracts" / "malformed_tracts.tsv")
    overlaps = overlapping_pairs(retained)
    write_tracts(overlaps, output / "qc" / "overlapping_tracts.tsv")
    _write_frames(summarize_tracts(retained), output / "qc")
    print(json.dumps({"retained": len(retained), "excluded": len(excluded), "overlaps": len(overlaps)}, indent=2))
    return 0


def _print_fit(fit: dict[str, Any], output: Path, name: str) -> int:
    target = output / "models" / f"{name}.json"
    atomic_write_json(target, fit)
    print(json.dumps(fit, indent=2))
    return 0


def command_fit_single(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts = _select(read_tracts(_tract_path(args, output)), args)
    bounds = tuple(config["dating"]["single_pulse"]["bounds_generations"])
    fit = fit_single_pulse(tracts["length_cm"], bounds_generations=bounds, **_fit_kwargs(config))
    return _print_fit(fit, output, "single_pulse")


def command_fit_two(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts = _select(read_tracts(_tract_path(args, output)), args)
    fit = fit_two_pulse(
        tracts["length_cm"],
        minimum_separation_generations=config["dating"]["two_pulse"]["minimum_separation_generations"],
        **_fit_kwargs(config),
    )
    return _print_fit(fit, output, "two_pulse")


def command_fit_continuous(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts = _select(read_tracts(_tract_path(args, output)), args)
    fit = fit_continuous_flow(tracts["length_cm"], **_fit_kwargs(config))
    return _print_fit(fit, output, "continuous_flow")


def command_simulate(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    output.mkdir(parents=True, exist_ok=True)
    replicates = args.replicates or int(config["simulation"]["final_replicates_per_model"])
    summary, paused = simulate_model_set(
        config,
        output,
        replicates=replicates,
        resume=args.resume,
        deadline=_deadline(config),
        config_digest=config_hash(config),
    )
    print(json.dumps({"completed_summaries": len(summary), "paused": paused}, indent=2))
    return RESUMABLE_EXIT_CODE if paused else 0


def command_compare(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts = _select(read_tracts(_tract_path(args, output)), args)
    table, fits = compare_models(
        tracts["length_cm"],
        single_bounds=tuple(config["dating"]["single_pulse"]["bounds_generations"]),
        two_minimum_separation_generations=config["dating"]["two_pulse"]["minimum_separation_generations"],
        **_fit_kwargs(config),
    )
    write_tracts(table, output / "models" / "model_comparison.tsv")
    atomic_write_json(output / "models" / "fits.json", fits)
    print(table.drop(columns=["fit_json"]).to_string(index=False))
    return 0


def command_sensitivity(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts = _select(read_tracts(_tract_path(args, output)), args)
    table = run_sensitivity(tracts, config)
    write_tracts(table, output / "models" / "sensitivity.tsv")
    print(table.to_string(index=False))
    return 0


def _load_report_inputs(output: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, Any]]:
    tracts = read_tracts(output / "tracts" / "filtered_tracts.tsv")
    summaries = {}
    for name in (
        "overall",
        "individual",
        "population",
        "chromosome",
        "sample_missingness",
        "caller_batch",
        "source_class",
    ):
        path = output / "qc" / f"{name}.tsv"
        if path.exists():
            summaries[name] = read_tracts(path)
    model_table = read_tracts(output / "models" / "model_comparison.tsv")
    with (output / "models" / "fits.json").open("r", encoding="utf-8") as handle:
        fits = json.load(handle)
    return tracts, summaries, model_table, fits


def command_report(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    tracts, summaries, model_table, fits = _load_report_inputs(output)
    bootstrap_path = output / "models" / "bootstrap_single.tsv"
    bootstrap = read_tracts(bootstrap_path) if bootstrap_path.exists() else None
    figure = make_qc_figure(tracts, model_table, output / "report" / "qc_model_comparison.png")
    simulation_path = output / "simulations" / "simulation_summaries.tsv"
    simulations = read_tracts(simulation_path) if simulation_path.exists() else pd.DataFrame()
    sensitivity_path = output / "models" / "sensitivity.tsv"
    sensitivity = read_tracts(sensitivity_path) if sensitivity_path.exists() else pd.DataFrame()
    calibration_path = output / "simulations" / "parameter_recovery.tsv"
    calibration = read_tracts(calibration_path) if calibration_path.exists() else pd.DataFrame()
    confusion_path = output / "simulations" / "model_confusion.tsv"
    confusion = read_tracts(confusion_path) if confusion_path.exists() else pd.DataFrame()
    uncertainty_path = output / "models" / "bootstrap_intervals.tsv"
    uncertainty = read_tracts(uncertainty_path) if uncertainty_path.exists() else pd.DataFrame()
    warnings = collect_warnings(tracts, model_table, bootstrap, calibration)
    _, html_path = write_report(
        output / "report",
        config=config,
        tracts=tracts,
        summaries=summaries,
        model_table=model_table,
        fits=fits,
        warnings=warnings,
        command=" ".join(sys.argv),
        repository=REPOSITORY_ROOT,
        figure=figure,
        simulation_summary=simulations,
        sensitivity_table=sensitivity,
        calibration_table=calibration,
        confusion_table=confusion,
        uncertainty_table=uncertainty,
    )
    print(html_path)
    return 0


def _output_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and ".part" not in path.name:
            fingerprint = file_fingerprint(path)
            fingerprint["relative_path"] = path.relative_to(output).as_posix()
            rows.append(fingerprint)
    atomic_write_json(output / "output_manifest.json", {"files": rows})


def command_run_all(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    output.mkdir(parents=True, exist_ok=True)
    save_snapshot(config, output / "config.snapshot.yaml")
    atomic_write_text(output / "command.txt", " ".join(sys.argv) + "\n")
    print(
        json.dumps(
            {
                "output": str(output),
                "profile": config["runtime"]["profile"],
                "memory_ceiling_gb": config["runtime"]["max_memory_gb"],
                "threads": config["runtime"]["max_threads"],
                "simulation_units": len(config["models"])
                * int(config["simulation"]["final_replicates_per_model"]),
                "checkpoint_directory": str(output / "checkpoints"),
                "maximum_block_minutes": config["runtime"]["max_block_minutes"],
            },
            indent=2,
        )
    )
    digest = config_hash(config)
    checkpoint = CheckpointStore(
        output / "checkpoints" / "run_all.json",
        "run_all",
        digest,
        config["project"]["random_seed"],
    )
    state = checkpoint.load(resume=args.resume)
    deadline = _deadline(config)

    filtered_path = output / "tracts" / "filtered_tracts.tsv"
    if not (args.resume and checkpoint.unit_valid(state, "tracts")):
        if args.input:
            imported, malformed = import_tracts(
                args.input,
                caller=args.caller,
                column_map=parse_column_map(args.map_column),
                population=args.population,
            )
            atomic_write_json(
                output / "input_manifest.json",
                {
                    "input": file_fingerprint(args.input),
                    "caller": args.caller,
                    "population_override": args.population,
                    "access_authorization": (
                        "user-supplied path; operator is responsible for documented permission"
                    ),
                },
            )
        elif config["runtime"]["profile"] == "smoke":
            imported = simulate_tracts(
                "SMOKE_OBS",
                config["models"]["M5"],
                n_tracts=int(config["simulation"]["observed_tract_count"]),
                generation_time_years=config["project"]["generation_time_years"],
                minimum_length_cm=config["tracts"]["minimum_length_cm"],
                seed=config["project"]["random_seed"],
            )
            malformed = pd.DataFrame()
            atomic_write_json(
                output / "input_manifest.json",
                {
                    "input_type": "synthetic",
                    "model": "M5",
                    "seed": config["project"]["random_seed"],
                    "contains_human_genomic_data": False,
                },
            )
        else:
            raise ValueError("run-all requires --input outside the smoke profile")
        imported = _select(imported, args)
        tract_settings = config["tracts"]
        retained, excluded = filter_tracts(
            imported,
            minimum_length_cm=tract_settings["minimum_length_cm"],
            minimum_confidence=tract_settings["minimum_confidence"],
            minimum_callable_fraction=tract_settings["minimum_callable_fraction"],
            masks=_configured_masks(config),
        )
        if retained.empty:
            raise ValueError("No tracts remain after QC; dating is not estimable")
        write_tracts(imported, output / "tracts" / "imported_tracts.tsv")
        write_tracts(retained, filtered_path)
        write_tracts(excluded, output / "tracts" / "qc_excluded_tracts.tsv")
        if len(malformed):
            write_tracts(malformed, output / "tracts" / "malformed_tracts.tsv")
        checkpoint.mark_completed(state, "tracts", [filtered_path])
    tracts = read_tracts(filtered_path)

    summary_names = (
        "overall",
        "individual",
        "population",
        "chromosome",
        "sample_missingness",
        "caller_batch",
        "source_class",
    )
    summary_paths = [output / "qc" / f"{name}.tsv" for name in summary_names]
    if not (args.resume and checkpoint.unit_valid(state, "qc")):
        summaries = summarize_tracts(tracts)
        _write_frames(summaries, output / "qc")
        write_tracts(overlapping_pairs(tracts), output / "qc" / "overlapping_tracts.tsv")
        checkpoint.mark_completed(state, "qc", summary_paths)
    summaries = {
        name: read_tracts(output / "qc" / f"{name}.tsv")
        for name in summary_names
    }

    comparison_path = output / "models" / "model_comparison.tsv"
    fits_path = output / "models" / "fits.json"
    if not (args.resume and checkpoint.unit_valid(state, "models")):
        table, fits = compare_models(
            tracts["length_cm"],
            single_bounds=tuple(config["dating"]["single_pulse"]["bounds_generations"]),
            two_minimum_separation_generations=config["dating"]["two_pulse"]["minimum_separation_generations"],
            **_fit_kwargs(config),
        )
        write_tracts(table, comparison_path)
        atomic_write_json(fits_path, fits)
        checkpoint.mark_completed(state, "models", [comparison_path, fits_path])
    model_table = read_tracts(comparison_path)
    with fits_path.open("r", encoding="utf-8") as handle:
        fits = json.load(handle)

    bootstrap_path = output / "models" / "bootstrap_single.tsv"
    sample_bootstrap_path = output / "models" / "bootstrap_sample.tsv"
    two_bootstrap_path = output / "models" / "bootstrap_two_pulse.tsv"
    continuous_bootstrap_path = output / "models" / "bootstrap_continuous.tsv"
    uncertainty_path = output / "models" / "bootstrap_intervals.tsv"
    if not (args.resume and checkpoint.unit_valid(state, "bootstrap")):
        bootstrap = bootstrap_fit(
            tracts,
            fit_single_pulse,
            replicates=int(config["dating"]["bootstrap_replicates"]),
            seed=int(config["project"]["random_seed"]) + 700,
            group_column="chromosome",
            fitter_kwargs={
                **_fit_kwargs(config),
                "bounds_generations": tuple(config["dating"]["single_pulse"]["bounds_generations"]),
            },
        )
        write_tracts(bootstrap, bootstrap_path)
        sample_bootstrap = bootstrap_fit(
            tracts,
            fit_single_pulse,
            replicates=int(config["dating"]["bootstrap_replicates"]),
            seed=int(config["project"]["random_seed"]) + 701,
            group_column="sample_id",
            fitter_kwargs={
                **_fit_kwargs(config),
                "bounds_generations": tuple(config["dating"]["single_pulse"]["bounds_generations"]),
            },
        )
        write_tracts(sample_bootstrap, sample_bootstrap_path)
        two_bootstrap = bootstrap_fit(
            tracts,
            fit_two_pulse,
            replicates=int(config["dating"]["bootstrap_replicates"]),
            seed=int(config["project"]["random_seed"]) + 702,
            group_column="chromosome",
            fitter_kwargs={
                **_fit_kwargs(config),
                "minimum_separation_generations": config["dating"]["two_pulse"][
                    "minimum_separation_generations"
                ],
            },
        )
        write_tracts(two_bootstrap, two_bootstrap_path)
        continuous_bootstrap = bootstrap_fit(
            tracts,
            fit_continuous_flow,
            replicates=int(config["dating"]["bootstrap_replicates"]),
            seed=int(config["project"]["random_seed"]) + 703,
            group_column="chromosome",
            fitter_kwargs=_fit_kwargs(config),
        )
        write_tracts(continuous_bootstrap, continuous_bootstrap_path)
        interval_rows = []
        for model_id, table, columns in (
            ("single_pulse", bootstrap, ["kya", "generations"]),
            ("two_pulse", two_bootstrap, ["older_kya", "younger_kya", "weight_older"]),
            (
                "continuous_flow",
                continuous_bootstrap,
                ["older_kya", "younger_kya", "duration_generations"],
            ),
        ):
            for parameter in columns:
                bounds = interval(table, parameter)
                interval_rows.append(
                    {"model_id": model_id, "parameter": parameter, **bounds}
                )
        uncertainty = pd.DataFrame(interval_rows)
        write_tracts(uncertainty, uncertainty_path)
        checkpoint.mark_completed(
            state,
            "bootstrap",
            [
                bootstrap_path,
                sample_bootstrap_path,
                two_bootstrap_path,
                continuous_bootstrap_path,
                uncertainty_path,
            ],
        )
    bootstrap = read_tracts(bootstrap_path)
    uncertainty = read_tracts(uncertainty_path)

    sensitivity_path = output / "models" / "sensitivity.tsv"
    if not (args.resume and checkpoint.unit_valid(state, "sensitivity")):
        sensitivity = run_sensitivity(tracts, config)
        write_tracts(sensitivity, sensitivity_path)
        checkpoint.mark_completed(state, "sensitivity", [sensitivity_path])
    sensitivity = read_tracts(sensitivity_path)

    replicate_count = int(config["simulation"]["final_replicates_per_model"])
    simulation_summary, paused = simulate_model_set(
        config,
        output,
        replicates=replicate_count,
        resume=args.resume,
        deadline=deadline,
        config_digest=digest,
    )
    if paused:
        state["state"] = "paused"
        checkpoint.save(state)
        print(f"Paused safely before deadline. Resume with the same command. Output: {output}")
        return RESUMABLE_EXIT_CODE

    recovery_path = output / "simulations" / "parameter_recovery.tsv"
    confusion_path = output / "simulations" / "model_confusion.tsv"
    if not (args.resume and checkpoint.unit_valid(state, "calibration")):
        recovery, confusion = calibrate_simulations(output, config)
        write_tracts(recovery, recovery_path)
        write_tracts(confusion, confusion_path)
        checkpoint.mark_completed(state, "calibration", [recovery_path, confusion_path])
    recovery = read_tracts(recovery_path)
    confusion = read_tracts(confusion_path)
    valid_recovery = recovery.loc[recovery["status"] == "ok"].copy()
    if not valid_recovery.empty:
        if valid_recovery["correct"].dtype != bool:
            valid_recovery["correct"] = (
                valid_recovery["correct"]
                .astype(str)
                .str.lower()
                .map({"true": True, "false": False})
            )
        accuracies = valid_recovery.groupby("expected_family")["correct"].mean()
        errors = valid_recovery.groupby("expected_family")["relative_parameter_error"].median()
        model_table["simulation_classification_accuracy"] = model_table["model_id"].map(accuracies)
        model_table["parameter_recovery_quality"] = model_table["model_id"].map(
            lambda model_id: (
                "good"
                if model_id in errors
                and pd.notna(errors[model_id])
                and errors[model_id] <= 0.20
                and accuracies.get(model_id, 0) >= 0.70
                else "poor"
                if model_id in errors and pd.notna(errors[model_id])
                else "not_applicable"
            )
        )
        write_tracts(model_table, comparison_path)
        checkpoint.mark_completed(state, "models", [comparison_path, fits_path])

    if deadline.should_stop(estimated_next_seconds=5):
        state["state"] = "paused"
        checkpoint.save(state)
        print(f"Paused safely before report stage. Resume with the same command. Output: {output}")
        return RESUMABLE_EXIT_CODE

    warnings = collect_warnings(tracts, model_table, bootstrap, recovery)
    figure = make_qc_figure(tracts, model_table, output / "report" / "qc_model_comparison.png")
    markdown_path, html_path = write_report(
        output / "report",
        config=config,
        tracts=tracts,
        summaries=summaries,
        model_table=model_table,
        fits=fits,
        warnings=warnings,
        command=" ".join(sys.argv),
        repository=REPOSITORY_ROOT,
        figure=figure,
        simulation_summary=simulation_summary,
        sensitivity_table=sensitivity,
        calibration_table=recovery,
        confusion_table=confusion,
        uncertainty_table=uncertainty,
    )
    checkpoint.mark_completed(state, "report", [markdown_path, html_path])
    state["state"] = "complete"
    checkpoint.save(state)
    _output_manifest(output)
    print(html_path)
    return 0


def command_status(args: argparse.Namespace) -> int:
    config = _resolved(args)
    output = _run_dir(args, config)
    paths = {
        "run": output / "checkpoints" / "run_all.json",
        "simulations": output / "checkpoints" / "simulations.json",
        "report": output / "report" / "report.html",
    }
    result = {"output": str(output), "exists": output.exists()}
    for name, path in paths.items():
        if path.suffix == ".json" and path.exists():
            result[name] = json.loads(path.read_text(encoding="utf-8"))
        else:
            result[name] = {"path": str(path), "exists": path.exists()}
    print(json.dumps(result, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="archaic-admixture-dating",
        description="Restartable, uncertainty-aware Papuan Denisovan tract dating",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = {
        "init-project": command_init,
        "inspect-data": command_inspect,
        "estimate-storage": command_estimate_storage,
        "download-status": command_download_status,
        "verify-downloads": command_verify_downloads,
        "qc": command_qc,
        "fit-single-pulse": command_fit_single,
        "fit-two-pulse": command_fit_two,
        "fit-continuous": command_fit_continuous,
        "compare-models": command_compare,
        "sensitivity": command_sensitivity,
        "report": command_report,
        "status": command_status,
    }
    for name, function in commands.items():
        sub = subparsers.add_parser(name)
        _common(sub)
        sub.set_defaults(function=function)
        if name in {"inspect-data", "qc", "fit-single-pulse", "fit-two-pulse", "fit-continuous", "compare-models", "sensitivity"}:
            sub.add_argument("--input")

    download_parser = subparsers.add_parser("download")
    _common(download_parser)
    download_parser.add_argument("--dataset", required=True)
    download_parser.set_defaults(function=command_download)

    import_parser = subparsers.add_parser("import-tracts")
    _common(import_parser)
    import_parser.add_argument("--input", required=True)
    import_parser.add_argument("--caller", default="generic", choices=["generic", "ibdmix", "admixfrog", "hmmix", "archaicseeker2"])
    import_parser.add_argument("--population")
    import_parser.add_argument("--map-column", action="append")
    import_parser.set_defaults(function=command_import)

    simulate_parser = subparsers.add_parser("simulate")
    _common(simulate_parser)
    simulate_parser.add_argument("--replicates", type=int)
    simulate_parser.set_defaults(function=command_simulate)

    run_parser = subparsers.add_parser("run-all")
    _common(run_parser)
    run_parser.add_argument("--input")
    run_parser.add_argument("--caller", default="generic", choices=["generic", "ibdmix", "admixfrog", "hmmix", "archaicseeker2"])
    run_parser.add_argument("--population")
    run_parser.add_argument("--map-column", action="append")
    run_parser.set_defaults(function=command_run_all)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.verbose:
        os.environ["ARCHAIC_LOG_LEVEL"] = "DEBUG"
    try:
        return int(args.function(args))
    except (DownloadError, ValueError, RuntimeError, OSError) as error:
        LOG.error("%s", error)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
