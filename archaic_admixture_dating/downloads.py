"""Resumable, checksum-verified, disk-guarded HTTP downloads."""

from __future__ import annotations

import hashlib
import os
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

from .checkpointing import Deadline, sha256_file
from .manifests import DownloadManifest, DownloadRecord


class DownloadError(RuntimeError):
    pass


class AccessRestrictedError(DownloadError):
    pass


def free_bytes(path: str | Path) -> int:
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return int(shutil.disk_usage(target).free)


def required_storage_bytes(expected_size: int, existing_partial: int = 0, headroom: float = 1.10) -> int:
    remaining = max(0, int(expected_size) - int(existing_partial))
    return int(remaining * headroom)


def ensure_storage(destination_dir: str | Path, expected_size: int, existing_partial: int = 0) -> None:
    required = required_storage_bytes(expected_size, existing_partial)
    available = free_bytes(destination_dir)
    if available < required:
        raise DownloadError(
            f"Insufficient disk space: need {required:,} bytes including headroom; "
            f"{available:,} bytes available"
        )


def _checksum(path: Path, algorithm: str) -> str:
    if algorithm.lower() == "sha256":
        return sha256_file(path)
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def probe_http(url: str, timeout: float = 30.0) -> dict[str, Any]:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            headers = response.headers
            length = headers.get("Content-Length")
            return {
                "url": response.geturl(),
                "expected_size": int(length) if length else None,
                "accept_ranges": "bytes" in (headers.get("Accept-Ranges") or "").lower(),
                "etag": headers.get("ETag"),
                "last_modified": headers.get("Last-Modified"),
            }
    except urllib.error.HTTPError as error:
        if error.code in {401, 403}:
            raise AccessRestrictedError(
                f"Source requires authorization ({error.code}); access controls will not be bypassed"
            ) from error
        raise DownloadError(f"Unable to inspect source: HTTP {error.code}") from error


def download(
    record: DownloadRecord,
    manifest: DownloadManifest,
    *,
    resume: bool = True,
    dry_run: bool = False,
    chunk_size_mb: int = 64,
    retries: int = 8,
    retry_backoff_seconds: float = 30.0,
    bandwidth_limit_mbps: float | None = None,
    deadline: Deadline | None = None,
    opener=urllib.request.urlopen,
) -> DownloadRecord:
    if record.access != "public":
        raise AccessRestrictedError(
            f"Dataset {record.dataset_id!r} is marked {record.access!r}; automated download refused"
        )
    destination = Path(record.destination)
    partial = destination.with_name(destination.name + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)
    existing = partial.stat().st_size if resume and partial.exists() else 0
    if partial.exists() and not resume:
        raise DownloadError(f"Partial file exists at {partial}; use --resume to preserve it")
    if record.expected_size is None:
        probe = probe_http(record.source_url)
        record.expected_size = probe["expected_size"]
        record.etag = probe["etag"]
        record.last_modified = probe["last_modified"]
        if existing and not probe["accept_ranges"]:
            raise DownloadError("Source does not advertise byte-range support; manual/component download required")
    if record.expected_size is None:
        raise DownloadError("Source size is unknown; refusing an unbounded download")
    ensure_storage(destination.parent, record.expected_size, existing)
    record.bytes_completed = existing
    record.completion_state = "dry-run" if dry_run else "downloading"
    record.touch()
    manifest.save_record(record)
    if dry_run:
        return record

    chunk_size = max(1, int(chunk_size_mb)) * 1024 * 1024
    mode = "ab" if existing else "wb"
    for attempt in range(retries + 1):
        if deadline and deadline.should_stop():
            record.completion_state = "paused"
            record.touch()
            manifest.save_record(record)
            return record
        headers = {"User-Agent": "archaic-admixture-dating/0.1"}
        if existing:
            headers["Range"] = f"bytes={existing}-"
        request = urllib.request.Request(record.source_url, headers=headers)
        try:
            with opener(request, timeout=60) as response, partial.open(mode) as handle:
                status = getattr(response, "status", response.getcode())
                if existing and status != 206:
                    raise DownloadError(
                        "Server did not honor the range request; partial file preserved"
                    )
                while True:
                    if deadline and deadline.should_stop():
                        record.completion_state = "paused"
                        record.bytes_completed = partial.stat().st_size
                        record.touch()
                        manifest.save_record(record)
                        return record
                    read_started = time.monotonic()
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    handle.write(chunk)
                    handle.flush()
                    os.fsync(handle.fileno())
                    record.bytes_completed = handle.tell()
                    record.touch()
                    manifest.save_record(record)
                    if bandwidth_limit_mbps:
                        target_seconds = len(chunk) * 8.0 / (float(bandwidth_limit_mbps) * 1_000_000)
                        delay = target_seconds - (time.monotonic() - read_started)
                        if delay > 0:
                            time.sleep(min(delay, 55.0))
            break
        except (OSError, urllib.error.URLError, DownloadError) as error:
            record.retry_count += 1
            record.error = str(error)
            record.completion_state = "retrying" if attempt < retries else "failed"
            record.bytes_completed = partial.stat().st_size if partial.exists() else 0
            record.touch()
            manifest.save_record(record)
            if attempt >= retries or isinstance(error, DownloadError):
                raise
            time.sleep(min(55.0, retry_backoff_seconds * (attempt + 1)))
            existing = record.bytes_completed
            mode = "ab" if existing else "wb"

    size = partial.stat().st_size
    if size != record.expected_size:
        record.completion_state = "failed"
        record.error = f"Size mismatch: expected {record.expected_size}, received {size}"
        manifest.save_record(record)
        raise DownloadError(record.error)
    if record.checksum and record.checksum_algorithm:
        observed = _checksum(partial, record.checksum_algorithm)
        if observed.lower() != record.checksum.lower():
            record.completion_state = "failed"
            record.error = "Checksum mismatch; partial file preserved"
            manifest.save_record(record)
            raise DownloadError(record.error)
    os.replace(partial, destination)
    record.bytes_completed = destination.stat().st_size
    record.completion_state = "complete"
    record.error = None
    record.touch()
    manifest.save_record(record)
    return record


def verify_record(record: DownloadRecord) -> bool:
    path = Path(record.destination)
    if not path.exists() or record.expected_size is None:
        return False
    if path.stat().st_size != record.expected_size:
        return False
    if record.checksum and record.checksum_algorithm:
        return _checksum(path, record.checksum_algorithm).lower() == record.checksum.lower()
    return True
