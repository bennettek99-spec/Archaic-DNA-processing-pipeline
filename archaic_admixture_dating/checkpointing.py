"""Atomic checkpoints, fingerprints, and voluntary deadline handling."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__

RESUMABLE_EXIT_CODE = 75


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_text(path: str | Path, text: str, encoding: str = "utf-8") -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent))
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_write_json(path: str | Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: str | Path, hash_limit_bytes: int = 128 * 1024 * 1024) -> dict[str, Any]:
    target = Path(path)
    stat = target.stat()
    value: dict[str, Any] = {
        "path": str(target.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if stat.st_size <= hash_limit_bytes:
        value["sha256"] = sha256_file(target)
    return value


class Deadline:
    """Tracks a block deadline and stops before another atomic unit begins."""

    def __init__(self, max_minutes: float, stop_before_minutes: float = 5.0):
        if max_minutes <= stop_before_minutes:
            raise ValueError("max_minutes must exceed stop_before_minutes")
        self.started = time.monotonic()
        self.deadline = self.started + (max_minutes * 60.0)
        self.buffer_seconds = stop_before_minutes * 60.0

    @property
    def remaining_seconds(self) -> float:
        return max(0.0, self.deadline - time.monotonic())

    def should_stop(self, estimated_next_seconds: float = 0.0) -> bool:
        return self.remaining_seconds <= self.buffer_seconds + max(0.0, estimated_next_seconds)


class CheckpointStore:
    """One atomic JSON checkpoint per restartable task."""

    def __init__(
        self,
        path: str | Path,
        task: str,
        config_digest: str,
        seed: int,
        input_fingerprints: list[dict[str, Any]] | None = None,
    ):
        self.path = Path(path)
        self.task = task
        self.config_digest = config_digest
        self.seed = int(seed)
        self.input_fingerprints = input_fingerprints or []

    def _new_state(self) -> dict[str, Any]:
        return {
            "task": self.task,
            "config_hash": self.config_digest,
            "software_version": __version__,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "random_seed": self.seed,
            "input_fingerprints": self.input_fingerprints,
            "completed_units": {},
            "failed_units": {},
            "created_at": utc_now(),
            "updated_at": utc_now(),
            "state": "running",
        }

    def load(self, resume: bool = True) -> dict[str, Any]:
        if not resume or not self.path.exists():
            return self._new_state()
        with self.path.open("r", encoding="utf-8") as handle:
            state = json.load(handle)
        if state.get("task") != self.task:
            raise ValueError(f"Checkpoint task mismatch at {self.path}")
        if state.get("config_hash") != self.config_digest:
            raise ValueError(
                f"Checkpoint configuration mismatch at {self.path}; use a new run directory"
            )
        if state.get("input_fingerprints", []) != self.input_fingerprints:
            raise ValueError(f"Checkpoint inputs changed at {self.path}")
        return state

    def save(self, state: dict[str, Any]) -> None:
        state["updated_at"] = utc_now()
        atomic_write_json(self.path, state)

    def mark_completed(self, state: dict[str, Any], unit: str, outputs: list[str | Path]) -> None:
        fingerprints = [file_fingerprint(path) for path in outputs]
        state["completed_units"][str(unit)] = {
            "completed_at": utc_now(),
            "outputs": fingerprints,
        }
        state["failed_units"].pop(str(unit), None)
        self.save(state)

    def mark_failed(self, state: dict[str, Any], unit: str, error: Exception | str) -> None:
        state["failed_units"][str(unit)] = {
            "failed_at": utc_now(),
            "error": str(error),
        }
        self.save(state)

    @staticmethod
    def unit_valid(state: dict[str, Any], unit: str) -> bool:
        record = state.get("completed_units", {}).get(str(unit))
        if not record:
            return False
        for expected in record.get("outputs", []):
            path = Path(expected["path"])
            if not path.exists():
                return False
            current = file_fingerprint(path)
            for key in ("size", "sha256"):
                if key in expected and current.get(key) != expected[key]:
                    return False
        return True
