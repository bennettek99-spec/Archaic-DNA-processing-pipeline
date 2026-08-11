"""Machine-readable source and download manifest records."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .checkpointing import atomic_write_json, utc_now


@dataclass
class DownloadRecord:
    dataset_id: str
    source_url: str
    destination: str
    expected_size: int | None = None
    bytes_completed: int = 0
    checksum_algorithm: str | None = None
    checksum: str | None = None
    last_successful_update: str | None = None
    retry_count: int = 0
    completion_state: str = "pending"
    access: str = "public"
    etag: str | None = None
    last_modified: str | None = None
    error: str | None = None

    def touch(self) -> None:
        self.last_successful_update = utc_now()


class DownloadManifest:
    def __init__(self, path: str | Path):
        self.path = Path(path)

    def load(self) -> dict[str, Any]:
        if not self.path.exists():
            return {"schema_version": 1, "updated_at": utc_now(), "records": {}}
        with self.path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    def save_record(self, record: DownloadRecord) -> None:
        data = self.load()
        data["records"][record.dataset_id] = asdict(record)
        data["updated_at"] = utc_now()
        atomic_write_json(self.path, data)

    def get(self, dataset_id: str) -> DownloadRecord | None:
        data = self.load()
        value = data.get("records", {}).get(dataset_id)
        return DownloadRecord(**value) if value else None
