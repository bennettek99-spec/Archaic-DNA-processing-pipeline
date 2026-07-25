from __future__ import annotations

import hashlib
import io
import json
import urllib.error

import pytest

from archaic_admixture_dating import downloads
from archaic_admixture_dating.downloads import (
    AccessRestrictedError,
    DownloadError,
    download,
    ensure_storage,
)
from archaic_admixture_dating.manifests import DownloadManifest, DownloadRecord


class FakeResponse:
    def __init__(self, chunks, status=200):
        self.chunks = iter(chunks)
        self.status = status

    def read(self, _size):
        value = next(self.chunks, b"")
        if isinstance(value, Exception):
            raise value
        return value

    def getcode(self):
        return self.status

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False


def test_interrupted_download_resumes_and_preserves_partial(tmp_path):
    payload = b"abcdef"
    destination = tmp_path / "payload.bin"
    manifest = DownloadManifest(tmp_path / "manifest.json")
    record = DownloadRecord(
        dataset_id="fixture",
        source_url="https://example.invalid/payload",
        destination=str(destination),
        expected_size=len(payload),
        checksum_algorithm="sha256",
        checksum=hashlib.sha256(payload).hexdigest(),
    )

    def interrupted(_request, timeout=60):
        return FakeResponse([b"abc", urllib.error.URLError("connection lost")], status=200)

    with pytest.raises(urllib.error.URLError):
        download(record, manifest, retries=0, retry_backoff_seconds=0, opener=interrupted)
    partial = destination.with_name(destination.name + ".part")
    assert partial.read_bytes() == b"abc"

    def resumed(request, timeout=60):
        assert request.headers["Range"] == "bytes=3-"
        return FakeResponse([b"def", b""], status=206)

    completed = download(record, manifest, retries=0, opener=resumed)
    assert completed.completion_state == "complete"
    assert destination.read_bytes() == payload
    assert not partial.exists()


def test_checksum_mismatch_is_detected_and_partial_is_retained(tmp_path):
    destination = tmp_path / "payload.bin"
    manifest = DownloadManifest(tmp_path / "manifest.json")
    record = DownloadRecord(
        dataset_id="bad-checksum",
        source_url="https://example.invalid/payload",
        destination=str(destination),
        expected_size=3,
        checksum_algorithm="sha256",
        checksum="0" * 64,
    )

    with pytest.raises(DownloadError, match="Checksum mismatch"):
        download(record, manifest, retries=0, opener=lambda *_args, **_kwargs: FakeResponse([b"abc", b""]))
    assert destination.with_name(destination.name + ".part").read_bytes() == b"abc"


def test_insufficient_disk_space_blocks_download(monkeypatch, tmp_path):
    monkeypatch.setattr(downloads, "free_bytes", lambda _path: 10)
    with pytest.raises(DownloadError, match="Insufficient disk space"):
        ensure_storage(tmp_path, expected_size=1000)


def test_manifest_updates_are_valid_json(tmp_path):
    manifest = DownloadManifest(tmp_path / "manifest.json")
    record = DownloadRecord("fixture", "https://example.invalid", str(tmp_path / "x"), expected_size=1)
    manifest.save_record(record)
    data = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert data["records"]["fixture"]["expected_size"] == 1
    assert not list(tmp_path.glob("*.tmp"))


def test_controlled_access_cannot_be_forced_through_downloader(tmp_path):
    manifest = DownloadManifest(tmp_path / "manifest.json")
    record = DownloadRecord(
        "controlled",
        "https://example.invalid",
        str(tmp_path / "x"),
        expected_size=1,
        access="controlled",
    )
    with pytest.raises(AccessRestrictedError):
        download(record, manifest, opener=lambda *_args, **_kwargs: io.BytesIO(b"x"))
