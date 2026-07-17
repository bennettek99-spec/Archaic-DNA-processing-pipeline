#!/usr/bin/env python3
"""Check repository-facing documentation, version metadata, and data provenance."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

import yaml

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]
LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
PERSONAL_PATH_RE = re.compile(
    r"(?:C:[\\/]" + r"Users[\\/]" + r"benne"
    + r"|/Users/" + r"benne"
    + r"|/home/" + r"benne"
    + r"|OneDrive[\\/]Desktop[\\/]" + r"Claude Work)",
    re.IGNORECASE,
)


def tracked_files(pattern: str) -> list[Path]:
    output = subprocess.check_output(
        ["git", "ls-files", pattern], cwd=ROOT, text=True, encoding="utf-8"
    )
    return [ROOT / line for line in output.splitlines() if line]


def all_tracked_files() -> list[Path]:
    output = subprocess.check_output(
        ["git", "ls-files"], cwd=ROOT, text=True, encoding="utf-8"
    )
    return [ROOT / line for line in output.splitlines() if line]


def strip_code(markdown: str) -> str:
    markdown = re.sub(r"```.*?```", "", markdown, flags=re.DOTALL)
    return re.sub(r"`[^`]*`", "", markdown)


def check_links(errors: list[str]) -> None:
    for path in tracked_files("*.md"):
        text = strip_code(path.read_text(encoding="utf-8"))
        for match in LINK_RE.finditer(text):
            raw = match.group(1).strip()
            if raw.startswith("<") and raw.endswith(">"):
                raw = raw[1:-1]
            target = unquote(raw.split("#", 1)[0].strip())
            if not target or re.match(r"^(?:https?://|mailto:)", target):
                continue
            candidate = (path.parent / target).resolve()
            if not candidate.exists():
                rel = path.relative_to(ROOT)
                errors.append(f"broken local link: {rel} -> {target}")


def check_versions(errors: list[str]) -> None:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        project_version = tomllib.load(handle)["project"]["version"]

    init_text = (ROOT / "archaic" / "__init__.py").read_text(encoding="utf-8")
    init_match = re.search(r'__version__\s*=\s*"([^"]+)"', init_text)

    citation_text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    citation_match = re.search(r"^version:\s*([^\s]+)", citation_text, re.MULTILINE)

    citation = yaml.safe_load(citation_text)
    for field in ("cff-version", "title", "message", "type", "authors", "license"):
        if not citation.get(field):
            errors.append(f"CITATION.cff is missing required field {field!r}")

    zenodo_version = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))["version"]

    versions = {
        "pyproject.toml": project_version,
        "archaic/__init__.py": init_match.group(1) if init_match else None,
        "CITATION.cff": citation_match.group(1) if citation_match else None,
        ".zenodo.json": zenodo_version,
    }
    for source, version in versions.items():
        if version != project_version:
            errors.append(
                f"version mismatch: {source} has {version!r}, expected {project_version!r}"
            )

    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    if f"## [{project_version}]" not in changelog:
        errors.append(f"CHANGELOG.md has no release heading for {project_version}")

    citation_date = str(citation.get("date-released", ""))
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    if citation_date != zenodo.get("publication_date"):
        errors.append(
            "release-date mismatch: CITATION.cff and .zenodo.json must agree"
        )


def check_paths(errors: list[str]) -> None:
    binary_suffixes = {".bam", ".png", ".jpg", ".jpeg", ".gif", ".pdf"}
    for path in all_tracked_files():
        if path.suffix.lower() in binary_suffixes:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        match = PERSONAL_PATH_RE.search(text)
        if match:
            errors.append(
                f"personal machine path in {path.relative_to(ROOT)}: {match.group(0)}"
            )


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def check_tracked_bams(errors: list[str]) -> None:
    data_doc = (ROOT / "DATA.md").read_text(encoding="utf-8").lower()
    bams = tracked_files("*.bam")
    if not bams:
        errors.append("expected documented PRJEB10597 BAM exception is missing")
        return
    for path in bams:
        relative = path.relative_to(ROOT).as_posix()
        md5 = digest(path, "md5")
        sha256 = digest(path, "sha256")
        for value, label in ((path.name, "filename"), (md5, "MD5"), (sha256, "SHA-256")):
            if value.lower() not in data_doc:
                errors.append(f"DATA.md is missing {label} for tracked BAM {relative}")


def main() -> int:
    errors: list[str] = []
    check_links(errors)
    check_versions(errors)
    check_paths(errors)
    check_tracked_bams(errors)

    if errors:
        print("Repository documentation checks failed:")
        for error in errors:
            print(f"- {error}")
        return 1

    print("Repository documentation checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
