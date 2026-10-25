"""Strict file I/O and provenance checks for the bundled evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if not reader.fieldnames or not rows:
            raise ValueError(f"Empty table: {path.name}")
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f"Malformed table: {path.name}")
        return rows


def number(value: str | float) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Non-finite measurement")
    return result


def ratio(numerator: str | float, denominator: str | float) -> float:
    denominator = number(denominator)
    if denominator <= 0:
        raise ValueError("Uncertainty or scale must be positive")
    return number(numerator) / denominator


def verify_manifest(root: Path) -> int:
    manifest = json.loads((root / "data/manifest.json").read_text(encoding="utf-8"))
    seen = set()
    if not manifest.get("files"):
        raise ValueError("Empty provenance manifest")
    for entry in manifest["files"]:
        name = entry["path"]
        if name in seen:
            raise ValueError(f"Duplicate manifest entry: {name}")
        seen.add(name)
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError("Manifest path escapes repository")
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Evidence checksum mismatch: {name}")
    return len(seen)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty reconciliation table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
