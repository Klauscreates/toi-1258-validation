"""Offline audit entry point. Does not download or modify source evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .analysis import analyze
from .common import verify_manifest, write_csv


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Reproducible evidence checks for the TOI-1258.01 statistical-validation analysis."
    )
    parser.add_argument(
        "--root", type=Path, default=Path.cwd(), help="Repository root; defaults to current directory."
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results"), help="Directory for regenerated audit outputs."
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    out = args.output.resolve()
    if (
        out == root
        or out.is_relative_to(root / "data")
        or out.is_relative_to(root / "src")
        or out.is_relative_to(root / "research")
    ):
        parser.error("Output must not overwrite repository inputs or source")
    try:
        verified = verify_manifest(root)
        summary, rows = analyze(root / "data")
        summary["verified_provenance_files"] = verified
        out.mkdir(parents=True, exist_ok=True)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        write_csv(out / "audit.csv", rows)
        (out / "REPORT.md").write_text(
            "# TOI-1258.01 | Transit Validation\n\n"
            + summary["interpretation"]
            + "\n\nThis is an offline reanalysis of bundled historical tables, not a fresh raw-data reduction.\n\n"
            + "```json\n"
            + json.dumps(summary, indent=2, allow_nan=False)
            + "\n```\n",
            encoding="utf-8",
        )
        print(json.dumps(summary, indent=2, allow_nan=False))
        return 0
    except (OSError, ValueError, KeyError) as exc:
        print(f"Audit failed: {exc}", file=sys.stderr)
        return 1
