#!/usr/bin/env python
"""Create a new mission scaffold for the Astro Candidate Pipeline.

Example:
    python scripts/astro_mission_pipeline_scaffold.py ^
      --mission sirius_followup ^
      --target IRPCA-002 ^
      --primary "Sirius A" ^
      --ra 101.28715533 ^
      --dec -16.71611586 ^
      --sep 1.185 ^
      --pa 129.5
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path
from typing import Any


WORKSPACE = Path(r".")
PIPELINE = WORKSPACE / "pipelines" / "astro_candidate_pipeline"
TEMPLATE = PIPELINE / "TARGET_MANIFEST_TEMPLATE.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mission", required=True, help="Mission/output slug, e.g. alpha_cen_candidate")
    parser.add_argument("--target", required=True, help="Candidate label")
    parser.add_argument("--primary", required=True, help="Primary object name")
    parser.add_argument("--ra", type=float, required=True, help="Primary RA in decimal degrees")
    parser.add_argument("--dec", type=float, required=True, help="Primary Dec in decimal degrees")
    parser.add_argument("--sep", type=float, default=None, help="Candidate separation in arcsec")
    parser.add_argument("--pa", type=float, default=None, help="Candidate PA in degrees east of north")
    parser.add_argument("--inner", type=float, default=None, help="Search annulus inner radius in arcsec")
    parser.add_argument("--outer", type=float, default=None, help="Search annulus outer radius in arcsec")
    parser.add_argument("--hypothesis", default="unclassified high-contrast candidate")
    parser.add_argument("--motivation", default="")
    args = parser.parse_args()

    manifest = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    today = date.today().isoformat()
    package_name = f"{safe_slug(args.mission).upper()}_CANDIDATE_PAPER_PACKAGE_{today}"
    out_dir = WORKSPACE / "outputs" / package_name

    east = north = target_ra = target_dec = None
    if args.sep is not None and args.pa is not None:
        east, north = offset_from_sep_pa(args.sep, args.pa)
        target_ra, target_dec = absolute_from_offset(args.ra, args.dec, east, north)

    manifest["mission_name"] = args.mission
    manifest["target_label"] = args.target
    manifest["primary"]["name"] = args.primary
    manifest["primary"]["ra_deg"] = args.ra
    manifest["primary"]["dec_deg"] = args.dec
    manifest["search_geometry"]["separation_arcsec"] = args.sep
    manifest["search_geometry"]["pa_deg"] = args.pa
    manifest["search_geometry"]["annulus_inner_arcsec"] = args.inner
    manifest["search_geometry"]["annulus_outer_arcsec"] = args.outer
    manifest["search_geometry"]["east_arcsec"] = east
    manifest["search_geometry"]["north_arcsec"] = north
    manifest["search_geometry"]["target_ra_deg_estimate"] = target_ra
    manifest["search_geometry"]["target_dec_deg_estimate"] = target_dec
    manifest["hypothesis"]["working_interpretation"] = args.hypothesis
    manifest["hypothesis"]["motivation"] = args.motivation
    manifest["output"]["base_dir"] = str(out_dir)
    manifest["output"]["date"] = today

    create_package_dirs(out_dir)
    manifest_path = out_dir / "target_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    readme_path = out_dir / "README.md"
    readme_path.write_text(readme_text(args, manifest_path, out_dir), encoding="utf-8")

    print(json.dumps({
        "mission": args.mission,
        "target": args.target,
        "output_dir": str(out_dir),
        "manifest": str(manifest_path),
        "east_arcsec": east,
        "north_arcsec": north,
        "target_ra_deg_estimate": target_ra,
        "target_dec_deg_estimate": target_dec,
    }, indent=2))
    return 0


def safe_slug(value: str) -> str:
    keep = []
    for char in value.strip().lower():
        if char.isalnum():
            keep.append(char)
        elif char in {"-", "_", " ", "."}:
            keep.append("_")
    slug = "".join(keep).strip("_")
    return slug or "mission"


def offset_from_sep_pa(sep_arcsec: float, pa_deg: float) -> tuple[float, float]:
    theta = math.radians(pa_deg)
    east = sep_arcsec * math.sin(theta)
    north = sep_arcsec * math.cos(theta)
    return east, north


def absolute_from_offset(ra_deg: float, dec_deg: float, east_arcsec: float, north_arcsec: float) -> tuple[float, float]:
    target_ra = ra_deg + east_arcsec / (3600.0 * math.cos(math.radians(dec_deg)))
    target_dec = dec_deg + north_arcsec / 3600.0
    return target_ra, target_dec


def create_package_dirs(out_dir: Path) -> None:
    for rel in [
        "manuscript",
        "evidence/figures",
        "evidence/tables",
        "evidence/reports",
        "evidence/archive_gates",
        "evidence/scripts",
    ]:
        (out_dir / rel).mkdir(parents=True, exist_ok=True)


def readme_text(args: argparse.Namespace, manifest_path: Path, out_dir: Path) -> str:
    return f"""# {args.mission} Candidate Pipeline Package

Created by `scripts/astro_mission_pipeline_scaffold.py`.

Target: `{args.target}`

Primary: `{args.primary}`

Manifest:

`{manifest_path}`

Pipeline:

`{PIPELINE}`

Next steps:

1. Open `target_manifest.json` and fill any missing mission-specific fields.
2. Follow `pipelines/astro_candidate_pipeline/PLAYBOOK.md`.
3. Save every output inside this package:

`{out_dir}`
"""


if __name__ == "__main__":
    raise SystemExit(main())
