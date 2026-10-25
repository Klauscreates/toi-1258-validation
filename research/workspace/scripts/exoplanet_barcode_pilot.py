#!/usr/bin/env python
"""Build a first-pass exoplanet atmosphere-barcode pilot package.

This script intentionally starts with known transiting planets that already
have public atmosphere spectra. That gives us an end-to-end calibration target
before using the same workflow for new TESS/Kepler candidate discovery.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import shutil
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable


WORKSPACE = Path(r".")
TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EXOMAST_SPECTRA = "https://exo.mast.stsci.edu/api/v0.1/spectra"

CURATED_TARGETS = [
    "WASP-39 b",
    "HD 189733 b",
    "HD 209458 b",
    "HAT-P-11 b",
    "WASP-43 b",
    "WASP-96 b",
    "K2-18 b",
    "TRAPPIST-1 b",
]

JWST_CONTEXT_TARGETS = {
    "WASP-39 b": "NASA/JWST public atmosphere plot includes NIRISS, NIRCam, NIRSpec G395H, and NIRSpec PRISM spectral features.",
    "WASP-96 b": "JWST public early-release transmission spectrum target; Exo.MAST curated text file was not present in this run.",
    "K2-18 b": "Public JWST literature target, but Exo.MAST spectra endpoint did not provide a curated file in this run.",
}

MOLECULAR_BANDS = [
    ("Na", 0.52, 0.70),
    ("K", 0.72, 0.82),
    ("H2O", 0.90, 1.50),
    ("H2O", 1.72, 2.16),
    ("CO", 2.30, 2.40),
    ("H2O", 2.40, 3.25),
    ("CO2", 2.65, 2.90),
    ("SO2", 3.96, 4.10),
    ("CO2", 4.16, 4.58),
    ("CO", 4.60, 4.96),
    ("H2O", 5.00, 5.35),
]


@dataclass
class SpectrumPoint:
    wavelength_micron: float
    delta_wavelength_micron: float
    transit_depth: float
    transit_depth_unc: float


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output package directory. Defaults to outputs/EXOPLANET_BARCODE_PILOT_<date>.",
    )
    parser.add_argument(
        "--target",
        default="WASP-39 b",
        help="Primary pilot target for the mission manifest.",
    )
    args = parser.parse_args()

    out_dir = (
        Path(args.output_dir)
        if args.output_dir
        else WORKSPACE / "outputs" / f"EXOPLANET_BARCODE_PILOT_{date.today().isoformat()}"
    )
    dirs = create_dirs(out_dir)

    target_rows = fetch_candidate_parameters(CURATED_TARGETS)
    write_csv(dirs["tables"] / "nasa_exoplanet_archive_target_parameters.csv", target_rows)

    spectra_probe_rows = []
    feature_rows = []
    spectra_by_target: dict[str, list[SpectrumPoint]] = {}

    for target in CURATED_TARGETS:
        filelist = fetch_exomast_filelist(target)
        spectra_probe_rows.append(
            {
                "target": target,
                "exomast_file_count": len(filelist),
                "exomast_files": ";".join(filelist),
                "jwst_context_note": JWST_CONTEXT_TARGETS.get(target, ""),
            }
        )
        for filename in filelist:
            text = fetch_exomast_spectrum_file(target, filename)
            save_name = safe_slug(f"{target}_{filename}")
            raw_path = dirs["raw_spectra"] / save_name
            raw_path.write_text(text, encoding="utf-8")
            points = parse_spectrum_text(text)
            spectra_by_target.setdefault(target, []).extend(points)
            for row in summarize_molecular_bands(target, filename, points):
                feature_rows.append(row)

    write_csv(dirs["tables"] / "exomast_spectra_probe.csv", spectra_probe_rows)
    write_csv(dirs["tables"] / "molecular_band_summary.csv", feature_rows)

    shortlist = rank_targets(target_rows, spectra_probe_rows, spectra_by_target)
    write_csv(dirs["tables"] / "exoplanet_barcode_target_shortlist.csv", shortlist)

    plot_paths = []
    for target, points in spectra_by_target.items():
        if not points:
            continue
        plot_paths.append(plot_spectrum(target, points, dirs["figures"]))

    manifest = build_manifest(args.target, target_rows, spectra_probe_rows, shortlist, out_dir)
    (out_dir / "target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report = render_report(args.target, target_rows, spectra_probe_rows, shortlist, feature_rows, plot_paths, out_dir)
    (dirs["reports"] / "EXOPLANET_BARCODE_PILOT_REPORT.md").write_text(report, encoding="utf-8")

    shutil.copy2(Path(__file__), dirs["scripts"] / Path(__file__).name)

    print(json.dumps({
        "output_dir": str(out_dir),
        "primary_target": args.target,
        "targets_ranked": len(shortlist),
        "spectra_targets": sorted(spectra_by_target),
        "plots": [str(p) for p in plot_paths],
        "report": str(dirs["reports"] / "EXOPLANET_BARCODE_PILOT_REPORT.md"),
    }, indent=2))
    return 0


def create_dirs(out_dir: Path) -> dict[str, Path]:
    paths = {
        "base": out_dir,
        "tables": out_dir / "evidence" / "tables",
        "figures": out_dir / "evidence" / "figures",
        "reports": out_dir / "evidence" / "reports",
        "raw_spectra": out_dir / "evidence" / "raw_spectra",
        "scripts": out_dir / "evidence" / "scripts",
        "archive_gates": out_dir / "evidence" / "archive_gates",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def fetch_candidate_parameters(targets: list[str]) -> list[dict[str, str]]:
    names = ",".join("'" + name.replace("'", "''") + "'" for name in targets)
    query = f"""
SELECT pl_name,hostname,sy_dist,discoverymethod,disc_facility,tran_flag,
       pl_orbper,pl_rade,pl_bmasse,pl_eqt,st_teff,st_rad,st_mass,
       sy_vmag,sy_jmag,sy_kmag,ra,dec
FROM pscomppars
WHERE pl_name IN ({names})
"""
    return tap_csv(query)


def tap_csv(query: str) -> list[dict[str, str]]:
    url = TAP_SYNC + "?" + urllib.parse.urlencode({"query": one_line(query), "format": "csv"})
    text = fetch_text(url, timeout=180)
    return list(csv.DictReader(io.StringIO(text)))


def fetch_exomast_filelist(target: str) -> list[str]:
    url = f"{EXOMAST_SPECTRA}/{urllib.parse.quote(target, safe='')}/filelist/"
    try:
        payload = json.loads(fetch_text(url, timeout=60))
    except Exception:
        return []
    names = payload.get("filenames") or payload.get("filename") or []
    return [str(name) for name in names]


def fetch_exomast_spectrum_file(target: str, filename: str) -> str:
    url = (
        f"{EXOMAST_SPECTRA}/"
        f"{urllib.parse.quote(target, safe='')}/file/"
        f"{urllib.parse.quote(filename, safe='')}/"
    )
    return fetch_text(url, timeout=60)


def fetch_text(url: str, timeout: int = 60, tries: int = 3) -> str:
    last_error: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "codex-exoplanet-barcode-pilot/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except Exception as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def parse_spectrum_text(text: str) -> list[SpectrumPoint]:
    points = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            values = [float(part) for part in parts[:4]]
        except ValueError:
            continue
        points.append(SpectrumPoint(*values))
    return points


def summarize_molecular_bands(target: str, filename: str, points: list[SpectrumPoint]) -> list[dict[str, object]]:
    if not points:
        return []
    depths = [p.transit_depth for p in points]
    global_median = median(depths)
    rows = []
    for label, lo, hi in MOLECULAR_BANDS:
        band = [p for p in points if lo <= p.wavelength_micron <= hi]
        if not band:
            continue
        band_depth = weighted_mean(band)
        band_unc = math.sqrt(sum(p.transit_depth_unc ** 2 for p in band)) / max(len(band), 1)
        contrast = band_depth - global_median
        z = contrast / band_unc if band_unc > 0 else float("nan")
        rows.append(
            {
                "target": target,
                "filename": filename,
                "feature": label,
                "band_min_micron": f"{lo:.3f}",
                "band_max_micron": f"{hi:.3f}",
                "n_points": len(band),
                "band_weighted_depth": f"{band_depth:.8f}",
                "global_median_depth": f"{global_median:.8f}",
                "band_minus_global_median": f"{contrast:.8f}",
                "approx_z_vs_global_median": f"{z:.3f}",
            }
        )
    return rows


def weighted_mean(points: list[SpectrumPoint]) -> float:
    weights = []
    values = []
    for point in points:
        if point.transit_depth_unc > 0:
            weights.append(1.0 / point.transit_depth_unc**2)
            values.append(point.transit_depth)
    if not weights:
        return sum(p.transit_depth for p in points) / len(points)
    return sum(v * w for v, w in zip(values, weights)) / sum(weights)


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def rank_targets(
    target_rows: list[dict[str, str]],
    spectra_probe_rows: list[dict[str, str]],
    spectra_by_target: dict[str, list[SpectrumPoint]],
) -> list[dict[str, object]]:
    probe = {row["target"]: row for row in spectra_probe_rows}
    ranked = []
    for row in target_rows:
        name = row["pl_name"]
        points = spectra_by_target.get(name, [])
        file_count = int(probe.get(name, {}).get("exomast_file_count", 0) or 0)
        score = 0.0
        reasons = []

        if file_count:
            score += 60
            reasons.append("downloadable Exo.MAST transmission spectrum")
        if name in JWST_CONTEXT_TARGETS:
            score += 20
            reasons.append("public JWST atmosphere context")
        if to_float(row.get("pl_rade")) and to_float(row.get("pl_rade")) >= 8:
            score += 10
            reasons.append("large planet radius / high atmospheric signal")
        if to_float(row.get("sy_jmag")) and to_float(row.get("sy_jmag")) <= 8:
            score += 8
            reasons.append("bright near-IR host")
        if to_float(row.get("sy_dist")) and to_float(row.get("sy_dist")) <= 60:
            score += 6
            reasons.append("nearby host")
        if points:
            wavelength_span = max(p.wavelength_micron for p in points) - min(p.wavelength_micron for p in points)
            score += min(15, wavelength_span * 4)
            reasons.append(f"{len(points)} spectral points, {wavelength_span:.2f} micron span")

        ranked.append(
            {
                "rank_score": f"{score:.2f}",
                "target": name,
                "host": row.get("hostname", ""),
                "ra_deg": row.get("ra", ""),
                "dec_deg": row.get("dec", ""),
                "distance_pc": row.get("sy_dist", ""),
                "planet_radius_rearth": row.get("pl_rade", ""),
                "equilibrium_temp_k": row.get("pl_eqt", ""),
                "host_jmag": row.get("sy_jmag", ""),
                "host_vmag": row.get("sy_vmag", ""),
                "exomast_file_count": file_count,
                "spectral_points_downloaded": len(points),
                "selection_reasons": "; ".join(reasons),
            }
        )
    ranked.sort(key=lambda r: float(r["rank_score"]), reverse=True)
    return ranked


def plot_spectrum(target: str, points: list[SpectrumPoint], fig_dir: Path) -> Path:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError("matplotlib is required for spectrum plotting") from exc

    points = sorted(points, key=lambda p: p.wavelength_micron)
    x = [p.wavelength_micron for p in points]
    y = [100.0 * p.transit_depth for p in points]
    yerr = [100.0 * p.transit_depth_unc for p in points]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.errorbar(x, y, yerr=yerr, fmt="o", ms=3, lw=0.7, alpha=0.78, color="#1f5a8a", ecolor="#9ab7cc")
    colors = {
        "Na": "#7aa6ff",
        "K": "#a0a0a0",
        "H2O": "#78c8ff",
        "CO": "#df6b57",
        "CO2": "#d8bf33",
        "SO2": "#6dbb6d",
    }
    used = set()
    for label, lo, hi in MOLECULAR_BANDS:
        if max(x) < lo or min(x) > hi:
            continue
        ax.axvspan(lo, hi, alpha=0.18, color=colors.get(label, "#cccccc"), label=label if label not in used else None)
        used.add(label)
    ax.set_title(f"{target} public transmission spectrum")
    ax.set_xlabel("Wavelength (microns)")
    ax.set_ylabel("Transit depth (%)")
    ax.grid(alpha=0.25)
    if used:
        ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    out_path = fig_dir / f"{safe_slug(target)}_transmission_spectrum.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def build_manifest(
    primary_target: str,
    target_rows: list[dict[str, str]],
    spectra_probe_rows: list[dict[str, str]],
    shortlist: list[dict[str, object]],
    out_dir: Path,
) -> dict[str, object]:
    target_row = next((row for row in target_rows if row.get("pl_name") == primary_target), target_rows[0])
    return {
        "mission_name": "exoplanet_barcode_pilot",
        "date": date.today().isoformat(),
        "primary_target": primary_target,
        "status": "known-planet atmosphere-barcode validation package",
        "claim_boundary": "This package validates the barcode extraction workflow on known transiting exoplanets; it is not a new-planet discovery claim.",
        "target": target_row,
        "archive_sources": {
            "nasa_exoplanet_archive_tap": TAP_SYNC,
            "exomast_spectra_api": EXOMAST_SPECTRA,
            "nasa_wasp39b_jwst_context": "https://science.nasa.gov/asset/webb/exoplanet-wasp-39-b-transmission-spectra/",
        },
        "required_next_gates_for_new_discovery": [
            "TESS/Kepler light-curve search",
            "transit model fit and false-positive vetting",
            "Gaia neighbor/contamination audit",
            "ephemeris prediction",
            "spectroscopy availability or follow-up request",
        ],
        "spectra_probe": spectra_probe_rows,
        "shortlist_top": shortlist[:5],
        "output_dir": str(out_dir),
    }


def render_report(
    primary_target: str,
    target_rows: list[dict[str, str]],
    spectra_probe_rows: list[dict[str, str]],
    shortlist: list[dict[str, object]],
    feature_rows: list[dict[str, object]],
    plot_paths: list[Path],
    out_dir: Path,
) -> str:
    by_target = {row["pl_name"]: row for row in target_rows}
    primary = by_target.get(primary_target, target_rows[0])
    lines = [
        "# Exoplanet Barcode Pilot Report",
        "",
        f"Created: {date.today().isoformat()}",
        "",
        "## Result",
        "",
        f"Primary pilot target: `{primary_target}`.",
        "",
        "This is the correct first step for a real exoplanet discovery program: validate the full atmosphere-barcode extraction workflow on known public data before claiming a new planet. The strongest pilot target is WASP-39 b because it has downloadable public transmission-spectrum data and public JWST spectral-context products.",
        "",
        "Claim boundary: this package is a workflow validation and target-selection result. It does not claim a newly discovered planet yet.",
        "",
        "## Primary Target Parameters",
        "",
        f"- Host: `{primary.get('hostname', '')}`",
        f"- RA/Dec: `{primary.get('ra', '')}`, `{primary.get('dec', '')}` deg",
        f"- Distance: `{primary.get('sy_dist', '')}` pc",
        f"- Planet radius: `{primary.get('pl_rade', '')}` Earth radii",
        f"- Orbital period: `{primary.get('pl_orbper', '')}` days",
        f"- Equilibrium temperature: `{primary.get('pl_eqt', '')}` K",
        f"- Host J/V mag: `{primary.get('sy_jmag', '')}` / `{primary.get('sy_vmag', '')}`",
        "",
        "## Ranked Barcode Targets",
        "",
        "| Rank | Target | Score | Exo.MAST files | Spectral points | Why |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for idx, row in enumerate(shortlist, start=1):
        lines.append(
            f"| {idx} | {row['target']} | {row['rank_score']} | {row['exomast_file_count']} | "
            f"{row['spectral_points_downloaded']} | {row['selection_reasons']} |"
        )

    lines.extend([
        "",
        "## Spectrum Files Found",
        "",
        "| Target | File count | Files | JWST context |",
        "|---|---:|---|---|",
    ])
    for row in spectra_probe_rows:
        lines.append(
            f"| {row['target']} | {row['exomast_file_count']} | `{row['exomast_files']}` | {row['jwst_context_note']} |"
        )

    lines.extend([
        "",
        "## Molecular Band Summary",
        "",
        "The band table is a first-pass diagnostic, not a full atmospheric retrieval. It measures weighted transit depth inside common feature windows against each spectrum's global median depth.",
        "",
        f"Full table: `{out_dir / 'evidence' / 'tables' / 'molecular_band_summary.csv'}`",
        "",
        "Top absolute band deviations:",
        "",
        "| Target | Feature | Band (micron) | N | Approx z |",
        "|---|---|---|---:|---:|",
    ])
    top_features = sorted(
        feature_rows,
        key=lambda r: abs(float(r["approx_z_vs_global_median"])) if is_number(r["approx_z_vs_global_median"]) else -1,
        reverse=True,
    )[:12]
    for row in top_features:
        lines.append(
            f"| {row['target']} | {row['feature']} | {row['band_min_micron']}-{row['band_max_micron']} | "
            f"{row['n_points']} | {row['approx_z_vs_global_median']} |"
        )

    lines.extend([
        "",
        "## Figures",
        "",
    ])
    for path in plot_paths:
        lines.append(f"- `{path}`")

    lines.extend([
        "",
        "## Discovery-Search Expansion",
        "",
        "After this validation pass, the next executable stage is not spectroscopy first. It is transit discovery first:",
        "",
        "1. Pull TESS/Kepler light curves for bright G/K/M stars with low crowding.",
        "2. Run BLS/TLS-style transit search and reject common false positives.",
        "3. Use Gaia to check nearby contaminating stars and stellar parameters.",
        "4. Promote only surviving candidates to spectroscopy/archive checks.",
        "5. If spectra exist, run this same barcode extraction; if not, produce an ephemeris/follow-up package.",
        "",
        "## Sources Used By The Pipeline",
        "",
        "- NASA Exoplanet Archive TAP: https://exoplanetarchive.ipac.caltech.edu/TAP/sync",
        "- Exo.MAST spectra API: https://exo.mast.stsci.edu/api/v0.1/spectra",
        "- NASA WASP-39 b JWST transmission-spectrum context: https://science.nasa.gov/asset/webb/exoplanet-wasp-39-b-transmission-spectra/",
        "",
    ])
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def one_line(query: str) -> str:
    return " ".join(query.split())


def safe_slug(value: str) -> str:
    chars = []
    for char in value:
        if char.isalnum():
            chars.append(char.lower())
        elif char in {" ", "-", "_", ".", "/"}:
            chars.append("_")
    return "".join(chars).strip("_") or "target"


def to_float(value: object) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def is_number(value: object) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    raise SystemExit(main())
