#!/usr/bin/env python
"""Deep-vet one TOI lead from the TESS discovery pilot."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
from astroquery.gaia import Gaia


WORKSPACE = Path(r".")
TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EXOMAST_SPECTRA = "https://exo.mast.stsci.edu/api/v0.1/spectra"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-dir", required=True)
    parser.add_argument("--toi", default="TOI-1258.01")
    parser.add_argument(
        "--output-dir",
        default=str(WORKSPACE / "outputs" / f"EXOPLANET_TOI1258_DEEP_VET_{date.today().isoformat()}"),
    )
    args = parser.parse_args()

    pilot_dir = Path(args.pilot_dir)
    out_dir = Path(args.output_dir)
    dirs = create_dirs(out_dir)

    selected = read_csv(pilot_dir / "evidence" / "tables" / "selected_toi_candidates.csv")
    summary = read_csv(pilot_dir / "evidence" / "tables" / "tess_bls_candidate_summary.csv")
    toi_row = first(row for row in selected if row["toidisplay"] == args.toi)
    bls_row = first(row for row in summary if row["toi"] == args.toi)

    host = args.toi.split(".")[0]
    confirmed_matches = confirmed_crossmatch(host, float(toi_row["ra"]), float(toi_row["dec"]))
    write_csv(dirs["tables"] / "confirmed_planet_crossmatch.csv", confirmed_matches)

    gaia_rows = gaia_neighbor_query(float(toi_row["ra"]), float(toi_row["dec"]), radius_arcsec=120)
    contamination = contamination_summary(gaia_rows)
    write_csv(dirs["tables"] / "gaia_dr3_neighbors_120arcsec.csv", gaia_rows)
    write_csv(dirs["tables"] / "gaia_contamination_summary.csv", contamination)
    plot_gaia_neighbors(gaia_rows, dirs["figures"] / "toi1258_gaia_neighbor_map.png")

    spectra_probe = probe_spectra(args.toi, host)
    write_csv(dirs["tables"] / "spectra_archive_probe.csv", spectra_probe)

    verdict = classify(toi_row, bls_row, confirmed_matches, contamination)
    report = render_report(args.toi, toi_row, bls_row, confirmed_matches, contamination, spectra_probe, verdict, out_dir)
    (dirs["reports"] / "TOI1258_DEEP_VETTING_REPORT.md").write_text(report, encoding="utf-8")

    manifest = {
        "date": date.today().isoformat(),
        "toi": args.toi,
        "input_pilot_dir": str(pilot_dir),
        "output_dir": str(out_dir),
        "verdict": verdict,
        "claim_boundary": "Deep vetting supports candidate prioritization only, not planet confirmation.",
    }
    (out_dir / "target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(json.dumps({
        "output_dir": str(out_dir),
        "toi": args.toi,
        "verdict": verdict,
        "report": str(dirs["reports"] / "TOI1258_DEEP_VETTING_REPORT.md"),
        "confirmed_matches": len(confirmed_matches),
        "gaia_neighbors": len(gaia_rows),
    }, indent=2))
    return 0


def create_dirs(out_dir: Path) -> dict[str, Path]:
    paths = {
        "base": out_dir,
        "tables": out_dir / "evidence" / "tables",
        "figures": out_dir / "evidence" / "figures",
        "reports": out_dir / "evidence" / "reports",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def confirmed_crossmatch(host: str, ra: float, dec: float) -> list[dict[str, str]]:
    q1 = f"""
SELECT pl_name, hostname, discoverymethod, disc_year, sy_dist, ra, dec, pl_orbper, pl_rade
FROM pscomppars
WHERE hostname='{host}'
"""
    rows = tap_csv(q1)
    if rows:
        return rows
    radius = 0.01
    q2 = f"""
SELECT pl_name, hostname, discoverymethod, disc_year, sy_dist, ra, dec, pl_orbper, pl_rade
FROM pscomppars
WHERE ra BETWEEN {ra-radius} AND {ra+radius}
  AND dec BETWEEN {dec-radius} AND {dec+radius}
"""
    return tap_csv(q2)


def gaia_neighbor_query(ra: float, dec: float, radius_arcsec: float) -> list[dict[str, object]]:
    radius_deg = radius_arcsec / 3600.0
    q = f"""
SELECT TOP 200
  source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
  parallax, pmra, pmdec,
  DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec}))*3600.0 AS sep_arcsec
FROM gaiadr3.gaia_source
WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))=1
ORDER BY sep_arcsec ASC
"""
    table = Gaia.launch_job(q).get_results()
    rows = []
    target_mag = float(table[0]["phot_g_mean_mag"]) if len(table) else np.nan
    for row in table:
        gmag = to_float(row["phot_g_mean_mag"])
        flux_ratio = 10 ** (-0.4 * (gmag - target_mag)) if np.isfinite(gmag) and np.isfinite(target_mag) else np.nan
        rows.append(
            {
                "source_id": str(row["source_id"]),
                "ra": f"{float(row['ra']):.10f}",
                "dec": f"{float(row['dec']):.10f}",
                "sep_arcsec": f"{float(row['sep_arcsec']):.4f}",
                "phot_g_mean_mag": "" if not np.isfinite(gmag) else f"{gmag:.4f}",
                "phot_bp_mean_mag": safe_float_str(row["phot_bp_mean_mag"]),
                "phot_rp_mean_mag": safe_float_str(row["phot_rp_mean_mag"]),
                "parallax_mas": safe_float_str(row["parallax"]),
                "pmra_masyr": safe_float_str(row["pmra"]),
                "pmdec_masyr": safe_float_str(row["pmdec"]),
                "flux_ratio_vs_target_g": "" if not np.isfinite(flux_ratio) else f"{flux_ratio:.8f}",
            }
        )
    return rows


def contamination_summary(gaia_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for aperture in [21, 42, 63, 84, 120]:
        neighbors = [
            row for row in gaia_rows[1:]
            if to_float(row["sep_arcsec"]) is not None and to_float(row["sep_arcsec"]) <= aperture
        ]
        total_flux_ratio = sum(to_float(row["flux_ratio_vs_target_g"]) or 0.0 for row in neighbors)
        brightest = max((to_float(row["flux_ratio_vs_target_g"]) or 0.0 for row in neighbors), default=0.0)
        rows.append(
            {
                "aperture_radius_arcsec": aperture,
                "neighbor_count_excluding_target": len(neighbors),
                "summed_neighbor_flux_ratio_g": f"{total_flux_ratio:.6f}",
                "brightest_neighbor_flux_ratio_g": f"{brightest:.6f}",
                "estimated_depth_dilution_percent": f"{100.0 * total_flux_ratio / (1 + total_flux_ratio):.3f}",
            }
        )
    return rows


def probe_spectra(toi: str, host: str) -> list[dict[str, object]]:
    names = [toi, host, f"{host} b"]
    rows = []
    for name in names:
        url = f"{EXOMAST_SPECTRA}/{urllib.parse.quote(name, safe='')}/filelist/"
        try:
            payload = json.loads(fetch_text(url, timeout=60))
            files = payload.get("filenames") or []
            rows.append({"name": name, "status": "OK", "file_count": len(files), "files": ";".join(files)})
        except Exception as exc:
            rows.append({"name": name, "status": "ERROR", "file_count": 0, "files": repr(exc)})
    return rows


def classify(toi_row, bls_row, confirmed_matches, contamination) -> str:
    if confirmed_matches:
        return "KNOWN_CONFIRMED_PLANET"
    period_error = float(bls_row["period_error_percent"])
    primary_snr = float(bls_row["primary_depth_snr"])
    secondary_snr = abs(float(bls_row["secondary_depth_snr"]))
    odd_even_delta = abs(float(bls_row["odd_even_depth_delta_ppm"]))
    depth = float(bls_row["catalog_depth_ppm"])
    aperture42 = next(row for row in contamination if int(row["aperture_radius_arcsec"]) == 42)
    dilution = float(aperture42["estimated_depth_dilution_percent"])
    if period_error < 0.1 and primary_snr > 10 and secondary_snr < 3 and odd_even_delta < max(0.25 * depth, 150) and dilution < 20:
        return "HIGH_PRIORITY_UNCONFIRMED_TOI"
    if period_error < 1.0 and primary_snr > 8 and secondary_snr < 5:
        return "PRIORITY_UNCONFIRMED_TOI"
    return "NEEDS_ADDITIONAL_VETTING"


def render_report(toi, toi_row, bls_row, confirmed_matches, contamination, spectra_probe, verdict, out_dir: Path) -> str:
    lines = [
        "# TOI-1258.01 Deep Vetting Report",
        "",
        f"Created: {date.today().isoformat()}",
        "",
        "## Verdict",
        "",
        f"`{verdict}`",
        "",
        "Claim boundary: this is a candidate-prioritization report. It does not confirm a new planet.",
        "",
        "## Target",
        "",
        f"- TOI: `{toi}`",
        f"- TIC: `{toi_row['tid']}`",
        f"- RA/Dec: `{toi_row['ra']}`, `{toi_row['dec']}` deg",
        f"- Tmag: `{toi_row['st_tmag']}`",
        f"- TOI disposition: `{toi_row['tfopwg_disp']}`",
        f"- Catalog period: `{toi_row['pl_orbper']}` d",
        f"- Catalog radius: `{toi_row['pl_rade']}` Earth radii",
        f"- Catalog depth: `{toi_row['pl_trandep']}` ppm",
        "",
        "## Light-Curve Gate",
        "",
        f"- BLS period: `{bls_row['best_bls_period_days']}` d",
        f"- Period error: `{bls_row['period_error_percent']}` percent",
        f"- Primary SNR: `{bls_row['primary_depth_snr']}`",
        f"- Secondary SNR: `{bls_row['secondary_depth_snr']}`",
        f"- Odd/even depths: `{bls_row['odd_depth_ppm']}` / `{bls_row['even_depth_ppm']}` ppm",
        f"- Odd/even delta: `{bls_row['odd_even_depth_delta_ppm']}` ppm",
        "",
        "## Confirmed-Planet Crossmatch",
        "",
    ]
    if confirmed_matches:
        lines.append("Confirmed-planet match found:")
        for row in confirmed_matches:
            lines.append(f"- `{row.get('pl_name')}` around `{row.get('hostname')}`, period `{row.get('pl_orbper')}` d")
    else:
        lines.append("No matching confirmed planet was found in `pscomppars` by TOI hostname or tight coordinate box.")

    lines.extend([
        "",
        "## Gaia Contamination Gate",
        "",
        "| Aperture radius | Neighbors | Summed G-flux ratio | Brightest neighbor ratio | Estimated dilution |",
        "|---:|---:|---:|---:|---:|",
    ])
    for row in contamination:
        lines.append(
            f"| {row['aperture_radius_arcsec']} arcsec | {row['neighbor_count_excluding_target']} | "
            f"{row['summed_neighbor_flux_ratio_g']} | {row['brightest_neighbor_flux_ratio_g']} | "
            f"{row['estimated_depth_dilution_percent']}% |"
        )

    lines.extend([
        "",
        "## Spectra/Barcode Archive Gate",
        "",
        "| Name queried | Status | File count | Files / error |",
        "|---|---|---:|---|",
    ])
    for row in spectra_probe:
        lines.append(f"| {row['name']} | {row['status']} | {row['file_count']} | `{row['files']}` |")

    lines.extend([
        "",
        "## Files",
        "",
        f"- Gaia neighbors: `{out_dir / 'evidence' / 'tables' / 'gaia_dr3_neighbors_120arcsec.csv'}`",
        f"- Contamination summary: `{out_dir / 'evidence' / 'tables' / 'gaia_contamination_summary.csv'}`",
        f"- Confirmed crossmatch: `{out_dir / 'evidence' / 'tables' / 'confirmed_planet_crossmatch.csv'}`",
        f"- Gaia neighbor map: `{out_dir / 'evidence' / 'figures' / 'toi1258_gaia_neighbor_map.png'}`",
        "",
        "## Next Required Tests",
        "",
        "1. Centroid motion during transit using TESS target-pixel files.",
        "2. Full-sector/all-sector light-curve stitch, not just two products.",
        "3. Gaia-based aperture-level dilution correction in the transit fit.",
        "4. ExoFOP/manual literature check before public discovery language.",
        "5. If it survives, prepare ephemeris and follow-up package.",
        "",
    ])
    return "\n".join(lines)


def plot_gaia_neighbors(rows: list[dict[str, object]], path: Path) -> None:
    import matplotlib.pyplot as plt

    if not rows:
        return
    target_ra = float(rows[0]["ra"])
    target_dec = float(rows[0]["dec"])
    xs = []
    ys = []
    sizes = []
    labels = []
    for idx, row in enumerate(rows):
        ra = float(row["ra"])
        dec = float(row["dec"])
        east = (ra - target_ra) * math.cos(math.radians(target_dec)) * 3600
        north = (dec - target_dec) * 3600
        gmag = to_float(row["phot_g_mean_mag"]) or 20
        xs.append(east)
        ys.append(north)
        sizes.append(max(15, 300 * 10 ** (-0.25 * (gmag - 10))))
        labels.append(idx)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(xs, ys, s=sizes, alpha=0.65, color="#1f5a8a", edgecolors="white", linewidths=0.5)
    for radius in [21, 42, 63]:
        ax.add_patch(plt.Circle((0, 0), radius, fill=False, ls="--", lw=1, alpha=0.5))
    ax.scatter([0], [0], marker="+", s=180, color="#d4472f", label="target")
    ax.set_xlabel("arcsec east of target")
    ax.set_ylabel("arcsec north of target")
    ax.set_title("TOI-1258.01 Gaia DR3 neighbor field")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def tap_csv(query: str) -> list[dict[str, str]]:
    url = TAP_SYNC + "?" + urllib.parse.urlencode({"query": " ".join(query.split()), "format": "csv"})
    with urllib.request.urlopen(url, timeout=180) as response:
        text = response.read().decode("utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(text)))


def fetch_text(url: str, timeout: int = 60) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def read_csv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(newline="", encoding="utf-8")))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def first(items):
    for item in items:
        return item
    raise RuntimeError("required row not found")


def to_float(value) -> float | None:
    try:
        if value is None or str(value) in {"", "nan", "NaN", "--"}:
            return None
        out = float(value)
        if not np.isfinite(out):
            return None
        return out
    except Exception:
        return None


def safe_float_str(value) -> str:
    out = to_float(value)
    return "" if out is None else f"{out:.6f}"


if __name__ == "__main__":
    raise SystemExit(main())
