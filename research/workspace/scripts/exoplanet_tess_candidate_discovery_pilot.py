#!/usr/bin/env python
"""Stage 3 exoplanet candidate discovery pilot.

The script queries live TESS Objects of Interest (TOIs), chooses bright
unconfirmed planet-candidate rows, downloads public TESS light-curve FITS files
from MAST, runs a local Box Least Squares search, and compares the recovered
period against the TOI catalog period.

This is a candidate-validation/discovery-pilot layer, not a confirmation claim.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import shutil
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.timeseries import BoxLeastSquares
from astroquery.mast import Observations


WORKSPACE = Path(r".")
TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


@dataclass
class Candidate:
    tid: str
    toi: str
    toidisplay: str
    disposition: str
    tmag: float
    ra: float
    dec: float
    period_days: float
    duration_hours: float
    depth_ppm: float
    radius_rearth: float | None
    eq_temp_k: float | None
    distance_pc: float | None
    teff_k: float | None
    stellar_radius_rsun: float | None
    row: dict[str, str]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=str(WORKSPACE / "outputs" / f"EXOPLANET_TESS_DISCOVERY_PILOT_{date.today().isoformat()}"),
    )
    parser.add_argument("--max-targets", type=int, default=5)
    parser.add_argument("--max-lc-products-per-target", type=int, default=3)
    parser.add_argument("--min-period", type=float, default=0.5)
    parser.add_argument("--max-period", type=float, default=20.0)
    parser.add_argument("--n-periods", type=int, default=5000)
    parser.add_argument("--exclude-confirmed", action="store_true", help="Exclude TOI hosts already present in pscomppars.")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    dirs = create_dirs(out_dir)

    candidates = fetch_toi_candidates(limit=max(args.max_targets * 8, 40))
    if args.exclude_confirmed:
        confirmed_hosts = fetch_confirmed_toi_hosts()
        candidates = [candidate for candidate in candidates if toi_hostname(candidate) not in confirmed_hosts]
    ranked = rank_candidates(candidates)
    selected = ranked[: args.max_targets]
    write_csv(dirs["tables"] / "live_toi_ranked_candidates.csv", [candidate.row for candidate in ranked])
    write_csv(dirs["tables"] / "selected_toi_candidates.csv", [candidate.row for candidate in selected])

    summaries = []
    for candidate in selected:
        try:
            result = process_candidate(candidate, dirs, args)
            summaries.append(result)
        except Exception as exc:
            summaries.append(
                {
                    "tid": candidate.tid,
                    "toi": candidate.toidisplay,
                    "status": "ERROR",
                    "error": repr(exc),
                }
            )

    write_csv(dirs["tables"] / "tess_bls_candidate_summary.csv", summaries)
    report = render_report(selected, summaries, out_dir)
    (dirs["reports"] / "TESS_CANDIDATE_DISCOVERY_PILOT_REPORT.md").write_text(report, encoding="utf-8")
    shutil.copy2(Path(__file__), dirs["scripts"] / Path(__file__).name)

    print(
        json.dumps(
            {
                "output_dir": str(out_dir),
                "selected": [c.toidisplay for c in selected],
                "summary": summaries,
                "report": str(dirs["reports"] / "TESS_CANDIDATE_DISCOVERY_PILOT_REPORT.md"),
            },
            indent=2,
        )
    )
    return 0


def create_dirs(out_dir: Path) -> dict[str, Path]:
    paths = {
        "base": out_dir,
        "archive": out_dir / "evidence" / "archive_gates",
        "tables": out_dir / "evidence" / "tables",
        "figures": out_dir / "evidence" / "figures",
        "reports": out_dir / "evidence" / "reports",
        "scripts": out_dir / "evidence" / "scripts",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def fetch_toi_candidates(limit: int) -> list[Candidate]:
    query = f"""
SELECT TOP {limit}
       tid, toi, toidisplay, tfopwg_disp, st_tmag, ra, dec, pl_orbper,
       pl_trandurh, pl_trandep, pl_rade, pl_eqt, st_dist, st_teff, st_rad,
       sectors, rowupdate, release_date
FROM toi
WHERE tfopwg_disp='PC'
  AND st_tmag < 10.8
  AND pl_orbper IS NOT NULL
  AND pl_trandurh IS NOT NULL
  AND pl_trandep IS NOT NULL
  AND pl_orbper BETWEEN 0.5 AND 20
  AND pl_trandep BETWEEN 100 AND 20000
ORDER BY st_tmag ASC
"""
    rows = tap_csv(query)
    return [candidate_from_row(row) for row in rows]


def fetch_confirmed_toi_hosts() -> set[str]:
    query = """
SELECT DISTINCT hostname
FROM pscomppars
WHERE hostname LIKE 'TOI-%'
"""
    return {row["hostname"] for row in tap_csv(query) if row.get("hostname")}


def toi_hostname(candidate: Candidate) -> str:
    prefix = candidate.toidisplay.split(".")[0]
    return prefix


def candidate_from_row(row: dict[str, str]) -> Candidate:
    return Candidate(
        tid=row["tid"],
        toi=row["toi"],
        toidisplay=row["toidisplay"],
        disposition=row["tfopwg_disp"],
        tmag=float(row["st_tmag"]),
        ra=float(row["ra"]),
        dec=float(row["dec"]),
        period_days=float(row["pl_orbper"]),
        duration_hours=float(row["pl_trandurh"]),
        depth_ppm=float(row["pl_trandep"]),
        radius_rearth=to_float(row.get("pl_rade")),
        eq_temp_k=to_float(row.get("pl_eqt")),
        distance_pc=to_float(row.get("st_dist")),
        teff_k=to_float(row.get("st_teff")),
        stellar_radius_rsun=to_float(row.get("st_rad")),
        row=row,
    )


def rank_candidates(candidates: list[Candidate]) -> list[Candidate]:
    def score(candidate: Candidate) -> float:
        score_value = 0.0
        score_value += max(0.0, 12.0 - candidate.tmag) * 8.0
        score_value += min(20.0, candidate.depth_ppm / 500.0)
        if candidate.distance_pc:
            score_value += max(0.0, 200.0 - candidate.distance_pc) / 20.0
        if candidate.radius_rearth and 1.0 <= candidate.radius_rearth <= 4.0:
            score_value += 8.0
        if candidate.teff_k and 3500 <= candidate.teff_k <= 6200:
            score_value += 5.0
        return score_value

    return sorted(candidates, key=score, reverse=True)


def process_candidate(candidate: Candidate, dirs: dict[str, Path], args: argparse.Namespace) -> dict[str, object]:
    target_dir = dirs["archive"] / safe_slug(candidate.toidisplay)
    target_dir.mkdir(parents=True, exist_ok=True)

    products = find_lightcurve_products(candidate)
    selected_products = products[: args.max_lc_products_per_target]
    if not selected_products:
        raise RuntimeError("no light-curve FITS products found")

    downloaded = []
    time_parts = []
    flux_parts = []
    err_parts = []
    for product in selected_products:
        manifest = Observations.download_products(product, download_dir=str(target_dir), cache=True)
        local_path = Path(str(manifest["Local Path"][0]))
        downloaded.append(local_path)
        time, flux, err = read_tess_lightcurve(local_path)
        if len(time):
            time_parts.append(time)
            flux_parts.append(flux)
            err_parts.append(err)

    if not time_parts:
        raise RuntimeError("no usable finite light-curve samples")

    time = np.concatenate(time_parts)
    flux = np.concatenate(flux_parts)
    err = np.concatenate(err_parts)
    order = np.argsort(time)
    time = time[order]
    flux = flux[order]
    err = err[order]

    bls_result = run_bls(time, flux, err, args.min_period, args.max_period, args.n_periods)
    catalog_period = candidate.period_days
    best_period = bls_result["best_period_days"]
    period_error = period_match_percent(best_period, catalog_period)
    harmonic = best_harmonic(best_period, catalog_period)
    catalog_power = evaluate_period(time, flux, err, catalog_period, candidate.duration_hours / 24.0)
    quick_vet = quick_vetting_metrics(
        time,
        flux,
        bls_result["best_period_days"],
        bls_result["best_t0"],
        bls_result["best_duration_days"],
    )

    verdict = "BLS_MATCH" if period_error <= 1.0 or harmonic["error_percent"] <= 1.0 else "BLS_MISMATCH_CHECK"
    folded_path = dirs["figures"] / f"{safe_slug(candidate.toidisplay)}_folded_bls.png"
    periodogram_path = dirs["figures"] / f"{safe_slug(candidate.toidisplay)}_bls_periodogram.png"
    plot_periodogram(candidate, bls_result, periodogram_path)
    plot_folded(candidate, time, flux, best_period, bls_result["best_t0"], folded_path)

    return {
        "tid": candidate.tid,
        "toi": candidate.toidisplay,
        "status": "OK",
        "tfopwg_disp": candidate.disposition,
        "tmag": f"{candidate.tmag:.4f}",
        "catalog_period_days": f"{catalog_period:.8f}",
        "catalog_duration_hours": f"{candidate.duration_hours:.4f}",
        "catalog_depth_ppm": f"{candidate.depth_ppm:.2f}",
        "best_bls_period_days": f"{best_period:.8f}",
        "best_bls_power": f"{bls_result['best_power']:.4g}",
        "best_bls_t0": f"{bls_result['best_t0']:.8f}",
        "best_bls_duration_days": f"{bls_result['best_duration_days']:.8f}",
        "best_bls_depth_fraction": f"{bls_result['best_depth']:.8f}",
        "period_error_percent": f"{period_error:.4f}",
        "best_harmonic_ratio": harmonic["ratio"],
        "best_harmonic_error_percent": f"{harmonic['error_percent']:.4f}",
        "catalog_period_power": f"{catalog_power:.4g}",
        "odd_depth_ppm": f"{quick_vet['odd_depth_ppm']:.2f}",
        "even_depth_ppm": f"{quick_vet['even_depth_ppm']:.2f}",
        "odd_even_depth_delta_ppm": f"{quick_vet['odd_even_depth_delta_ppm']:.2f}",
        "secondary_depth_ppm": f"{quick_vet['secondary_depth_ppm']:.2f}",
        "primary_depth_snr": f"{quick_vet['primary_depth_snr']:.3f}",
        "secondary_depth_snr": f"{quick_vet['secondary_depth_snr']:.3f}",
        "n_points": len(time),
        "n_products": len(downloaded),
        "downloaded_files": ";".join(str(path) for path in downloaded),
        "folded_plot": str(folded_path),
        "periodogram_plot": str(periodogram_path),
        "verdict": verdict,
    }


def find_lightcurve_products(candidate: Candidate):
    obs = Observations.query_criteria(
        coordinates=f"{candidate.ra} {candidate.dec}",
        radius="0.02 deg",
        dataproduct_type="timeseries",
        obs_collection="TESS",
    )
    if len(obs) == 0:
        return []
    products = Observations.get_product_list(obs)
    products = Observations.filter_products(products, extension="fits")
    lc_products = []
    for product in products:
        filename = str(product["productFilename"])
        subgroup = str(product.get("productSubGroupDescription", ""))
        if filename.endswith("_lc.fits") and "-tp.fits" not in filename:
            lc_products.append(product)
        elif subgroup.upper() == "LC" and filename.endswith(".fits"):
            lc_products.append(product)

    def product_key(product):
        filename = str(product["productFilename"])
        size = int(product["size"]) if str(product["size"]).isdigit() else 10**12
        cadence_rank = 0 if "-s_lc.fits" in filename else 1
        return (cadence_rank, size)

    return sorted(lc_products, key=product_key)


def read_tess_lightcurve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(path) as hdul:
        data = hdul[1].data
        time = np.array(data["TIME"], dtype=float)
        if "PDCSAP_FLUX" in data.columns.names:
            flux = np.array(data["PDCSAP_FLUX"], dtype=float)
            err = np.array(data["PDCSAP_FLUX_ERR"], dtype=float)
        else:
            flux = np.array(data["SAP_FLUX"], dtype=float)
            err = np.array(data["SAP_FLUX_ERR"], dtype=float)
        quality = np.array(data["QUALITY"], dtype=int) if "QUALITY" in data.columns.names else np.zeros(len(time), dtype=int)

    mask = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) & (err > 0) & (quality == 0)
    time = time[mask]
    flux = flux[mask]
    err = err[mask]
    if len(time) == 0:
        return time, flux, err

    median_flux = np.nanmedian(flux)
    flux = flux / median_flux
    err = err / median_flux
    clipped = sigma_clip(flux, sigma=6, maxiters=3)
    keep = ~np.array(clipped.mask, dtype=bool)
    return time[keep], np.array(clipped.data[keep], dtype=float), err[keep]


def run_bls(time: np.ndarray, flux: np.ndarray, err: np.ndarray, min_period: float, max_period: float, n_periods: int) -> dict[str, object]:
    periods = np.linspace(min_period, max_period, n_periods)
    durations = np.linspace(0.04, 0.35, 10)
    bls = BoxLeastSquares(time, flux, dy=err)
    power = bls.power(periods, durations)
    idx = int(np.argmax(power.power))
    stats = bls.compute_stats(power.period[idx], power.duration[idx], power.transit_time[idx])
    return {
        "periods": np.array(power.period, dtype=float),
        "power": np.array(power.power, dtype=float),
        "best_period_days": float(power.period[idx]),
        "best_duration_days": float(power.duration[idx]),
        "best_t0": float(power.transit_time[idx]),
        "best_depth": float(power.depth[idx]),
        "best_power": float(power.power[idx]),
        "depth_snr": float(stats.get("depth_snr", np.nan)),
    }


def evaluate_period(time: np.ndarray, flux: np.ndarray, err: np.ndarray, period: float, duration: float) -> float:
    duration = max(0.04, min(0.35, duration))
    bls = BoxLeastSquares(time, flux, dy=err)
    power = bls.power(np.array([period], dtype=float), np.array([duration], dtype=float))
    return float(power.power[0])


def quick_vetting_metrics(time: np.ndarray, flux: np.ndarray, period: float, t0: float, duration: float) -> dict[str, float]:
    epoch = np.rint((time - t0) / period).astype(int)
    primary_phase_days = time - (t0 + epoch * period)
    secondary_phase_days = time - (t0 + (epoch + 0.5) * period)
    half_window = max(duration / 2.0, 0.03)
    oot = (np.abs(primary_phase_days) > 2.5 * half_window) & (np.abs(secondary_phase_days) > 2.5 * half_window)
    primary = np.abs(primary_phase_days) <= half_window
    secondary = np.abs(secondary_phase_days) <= half_window
    odd = primary & ((epoch % 2) != 0)
    even = primary & ((epoch % 2) == 0)

    baseline = robust_median(flux[oot]) if np.any(oot) else robust_median(flux)
    odd_depth = baseline - robust_median(flux[odd]) if np.any(odd) else np.nan
    even_depth = baseline - robust_median(flux[even]) if np.any(even) else np.nan
    primary_depth = baseline - robust_median(flux[primary]) if np.any(primary) else np.nan
    secondary_depth = baseline - robust_median(flux[secondary]) if np.any(secondary) else np.nan
    oot_std = robust_std(flux[oot]) if np.any(oot) else robust_std(flux)

    return {
        "odd_depth_ppm": 1e6 * odd_depth,
        "even_depth_ppm": 1e6 * even_depth,
        "odd_even_depth_delta_ppm": 1e6 * abs(odd_depth - even_depth),
        "secondary_depth_ppm": 1e6 * secondary_depth,
        "primary_depth_snr": primary_depth / (oot_std / math.sqrt(max(np.count_nonzero(primary), 1))) if oot_std > 0 else np.nan,
        "secondary_depth_snr": secondary_depth / (oot_std / math.sqrt(max(np.count_nonzero(secondary), 1))) if oot_std > 0 else np.nan,
    }


def robust_median(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    return float(np.nanmedian(values))


def robust_std(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    if mad > 0:
        return float(1.4826 * mad)
    return float(np.nanstd(values))


def period_match_percent(recovered: float, catalog: float) -> float:
    return 100.0 * abs(recovered - catalog) / catalog


def best_harmonic(recovered: float, catalog: float) -> dict[str, object]:
    ratios = [0.25, 1 / 3, 0.5, 1.0, 2.0, 3.0, 4.0]
    best = {"ratio": "none", "error_percent": float("inf")}
    for ratio in ratios:
        expected = catalog * ratio
        error = 100.0 * abs(recovered - expected) / expected
        if error < best["error_percent"]:
            best = {"ratio": f"{ratio:g}", "error_percent": error}
    return best


def plot_periodogram(candidate: Candidate, result: dict[str, object], path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(result["periods"], result["power"], color="#174a7c", lw=1.0)
    ax.axvline(candidate.period_days, color="#d4472f", lw=1.1, label="TOI catalog period")
    ax.axvline(result["best_period_days"], color="#222222", ls="--", lw=1.1, label="best BLS")
    ax.set_title(f"{candidate.toidisplay} BLS period search")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("BLS power")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_folded(candidate: Candidate, time: np.ndarray, flux: np.ndarray, period: float, t0: float, path: Path) -> None:
    import matplotlib.pyplot as plt

    phase = ((time - t0 + 0.5 * period) % period) / period - 0.5
    order = np.argsort(phase)
    phase = phase[order]
    flux = flux[order]
    keep = np.abs(phase) < max(0.08, 3 * candidate.duration_hours / 24.0 / period)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(phase[keep], flux[keep], s=5, alpha=0.35, color="#1f5a8a", edgecolors="none")
    ax.set_title(f"{candidate.toidisplay} folded on recovered period")
    ax.set_xlabel("Orbital phase around transit")
    ax.set_ylabel("Normalized PDCSAP flux")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_report(candidates: list[Candidate], summaries: list[dict[str, object]], out_dir: Path) -> str:
    lines = [
        "# TESS Candidate Discovery Pilot Report",
        "",
        f"Created: {date.today().isoformat()}",
        "",
        "## Result",
        "",
        "This run queried live NASA Exoplanet Archive TOI records, selected bright planet-candidate rows, downloaded public TESS light-curve FITS products from MAST, and ran a local Box Least Squares search.",
        "",
        "Claim boundary: these are not confirmed new planets. They are candidate-validation outputs. A candidate needs false-positive vetting, centroid checks, Gaia contamination checks, and independent follow-up before confirmation language is justified.",
        "",
        "## Selected Targets",
        "",
        "| TOI | TIC | Tmag | Period d | Depth ppm | Radius Re | Disposition |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for candidate in candidates:
        radius = "" if candidate.radius_rearth is None else f"{candidate.radius_rearth:.3f}"
        lines.append(
            f"| {candidate.toidisplay} | {candidate.tid} | {candidate.tmag:.3f} | "
            f"{candidate.period_days:.6f} | {candidate.depth_ppm:.1f} | {radius} | {candidate.disposition} |"
        )

    lines.extend([
        "",
        "## BLS Results",
        "",
        "| TOI | Status | Catalog period d | BLS period d | Error % | Harmonic | Primary SNR | Secondary SNR | Points | Verdict |",
        "|---|---|---:|---:|---:|---|---:|---:|---:|---|",
    ])
    for summary in summaries:
        if summary.get("status") != "OK":
            lines.append(f"| {summary.get('toi')} | ERROR | | | | | | | | `{summary.get('error')}` |")
            continue
        lines.append(
            f"| {summary['toi']} | OK | {summary['catalog_period_days']} | {summary['best_bls_period_days']} | "
            f"{summary['period_error_percent']} | {summary['best_harmonic_ratio']} ({summary['best_harmonic_error_percent']}%) | "
            f"{summary['primary_depth_snr']} | {summary['secondary_depth_snr']} | {summary['n_points']} | {summary['verdict']} |"
        )

    lines.extend([
        "",
        "## Files",
        "",
        f"- Ranked candidates: `{out_dir / 'evidence' / 'tables' / 'live_toi_ranked_candidates.csv'}`",
        f"- Selected candidates: `{out_dir / 'evidence' / 'tables' / 'selected_toi_candidates.csv'}`",
        f"- BLS summary: `{out_dir / 'evidence' / 'tables' / 'tess_bls_candidate_summary.csv'}`",
        f"- Figures: `{out_dir / 'evidence' / 'figures'}`",
        "",
        "## Next Validation Gates",
        "",
        "1. Run odd/even transit-depth comparison.",
        "2. Check secondary eclipse at phase 0.5.",
        "3. Pull Gaia DR3 neighbors inside the TESS aperture and estimate contamination.",
        "4. Inspect centroid motion during transit.",
        "5. Crossmatch against confirmed planets and known false positives before any discovery claim.",
        "6. Only then promote survivors to the atmosphere-barcode archive check.",
        "",
        "## Sources",
        "",
        "- NASA Exoplanet Archive TOI table via TAP.",
        "- MAST TESS light-curve products via astroquery.mast Observations.",
        "- Astropy BoxLeastSquares for period search.",
        "",
    ])
    return "\n".join(lines)


def tap_csv(query: str) -> list[dict[str, str]]:
    url = TAP_SYNC + "?" + urllib.parse.urlencode({"query": " ".join(query.split()), "format": "csv"})
    with urllib.request.urlopen(url, timeout=180) as response:
        text = response.read().decode("utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(text)))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_slug(value: str) -> str:
    keep = []
    for char in value.lower():
        if char.isalnum():
            keep.append(char)
        elif char in {" ", "-", "_", "."}:
            keep.append("_")
    return "".join(keep).strip("_") or "target"


def to_float(value: str | None) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except ValueError:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
