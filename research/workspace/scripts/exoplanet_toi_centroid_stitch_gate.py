#!/usr/bin/env python
"""Centroid-motion and stitched-light-curve gate for one TOI.

This script is intentionally evidence-packaging oriented: it downloads public
TESS light curves and target-pixel files, stitches all usable sectors, recovers
the transit period with BLS, and checks whether in-transit target-pixel
centroids move in a way that suggests the dip is coming from another source.

It is a validation gate, not a planet-confirmation claim.
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

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.timeseries import BoxLeastSquares
from astroquery.mast import Observations


WORKSPACE = Path(r".")
TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
BTJD_OFFSET = 2457000.0
TESS_PIXEL_SCALE_ARCSEC = 21.0


@dataclass
class Target:
    toi: str
    tic: str
    ra: float
    dec: float
    period_days: float
    t0_btjd: float
    duration_days: float
    depth_ppm: float
    tmag: float | None
    disposition: str
    source_row: dict[str, str]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toi", default="TOI-1258.01")
    parser.add_argument("--tic", default="288132261")
    parser.add_argument("--ra", type=float, default=209.4843130)
    parser.add_argument("--dec", type=float, default=79.5833220)
    parser.add_argument(
        "--output-dir",
        default=str(WORKSPACE / "outputs" / f"EXOPLANET_TOI1258_CENTROID_STITCH_GATE_{date.today().isoformat()}"),
    )
    parser.add_argument("--period-min", type=float, default=0.5)
    parser.add_argument("--period-max", type=float, default=20.0)
    parser.add_argument("--n-periods", type=int, default=6000)
    parser.add_argument("--max-tpf-products", type=int, default=0, help="0 means use all standard TPF products.")
    parser.add_argument("--download-cache", action="store_true", help="Allow astroquery to reuse cached MAST files.")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    dirs = create_dirs(out_dir)

    target = fetch_target(args.toi, args.tic, args.ra, args.dec)
    products = find_tess_products(target)
    write_csv(dirs["tables"] / "mast_tess_product_inventory.csv", products["inventory"])

    lc_rows = products["lc"]
    tp_rows = products["tp"]
    if args.max_tpf_products and args.max_tpf_products > 0:
        tp_rows = tp_rows[: args.max_tpf_products]

    lc_paths = download_product_rows(lc_rows, dirs["archive"] / "lightcurves", cache=args.download_cache)
    tp_paths = download_product_rows(tp_rows, dirs["archive"] / "target_pixel_files", cache=args.download_cache)

    stitched = build_stitched_lightcurve(lc_paths)
    write_csv(dirs["tables"] / "stitched_lightcurve_summary_by_sector.csv", stitched["sector_summary"])
    write_csv(dirs["tables"] / "stitched_lightcurve_downsample.csv", downsample_rows(stitched["rows"], max_rows=25000))

    bls_summary = run_bls_and_plots(
        target=target,
        stitched=stitched,
        dirs=dirs,
        period_min=args.period_min,
        period_max=args.period_max,
        n_periods=args.n_periods,
    )
    write_csv(dirs["tables"] / "stitched_bls_summary.csv", [bls_summary])

    transit_metrics = compute_transit_metrics(
        stitched["time"],
        stitched["flux_rel"],
        target.period_days,
        target.t0_btjd,
        target.duration_days,
        sector=stitched["sector"],
    )
    write_csv(dirs["tables"] / "stitched_transit_metrics.csv", [transit_metrics["summary"]])
    write_csv(dirs["tables"] / "stitched_transit_depth_by_sector.csv", transit_metrics["sector_rows"])

    centroid_result = run_centroid_gate(target, tp_paths, dirs)
    write_csv(dirs["tables"] / "centroid_shift_summary_by_sector.csv", centroid_result["sector_rows"])
    write_csv(dirs["tables"] / "centroid_gate_summary.csv", [centroid_result["summary"]])

    report = render_report(
        target=target,
        products=products,
        lc_paths=lc_paths,
        tp_paths=tp_paths,
        bls_summary=bls_summary,
        transit_metrics=transit_metrics,
        centroid_result=centroid_result,
        out_dir=out_dir,
    )
    report_path = dirs["reports"] / "TOI1258_CENTROID_STITCH_GATE_REPORT.md"
    report_path.write_text(report, encoding="utf-8")
    shutil.copy2(Path(__file__), dirs["scripts"] / Path(__file__).name)

    manifest = {
        "created": date.today().isoformat(),
        "target": {
            "toi": target.toi,
            "tic": target.tic,
            "ra_deg": target.ra,
            "dec_deg": target.dec,
            "period_days": target.period_days,
            "t0_btjd": target.t0_btjd,
            "duration_days": target.duration_days,
            "depth_ppm": target.depth_ppm,
        },
        "output_dir": str(out_dir),
        "report": str(report_path),
        "claim_boundary": "Centroid and stitch tests can reject common false-positive behavior, but do not by themselves confirm a planet.",
    }
    (out_dir / "target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "output_dir": str(out_dir),
                "report": str(report_path),
                "lc_products": len(lc_paths),
                "tpf_products": len(tp_paths),
                "bls_period_days": bls_summary["best_bls_period_days"],
                "period_error_percent": bls_summary["period_error_percent"],
                "median_abs_diff_centroid_offset_arcsec": centroid_result["summary"][
                    "median_abs_difference_image_offset_arcsec"
                ],
                "centroid_verdict": centroid_result["summary"]["centroid_gate_verdict"],
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


def fetch_target(toi: str, fallback_tic: str, fallback_ra: float, fallback_dec: float) -> Target:
    query = f"""
SELECT TOP 1
       tid, toidisplay, tfopwg_disp, st_tmag, ra, dec, pl_orbper,
       pl_tranmid, pl_trandurh, pl_trandep
FROM toi
WHERE toidisplay='{toi}'
"""
    rows = tap_csv(query)
    if rows:
        row = rows[0]
        tic = row.get("tid") or fallback_tic
        ra = float(row.get("ra") or fallback_ra)
        dec = float(row.get("dec") or fallback_dec)
        period_days = float(row["pl_orbper"])
        t0_btjd = float(row["pl_tranmid"]) - BTJD_OFFSET
        duration_days = float(row["pl_trandurh"]) / 24.0
        depth_ppm = float(row["pl_trandep"])
        tmag = to_float(row.get("st_tmag"))
        disposition = row.get("tfopwg_disp", "")
        return Target(toi, tic, ra, dec, period_days, t0_btjd, duration_days, depth_ppm, tmag, disposition, row)

    return Target(
        toi=toi,
        tic=fallback_tic,
        ra=fallback_ra,
        dec=fallback_dec,
        period_days=10.2613397,
        t0_btjd=3333.6487860,
        duration_days=4.855 / 24.0,
        depth_ppm=726.0,
        tmag=None,
        disposition="",
        source_row={},
    )


def find_tess_products(target: Target) -> dict[str, object]:
    obs = Observations.query_criteria(
        coordinates=f"{target.ra} {target.dec}",
        radius="0.02 deg",
        dataproduct_type="timeseries",
        obs_collection="TESS",
    )
    products = Observations.get_product_list(obs)
    inventory = []
    lc_indices = []
    tp_indices = []
    tic_token = f"{int(target.tic):016d}" if str(target.tic).isdigit() else str(target.tic)
    for idx, product in enumerate(products):
        filename = str(product["productFilename"])
        subgroup = str(product.get("productSubGroupDescription", ""))
        size = int(product.get("size", 0) or 0)
        obs_id = str(product.get("obsID", ""))
        is_target = tic_token in filename
        is_standard_lc = is_target and filename.endswith("-s_lc.fits") and subgroup == "LC"
        is_standard_tp = is_target and filename.endswith("-s_tp.fits") and subgroup == "TP"
        if is_standard_lc:
            lc_indices.append(idx)
        if is_standard_tp:
            tp_indices.append(idx)
        if is_target:
            inventory.append(
                {
                    "productFilename": filename,
                    "productSubGroupDescription": subgroup,
                    "size_bytes": size,
                    "size_mb": f"{size / 1024 / 1024:.3f}",
                    "obsID": obs_id,
                    "selected_for_lc_stitch": int(is_standard_lc),
                    "selected_for_centroid_gate": int(is_standard_tp),
                }
            )

    lc = products[lc_indices]
    tp = products[tp_indices]
    lc.sort("productFilename")
    tp.sort("productFilename")
    inventory = sorted(inventory, key=lambda row: row["productFilename"])
    return {"inventory": inventory, "lc": lc, "tp": tp}


def download_product_rows(product_rows, download_dir: Path, cache: bool) -> list[Path]:
    if len(product_rows) == 0:
        return []
    download_dir.mkdir(parents=True, exist_ok=True)
    manifest = Observations.download_products(product_rows, download_dir=str(download_dir), cache=cache)
    paths = [Path(str(path)) for path in manifest["Local Path"] if str(path)]
    return sorted(paths, key=lambda p: p.name)


def build_stitched_lightcurve(lc_paths: list[Path]) -> dict[str, object]:
    time_parts = []
    flux_parts = []
    err_parts = []
    sector_parts = []
    rows = []
    sector_summary = []

    for path in lc_paths:
        with fits.open(path, memmap=True) as hdul:
            data = hdul[1].data
            header = hdul[0].header
            sector = int(header.get("SECTOR", -1))
            time = np.asarray(data["TIME"], dtype=float)
            flux = np.asarray(data["PDCSAP_FLUX"], dtype=float)
            err = np.asarray(data["PDCSAP_FLUX_ERR"], dtype=float)
            quality_name = "SAP_QUALITY" if "SAP_QUALITY" in data.names else "QUALITY"
            quality = np.asarray(data[quality_name], dtype=int)

        mask = np.isfinite(time) & np.isfinite(flux) & (quality == 0)
        time = time[mask]
        flux = flux[mask]
        err = err[mask]
        if len(time) < 10:
            continue

        clipped = sigma_clip(flux, sigma=5.0, maxiters=5, masked=True)
        med = float(np.ma.median(clipped))
        if not np.isfinite(med) or med <= 0:
            continue
        rel = flux / med
        rel_err = err / med
        good = np.isfinite(rel) & np.isfinite(rel_err)
        time = time[good]
        rel = rel[good]
        rel_err = rel_err[good]
        sector_arr = np.full(len(time), sector, dtype=int)

        time_parts.append(time)
        flux_parts.append(rel)
        err_parts.append(rel_err)
        sector_parts.append(sector_arr)
        sector_summary.append(
            {
                "sector": sector,
                "lc_path": str(path),
                "cadence_count_quality0": len(time),
                "time_min_btjd": f"{float(np.min(time)):.6f}",
                "time_max_btjd": f"{float(np.max(time)):.6f}",
                "median_pdcsap_flux": f"{med:.8g}",
                "rms_ppm": f"{robust_std((rel - 1.0) * 1e6):.3f}",
            }
        )
        for t, f, e in zip(time, rel, rel_err):
            rows.append({"sector": sector, "time_btjd": f"{t:.8f}", "flux_rel": f"{f:.10f}", "flux_err_rel": f"{e:.10f}"})

    if not time_parts:
        raise RuntimeError("No usable light-curve cadences found.")

    time_all = np.concatenate(time_parts)
    flux_all = np.concatenate(flux_parts)
    err_all = np.concatenate(err_parts)
    sector_all = np.concatenate(sector_parts)
    order = np.argsort(time_all)
    return {
        "time": time_all[order],
        "flux_rel": flux_all[order],
        "err_rel": err_all[order],
        "sector": sector_all[order],
        "sector_summary": sorted(sector_summary, key=lambda row: row["sector"]),
        "rows": sorted(rows, key=lambda row: (int(row["sector"]), float(row["time_btjd"]))),
    }


def run_bls_and_plots(
    target: Target,
    stitched: dict[str, object],
    dirs: dict[str, Path],
    period_min: float,
    period_max: float,
    n_periods: int,
) -> dict[str, object]:
    time = stitched["time"]
    flux_rel = stitched["flux_rel"]
    err_rel = stitched["err_rel"]
    y = flux_rel - 1.0
    finite = np.isfinite(time) & np.isfinite(y) & np.isfinite(err_rel)
    time = time[finite]
    y = y[finite]
    err_rel = np.clip(err_rel[finite], 1e-7, np.inf)

    periods = np.linspace(period_min, period_max, n_periods)
    duration_min = max(0.03, target.duration_days * 0.35)
    duration_max = min(0.75, target.duration_days * 2.25)
    durations = np.linspace(duration_min, duration_max, 14)
    bls = BoxLeastSquares(time, y, dy=err_rel)
    result = bls.power(periods, durations, objective="snr")
    best_idx = int(np.nanargmax(result.power))
    best_period = float(result.period[best_idx])
    best_duration = float(result.duration[best_idx])
    best_t0 = float(result.transit_time[best_idx])
    best_depth = float(result.depth[best_idx])
    best_snr = float(result.depth_snr[best_idx])
    period_error_percent = abs(best_period - target.period_days) / target.period_days * 100.0

    plot_periodogram(result.period, result.power, target.period_days, dirs["figures"] / "toi1258_allsector_bls_periodogram.png")
    plot_folded_lightcurve(
        time,
        flux_rel[finite],
        target.period_days,
        target.t0_btjd,
        target.duration_days,
        dirs["figures"] / "toi1258_allsector_folded_catalog_ephemeris.png",
        title="TOI-1258.01 all-sector stitch folded on TOI catalog ephemeris",
    )
    plot_folded_lightcurve(
        time,
        flux_rel[finite],
        best_period,
        best_t0,
        best_duration,
        dirs["figures"] / "toi1258_allsector_folded_bls_recovery.png",
        title="TOI-1258.01 all-sector stitch folded on recovered BLS ephemeris",
    )
    plot_time_series(
        stitched["time"],
        stitched["flux_rel"],
        stitched["sector"],
        dirs["figures"] / "toi1258_allsector_stitched_time_series.png",
    )

    return {
        "catalog_period_days": f"{target.period_days:.9f}",
        "best_bls_period_days": f"{best_period:.9f}",
        "best_bls_t0_btjd": f"{best_t0:.6f}",
        "best_bls_duration_days": f"{best_duration:.6f}",
        "best_bls_depth_ppm": f"{best_depth * 1e6:.3f}",
        "best_bls_depth_snr": f"{best_snr:.3f}",
        "period_error_percent": f"{period_error_percent:.6f}",
        "n_cadences": len(time),
        "n_sectors": len(set(map(int, stitched["sector"]))),
        "period_grid_min_days": period_min,
        "period_grid_max_days": period_max,
        "n_periods": n_periods,
    }


def compute_transit_metrics(
    time: np.ndarray,
    flux_rel: np.ndarray,
    period: float,
    t0: float,
    duration: float,
    sector: np.ndarray,
) -> dict[str, object]:
    primary = transit_mask(time, period, t0, duration * 1.20)
    out = out_of_transit_mask(time, period, t0, duration)
    secondary_t0 = t0 + 0.5 * period
    secondary = transit_mask(time, period, secondary_t0, duration * 1.20)

    summary = depth_summary(time, flux_rel, primary, out, prefix="primary")
    summary.update(depth_summary(time, flux_rel, secondary, out, prefix="secondary"))
    transit_number = np.rint((time - t0) / period).astype(int)
    odd_mask = primary & (transit_number % 2 != 0)
    even_mask = primary & (transit_number % 2 == 0)
    odd = depth_summary(time, flux_rel, odd_mask, out, prefix="odd")
    even = depth_summary(time, flux_rel, even_mask, out, prefix="even")
    summary.update(odd)
    summary.update(even)
    odd_depth = to_float(summary.get("odd_depth_ppm")) or np.nan
    even_depth = to_float(summary.get("even_depth_ppm")) or np.nan
    summary["odd_even_depth_delta_ppm"] = f"{abs(odd_depth - even_depth):.3f}" if np.isfinite(odd_depth + even_depth) else ""
    summary["catalog_period_days"] = f"{period:.9f}"
    summary["catalog_t0_btjd"] = f"{t0:.6f}"
    summary["catalog_duration_days"] = f"{duration:.6f}"

    sector_rows = []
    for sec in sorted(set(map(int, sector))):
        sec_mask = sector == sec
        row = {"sector": sec}
        row.update(depth_summary(time[sec_mask], flux_rel[sec_mask], primary[sec_mask], out[sec_mask], prefix="primary"))
        row["transit_cadence_count"] = int(np.sum(primary[sec_mask]))
        row["oot_cadence_count"] = int(np.sum(out[sec_mask]))
        sector_rows.append(row)
    return {"summary": summary, "sector_rows": sector_rows}


def run_centroid_gate(target: Target, tp_paths: list[Path], dirs: dict[str, Path]) -> dict[str, object]:
    sector_rows = []
    image_panels = []

    for path in tp_paths:
        with fits.open(path, memmap=True) as hdul:
            primary_header = hdul[0].header
            data = hdul[1].data
            aperture = np.asarray(hdul[2].data, dtype=float) if len(hdul) > 2 and hdul[2].data is not None else None
            sector = int(primary_header.get("SECTOR", -1))
            time = np.asarray(data["TIME"], dtype=float)
            flux_cube = np.asarray(data["FLUX"], dtype=float)
            quality_name = "QUALITY" if "QUALITY" in data.names else "SAP_QUALITY"
            quality = np.asarray(data[quality_name], dtype=int)

        valid = np.isfinite(time) & (quality == 0) & np.all(np.isfinite(flux_cube), axis=(1, 2))
        if np.sum(valid) < 30:
            continue
        time = time[valid]
        cube = flux_cube[valid]

        in_mask = transit_mask(time, target.period_days, target.t0_btjd, target.duration_days * 1.25)
        out_mask = out_of_transit_mask(time, target.period_days, target.t0_btjd, target.duration_days)
        if np.sum(in_mask) < 8 or np.sum(out_mask) < 30:
            sector_rows.append(
                {
                    "sector": sector,
                    "tpf_path": str(path),
                    "status": "SKIPPED_INSUFFICIENT_TRANSIT_CADENCES",
                    "in_cadences": int(np.sum(in_mask)),
                    "out_cadences": int(np.sum(out_mask)),
                }
            )
            continue

        aperture_mask = aperture_to_bool(aperture, cube.shape[1:])
        xgrid, ygrid = xy_grids(cube.shape[1:])

        all_centroids = frame_centroids(cube, aperture_mask, xgrid, ygrid)
        in_x = all_centroids["x"][in_mask]
        in_y = all_centroids["y"][in_mask]
        out_x = all_centroids["x"][out_mask]
        out_y = all_centroids["y"][out_mask]
        in_x_mean, in_x_sem = mean_sem(in_x)
        in_y_mean, in_y_sem = mean_sem(in_y)
        out_x_mean, out_x_sem = mean_sem(out_x)
        out_y_mean, out_y_sem = mean_sem(out_y)
        dx = in_x_mean - out_x_mean
        dy = in_y_mean - out_y_mean
        shift_pix = math.hypot(dx, dy)
        shift_err_pix = math.sqrt(in_x_sem**2 + in_y_sem**2 + out_x_sem**2 + out_y_sem**2)
        shift_sig = shift_pix / shift_err_pix if shift_err_pix > 0 else np.nan

        in_img = np.nanmedian(cube[in_mask], axis=0)
        out_img = np.nanmedian(cube[out_mask], axis=0)
        diff_img = out_img - in_img
        target_centroid = image_centroid(out_img, aperture_mask, xgrid, ygrid, positive=True)
        diff_centroid = image_centroid(diff_img, np.ones_like(aperture_mask, dtype=bool), xgrid, ygrid, positive=True)
        diff_dx = diff_centroid[0] - target_centroid[0]
        diff_dy = diff_centroid[1] - target_centroid[1]
        diff_offset_pix = math.hypot(diff_dx, diff_dy)
        diff_peak_snr = image_peak_snr(diff_img)

        row = {
            "sector": sector,
            "tpf_path": str(path),
            "status": "OK",
            "in_cadences": int(np.sum(in_mask)),
            "out_cadences": int(np.sum(out_mask)),
            "cadence_centroid_dx_pix": f"{dx:.6f}",
            "cadence_centroid_dy_pix": f"{dy:.6f}",
            "cadence_centroid_shift_pix": f"{shift_pix:.6f}",
            "cadence_centroid_shift_arcsec": f"{shift_pix * TESS_PIXEL_SCALE_ARCSEC:.4f}",
            "cadence_centroid_shift_sigma": f"{shift_sig:.3f}" if np.isfinite(shift_sig) else "",
            "difference_image_target_x_pix": f"{target_centroid[0]:.6f}",
            "difference_image_target_y_pix": f"{target_centroid[1]:.6f}",
            "difference_image_deficit_x_pix": f"{diff_centroid[0]:.6f}",
            "difference_image_deficit_y_pix": f"{diff_centroid[1]:.6f}",
            "difference_image_offset_pix": f"{diff_offset_pix:.6f}",
            "difference_image_offset_arcsec": f"{diff_offset_pix * TESS_PIXEL_SCALE_ARCSEC:.4f}",
            "difference_image_peak_snr": f"{diff_peak_snr:.3f}",
        }
        sector_rows.append(row)
        image_panels.append((sector, out_img, diff_img, aperture_mask, target_centroid, diff_centroid, row))

    ok_rows = [row for row in sector_rows if row.get("status") == "OK"]
    summary = summarize_centroids(ok_rows)
    if ok_rows:
        plot_centroid_vectors(ok_rows, dirs["figures"] / "toi1258_centroid_shift_vectors.png")
        plot_difference_offsets(ok_rows, dirs["figures"] / "toi1258_difference_image_offsets.png")
    if image_panels:
        best_panel = max(image_panels, key=lambda item: to_float(item[6].get("difference_image_peak_snr")) or -np.inf)
        plot_tpf_difference_panel(best_panel, dirs["figures"] / "toi1258_representative_tpf_difference_image.png")

    return {"sector_rows": sector_rows, "summary": summary}


def summarize_centroids(rows: list[dict[str, object]]) -> dict[str, object]:
    if not rows:
        return {
            "sectors_used": 0,
            "centroid_gate_verdict": "NO_USABLE_TPF_CENTROID_RESULT",
            "median_abs_cadence_centroid_shift_arcsec": "",
            "median_abs_difference_image_offset_arcsec": "",
            "max_difference_image_offset_arcsec": "",
        }
    cadence_arcsec = np.array([to_float(row["cadence_centroid_shift_arcsec"]) for row in rows], dtype=float)
    diff_arcsec = np.array([to_float(row["difference_image_offset_arcsec"]) for row in rows], dtype=float)
    diff_snr = np.array([to_float(row["difference_image_peak_snr"]) for row in rows], dtype=float)
    med_cad = float(np.nanmedian(np.abs(cadence_arcsec)))
    med_diff = float(np.nanmedian(np.abs(diff_arcsec)))
    max_diff = float(np.nanmax(np.abs(diff_arcsec)))

    if med_diff <= 7.0 and max_diff <= 21.0:
        verdict = "PASSES_TESS_PIXEL_CENTROID_GATE"
    elif med_diff <= 14.0:
        verdict = "BORDERLINE_TESS_PIXEL_CENTROID_GATE"
    else:
        verdict = "CENTROID_OFFSET_REQUIRES_FOLLOWUP"

    return {
        "sectors_used": len(rows),
        "centroid_gate_verdict": verdict,
        "median_abs_cadence_centroid_shift_arcsec": f"{med_cad:.4f}",
        "median_abs_difference_image_offset_arcsec": f"{med_diff:.4f}",
        "max_difference_image_offset_arcsec": f"{max_diff:.4f}",
        "median_difference_image_peak_snr": f"{float(np.nanmedian(diff_snr)):.3f}",
        "tess_pixel_scale_arcsec_per_pix": TESS_PIXEL_SCALE_ARCSEC,
        "interpretation_boundary": "TESS centroid consistency reduces background-source risk but cannot localize below the TESS pixel scale.",
    }


def transit_mask(time: np.ndarray, period: float, t0: float, duration: float) -> np.ndarray:
    phase = phase_days(time, period, t0)
    return np.abs(phase) <= 0.5 * duration


def out_of_transit_mask(time: np.ndarray, period: float, t0: float, duration: float) -> np.ndarray:
    phase = np.abs(phase_days(time, period, t0))
    return (phase >= 2.5 * duration) & (phase <= 0.45 * period)


def phase_days(time: np.ndarray, period: float, t0: float) -> np.ndarray:
    return ((time - t0 + 0.5 * period) % period) - 0.5 * period


def depth_summary(time: np.ndarray, flux_rel: np.ndarray, in_mask: np.ndarray, out_mask: np.ndarray, prefix: str) -> dict[str, object]:
    row = {
        f"{prefix}_in_cadences": int(np.sum(in_mask)),
        f"{prefix}_oot_cadences": int(np.sum(out_mask)),
    }
    if np.sum(in_mask) < 5 or np.sum(out_mask) < 20:
        row.update({f"{prefix}_depth_ppm": "", f"{prefix}_depth_snr": ""})
        return row
    in_flux = flux_rel[in_mask]
    out_flux = flux_rel[out_mask]
    depth = (np.nanmedian(out_flux) - np.nanmedian(in_flux)) * 1e6
    err = math.sqrt((robust_std(in_flux * 1e6) ** 2) / len(in_flux) + (robust_std(out_flux * 1e6) ** 2) / len(out_flux))
    snr = depth / err if err > 0 else np.nan
    row[f"{prefix}_depth_ppm"] = f"{depth:.3f}"
    row[f"{prefix}_depth_snr"] = f"{snr:.3f}" if np.isfinite(snr) else ""
    return row


def frame_centroids(cube: np.ndarray, mask: np.ndarray, xgrid: np.ndarray, ygrid: np.ndarray) -> dict[str, np.ndarray]:
    bg = np.nanmedian(cube[:, border_mask(cube.shape[1:])], axis=1)
    net = cube - bg[:, None, None]
    net = np.where(mask[None, :, :], net, 0.0)
    net = np.where(net > 0, net, 0.0)
    total = np.nansum(net, axis=(1, 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.nansum(net * xgrid[None, :, :], axis=(1, 2)) / total
        y = np.nansum(net * ygrid[None, :, :], axis=(1, 2)) / total
    return {"x": x, "y": y, "total": total}


def image_centroid(image: np.ndarray, mask: np.ndarray, xgrid: np.ndarray, ygrid: np.ndarray, positive: bool = True) -> tuple[float, float]:
    bg = np.nanmedian(image[border_mask(image.shape)])
    net = image - bg
    net = np.where(mask, net, 0.0)
    if positive:
        net = np.where(net > 0, net, 0.0)
    total = float(np.nansum(net))
    if not np.isfinite(total) or total <= 0:
        return (float("nan"), float("nan"))
    x = float(np.nansum(net * xgrid) / total)
    y = float(np.nansum(net * ygrid) / total)
    return (x, y)


def image_peak_snr(image: np.ndarray) -> float:
    bg = np.nanmedian(image[border_mask(image.shape)])
    net = image - bg
    sigma = robust_std(net[border_mask(image.shape)])
    if sigma <= 0 or not np.isfinite(sigma):
        sigma = robust_std(net.ravel())
    return float(np.nanmax(net) / sigma) if sigma > 0 else float("nan")


def aperture_to_bool(aperture: np.ndarray | None, shape: tuple[int, int]) -> np.ndarray:
    if aperture is None or aperture.shape != shape:
        return np.ones(shape, dtype=bool)
    mask = np.asarray(aperture) > 0
    if np.sum(mask) < 3:
        return np.ones(shape, dtype=bool)
    return mask


def xy_grids(shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    y, x = np.indices(shape)
    return x.astype(float), y.astype(float)


def border_mask(shape: tuple[int, int]) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    mask[0, :] = True
    mask[-1, :] = True
    mask[:, 0] = True
    mask[:, -1] = True
    return mask


def mean_sem(values: np.ndarray) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan"), float("nan")
    mean = float(np.nanmean(values))
    if len(values) < 2:
        return mean, float("nan")
    sem = float(np.nanstd(values, ddof=1) / math.sqrt(len(values)))
    return mean, sem


def plot_periodogram(periods: np.ndarray, power: np.ndarray, catalog_period: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(periods, power, color="#1f77b4", lw=1.0)
    ax.axvline(catalog_period, color="#d62728", lw=1.5, ls="--", label=f"Catalog {catalog_period:.6f} d")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("BLS power / SNR objective")
    ax.set_title("All-sector stitched BLS periodogram")
    ax.legend(loc="best")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_folded_lightcurve(time: np.ndarray, flux_rel: np.ndarray, period: float, t0: float, duration: float, path: Path, title: str) -> None:
    phase = phase_days(time, period, t0)
    phase_hours = phase * 24.0
    flux_ppm = (flux_rel - 1.0) * 1e6
    order = np.argsort(phase_hours)
    phase_hours = phase_hours[order]
    flux_ppm = flux_ppm[order]
    bins = np.linspace(-0.65 * period * 24.0, 0.65 * period * 24.0, 260)
    centers = 0.5 * (bins[:-1] + bins[1:])
    binned = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (phase_hours >= lo) & (phase_hours < hi)
        binned.append(np.nanmedian(flux_ppm[m]) if np.sum(m) else np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    sample = deterministic_sample(len(phase_hours), 35000)
    ax.scatter(phase_hours[sample], flux_ppm[sample], s=2, alpha=0.10, color="#555555", rasterized=True)
    ax.plot(centers, binned, color="#d62728", lw=1.8, label="Phase-bin median")
    ax.axvspan(-0.5 * duration * 24.0, 0.5 * duration * 24.0, color="#1f77b4", alpha=0.16, label="Transit window")
    ax.set_xlim(-2.5 * duration * 24.0, 2.5 * duration * 24.0)
    ax.set_xlabel("Phase from transit center (hours)")
    ax.set_ylabel("Relative flux (ppm)")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_time_series(time: np.ndarray, flux_rel: np.ndarray, sector: np.ndarray, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    sample = deterministic_sample(len(time), 60000)
    sc = ax.scatter(time[sample], (flux_rel[sample] - 1.0) * 1e6, c=sector[sample], s=2, alpha=0.28, cmap="viridis", rasterized=True)
    fig.colorbar(sc, ax=ax, label="TESS sector")
    ax.set_xlabel("Time (BTJD)")
    ax.set_ylabel("Relative flux (ppm)")
    ax.set_title("TOI-1258.01 stitched TESS light curve")
    ax.grid(alpha=0.20)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_centroid_vectors(rows: list[dict[str, object]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    for row in rows:
        dx = to_float(row["cadence_centroid_dx_pix"]) or 0.0
        dy = to_float(row["cadence_centroid_dy_pix"]) or 0.0
        ax.arrow(0, 0, dx * TESS_PIXEL_SCALE_ARCSEC, dy * TESS_PIXEL_SCALE_ARCSEC, head_width=0.08, length_includes_head=True, alpha=0.75)
        ax.text(dx * TESS_PIXEL_SCALE_ARCSEC, dy * TESS_PIXEL_SCALE_ARCSEC, str(row["sector"]), fontsize=8)
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("In-transit minus out-of-transit centroid dx (arcsec)")
    ax.set_ylabel("In-transit minus out-of-transit centroid dy (arcsec)")
    ax.set_title("TESS cadence-level centroid shifts by sector")
    ax.grid(alpha=0.25)
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_difference_offsets(rows: list[dict[str, object]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.scatter([0], [0], marker="*", s=160, color="#111111", label="Out-of-transit flux centroid")
    for row in rows:
        dx = (to_float(row["difference_image_deficit_x_pix"]) or np.nan) - (to_float(row["difference_image_target_x_pix"]) or np.nan)
        dy = (to_float(row["difference_image_deficit_y_pix"]) or np.nan) - (to_float(row["difference_image_target_y_pix"]) or np.nan)
        ax.scatter(dx * TESS_PIXEL_SCALE_ARCSEC, dy * TESS_PIXEL_SCALE_ARCSEC, s=36)
        ax.text(dx * TESS_PIXEL_SCALE_ARCSEC, dy * TESS_PIXEL_SCALE_ARCSEC, str(row["sector"]), fontsize=8)
    circle = plt.Circle((0, 0), TESS_PIXEL_SCALE_ARCSEC, fill=False, color="#d62728", ls="--", lw=1.1, label="1 TESS pixel")
    ax.add_patch(circle)
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Transit-deficit centroid offset dx (arcsec)")
    ax.set_ylabel("Transit-deficit centroid offset dy (arcsec)")
    ax.set_title("Difference-image centroid offsets by sector")
    ax.legend(loc="best")
    ax.grid(alpha=0.25)
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_tpf_difference_panel(panel, path: Path) -> None:
    sector, out_img, diff_img, aperture_mask, target_centroid, diff_centroid, row = panel
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    im0 = axes[0].imshow(out_img, origin="lower", cmap="viridis")
    axes[0].contour(aperture_mask, levels=[0.5], colors="white", linewidths=0.8)
    axes[0].scatter([target_centroid[0]], [target_centroid[1]], marker="+", s=140, color="red")
    axes[0].set_title(f"Sector {sector} out-of-transit image")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(diff_img, origin="lower", cmap="magma")
    axes[1].scatter([target_centroid[0]], [target_centroid[1]], marker="+", s=140, color="cyan", label="Target centroid")
    axes[1].scatter([diff_centroid[0]], [diff_centroid[1]], marker="x", s=110, color="white", label="Transit deficit")
    axes[1].set_title("Out minus in-transit image")
    axes[1].legend(loc="best", fontsize=8)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    resid = diff_img - np.nanmedian(diff_img[border_mask(diff_img.shape)])
    im2 = axes[2].imshow(resid, origin="lower", cmap="coolwarm")
    axes[2].scatter([target_centroid[0]], [target_centroid[1]], marker="+", s=140, color="black")
    axes[2].scatter([diff_centroid[0]], [diff_centroid[1]], marker="x", s=110, color="white")
    axes[2].set_title(f"Deficit peak SNR {row.get('difference_image_peak_snr')}")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        ax.set_xlabel("TPF column index")
        ax.set_ylabel("TPF row index")
    fig.suptitle("Representative TESS target-pixel centroid gate")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def deterministic_sample(n: int, max_n: int) -> np.ndarray:
    if n <= max_n:
        return np.arange(n)
    step = max(1, int(math.ceil(n / max_n)))
    return np.arange(0, n, step)[:max_n]


def downsample_rows(rows: list[dict[str, object]], max_rows: int) -> list[dict[str, object]]:
    if len(rows) <= max_rows:
        return rows
    idx = deterministic_sample(len(rows), max_rows)
    return [rows[int(i)] for i in idx]


def render_report(
    target: Target,
    products: dict[str, object],
    lc_paths: list[Path],
    tp_paths: list[Path],
    bls_summary: dict[str, object],
    transit_metrics: dict[str, object],
    centroid_result: dict[str, object],
    out_dir: Path,
) -> str:
    centroid_summary = centroid_result["summary"]
    transit_summary = transit_metrics["summary"]
    lines = [
        "# TOI-1258.01 Centroid and Stitch Gate",
        "",
        f"Created: {date.today().isoformat()}",
        "",
        "## Target",
        "",
        f"- TOI: `{target.toi}`",
        f"- TIC: `{target.tic}`",
        f"- RA/Dec: `{target.ra:.7f}`, `{target.dec:.7f}` deg",
        f"- TOI disposition: `{target.disposition}`",
        f"- Catalog period: `{target.period_days:.9f}` d",
        f"- Catalog epoch: `{target.t0_btjd:.6f}` BTJD",
        f"- Catalog duration: `{target.duration_days:.6f}` d",
        f"- Catalog depth: `{target.depth_ppm:.3f}` ppm",
        "",
        "## Product Inventory",
        "",
        f"- Standard TESS light-curve products used: `{len(lc_paths)}`",
        f"- Standard TESS target-pixel products used: `{len(tp_paths)}`",
        f"- Product inventory table: `evidence/tables/mast_tess_product_inventory.csv`",
        "",
        "## Stitch Gate",
        "",
        f"- All-sector cadences used in BLS: `{bls_summary['n_cadences']}`",
        f"- TESS sectors stitched: `{bls_summary['n_sectors']}`",
        f"- Recovered BLS period: `{bls_summary['best_bls_period_days']}` d",
        f"- Catalog period: `{bls_summary['catalog_period_days']}` d",
        f"- Period error: `{bls_summary['period_error_percent']}` percent",
        f"- Recovered BLS depth: `{bls_summary['best_bls_depth_ppm']}` ppm",
        f"- Recovered BLS depth SNR: `{bls_summary['best_bls_depth_snr']}`",
        "",
        "Direct fold on the catalog ephemeris:",
        "",
        f"- Primary depth: `{transit_summary.get('primary_depth_ppm', '')}` ppm",
        f"- Primary depth SNR: `{transit_summary.get('primary_depth_snr', '')}`",
        f"- Secondary depth: `{transit_summary.get('secondary_depth_ppm', '')}` ppm",
        f"- Secondary depth SNR: `{transit_summary.get('secondary_depth_snr', '')}`",
        f"- Odd/even depth delta: `{transit_summary.get('odd_even_depth_delta_ppm', '')}` ppm",
        "",
        "Interpretation: a coherent signal at the same period across the stitched sectors supports a persistent transit-like event on this TOI ephemeris. It does not by itself confirm a planet.",
        "",
        "## Centroid Gate",
        "",
        f"- Sectors with usable TPF centroid measurements: `{centroid_summary.get('sectors_used')}`",
        f"- Median cadence-level centroid shift: `{centroid_summary.get('median_abs_cadence_centroid_shift_arcsec')}` arcsec",
        f"- Median difference-image source offset: `{centroid_summary.get('median_abs_difference_image_offset_arcsec')}` arcsec",
        f"- Maximum difference-image source offset: `{centroid_summary.get('max_difference_image_offset_arcsec')}` arcsec",
        f"- Median difference-image peak SNR: `{centroid_summary.get('median_difference_image_peak_snr')}`",
        f"- Centroid verdict: `{centroid_summary.get('centroid_gate_verdict')}`",
        "",
        "Interpretation boundary: TESS pixels are about 21 arcsec wide. A centroid pass reduces the risk that the dip is a distant background eclipsing binary, but it cannot prove the source is exactly on the target at sub-pixel or high-resolution-imaging scale.",
        "",
        "## Figures",
        "",
        "- `evidence/figures/toi1258_allsector_bls_periodogram.png`",
        "- `evidence/figures/toi1258_allsector_folded_catalog_ephemeris.png`",
        "- `evidence/figures/toi1258_allsector_folded_bls_recovery.png`",
        "- `evidence/figures/toi1258_allsector_stitched_time_series.png`",
        "- `evidence/figures/toi1258_centroid_shift_vectors.png`",
        "- `evidence/figures/toi1258_difference_image_offsets.png`",
        "- `evidence/figures/toi1258_representative_tpf_difference_image.png`",
        "",
        "## Evidence Boundary",
        "",
        "This package tests two major false-positive gates: whether the signal survives all-sector stitching, and whether TESS target-pixel centroids shift during transit. A clean result here upgrades TOI-1258.01 as a candidate worth prioritizing. Confirmation still requires independent follow-up such as high-resolution imaging, radial velocities, or official validation/confirmation in the exoplanet literature.",
        "",
        "## Sources",
        "",
        "- NASA Exoplanet Archive TOI table via TAP: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html",
        "- MAST TESS products via astroquery.mast: https://astroquery.readthedocs.io/en/latest/mast/mast.html",
        "- TESS mission data products at MAST: https://archive.stsci.edu/tess/",
        "",
        f"Package root: `{out_dir}`",
    ]
    return "\n".join(lines) + "\n"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def tap_csv(query: str) -> list[dict[str, str]]:
    url = TAP_SYNC + "?" + urllib.parse.urlencode({"query": query, "format": "csv"})
    with urllib.request.urlopen(url, timeout=120) as response:
        text = response.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text)))


def to_float(value) -> float | None:
    try:
        if value is None or value == "":
            return None
        out = float(value)
        return out if np.isfinite(out) else None
    except Exception:
        return None


def robust_std(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    med = np.nanmedian(arr)
    mad = np.nanmedian(np.abs(arr - med))
    if mad > 0 and np.isfinite(mad):
        return float(1.4826 * mad)
    return float(np.nanstd(arr))


if __name__ == "__main__":
    raise SystemExit(main())
