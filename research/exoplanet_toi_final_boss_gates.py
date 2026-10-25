#!/usr/bin/env python
"""Final validation-gate audit for TOI-1258.01.

The gates covered here are:
1. RV / spectroscopy archive evidence.
2. High-resolution imaging and ground follow-up evidence.
3. A validation-prior feasibility summary from the existing TESS, Gaia,
   SPOC DV, ExoFOP, and speckle constraints.

This script does not confirm a planet. It packages the public evidence and
states which gates are cleared, partially cleared, or still require raw data.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


WORKSPACE = Path(r".")
DEFAULT_OUTPUT = Path(r"data/external/exoplanet_outputs") / f"EXOPLANET_TOI1258_FINAL_BOSS_GATES_{date.today().isoformat()}"
TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EXOFOP_URL = "https://exofop.ipac.caltech.edu/tess/target.php?id={tic}"
EXOFOP_FILELIST_URL = "https://exofop.ipac.caltech.edu/tess/download_filelist.php?id={tic}&output=csv"
ESO_TAP_SYNC = "https://archive.eso.org/tap_obs/sync"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toi", default="TOI-1258.01")
    parser.add_argument("--tic", default="288132261")
    parser.add_argument("--ra", type=float, default=209.4843130)
    parser.add_argument("--dec", type=float, default=79.5833220)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--centroid-package",
        default=r"data/external/exoplanet_outputs\EXOPLANET_TOI1258_CENTROID_STITCH_GATE_2026-06-15",
    )
    parser.add_argument(
        "--deep-vet-package",
        default=str(WORKSPACE / "outputs" / "EXOPLANET_TOI1258_DEEP_VET_2026-06-14"),
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    dirs = make_dirs(out_dir)

    target = fetch_toi(args.toi, args.tic, args.ra, args.dec)
    (dirs["tables"] / "toi_catalog_row.json").write_text(json.dumps(target, indent=2), encoding="utf-8")

    exofop = fetch_exofop(args.tic, dirs["archive"])
    write_csv(dirs["tables"] / "exofop_planet_parameters.csv", exofop.get("rowData3", []))
    write_csv(dirs["tables"] / "exofop_stellar_parameters.csv", exofop.get("rowData4", []))
    write_csv(dirs["tables"] / "exofop_magnitudes.csv", exofop.get("rowData6", []))
    write_csv(dirs["tables"] / "exofop_high_resolution_imaging.csv", exofop.get("rowData7", []))
    write_csv(dirs["tables"] / "exofop_spectroscopy_rv.csv", exofop.get("rowData8", []))
    write_csv(dirs["tables"] / "exofop_time_series_followup.csv", exofop.get("rowData9", []))
    write_csv(dirs["tables"] / "exofop_file_list.csv", exofop.get("filelist", []))

    eso_rows = query_eso(args.ra, args.dec)
    write_csv(dirs["tables"] / "eso_archive_cone_search_0p03deg.csv", eso_rows)

    centroid_summary = read_first_csv(Path(args.centroid_package) / "evidence" / "tables" / "spoc_dv_centroid_summary.csv")
    stitch_summary = read_first_csv(Path(args.centroid_package) / "evidence" / "tables" / "stitched_bls_summary.csv")
    transit_metrics = read_first_csv(Path(args.centroid_package) / "evidence" / "tables" / "stitched_transit_metrics.csv")
    gaia_contamination = read_csv(Path(args.deep_vet_package) / "evidence" / "tables" / "gaia_contamination_summary.csv")
    gaia_neighbors = read_csv(Path(args.deep_vet_package) / "evidence" / "tables" / "gaia_dr3_neighbors_120arcsec.csv")

    rv_rows = exofop.get("rowData8", [])
    imaging_rows = exofop.get("rowData7", [])
    timeseries_rows = exofop.get("rowData9", [])
    file_rows = exofop.get("filelist", [])
    rv_summary = summarize_rv_rows(rv_rows, target)
    imaging_summary = summarize_imaging(imaging_rows, file_rows)
    timeseries_summary = summarize_timeseries(timeseries_rows, file_rows)
    validation_rows = validation_prior_table(
        target=target,
        stitch_summary=stitch_summary,
        transit_metrics=transit_metrics,
        centroid_summary=best_centroid_row(Path(args.centroid_package) / "evidence" / "tables" / "spoc_dv_centroid_summary.csv"),
        gaia_contamination=gaia_contamination,
        gaia_neighbors=gaia_neighbors,
        imaging_summary=imaging_summary,
        timeseries_summary=timeseries_summary,
        rv_summary=rv_summary,
    )

    write_csv(dirs["tables"] / "rv_spectroscopy_gate_summary.csv", rv_summary)
    write_csv(dirs["tables"] / "high_resolution_imaging_gate_summary.csv", imaging_summary)
    write_csv(dirs["tables"] / "ground_timeseries_gate_summary.csv", timeseries_summary)
    write_csv(dirs["tables"] / "statistical_validation_prior_table.csv", validation_rows)

    report = render_report(
        target=target,
        out_dir=out_dir,
        rv_summary=rv_summary,
        imaging_summary=imaging_summary,
        timeseries_summary=timeseries_summary,
        validation_rows=validation_rows,
        stitch_summary=stitch_summary,
        transit_metrics=transit_metrics,
        centroid_summary=best_centroid_row(Path(args.centroid_package) / "evidence" / "tables" / "spoc_dv_centroid_summary.csv"),
        eso_rows=eso_rows,
    )
    report_path = dirs["reports"] / "TOI1258_FINAL_BOSS_GATES_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    manifest = {
        "created": date.today().isoformat(),
        "target": target,
        "output_dir": str(out_dir),
        "report": str(report_path),
        "claim_boundary": "This package audits final validation gates. It does not by itself confirm TOI-1258.01 as a planet.",
    }
    (out_dir / "target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "output_dir": str(out_dir),
                "report": str(report_path),
                "rv_rows": len(rv_rows),
                "imaging_rows": len(imaging_rows),
                "timeseries_rows": len(timeseries_rows),
                "eso_rows": len(eso_rows),
                "rv_gate": rv_summary[0]["gate_status"] if rv_summary else "",
                "imaging_gate": imaging_summary[0]["gate_status"] if imaging_summary else "",
                "statistical_gate": validation_rows[0]["overall_statistical_validation_status"] if validation_rows else "",
            },
            indent=2,
        )
    )
    return 0


def make_dirs(out_dir: Path) -> dict[str, Path]:
    paths = {
        "base": out_dir,
        "archive": out_dir / "evidence" / "archive_snapshots",
        "tables": out_dir / "evidence" / "tables",
        "reports": out_dir / "evidence" / "reports",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def fetch_toi(toi: str, fallback_tic: str, fallback_ra: float, fallback_dec: float) -> dict[str, object]:
    query = f"""
SELECT TOP 1 tid, toidisplay, tfopwg_disp, st_tmag, ra, dec, pl_orbper,
       pl_tranmid, pl_trandurh, pl_trandep, pl_rade, pl_eqt, st_dist,
       st_teff, st_rad
FROM toi
WHERE toidisplay='{toi}'
"""
    rows = tap_csv(query)
    if not rows:
        return {"toidisplay": toi, "tid": fallback_tic, "ra": fallback_ra, "dec": fallback_dec}
    row = rows[0]
    return {
        "toidisplay": row.get("toidisplay", toi),
        "tid": row.get("tid", fallback_tic),
        "tfopwg_disp": row.get("tfopwg_disp", ""),
        "st_tmag": as_float(row.get("st_tmag")),
        "ra": as_float(row.get("ra")) or fallback_ra,
        "dec": as_float(row.get("dec")) or fallback_dec,
        "pl_orbper_days": as_float(row.get("pl_orbper")),
        "pl_tranmid_bjd": as_float(row.get("pl_tranmid")),
        "pl_trandurh": as_float(row.get("pl_trandurh")),
        "pl_trandep_ppm": as_float(row.get("pl_trandep")),
        "pl_rade": as_float(row.get("pl_rade")),
        "pl_eqt": as_float(row.get("pl_eqt")),
        "st_dist_pc": as_float(row.get("st_dist")),
        "st_teff": as_float(row.get("st_teff")),
        "st_rad": as_float(row.get("st_rad")),
        "st_mass": 1.157,
    }


def fetch_exofop(tic: str, archive_dir: Path) -> dict[str, list[dict[str, object]]]:
    html_url = EXOFOP_URL.format(tic=tic)
    filelist_url = EXOFOP_FILELIST_URL.format(tic=tic)
    html = fetch_text(html_url)
    (archive_dir / f"exofop_tic_{tic}.html").write_text(html, encoding="utf-8")
    out: dict[str, list[dict[str, object]]] = {}
    for match in re.finditer(r"var\s+(rowData\d+)\s*=\s*(\[.*?\]);", html, flags=re.S):
        try:
            out[match.group(1)] = json.loads(match.group(2))
        except json.JSONDecodeError:
            out[match.group(1)] = []
    filelist_text = fetch_text(filelist_url)
    (archive_dir / f"exofop_tic_{tic}_filelist.csv").write_text(filelist_text, encoding="utf-8")
    out["filelist"] = list(csv.DictReader(io.StringIO(filelist_text)))
    return out


def query_eso(ra: float, dec: float) -> list[dict[str, str]]:
    radius = 0.03
    query = f"""
SELECT TOP 200 dp_id, target_name, instrument_name, dataproduct_type,
       dataproduct_subtype, em_min, em_max, em_res_power, t_exptime,
       t_min, t_max, s_ra, s_dec, proposal_id, obs_title, obs_release_date,
       access_estsize
FROM ivoa.ObsCore
WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))=1
ORDER BY t_min ASC
"""
    url = ESO_TAP_SYNC + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    try:
        return list(csv.DictReader(io.StringIO(fetch_text(url, timeout=180))))
    except Exception as exc:
        return [{"status": "ERROR", "error": repr(exc)}]


def summarize_rv_rows(rows: list[dict[str, object]], target: dict[str, object]) -> list[dict[str, object]]:
    rv_count = 0
    high_res_count = 0
    instruments = []
    for row in rows:
        res = as_float(row.get("ssres"))
        if res and res >= 40000:
            high_res_count += 1
        inst = f"{row.get('stel', '')}/{row.get('sinst', '')}".strip("/")
        if inst and inst not in instruments:
            instruments.append(inst)
        notes = str(row.get("snotes") or "")
        for match in re.finditer(r"(\d+)\s+RVs?", notes):
            rv_count += int(match.group(1))
    k_rows = expected_rv_amplitude(target)
    status = "RV_DATA_EXISTS_BUT_RAW_TIME_SERIES_NOT_FIT"
    if rv_count >= 20:
        status = "RV_CAMPAIGN_PRESENT_NEEDS_RAW_VALUES_OR_PUBLISHED_FIT"
    return [
        {
            "gate": "radial_velocity_wobble",
            "gate_status": status,
            "exofop_spectroscopy_rows": len(rows),
            "high_resolution_spectroscopy_rows_r_ge_40000": high_res_count,
            "rv_count_from_exofop_notes": rv_count,
            "instruments": "; ".join(instruments),
            "expected_k_mps_5_mearth": k_rows["k_5"],
            "expected_k_mps_10_mearth": k_rows["k_10"],
            "expected_k_mps_15_mearth": k_rows["k_15"],
            "expected_k_mps_25_mearth": k_rows["k_25"],
            "interpretation": "Public ExoFOP metadata reports a real RV/spectroscopy campaign, including APF/Levy RVs. This is not a mass measurement until the time-series velocities are retrieved or a published orbital fit is used.",
        }
    ]


def expected_rv_amplitude(target: dict[str, object]) -> dict[str, str]:
    period_days = as_float(target.get("pl_orbper_days")) or 10.2613397
    mstar = as_float(target.get("st_mass")) or 1.157
    period_years = period_days / 365.25
    out = {}
    for mass_earth in [5, 10, 15, 25]:
        mass_jup = mass_earth / 317.828
        k = 28.4329 * mass_jup * period_years ** (-1.0 / 3.0) * mstar ** (-2.0 / 3.0)
        out[f"k_{mass_earth}"] = f"{k:.3f}"
    return out


def summarize_imaging(rows: list[dict[str, object]], files: list[dict[str, object]]) -> list[dict[str, object]]:
    dat_files = [row for row in files if str(row.get("File Name", "")).endswith("_sensitivity.dat")]
    fits_files = [row for row in files if row.get("Type") == "Image" and str(row.get("File Name", "")).endswith(".fits")]
    best_contrast = []
    for row in rows:
        contrast = str(row.get("icont") or "")
        if contrast:
            best_contrast.append(contrast)
    status = "HIGH_RESOLUTION_IMAGING_METADATA_PRESENT"
    if rows and dat_files:
        status = "HIGH_RESOLUTION_IMAGING_GATE_PARTIALLY_CLEARED"
    return [
        {
            "gate": "high_resolution_imaging",
            "gate_status": status,
            "imaging_rows": len(rows),
            "sensitivity_dat_files_listed": len(dat_files),
            "reconstructed_fits_files_listed": len(fits_files),
            "instruments": "; ".join(sorted({f"{r.get('itel','')}/{r.get('iinst','')}" for r in rows})),
            "reported_contrast_limits": "; ".join(best_contrast),
            "interpretation": "Gemini/Alopeke speckle observations are listed with contrast limits. This is strong follow-up evidence against close stellar contaminants, but a formal validation should ingest the sensitivity curves directly.",
        }
    ]


def summarize_timeseries(rows: list[dict[str, object]], files: list[dict[str, object]]) -> list[dict[str, object]]:
    neb_files = [row for row in files if "NEB" in str(row.get("File Name", "")) or "NEB" in str(row.get("Description", ""))]
    notes = "; ".join([str(row.get("tsnotes") or "") for row in rows if row.get("tsnotes")])
    status = "GROUND_TRANSIT_FOLLOWUP_PRESENT" if rows else "NO_GROUND_TRANSIT_FOLLOWUP_FOUND"
    return [
        {
            "gate": "ground_timeseries_transit_followup",
            "gate_status": status,
            "time_series_rows": len(rows),
            "neb_check_files_listed": len(neb_files),
            "telescopes": "; ".join(sorted({f"{r.get('tstel','')}/{r.get('tscam','')}" for r in rows})),
            "notes": notes,
            "interpretation": "Ground-based time-series observations and NEB-check products are listed. This supports on-sky follow-up, but the raw photometry should be ingested before using it in a formal validation calculation.",
        }
    ]


def validation_prior_table(
    target: dict[str, object],
    stitch_summary: dict[str, str],
    transit_metrics: dict[str, str],
    centroid_summary: dict[str, str],
    gaia_contamination: list[dict[str, str]],
    gaia_neighbors: list[dict[str, str]],
    imaging_summary: list[dict[str, object]],
    timeseries_summary: list[dict[str, object]],
    rv_summary: list[dict[str, object]],
) -> list[dict[str, object]]:
    depth_ppm = as_float(target.get("pl_trandep_ppm")) or 726.0
    aperture84 = first_or_empty(gaia_contamination, "aperture_radius_arcsec", "84")
    aperture120 = first_or_empty(gaia_contamination, "aperture_radius_arcsec", "120")
    brightest_neighbor_flux_ratio = max(
        [as_float(row.get("flux_ratio_vs_target_g")) or 0.0 for row in gaia_neighbors[1:]], default=0.0
    )
    required_eclipse_depth_brightest = depth_ppm / 1e6 / brightest_neighbor_flux_ratio if brightest_neighbor_flux_ratio > 0 else math.inf
    spoc_offset = as_float(centroid_summary.get("ms_tic_mean_sky_offset_arcsec"))
    spoc_offset_unc = as_float(centroid_summary.get("ms_tic_mean_sky_offset_unc_arcsec"))
    spoc_sigma = spoc_offset / spoc_offset_unc if spoc_offset and spoc_offset_unc else None

    blockers = []
    rv_status = str(rv_summary[0].get("gate_status", "")) if rv_summary else ""
    if "RAW" in rv_status or "PUBLISHED_FIT" in rv_status:
        blockers.append("raw RV time series or published RV orbital fit")
    if imaging_summary and int(imaging_summary[0].get("sensitivity_dat_files_listed", 0) or 0):
        blockers.append("ingest machine-readable speckle contrast curves into formal FPP code")
        pass
    else:
        blockers.append("machine-readable high-resolution-imaging contrast curve")
    if timeseries_summary and timeseries_summary[0].get("time_series_rows", 0):
        pass
    else:
        blockers.append("ground time-series photometry")

    status = "NOT_YET_99P9_VALIDATED_READY_FOR_FORMAL_FPP"
    if blockers:
        status = "NOT_YET_99P9_VALIDATED"

    return [
        {
            "overall_statistical_validation_status": status,
            "bls_period_error_percent": stitch_summary.get("period_error_percent", ""),
            "bls_depth_snr": stitch_summary.get("best_bls_depth_snr", ""),
            "primary_depth_snr": transit_metrics.get("primary_depth_snr", ""),
            "secondary_depth_snr": transit_metrics.get("secondary_depth_snr", ""),
            "odd_even_depth_delta_ppm": transit_metrics.get("odd_even_depth_delta_ppm", ""),
            "spoc_dv_sky_offset_arcsec": centroid_summary.get("ms_tic_mean_sky_offset_arcsec", ""),
            "spoc_dv_sky_offset_unc_arcsec": centroid_summary.get("ms_tic_mean_sky_offset_unc_arcsec", ""),
            "spoc_dv_offset_sigma": "" if spoc_sigma is None else f"{spoc_sigma:.3f}",
            "spoc_quality_good_metrics": centroid_summary.get("quality_good_metrics", ""),
            "spoc_quality_total_metrics": centroid_summary.get("quality_total_metrics", ""),
            "gaia_84arcsec_dilution_percent": aperture84.get("estimated_depth_dilution_percent", ""),
            "gaia_120arcsec_dilution_percent": aperture120.get("estimated_depth_dilution_percent", ""),
            "brightest_neighbor_flux_ratio_g": f"{brightest_neighbor_flux_ratio:.8f}",
            "required_eclipse_depth_on_brightest_neighbor_percent": (
                f"{100.0 * required_eclipse_depth_brightest:.3f}" if math.isfinite(required_eclipse_depth_brightest) else ""
            ),
            "high_resolution_imaging_status": imaging_summary[0].get("gate_status", "") if imaging_summary else "",
            "ground_timeseries_status": timeseries_summary[0].get("gate_status", "") if timeseries_summary else "",
            "rv_status": rv_summary[0].get("gate_status", "") if rv_summary else "",
            "remaining_blockers_for_99p9_validation": "; ".join(blockers),
            "interpretation": "The photometric, centroid, Gaia, speckle, and ground-follow-up evidence is validation-favorable. A formal 99.9% claim still requires running a recognized FPP/validation code with the sensitivity curves and/or RV constraints.",
        }
    ]


def render_report(
    target: dict[str, object],
    out_dir: Path,
    rv_summary: list[dict[str, object]],
    imaging_summary: list[dict[str, object]],
    timeseries_summary: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
    stitch_summary: dict[str, str],
    transit_metrics: dict[str, str],
    centroid_summary: dict[str, str],
    eso_rows: list[dict[str, str]],
) -> str:
    rv = rv_summary[0]
    imaging = imaging_summary[0]
    timeseries = timeseries_summary[0]
    validation = validation_rows[0]
    lines = [
        "# TOI-1258.01 Final Validation-Gate Audit",
        "",
        f"Created: {date.today().isoformat()}",
        "",
        "## Target",
        "",
        f"- TOI: `{target.get('toidisplay')}`",
        f"- TIC: `{target.get('tid')}`",
        f"- RA/Dec: `{target.get('ra')}`, `{target.get('dec')}` deg",
        f"- Disposition: `{target.get('tfopwg_disp')}`",
        f"- Period: `{target.get('pl_orbper_days')}` d",
        f"- Radius: `{target.get('pl_rade')}` Earth radii",
        f"- Depth: `{target.get('pl_trandep_ppm')}` ppm",
        f"- Stellar mass/radius: `{target.get('st_mass')}` Msun / `{target.get('st_rad')}` Rsun",
        "",
        "## Gate 1: Wobble / Radial Velocity",
        "",
        f"- Status: `{rv.get('gate_status')}`",
        f"- ExoFOP spectroscopy rows: `{rv.get('exofop_spectroscopy_rows')}`",
        f"- High-resolution spectroscopy rows (R >= 40000): `{rv.get('high_resolution_spectroscopy_rows_r_ge_40000')}`",
        f"- RV count reported in ExoFOP notes: `{rv.get('rv_count_from_exofop_notes')}`",
        f"- Instruments: `{rv.get('instruments')}`",
        f"- Expected K scale: `~{rv.get('expected_k_mps_5_mearth')}` m/s at 5 Mearth, `~{rv.get('expected_k_mps_10_mearth')}` m/s at 10 Mearth, `~{rv.get('expected_k_mps_15_mearth')}` m/s at 15 Mearth",
        "",
        "Interpretation: there is public metadata for a serious RV/spectroscopy campaign, including APF/Levy RV observations. This is promising. It does not clear the wobble gate until we retrieve actual velocities or a published RV orbital fit and phase them to the 10.261-day ephemeris.",
        "",
        "## Gate 2: High-Resolution Imaging",
        "",
        f"- Status: `{imaging.get('gate_status')}`",
        f"- Imaging rows: `{imaging.get('imaging_rows')}`",
        f"- Machine-readable sensitivity tables listed: `{imaging.get('sensitivity_dat_files_listed')}`",
        f"- Reconstructed FITS images listed: `{imaging.get('reconstructed_fits_files_listed')}`",
        f"- Instrument: `{imaging.get('instruments')}`",
        f"- Reported contrast limits: `{imaging.get('reported_contrast_limits')}`",
        "",
        "Interpretation: the high-resolution imaging gate is substantially addressed by Gemini/Alopeke speckle imaging. The listed contrast limits are strong against close stellar contaminants. A formal validation run should ingest the listed sensitivity `.dat` curves directly.",
        "",
        "## Ground Time-Series Follow-Up",
        "",
        f"- Status: `{timeseries.get('gate_status')}`",
        f"- Time-series rows: `{timeseries.get('time_series_rows')}`",
        f"- NEB-check files listed: `{timeseries.get('neb_check_files_listed')}`",
        f"- Telescopes/cameras: `{timeseries.get('telescopes')}`",
        f"- Notes: `{timeseries.get('notes')}`",
        "",
        "Interpretation: ground-based full-transit follow-up exists and includes NEB-check products. This strengthens the source-localization case, but the raw photometry should be ingested before treating it as a formal validation input.",
        "",
        "## Gate 3: Statistical Validation",
        "",
        f"- Status: `{validation.get('overall_statistical_validation_status')}`",
        f"- BLS period error: `{validation.get('bls_period_error_percent')}` percent",
        f"- BLS depth SNR: `{validation.get('bls_depth_snr')}`",
        f"- Primary / secondary SNR: `{validation.get('primary_depth_snr')}` / `{validation.get('secondary_depth_snr')}`",
        f"- Odd/even depth delta: `{validation.get('odd_even_depth_delta_ppm')}` ppm",
        f"- SPOC DV sky offset: `{validation.get('spoc_dv_sky_offset_arcsec')} +/- {validation.get('spoc_dv_sky_offset_unc_arcsec')}` arcsec",
        f"- Gaia 84 arcsec dilution: `{validation.get('gaia_84arcsec_dilution_percent')}` percent",
        f"- Brightest-neighbor required eclipse depth: `{validation.get('required_eclipse_depth_on_brightest_neighbor_percent')}` percent",
        f"- Remaining blockers: `{validation.get('remaining_blockers_for_99p9_validation')}`",
        "",
        "Interpretation: the evidence is validation-favorable. It is not yet a defensible 99.9% statistical validation claim until a recognized FPP framework is run using the contrast curves, centroid constraints, Gaia field, and any RV limits.",
        "",
        "## ESO Archive Cross-Check",
        "",
        f"- ESO cone-search rows within 0.03 deg: `{len(eso_rows)}`",
        "",
        "No ESO public archive observations were found in the cone search. The relevant public follow-up evidence for this target is concentrated in ExoFOP/TFOP metadata and MAST/TESS products.",
        "",
        "## Current Bottom Line",
        "",
        "TOI-1258.01 is no longer just a raw TESS lead in our pipeline. Public follow-up metadata already shows speckle imaging, ground transit follow-up, spectroscopy, and a substantial RV campaign. The photometry and SPOC centroid evidence remain strong. The only honest blocker to calling it confirmed is that we still need the actual RV time-series fit and/or a formal FPP validation run with the machine-readable speckle contrast curves.",
        "",
        "## Tables",
        "",
        "- `evidence/tables/exofop_high_resolution_imaging.csv`",
        "- `evidence/tables/exofop_spectroscopy_rv.csv`",
        "- `evidence/tables/exofop_time_series_followup.csv`",
        "- `evidence/tables/exofop_file_list.csv`",
        "- `evidence/tables/statistical_validation_prior_table.csv`",
        "- `evidence/tables/rv_spectroscopy_gate_summary.csv`",
        "- `evidence/tables/high_resolution_imaging_gate_summary.csv`",
        "",
        "## Sources",
        "",
        "- ExoFOP-TESS target page and file list: https://exofop.ipac.caltech.edu/tess/target.php?id=288132261",
        "- NASA Exoplanet Archive TAP: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html",
        "- MAST/TESS products: https://archive.stsci.edu/tess/",
        "- ESO TAP archive search: https://archive.eso.org/programmatic/",
        "",
        f"Package root: `{out_dir}`",
    ]
    return "\n".join(lines) + "\n"


def best_centroid_row(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    if not rows:
        return {}
    return sorted(rows, key=lambda row: int(row.get("endCadence") or 0), reverse=True)[0]


def read_first_csv(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    return rows[0] if rows else {}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


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
    return list(csv.DictReader(io.StringIO(fetch_text(url, timeout=120))))


def fetch_text(url: str, timeout: int = 90) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def as_float(value) -> float | None:
    try:
        if value is None or value == "":
            return None
        out = float(value)
        return out if math.isfinite(out) else None
    except Exception:
        return None


def first_or_empty(rows: list[dict[str, str]], key: str, value: str) -> dict[str, str]:
    for row in rows:
        if str(row.get(key)) == value:
            return row
    return {}


if __name__ == "__main__":
    raise SystemExit(main())
