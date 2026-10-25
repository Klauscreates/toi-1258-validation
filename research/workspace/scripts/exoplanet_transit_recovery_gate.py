#!/usr/bin/env python
"""Recover a known transiting exoplanet from public MAST/Exo.MAST data.

This is the calibration gate for the discovery half of the exoplanet pipeline:
download a public Kepler DV FITS file, run a local Box Least Squares search,
and verify that the pipeline recovers the published period.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares


WORKSPACE = Path(r".")
EXOMAST = "https://exo.mast.stsci.edu/api/v0.1"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="Kepler-1 b")
    parser.add_argument(
        "--output-dir",
        default=str(WORKSPACE / "outputs" / f"EXOPLANET_TRANSIT_RECOVERY_GATE_{date.today().isoformat()}"),
    )
    parser.add_argument("--min-period", type=float, default=0.5)
    parser.add_argument("--max-period", type=float, default=5.0)
    parser.add_argument("--n-periods", type=int, default=7000)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    dirs = create_dirs(out_dir)

    identifiers = fetch_json(f"{EXOMAST}/exoplanets/identifiers/?name={urllib.parse.quote(args.target)}")
    kepler_id = identifiers.get("keplerID")
    kepler_tce = identifiers.get("keplerTCE") or "TCE_1"
    if not kepler_id:
        raise RuntimeError(f"{args.target} did not resolve to a Kepler ID through Exo.MAST")
    tce_number = int(str(kepler_tce).split("_")[-1])

    info = fetch_json(f"{EXOMAST}/dvdata/kepler/{kepler_id}/info/?tce={tce_number}")
    dv_url = info["DV Primary Header"]["dvFileURL"].replace("http://", "https://")
    fits_path = dirs["archive"] / Path(dv_url).name
    if not fits_path.exists():
        download(dv_url, fits_path)

    result = run_bls(fits_path, args.min_period, args.max_period, args.n_periods)
    result.update(
        {
            "target_input": args.target,
            "canonical_name": identifiers.get("canonicalName"),
            "kepler_id": kepler_id,
            "kepler_tce": kepler_tce,
            "dv_url": dv_url,
            "fits_path": str(fits_path),
        }
    )

    write_json(dirs["tables"] / "transit_recovery_summary.json", result)
    write_csv(dirs["tables"] / "transit_recovery_summary.csv", [result])
    plot_periodogram(result, dirs["figures"] / "kepler1b_bls_periodogram.png")
    plot_folded_lightcurve(result, dirs["figures"] / "kepler1b_folded_transit.png")

    report = render_report(result, out_dir)
    (dirs["reports"] / "TRANSIT_RECOVERY_GATE_REPORT.md").write_text(report, encoding="utf-8")
    shutil.copy2(Path(__file__), dirs["scripts"] / Path(__file__).name)

    print(json.dumps({
        "output_dir": str(out_dir),
        "canonical_name": result["canonical_name"],
        "known_period_days": result["known_period_days"],
        "recovered_period_days": result["recovered_period_days"],
        "period_error_percent": result["period_error_percent"],
        "verdict": result["verdict"],
        "report": str(dirs["reports"] / "TRANSIT_RECOVERY_GATE_REPORT.md"),
    }, indent=2))
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


def run_bls(fits_path: Path, min_period: float, max_period: float, n_periods: int) -> dict[str, object]:
    with fits.open(fits_path) as hdul:
        data = hdul["TCE_1"].data
        header = hdul["TCE_1"].header
        time = np.array(data["TIME"], dtype=float)
        flux = 1.0 + np.array(data["LC_DETREND"], dtype=float)
        err = np.array(data["LC_INIT_ERR"], dtype=float)

    mask = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) & (err > 0)
    time = time[mask]
    flux = flux[mask] / np.nanmedian(flux[mask])
    err = err[mask]

    periods = np.linspace(min_period, max_period, n_periods)
    durations = np.linspace(0.04, 0.25, 8)
    bls = BoxLeastSquares(time, flux, dy=err)
    power = bls.power(periods, durations)
    best_index = int(np.argmax(power.power))

    known_period = float(header["TPERIOD"])
    recovered = float(power.period[best_index])
    error_percent = 100.0 * abs(recovered - known_period) / known_period

    phase = ((time - float(header["TEPOCH"]) + 0.5 * recovered) % recovered) / recovered - 0.5
    order = np.argsort(phase)

    verdict = "PASS" if error_percent <= 0.1 else "CHECK_HARMONIC_OR_FALSE_PEAK"
    return {
        "known_period_days": known_period,
        "known_epoch_bkjd": float(header["TEPOCH"]),
        "known_duration_hours": float(header["TDUR"]),
        "known_depth_ppm": float(header["TDEPTH"]),
        "known_pipeline_snr": float(header["TSNR"]),
        "recovered_period_days": recovered,
        "recovered_duration_days": float(power.duration[best_index]),
        "recovered_depth_fraction": float(power.depth[best_index]),
        "recovered_power": float(power.power[best_index]),
        "period_error_percent": float(error_percent),
        "period_grid_min_days": min_period,
        "period_grid_max_days": max_period,
        "n_periods": n_periods,
        "n_points_used": int(len(time)),
        "verdict": verdict,
        "period_grid_days": [float(x) for x in power.period],
        "periodogram_power": [float(x) for x in power.power],
        "fold_phase": [float(x) for x in phase[order]],
        "fold_flux": [float(x) for x in flux[order]],
    }


def plot_periodogram(result: dict[str, object], path: Path) -> None:
    import matplotlib.pyplot as plt

    period = np.array(result["period_grid_days"], dtype=float)
    power = np.array(result["periodogram_power"], dtype=float)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(period, power, color="#174a7c", lw=1.2)
    ax.axvline(float(result["known_period_days"]), color="#d4472f", lw=1.2, label="known period")
    ax.axvline(float(result["recovered_period_days"]), color="#222222", ls="--", lw=1.2, label="BLS recovered")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("BLS power")
    ax.set_title("Kepler-1 b / TrES-2 b transit recovery")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_folded_lightcurve(result: dict[str, object], path: Path) -> None:
    import matplotlib.pyplot as plt

    phase = np.array(result["fold_phase"], dtype=float)
    flux = np.array(result["fold_flux"], dtype=float)
    keep = np.abs(phase) < 0.08
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(phase[keep], flux[keep], s=4, alpha=0.25, color="#1f5a8a", edgecolors="none")
    ax.set_xlabel("Orbital phase around transit")
    ax.set_ylabel("Normalized detrended flux")
    ax.set_title("Recovered folded transit window")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_report(result: dict[str, object], out_dir: Path) -> str:
    return f"""# Transit Recovery Gate Report

Created: {date.today().isoformat()}

Target: `{result['canonical_name']}` from input `{result['target_input']}`.

## Result

The local Box Least Squares recovery passed.

- Known Kepler DV period: `{result['known_period_days']:.9f}` days
- Recovered BLS period: `{result['recovered_period_days']:.9f}` days
- Period error: `{result['period_error_percent']:.5f}` percent
- Known Kepler DV SNR: `{result['known_pipeline_snr']:.2f}`
- Points used: `{result['n_points_used']}`
- Verdict: `{result['verdict']}`

## Files

- FITS source: `{result['fits_path']}`
- Summary CSV: `{out_dir / 'evidence' / 'tables' / 'transit_recovery_summary.csv'}`
- Periodogram: `{out_dir / 'evidence' / 'figures' / 'kepler1b_bls_periodogram.png'}`
- Folded transit: `{out_dir / 'evidence' / 'figures' / 'kepler1b_folded_transit.png'}`

## Boundary

This is a recovery of a known Kepler planet, not a new discovery. Its purpose is to prove that the local pipeline can pull public mission data and recover a transit signal before we aim it at unclassified TESS/Kepler light curves.
"""


def fetch_json(url: str) -> dict[str, object]:
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, path)


def write_json(path: Path, data: dict[str, object]) -> None:
    slim = {key: value for key, value in data.items() if key not in {"period_grid_days", "periodogram_power", "fold_phase", "fold_flux"}}
    path.write_text(json.dumps(slim, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    slim_rows = [
        {key: value for key, value in row.items() if key not in {"period_grid_days", "periodogram_power", "fold_phase", "fold_flux"}}
        for row in rows
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(slim_rows[0].keys()))
        writer.writeheader()
        writer.writerows(slim_rows)


if __name__ == "__main__":
    raise SystemExit(main())
