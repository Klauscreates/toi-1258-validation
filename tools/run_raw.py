"""Portable entry points around preserved historical source, with explicit inputs."""
from pathlib import Path
import argparse
import csv
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import time

from raw_inputs import safe_path, verify


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New directory; existing output is rejected")
    parser.add_argument("--local-map", type=Path, help="Optional private map for an existing verified cache")
    parser.add_argument("--dolphot-root", type=Path)
    parser.add_argument("--rows", help="PHAT one-based exposure rows, comma-separated; omit for all 32")
    parser.add_argument("--skip-centroids", action="store_true", help="TOI light-curve-only run")
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists() or out.is_relative_to(repo):
        raise ValueError("Choose a new output directory outside the repository")
    manifest = json.loads((repo / "reproduce/inputs.json").read_text())
    project = manifest["project"]
    mapping = json.loads(args.local_map.read_text())[project] if args.local_map else {}
    records = {r["path"]: r for r in manifest["inputs"]}

    def path(name):
        p = Path(mapping[name]) if name in mapping else safe_path(args.data_root, name)
        verify(p, records[name])
        return p

    def folder(prefix):
        paths = [path(n) for n in records if n.startswith(prefix + "/") and
                 "/" not in n[len(prefix) + 1:]]
        parents = {p.parent for p in paths}
        if len(parents) != 1:
            raise ValueError(f"Expected one complete directory for {prefix}")
        return parents.pop()

    sys.path.insert(0, str(repo / "research/workspace/scripts"))
    start = time.monotonic()
    consumed = []
    out.mkdir(parents=True)
    if project == "irpca-002":
        base = importlib.import_module("irpca_independent_reaudit")
        base.ROOT = out
        base.OUT = out / "outputs/IRPCA002_INDEPENDENT_REAUDIT_2026-10-01"
        base.OUT.mkdir(parents=True)
        base.RAW = folder("sirius/science_raw_20141207_fits")
        flat = path("sirius/reduction_noadi_20141207/products_flat/ird_flat_field.fits")
        path("sirius/reduction_noadi_20141207/products_sky/ird_sky_bg.fits")
        base.CAL = flat.parents[1]
        base.BETA = folder("sirius/beta_pic_flux_reference/decompressed")
        folder("sirius/calib_raw_20141206_08_fits")
        base.tests()
        base.templates(*base.calibrations())
        astrometry = importlib.import_module("irpca_astrometry_verify")
        astrometry.inspect()
        astrometry.fit_centers()
        corrected = importlib.import_module("irpca_original_corrected_test")
        corrected.main()
        consumed = [n for n in records if n.startswith("sirius/")]
        scope = "Independent centering, PSF construction, corrected original-coordinate extraction; not complete historical survey sweeps"
    elif project == "toi-1258-validation":
        code = importlib.import_module("exoplanet_toi_centroid_stitch_gate")
        t = json.loads((repo / "reproduce/target.json").read_text())
        target = code.Target(toi=t["toi"], tic=t["tic"], ra=t["ra_deg"], dec=t["dec_deg"],
                             period_days=t["period_days"], t0_btjd=t["t0_btjd"],
                             duration_days=t["duration_days"], depth_ppm=t["depth_ppm"],
                             tmag=None, disposition="historical PC", source_row={})
        names = sorted(n for n in records if n.startswith("tess/") and n.endswith("_lc.fits"))
        stitched = code.build_stitched_lightcurve([path(n) for n in names])
        dirs = code.create_dirs(out)
        lightcurve = dirs["tables"] / "stitched_lightcurve_downsample.csv"
        code.write_csv(lightcurve, code.downsample_rows(stitched["rows"], max_rows=25000))
        folding = importlib.import_module("exoplanet_toi_blocker_clearance")
        folding.prepare_folded_lightcurve(lightcurve, dirs["tables"], target.period_days,
                                         3333.648786, 0.202292)
        bls = code.run_bls_and_plots(target, stitched, dirs, .5, 20., 6000)
        metrics = code.compute_transit_metrics(stitched["time"], stitched["flux_rel"],
                  target.period_days, target.t0_btjd, target.duration_days, sector=stitched["sector"])
        code.write_csv(dirs["tables"] / "stitched_bls_summary.csv", [bls])
        code.write_csv(dirs["tables"] / "stitched_transit_metrics.csv", [metrics["summary"]])
        code.write_csv(dirs["tables"] / "stitched_transit_depth_by_sector.csv", metrics["sector_rows"])
        consumed = names
        if not args.skip_centroids:
            names = sorted(n for n in records if n.startswith("tess/") and n.endswith("_tp.fits"))
            centroids = code.run_centroid_gate(target, [path(n) for n in names], dirs)
            code.write_csv(dirs["tables"] / "centroid_gate_summary.csv", [centroids["summary"]])
            code.write_csv(dirs["tables"] / "centroid_shift_summary_by_sector.csv", centroids["sector_rows"])
            consumed += names
        scope = "Pinned 12-sector TESS LC/BLS/transit metrics and optional TPF centroids; does not rerun TRICERATOPS or measure RV mass"
    else:
        if args.dolphot_root is None:
            raise ValueError("PHAT requires --dolphot-root with compiled ACS/WFC3 modules and PSF/PAM files")
        code = importlib.import_module("candidate003_dolphot_module_forced_pass")
        code.WORKSPACE = out
        code.OUT = out
        code.RUN = out / "module_forced_catalog_v1"
        code.WORK = code.RUN / "work"
        code.TABLES = code.RUN / "tables"
        code.REPORTS = code.RUN / "reports"
        code.DOLPHOT_ROOT = args.dolphot_root.resolve()
        code.BIN = code.DOLPHOT_ROOT / "bin"
        for constant, name in (("DOLPHOT", "dolphot"), ("CALCSKY", "calcsky"),
                               ("SPLITGROUPS", "splitgroups"), ("ACSMASK", "acsmask"), ("WFC3MASK", "wfc3mask")):
            executable = code.BIN / (name + (".exe" if os.name == "nt" else ""))
            if not executable.is_file():
                raise FileNotFoundError(f"Missing DOLPHOT executable: {name}")
            setattr(code, constant, executable)
        with (repo / "data/epsf_measurements.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        selected = {int(v) for v in args.rows.split(",")} if args.rows else set(range(1, len(rows) + 1))
        if not selected or not selected <= set(range(1, len(rows) + 1)):
            raise ValueError("Invalid PHAT row selection")
        for i, row in enumerate(rows, 1):
            n = "hst/" + row["product_filename"]
            row["local_fits_path"] = str(path(n)) if i in selected else "not_selected"
            if i in selected:
                consumed.append(n)
        code.MEASUREMENTS = out / "input_measurements.csv"
        with code.MEASUREMENTS.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        os.environ["CANDIDATE003_DOLPHOT_ROWS"] = ",".join(map(str, sorted(selected)))
        code.main()
        summary = json.loads((code.TABLES / "candidate003_dolphot_module_forced_summary.json").read_text())
        with (code.TABLES / "candidate003_dolphot_module_forced_catalog.csv").open(newline="") as f:
            measurements = list(csv.DictReader(f))
        codes = ("mask_exit_code", "splitgroups_exit_code", "calcsky_exit_code", "dolphot_exit_code")
        if (summary["dolphot_success_frames"] != len(selected)
                or summary["target_catalog_rows_recovered"] != len(selected)
                or any(r[k] != "0" for r in measurements for k in codes)):
            raise RuntimeError("One or more DOLPHOT exposures failed; see run outputs")
        scope = "Instrument-specific DOLPHOT from exact FLC/FLT exposures at historical positions; separate from the recovered ePSF estimator"
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "astropy", "matplotlib")}
    record = dict(project=project, status="completed", scope=scope,
                  python=platform.python_version(), platform=platform.system(), versions=versions,
                  seconds=time.monotonic() - start, consumed_inputs=consumed)
    (out / "RUN.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
