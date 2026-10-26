"""Reduce the rerun outputs to the numbers quoted in rerun/README.md.

Usage: python summarize.py RESULTS_DIR SPOC_DV_CENTROID_CSV
Writes RESULTS_DIR/SUMMARY.json. Standard library only.
"""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

res = Path(sys.argv[1])
spoc_csv = Path(sys.argv[2])
# TIC 288132263 relative to the target (TRICERATOPS stars table: 17.029", PA 64.221 deg E of N)
T2 = (17.029 * math.sin(math.radians(64.221)), 17.029 * math.cos(math.radians(64.221)))
UNUSABLE_SAP = {60, 73}  # one or two transits each, strong SAP systematics


def wmean(vals, errs):
    w = [1 / e**2 for e in errs]
    m = sum(v * wi for v, wi in zip(vals, w)) / sum(w)
    chi2 = sum(((v - m) / e) ** 2 for v, e in zip(vals, errs))
    return m, sum(w) ** -0.5, chi2, len(vals) - 1


out = {}

# depth consistency and aperture trend
rows = json.loads((res / "depth_aperture.json").read_text())
for key in ("pdcsap", "sap", "r3.0"):
    m, e, chi2, dof = wmean([r[f"{key}_depth"] for r in rows], [r[f"{key}_err"] for r in rows])
    out[f"{key}_depth_ppm"] = dict(mean=round(m, 1), err=round(e, 1), chi2=round(chi2, 1), dof=dof)
trend = {r["sector"]: round(r["r4.0_depth"] - r["r2.0_depth"], 1) for r in rows}
good = [v for s, v in trend.items() if s not in UNUSABLE_SAP]
out["aperture_r2_to_r4_change_ppm"] = dict(
    per_sector=trend, usable_sectors=len(good), shallower=sum(v < 0 for v in good),
    median=round(statistics.median(good), 1), excluded_sectors=sorted(UNUSABLE_SAP))

# own difference-image centroids versus target and T2 (sectors with >= 2 transits)
cen = [r for r in json.loads((res / "diff_centroids.json").read_text()) if r["n_transits"] >= 2]
per_axis = [r["offset_err_arcsec"] / math.sqrt(2) for r in cen]
mra, se, _, _ = wmean([r["dRA_arcsec"] for r in cen], per_axis)
mde, _, _, _ = wmean([r["dDec_arcsec"] for r in cen], per_axis)
out["difference_image_centroid"] = dict(
    sectors=len(cen), mean_offset_arcsec=[round(mra, 2), round(mde, 2)], per_axis_err_arcsec=round(se, 2),
    t2_offset_arcsec=[round(T2[0], 2), round(T2[1], 2)],
    t2_exclusion_sigma=round(math.hypot(mra - T2[0], mde - T2[1]) / se, 1))

# SPOC DV centroids versus target and T2
spoc = []
for s in csv.DictReader(spoc_csv.open()):
    parts = s["filename"].split("-")
    if s["toiId_planet"] != "1258.01":
        continue
    ra, de = float(s["ms_tic_mean_ra_offset_arcsec"]), float(s["ms_tic_mean_dec_offset_arcsec"])
    er, ed = float(s["ms_tic_mean_ra_offset_unc_arcsec"]), float(s["ms_tic_mean_dec_offset_unc_arcsec"])
    spoc.append(dict(sectors=f"{parts[1]}-{parts[2]}", offset_arcsec=[ra, de], err_arcsec=[er, ed],
                     sigma_from_target=round(math.hypot(ra / er, de / ed), 1),
                     sigma_from_t2=round(math.hypot((ra - T2[0]) / er, (de - T2[1]) / ed), 1)))
out["spoc_dv_centroids"] = spoc

# TRICERATOPS runs
fpp = {}
for d in sorted(p for p in res.iterdir() if p.is_dir() and p.name.startswith("fpp_")):
    s = json.loads((d / "summary.json").read_text())
    scen = {}
    for f in sorted(d.glob("probs_seed*.csv")):
        # a scenario can appear once per nearby star; sum within each seed
        tot = {}
        for row in csv.DictReader(f.open()):
            tot[row["scenario"]] = tot.get(row["scenario"], 0.0) + float(row["prob"])
        for k, v in tot.items():
            scen.setdefault(k, []).append(v)
    fpp[d.name] = dict(
        curves=s["curves"], background=s["background"], seeds=len(s["runs"]), N=s["runs"][0]["N"],
        FPP_max=max(r["FPP"] for r in s["runs"]), NFPP_max=max(r["NFPP"] for r in s["runs"]),
        scenario_range={k: [min(v), max(v)] for k, v in scen.items() if max(v) > 1e-6})
out["triceratops"] = fpp

(res / "SUMMARY.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
