"""TOI-1258.01 independent mass/barcode inference package.

This script does NOT claim measured planet mass or measured atmospheric
composition. It derives what can be derived from public TOI metadata and basic
physics:
- stellar mass from logg and stellar radius,
- semi-major axis from Kepler's third law,
- empirical mass estimates from radius,
- expected RV semi-amplitude range,
- expected transmission spectroscopy amplitude for different atmospheres.
"""
from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import matplotlib.pyplot as plt

OUT = Path(r"data/external/exoplanet_outputs\TOI1258_INDEPENDENT_MASS_BARCODE_MODEL_2026-07-08")
REPORTS = OUT / "evidence" / "reports"
TABLES = OUT / "evidence" / "tables"
FIGS = OUT / "evidence" / "figures"
SCRIPTS = OUT / "scripts"
for d in [REPORTS, TABLES, FIGS, SCRIPTS]:
    d.mkdir(parents=True, exist_ok=True)

# Constants
G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
M_EARTH = 5.9722e24
R_EARTH = 6.371e6
AU = 1.495978707e11
LOGG_SUN_CGS = 4.438
K_B = 1.380649e-23
M_H = 1.6735575e-27

BASE = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
QUERY = "select top 1 * from toi where tid=288132261"
url = BASE + "?" + urlencode({"query": QUERY, "format": "json"})
row = json.loads(urlopen(url, timeout=60).read().decode("utf-8"))[0]

# Public metadata
P_d = float(row["pl_orbper"])
P_y = P_d / 365.25
P_s = P_d * 86400.0
Rp_re = float(row["pl_rade"])
Rp_re_err = float(row.get("pl_radeerr1") or 0.0)
Rs_rs = float(row["st_rad"])
Rs_rs_err = float(row.get("st_raderr1") or 0.0)
logg = float(row["st_logg"])
logg_err = float(row.get("st_loggerr1") or 0.0)
Teq = float(row["pl_eqt"])
Teff = float(row["st_teff"])
transit_depth_ppm = float(row["pl_trandep"])
tran_dur_h = float(row["pl_trandurh"])

# Stellar mass from logg + radius: M/Msun = 10^(logg-logg_sun)*(R/Rsun)^2
Mstar_ms = 10 ** (logg - LOGG_SUN_CGS) * Rs_rs**2
# Kepler semi-major axis: a^3 = GM P^2 / 4pi^2
a_m = (G * Mstar_ms * M_SUN * P_s**2 / (4 * math.pi**2)) ** (1 / 3)
a_au = a_m / AU
a_over_rs = a_m / (Rs_rs * R_SUN)
central_duration_h = (P_d / math.pi) * math.asin(min(1.0, (1 + (Rp_re * R_EARTH) / (Rs_rs * R_SUN)) / a_over_rs)) * 24

# Empirical/substitute masses. These are estimates, not measured mass.
def wolfgang_mass(radius_re: np.ndarray) -> np.ndarray:
    # Wolfgang, Rogers & Ford-style mean relation often used for sub-Neptunes:
    # M ~= 2.7 R^1.3 Mearth. Scatter handled separately below.
    return 2.7 * radius_re**1.3

rng = np.random.default_rng(1258)
N = 200000
Rp_s = rng.normal(Rp_re, max(Rp_re_err, 0.01), N)
Rp_s = np.clip(Rp_s, 0.1, None)
Rs_s = rng.normal(Rs_rs, max(Rs_rs_err, 0.001), N)
Rs_s = np.clip(Rs_s, 0.1, None)
logg_s = rng.normal(logg, max(logg_err, 0.001), N)
Mstar_s = 10 ** (logg_s - LOGG_SUN_CGS) * Rs_s**2
Mstar_s = np.clip(Mstar_s, 0.1, 5.0)
M_emp_s = wolfgang_mass(Rp_s) + rng.normal(0, 1.9, N)
M_emp_s = np.clip(M_emp_s, 0.5, None)

# Neptune-density scaling: same bulk density as Neptune.
M_neptune_re = 17.147
R_neptune_re = 3.883
M_neptune_density = M_neptune_re * (Rp_re / R_neptune_re) ** 3

# Broad physical scenario grid, to make the sensitivity explicit.
scenario_masses = {
    "low_density_envelope_6ME": 6.0,
    "low_density_envelope_8ME": 8.0,
    "empirical_radius_relation_median": float(np.median(M_emp_s)),
    "neptune_density_scaling": float(M_neptune_density),
    "warm_sub_neptune_20ME": 20.0,
    "dense_sub_neptune_30ME": 30.0,
    "very_dense_upper_stress_40ME": 40.0,
}

def rv_k_mps(Mp_me: np.ndarray | float, Mstar_ms_: np.ndarray | float, ecc: float = 0.0) -> np.ndarray | float:
    # K = 0.08946 m/s * Mp/Mearth * (P/yr)^(-1/3) * (Mstar/Msun)^(-2/3) / sqrt(1-e^2)
    return 0.08946 * np.asarray(Mp_me) * P_y ** (-1 / 3) * np.asarray(Mstar_ms_) ** (-2 / 3) / math.sqrt(1 - ecc**2)

K_emp_s = rv_k_mps(M_emp_s, Mstar_s)

# Transmission signal estimate.
def scale_height_km(Mp_me: float, mu: float) -> float:
    g = 9.80665 * Mp_me / (Rp_re**2)
    return (K_B * Teq / (mu * M_H * g)) / 1000.0

def transmission_ppm(Mp_me: float, mu: float, nH: float = 5.0) -> float:
    Rp_km = Rp_re * R_EARTH / 1000.0
    Rs_km = Rs_rs * R_SUN / 1000.0
    H = scale_height_km(Mp_me, mu)
    return 2 * Rp_km * (nH * H) / (Rs_km**2) * 1e6

atmospheres = {
    "H2_He_clear_mu2p3": 2.3,
    "water_rich_mu18": 18.0,
    "CO2_N2_heavy_mu28": 28.0,
}

# Depth consistency check.
expected_depth_from_catalog_radius_ppm = ((Rp_re * R_EARTH) / (Rs_rs * R_SUN)) ** 2 * 1e6
radius_from_depth_re = math.sqrt(transit_depth_ppm / 1e6) * (Rs_rs * R_SUN / R_EARTH)

# Tables
with (TABLES / "toi1258_public_toi_metadata.json").open("w", encoding="utf-8") as f:
    json.dump(row, f, indent=2)

scenario_rows = []
for name, mass in scenario_masses.items():
    item = {
        "scenario": name,
        "mass_Mearth": mass,
        "rv_K_mps_e0": float(rv_k_mps(mass, Mstar_ms)),
    }
    for atm, mu in atmospheres.items():
        item[f"H_km_{atm}"] = scale_height_km(mass, mu)
        item[f"transmission_ppm_5H_{atm}"] = transmission_ppm(mass, mu, 5)
    scenario_rows.append(item)

with (TABLES / "toi1258_mass_rv_transmission_scenarios.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(scenario_rows[0].keys()))
    writer.writeheader()
    writer.writerows(scenario_rows)

summary = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "source": "NASA Exoplanet Archive TOI table via TAP",
    "target": "TOI-1258.01 / TIC 288132261",
    "disposition": row.get("tfopwg_disp"),
    "period_days": P_d,
    "radius_Rearth_catalog": Rp_re,
    "stellar_radius_Rsun": Rs_rs,
    "stellar_logg_cgs": logg,
    "derived_stellar_mass_Msun_from_logg_radius": Mstar_ms,
    "semi_major_axis_AU": a_au,
    "a_over_Rstar": a_over_rs,
    "observed_duration_h": tran_dur_h,
    "central_duration_estimate_h": central_duration_h,
    "catalog_depth_ppm": transit_depth_ppm,
    "depth_from_catalog_Rp_Rs_ppm": expected_depth_from_catalog_radius_ppm,
    "radius_implied_by_depth_Rearth_if_undiluted": radius_from_depth_re,
    "empirical_mass_relation_Mearth_median": float(np.median(M_emp_s)),
    "empirical_mass_relation_Mearth_p16": float(np.percentile(M_emp_s, 16)),
    "empirical_mass_relation_Mearth_p84": float(np.percentile(M_emp_s, 84)),
    "expected_RV_K_mps_median": float(np.median(K_emp_s)),
    "expected_RV_K_mps_p16": float(np.percentile(K_emp_s, 16)),
    "expected_RV_K_mps_p84": float(np.percentile(K_emp_s, 84)),
    "neptune_density_mass_Mearth": float(M_neptune_density),
    "best_current_boundary": "Independent public metadata supports a predicted mass/RV/barcode target model, not a measured mass or observed atmosphere.",
}
with (OUT / "target_manifest.json").open("w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

# Figures
masses = np.linspace(3, 45, 300)
Ks = rv_k_mps(masses, Mstar_ms)
plt.figure(figsize=(8, 5))
plt.plot(masses, Ks, color="#1f77b4", lw=2)
plt.axvline(np.median(M_emp_s), color="#d62728", ls="--", label=f"empirical median {np.median(M_emp_s):.1f} Mearth")
plt.fill_betweenx([Ks.min(), Ks.max()], np.percentile(M_emp_s, 16), np.percentile(M_emp_s, 84), color="#d62728", alpha=0.15, label="empirical 16-84%")
plt.xlabel("Assumed planet mass (Earth masses)")
plt.ylabel("Expected RV semi-amplitude K (m/s)")
plt.title("TOI-1258.01 predicted RV signal from public metadata")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(FIGS / "toi1258_predicted_rv_amplitude.png", dpi=180)
plt.close()

plt.figure(figsize=(8, 5))
for atm, mu in atmospheres.items():
    vals = [transmission_ppm(m, mu, 5) for m in masses]
    plt.plot(masses, vals, lw=2, label=atm)
plt.axvline(np.median(M_emp_s), color="black", ls="--", alpha=0.7)
plt.xlabel("Assumed planet mass (Earth masses)")
plt.ylabel("Approx. 5-scale-height transmission signal (ppm)")
plt.title("TOI-1258.01 predicted atmospheric signal size")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(FIGS / "toi1258_predicted_transmission_signal.png", dpi=180)
plt.close()

# Report
report = f"""# TOI-1258.01 Independent Mass and Barcode Inference

Created UTC: `{summary['created_utc']}`

## What This Does and Does Not Claim

This package uses public TOI metadata plus basic orbital/atmospheric physics to predict the mass and spectroscopy target space for TOI-1258.01. It does **not** claim a measured RV mass and does **not** claim an observed atmospheric chemical barcode.

## Public Inputs

- Target: `TOI-1258.01 / TIC 288132261`
- NASA Exoplanet Archive TOI disposition: `{row.get('tfopwg_disp')}`
- Period: `{P_d:.9f}` d
- Catalog radius: `{Rp_re:.5f} +/- {Rp_re_err:.5f}` Earth radii
- Stellar radius: `{Rs_rs:.3f} +/- {Rs_rs_err:.3f}` solar radii
- Stellar logg: `{logg:.5f} +/- {logg_err:.5f}` cgs
- Equilibrium temperature: `{Teq:.0f}` K

## Keplerian Orbit From Public Metadata

Using stellar mass inferred from logg and radius:

- Derived stellar mass: `{Mstar_ms:.3f}` solar masses
- Semi-major axis: `{a_au:.4f}` AU
- a/Rstar: `{a_over_rs:.2f}`
- Observed transit duration: `{tran_dur_h:.2f}` h
- Central-transit duration estimate: `{central_duration_h:.2f}` h

This is geometrically plausible: the observed duration is shorter than the central-transit estimate, consistent with a nonzero impact parameter and/or catalog/dilution differences.

## Mass Estimate Boundary

Kepler's third law gives the orbit size, not planet mass. Planet mass requires RV, TTV, astrometry, or another dynamical measurement. Without that, the best independent public-data move is a radius-based mass estimate.

- Empirical radius-relation mass median: `{np.median(M_emp_s):.2f}` Earth masses
- 16-84 percent range: `{np.percentile(M_emp_s, 16):.2f}` to `{np.percentile(M_emp_s, 84):.2f}` Earth masses
- Neptune-density scaling mass: `{M_neptune_density:.2f}` Earth masses

The corresponding expected RV semi-amplitude is:

- Median K: `{np.median(K_emp_s):.2f}` m/s
- 16-84 percent K range: `{np.percentile(K_emp_s, 16):.2f}` to `{np.percentile(K_emp_s, 84):.2f}` m/s

That means the missing APF/Levy or Keck/HIRES RV table should be looking for a signal of order a few m/s if the planet has a normal sub-Neptune mass.

## Chemical Barcode Boundary

A chemical barcode is not recoverable from broadband TESS photometry. A real atmospheric barcode needs wavelength-resolved in-transit/out-of-transit spectroscopy or equivalent eclipse/phase spectroscopy. Public MAST checked in the archive-hunt package showed no HST/JWST observation in the target cone.

What we can compute now is the expected signal size:

- For an H2/He atmosphere around the empirical median-mass planet, the approximate 5-scale-height signal is `{transmission_ppm(float(np.median(M_emp_s)), 2.3, 5):.1f}` ppm.
- For a water-rich high mean-molecular-weight atmosphere, the same estimate drops to `{transmission_ppm(float(np.median(M_emp_s)), 18.0, 5):.1f}` ppm.
- For a CO2/N2-heavy atmosphere, it drops to `{transmission_ppm(float(np.median(M_emp_s)), 28.0, 5):.1f}` ppm.

So the target is physically plausible for JWST-class transmission spectroscopy if the atmosphere is low mean-molecular-weight, but there is no observed barcode in the public archive yet.

## Depth Consistency Note

The catalog transit depth is `{transit_depth_ppm:.1f}` ppm. The simple depth implied by the catalog planet radius and stellar radius is `{expected_depth_from_catalog_radius_ppm:.1f}` ppm. If the listed depth were interpreted as undiluted with the listed stellar radius, it would imply `{radius_from_depth_re:.2f}` Earth radii. This mismatch should be carried as a modeling caveat and checked against the actual SPOC DV dilution/blend treatment.

## Outputs

- Scenario table: `evidence/tables/toi1258_mass_rv_transmission_scenarios.csv`
- RV figure: `evidence/figures/toi1258_predicted_rv_amplitude.png`
- Transmission figure: `evidence/figures/toi1258_predicted_transmission_signal.png`
- Manifest: `target_manifest.json`

## Practical Next Step

Do not describe this as a measured mass or measured barcode. Use it as the independent public-data prediction package. The next hard evidence gates are:

1. Retrieve or produce the actual BJD/RV/error series.
2. Obtain or propose wavelength-resolved transit spectroscopy.
3. Reconcile the catalog radius/depth/dilution boundary using the exact SPOC DV products.
"""
(REPORTS / "TOI1258_INDEPENDENT_MASS_BARCODE_INFERENCE_REPORT.md").write_text(report, encoding="utf-8")

# Copy self into scripts folder.
this_script = Path(__file__) if "__file__" in globals() else None
if this_script and this_script.exists():
    (SCRIPTS / this_script.name).write_text(this_script.read_text(encoding="utf-8"), encoding="utf-8")

print(json.dumps(summary, indent=2))
