"""Three-panel figure for the TOI-1258 b Research Note.

(a) Phase-folded SPOC PDCSAP photometry with the median MCMC transit model.
(b) Gemini/'Alopeke 562 and 832 nm 5-sigma contrast curves.
(c) Transit-source centroid offsets relative to the target, with TIC 288132263 marked.

Usage: python make_figure.py DATA_DIR REPO_ROOT OUT_PNG
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.integrate

for _name, _value in (("int", int), ("trapz", getattr(np, "trapezoid", None)), ("NaN", np.nan), ("Inf", np.inf)):
    if not hasattr(np, _name):
        setattr(np, _name, _value)
if not hasattr(scipy.integrate, "trapz"):
    scipy.integrate.trapz = scipy.integrate.trapezoid

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from astropy.io import fits  # noqa: E402
from pytransit import QuadraticModel  # noqa: E402

data, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
fit = json.loads((repo / "rerun/outputs/transit_fit.json").read_text())
med = {k: fit[k]["median"] for k in ("T0_BTJD", "P_days", "k", "a_over_rstar", "inc_deg", "u1", "u2")}
P, T0 = med["P_days"], med["T0_BTJD"]
HALF, WIN = 0.2023 / 2, 0.5

ph_all, f_all = [], []
for p in sorted((data / "tess" / "lc").rglob("*lc.fits")):
    d = fits.getdata(p, 1)
    q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
    t, f = d["TIME"][q], d["PDCSAP_FLUX"][q]
    n = np.round((t - T0) / P)
    ph = t - (T0 + n * P)
    for k in np.unique(n):
        m = (n == k) & (np.abs(ph) < WIN)
        oot = m & (np.abs(ph) > HALF + 0.03)
        if oot.sum() < 60 or (oot & (ph < 0)).sum() < 20 or (oot & (ph > 0)).sum() < 20:
            continue
        c = np.polyfit(ph[oot], f[oot], 1)
        ph_all.append(ph[m])
        f_all.append(f[m] / np.polyval(c, ph[m]))
ph = np.concatenate(ph_all) * 24
fl = (np.concatenate(f_all) - 1) * 1e6

edges = np.arange(-12, 12.01, 20 / 60)
idx = np.digitize(ph, edges) - 1
bx = np.array([ph[idx == i].mean() for i in range(len(edges) - 1) if (idx == i).sum() > 10])
by = np.array([fl[idx == i].mean() for i in range(len(edges) - 1) if (idx == i).sum() > 10])
be = np.array([fl[idx == i].std() / math.sqrt((idx == i).sum()) for i in range(len(edges) - 1)
               if (idx == i).sum() > 10])

tg = np.linspace(-WIN, WIN, 2000)
tm = QuadraticModel()
tm.set_data(tg + T0)
model = tm.evaluate(med["k"], [med["u1"], med["u2"]], T0, P, med["a_over_rstar"], np.radians(med["inc_deg"]), 0, 0)

fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), gridspec_kw=dict(width_ratios=[1.45, 1, 1]))

a = ax[0]
a.plot(ph, fl, ".", ms=1, color="0.8", rasterized=True, zorder=1)
a.errorbar(bx, by, be, fmt="o", ms=3, color="k", elinewidth=0.8, zorder=2, label="20-min bins")
a.plot(tg * 24, (model - 1) * 1e6, color="C3", lw=1.8, zorder=3, label="median model")
a.set(xlim=(-12, 12), ylim=(-1600, 900), xlabel="Hours from mid-transit", ylabel="Relative flux (ppm)")
a.set_title(f"(a) TESS, {fit['n_transits']} transits, 12 sectors", fontsize=10)
a.legend(loc="lower left", fontsize=8, frameon=False)

a = ax[1]
for band, col, lab in (("562", "C0", "562 nm"), ("832", "C3", "832 nm")):
    cc = np.loadtxt(repo / f"rerun/inputs/cc_gemini_{band}.csv", delimiter=",")
    a.plot(cc[:, 0], cc[:, 1], color=col, lw=1.8, label=f"'Alopeke {lab}")
a.set(xlim=(0, 1.2), ylim=(9, 0), xlabel="Separation (arcsec)", ylabel=r"$\Delta$mag (5$\sigma$)")
a.set_title("(b) Gemini speckle, 2020 Feb 15", fontsize=10)
a.legend(loc="upper right", fontsize=8, frameon=False)

a = ax[2]
T2 = (17.029 * math.sin(math.radians(64.221)), 17.029 * math.cos(math.radians(64.221)))
for s in csv.DictReader((repo / "data/spoc_dv_centroid_summary.csv").open()):
    parts = s["filename"].split("-")
    if s["toiId_planet"] != "1258.01" or parts[1] != parts[2]:
        continue
    x, y = float(s["ms_tic_mean_ra_offset_arcsec"]), float(s["ms_tic_mean_dec_offset_arcsec"])
    ex, ey = float(s["ms_tic_mean_ra_offset_unc_arcsec"]), float(s["ms_tic_mean_dec_offset_unc_arcsec"])
    a.errorbar(x, y, ey, ex, fmt="o", ms=3, color="0.55", elinewidth=0.7, zorder=2)
    if parts[1] in ("s0026", "s0074"):
        a.annotate(f"S{int(parts[1][1:])}", (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7, color="0.4")
last = [s for s in csv.DictReader((repo / "data/spoc_dv_centroid_summary.csv").open()) if "s0078" in s["filename"]][0]
x, y = float(last["ms_tic_mean_ra_offset_arcsec"]), float(last["ms_tic_mean_dec_offset_arcsec"])
a.errorbar(x, y, float(last["ms_tic_mean_dec_offset_unc_arcsec"]), float(last["ms_tic_mean_ra_offset_unc_arcsec"]),
           fmt="s", ms=6, color="C0", zorder=4, label="SPOC multi-sector")
summ = json.loads((repo / "rerun/outputs/SUMMARY.json").read_text())["difference_image_centroid"]
mx, my = summ["mean_offset_arcsec"]
se = summ["per_axis_err_arcsec"]
a.errorbar(mx, my, se, se, fmt="D", ms=5, color="C2", zorder=4, label="this work (10 sectors)")
a.plot(0, 0, "*", ms=14, color="gold", mec="k", zorder=5, label="TOI-1258")
a.plot(*T2, "X", ms=10, color="C3", mec="k", zorder=5, label=r"TIC 288132263 ($\Delta T$=4.4)")
a.axhline(0, color="0.85", lw=0.6, zorder=0)
a.axvline(0, color="0.85", lw=0.6, zorder=0)
a.set(xlim=(-16, 22), ylim=(-16, 18), aspect="equal",
      xlabel=r"$\Delta\alpha\cos\delta$ (arcsec, E+)", ylabel=r"$\Delta\delta$ (arcsec, N+)")
a.set_title("(c) Transit-source location", fontsize=10)
a.legend(loc="upper left", fontsize=6.5, frameon=True, framealpha=0.9, edgecolor="0.85")

fig.tight_layout()
fig.savefig(out, dpi=200)
print("wrote", out)
