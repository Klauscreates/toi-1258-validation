"""Build TRICERATOPS inputs: folded PDCSAP light curve and contrast curves.

Also checks the historical contrast-curve conversion against the ExoFOP originals.
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

P = 10.2613397
T0 = 3333.648786
HALF = 0.2023 / 2
WIN = 0.5
BIN_MIN = 10.0

root, hist, out = (Path(a) for a in sys.argv[1:4])
out.mkdir(parents=True, exist_ok=True)
ex = root / "exofop"

# --- contrast curves --------------------------------------------------------
def gemini_curve(path):
    rows, fit = [], False
    for line in path.read_text().splitlines():
        if line.startswith("# fit of angular separation"):
            fit = True
            continue
        if fit and line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    return np.array(rows)


check = {}
for band in ("562", "832"):
    g = gemini_curve(ex / f"TOI1258I-ef20200215-{band}_sensitivity.dat")
    np.savetxt(out / f"cc_gemini_{band}.csv", g, delimiter=",", fmt="%.6f")
    h = np.loadtxt(hist / f"triceratops_contrast_curve_{band}nm.csv", delimiter=",")
    check[band] = dict(
        exofop_sep_range=[float(g[:, 0].min()), float(g[:, 0].max())],
        exofop_dmag_at_0p5=float(np.interp(0.5, g[:, 0], g[:, 1])),
        historical_col0_equals_col1=bool(np.allclose(h[:, 0], h[:, 1])),
        historical_col0_equals_exofop_dmag=bool(np.allclose(h[:, 0], g[:, 1])),
        historical_dmag_at_0p5=float(np.interp(0.5, h[:, 0], h[:, 1])),
    )
pal = []
for line in (ex / "TOI1258I-dc20260407-Kcont_plot.tbl").read_text().splitlines():
    parts = [p.strip() for p in line.split(",")]
    if len(parts) == 3:
        try:
            pal.append([float(parts[0]), float(parts[1])])
        except ValueError:
            pass
pal = np.array(pal)
np.savetxt(out / "cc_palomar_K.csv", pal, delimiter=",", fmt="%.6f")
check["palomar_K"] = dict(sep_range=[float(pal[:, 0].min()), float(pal[:, 0].max())],
                          dmag_at_0p5=float(np.interp(0.5, pal[:, 0], pal[:, 1])))

# --- folded PDCSAP light curve ---------------------------------------------
ts, fs = [], []
for p in sorted((root / "tess" / "lc").rglob("*lc.fits")):
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
        ts.append(ph[m])
        fs.append(f[m] / np.polyval(c, ph[m]))
t = np.concatenate(ts)
f = np.concatenate(fs)
edges = np.arange(-WIN, WIN + 1e-9, BIN_MIN / 1440)
idx = np.digitize(t, edges) - 1
bt, bf, bn = [], [], []
for i in range(len(edges) - 1):
    m = idx == i
    if m.sum() >= 10:
        bt.append(t[m].mean())
        bf.append(f[m].mean())
        bn.append(m.sum())
bt, bf = np.array(bt), np.array(bf)
oot = np.abs(bt) > HALF + 0.03
err = float(np.std(bf[oot]))
depth = float(1 - bf[np.abs(bt) < 0.8 * HALF].mean())
np.savetxt(out / "folded_pdcsap_10min.csv", np.c_[bt, bf], delimiter=",", header="time_from_mid_days,flux", fmt="%.8f")
meta = dict(n_points_unbinned=int(len(t)), n_bins=int(len(bt)), bin_minutes=BIN_MIN,
            flux_err_per_bin=err, central_depth_ppm=depth * 1e6, contrast_checks=check)
(out / "inputs_meta.json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta, indent=1))
