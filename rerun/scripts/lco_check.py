"""LCO 2022-05-08 zs light curve: on-target depth and a direct test of neighbour T2.

Depth: box-shaped transit with the TESS duration, free mid-time, fitted jointly
with a linear baseline in tot_C_cnts (the observer's detrend) and airmass.
T2: the same model applied to T2's own light curve, compared with the ~4.6%
eclipse T2 would need in order to produce the TESS signal.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

tbl = Path(sys.argv[1])
lines = tbl.read_text().splitlines()
cols = lines[0].split("\t")
def num(v):
    try:
        return float(v)
    except ValueError:
        return np.nan


data = np.array([[num(v) for v in ln.split("\t")] for ln in lines[1:]])
col = {c: data[:, i] for i, c in enumerate(cols)}

t = col["BJD_TDB"] - 2457000
DUR = 4.855 / 24  # ExoFOP duration
TC_PRED = 2459707.7034 - 2457000


def box(tt, tc, dur=DUR, ing=0.35 / 24):
    """Trapezoid: 1 inside, 0 outside, linear ingress/egress of width `ing`."""
    x = np.abs(tt - tc)
    return np.clip((dur / 2 - x) / ing + 0.5, 0, 1)


def fit(flux, err):
    good = np.isfinite(flux) & np.isfinite(col["tot_C_cnts"]) & np.isfinite(col["AIRMASS"])
    tt, f, e = t[good], flux[good] / np.nanmedian(flux[good]), err[good] / np.nanmedian(flux[good])
    c1 = (col["tot_C_cnts"][good] / np.nanmedian(col["tot_C_cnts"][good])) - 1
    am = col["AIRMASS"][good] - np.nanmean(col["AIRMASS"][good])

    def model(p, tc=None):
        d, tcp, a0, a1, a2 = p
        return a0 * (1 + a1 * c1 + a2 * am) * (1 - d * box(tt, tcp if tc is None else tc))

    r = least_squares(lambda p: (model(p) - f) / e, [1e-3, TC_PRED, 1, 0, 0])
    res = f - model(r.x)
    # inflate errors to reduced chi2 = 1, and add a time-correlated (beta) factor from 10-min binning
    scale = np.sqrt(np.sum((res / e) ** 2) / (len(f) - 5))
    nb = max(1, int(round(10 / (np.median(np.diff(tt)) * 1440))))
    m = len(res) // nb
    binned = res[: m * nb].reshape(m, nb).mean(1)
    beta = max(1.0, np.std(binned) / (np.std(res) / np.sqrt(nb)))
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * (scale * beta) ** 2
    return dict(depth_ppt=r.x[0] * 1e3, depth_err_ppt=float(np.sqrt(cov[0, 0]) * 1e3),
                tc_btjd=float(r.x[1]), tc_err_min=float(np.sqrt(cov[1, 1]) * 1440),
                rms_ppt=float(np.std(res) * 1e3), beta=float(beta), n=int(len(f)))


out = dict(T1=fit(col["rel_flux_T1"], col["rel_flux_err_T1"]))
# T2: fit with the mid-time fixed to the T1 solution (an NEB would eclipse at the same time)
tc1 = out["T1"]["tc_btjd"]
g = np.isfinite(col["rel_flux_T2"])
f2 = col["rel_flux_T2"] / np.nanmedian(col["rel_flux_T2"])
e2 = col["rel_flux_err_T2"] / np.nanmedian(col["rel_flux_T2"])
inside = box(t, tc1) > 0.99
outside = box(t, tc1) < 0.01
out["T2_in_minus_out_ppt"] = float((np.nanmean(f2[outside]) - np.nanmean(f2[inside])) * 1e3)
out["T2_point_rms_ppt"] = float(np.nanstd(f2[outside]) * 1e3)
out["T2_fit"] = fit(col["rel_flux_T2"], col["rel_flux_err_T2"])
# TRICERATOPS tdepth for TIC 288132263 (TESS band), from results/fpp_main/stars.csv
out["T2_required_depth_ppt_if_source"] = 45.8
print(json.dumps(out, indent=1))
Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
