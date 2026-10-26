"""Per-sector difference-image centroids for TOI-1258.01.

Difference image = mean out-of-transit - mean in-transit, built per transit
from a local window and averaged. Its centroid (PRF-like 2D Gaussian fit on a
7x7 box) is compared with the out-of-transit image centroid and with the
catalog target position. Uncertainties come from null "transits": the same
difference image built at shifted epochs away from the real transit, added to
the measured difference image, and the resulting centroid scatter. Sectors
with a single usable transit have unreliable uncertainties.
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

P = 10.2613397
T0 = 3333.648786
HALF = 0.2023 / 2
WIN = 0.6
PIX_ARCSEC = 21.0


def gauss_centroid(img, x0, y0, box=3):
    ny, nx = img.shape
    xi, yi = int(round(x0)), int(round(y0))
    ys, xs = slice(max(yi - box, 0), min(yi + box + 1, ny)), slice(max(xi - box, 0), min(xi + box + 1, nx))
    sub = img[ys, xs]
    yy, xx = np.indices(sub.shape)
    xx = xx + xs.start
    yy = yy + ys.start

    def model(p):
        a, cx, cy, s, b = p
        return a * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * s ** 2)) + b

    p0 = [sub.max() - np.median(sub), x0, y0, 0.8, np.median(sub)]
    r = least_squares(lambda p: (model(p) - sub).ravel(), p0,
                      bounds=([0, x0 - box, y0 - box, 0.3, -np.inf], [np.inf, x0 + box, y0 + box, 3, np.inf]))
    return r.x[1], r.x[2]


def per_transit_diffs(t, cube, shift=0.0):
    n = np.round((t - T0 - shift) / P)
    ph = t - (T0 + shift + n * P)
    diffs, oots = [], []
    for k in np.unique(n):
        m = (n == k) & (np.abs(ph) < WIN)
        oot = m & (np.abs(ph) > HALF + 0.03)
        it = m & (np.abs(ph) < 0.8 * HALF)
        if it.sum() < 30 or oot.sum() < 60:
            continue
        two_sided = (oot & (ph < 0)).sum() >= 20 and (oot & (ph > 0)).sum() >= 20
        # linear-in-time baseline per pixel when both sides exist, else constant
        A = np.vstack([np.ones(oot.sum()), ph[oot]]).T if two_sided else np.ones((oot.sum(), 1))
        coef, *_ = np.linalg.lstsq(A, cube[oot].reshape(oot.sum(), -1), rcond=None)
        base = coef[0].reshape(cube.shape[1:])
        intr = cube[it].mean(0)
        diffs.append(base - intr)
        oots.append(base)
    return np.array(diffs), np.array(oots)


def main(root, out):
    res = []
    for p in sorted(Path(root, "tpf").rglob("*tp.fits")):
        with fits.open(p) as h:
            sec = h[0].header["SECTOR"]
            d = h[1].data
            q = (d["QUALITY"] == 0) & np.isfinite(d["TIME"])
            t = d["TIME"][q]
            cube = np.nan_to_num(d["FLUX"][q])
            w = WCS(h[2].header)
            tx, ty = w.all_world2pix([[h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"]]], 0)[0]
        diffs, oots = per_transit_diffs(t, cube)
        n = len(diffs)
        D, oot_img = diffs.mean(0), oots.mean(0)
        dx, dy = gauss_centroid(D, tx, ty)
        ox, oy = gauss_centroid(oot_img, tx, ty)
        # empirical noise: null "transits" at shifted epochs away from the real one
        nulls = []
        for sh in np.arange(0.7, P - 0.7, 0.35):
            nd, _ = per_transit_diffs(t, cube, sh)
            if len(nd) >= max(1, n - 1):
                nulls.append(nd.mean(0))
        trials = []
        for N in nulls:
            jx, jy = gauss_centroid(D + N, tx, ty)
            trials.append((jx - dx, jy - dy))
        trials = np.array(trials)
        err = trials.std(0) if len(trials) > 3 else np.array([np.nan, np.nan])
        off = np.array([dx - ox, dy - oy])
        # sky offset: convert pixel offset with the WCS
        ra0, de0 = w.all_pix2world([[ox, oy]], 0)[0]
        ra1, de1 = w.all_pix2world([[dx, dy]], 0)[0]
        dra = (ra1 - ra0) * np.cos(np.deg2rad(de0)) * 3600
        dde = (de1 - de0) * 3600
        sep = float(np.hypot(dra, dde))
        sep_err = float(np.hypot(*err) * PIX_ARCSEC)
        row = dict(sector=sec, n_transits=n, n_null=len(nulls), diff_minus_oot_px=off.tolist(), err_px=err.tolist(),
                   dRA_arcsec=float(dra), dDec_arcsec=float(dde), offset_arcsec=sep, offset_err_arcsec=sep_err,
                   oot_minus_catalog_px=[float(ox - tx), float(oy - ty)],
                   diff_peak_to_noise=float(D.max() / np.std(D[D < np.percentile(D, 70)])))
        res.append(row)
        print(f"S{sec:2} n={n:2}  diff-OOT offset = {sep:5.2f} +/- {sep_err:4.2f} arcsec "
              f"(dRA {dra:+.2f}, dDec {dde:+.2f})  OOT-catalog = ({ox - tx:+.2f},{oy - ty:+.2f}) px", flush=True)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
