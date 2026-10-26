"""Per-sector transit depth vs photometric aperture for TOI-1258.01.

For every SPOC sector: PDCSAP depth (crowding corrected), SAP depth from the
SPOC aperture, and SAP depths from target-centred apertures of growing radius
built from the TPF. Depths use a fixed ephemeris and a local linear baseline
around each transit, so no global detrending choice leaks between apertures.
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

P = 10.2613397
T0 = 3333.648786  # BTJD
DUR = 0.2023
HALF = DUR / 2
WIN = 0.6  # days either side of mid-transit used for the local baseline
RADII = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]


def depth_local(t, f, inner=0.8):
    """Mean in-transit depth (central `inner` fraction) against per-transit linear baselines."""
    good = np.isfinite(t) & np.isfinite(f)
    t, f = t[good], f[good]
    n = np.round((t - T0) / P)
    ph = t - (T0 + n * P)
    ins, outs = [], []
    per_tr = []
    for k in np.unique(n):
        m = (n == k) & (np.abs(ph) < WIN)
        oot = m & (np.abs(ph) > HALF + 0.03)
        it = m & (np.abs(ph) < inner * HALF)
        if oot.sum() < 60 or it.sum() < 30:
            continue
        # both sides of the transit must be present
        if (ph[oot] < 0).sum() < 20 or (ph[oot] > 0).sum() < 20:
            continue
        c = np.polyfit(ph[oot], f[oot], 1)
        r = f[m] / np.polyval(c, ph[m]) - 1
        rin = r[np.abs(ph[m]) < inner * HALF]
        rout = r[np.abs(ph[m]) > HALF + 0.03]
        ins.append(rin)
        outs.append(rout)
        per_tr.append(-np.mean(rin) * 1e6)
    if not ins:
        return np.nan, np.nan, 0, []
    rin = np.concatenate(ins)
    rout = np.concatenate(outs)
    sig = np.std(rout)
    d = -np.mean(rin) * 1e6
    e = sig / np.sqrt(len(rin)) * 1e6
    return d, e, len(ins), per_tr


def circ_mask(shape, cx, cy, r):
    yy, xx = np.indices(shape)
    return (xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2


def main(root, out):
    rows = []
    lcs = sorted(Path(root, "lc").rglob("*lc.fits"))
    tpfs = {fits.getheader(p)["SECTOR"]: p for p in Path(root, "tpf").rglob("*tp.fits")}
    for lcp in lcs:
        with fits.open(lcp) as h:
            hd0, hd1 = h[0].header, h[1].header
            sec = hd0["SECTOR"]
            d = h[1].data
            q = d["QUALITY"] == 0
            t = d["TIME"][q]
            pdc = d["PDCSAP_FLUX"][q]
            sap = d["SAP_FLUX"][q]
            crowd, flfr = hd1.get("CROWDSAP"), hd1.get("FLFRCSAP")
        row = dict(sector=sec, camera=hd0["CAMERA"], ccd=hd0["CCD"], crowdsap=crowd, flfrcsap=flfr)
        row["pdcsap_depth"], row["pdcsap_err"], row["n_transits"], _ = depth_local(t, pdc)
        row["sap_depth"], row["sap_err"], _, _ = depth_local(t, sap)
        with fits.open(tpfs[sec]) as h:
            d = h[1].data
            ap = h[2].data
            spoc_ap = (ap & 2) > 0
            q = (d["QUALITY"] == 0) & np.isfinite(d["TIME"])
            tt = d["TIME"][q]
            cube = d["FLUX"][q]
            cube = np.nan_to_num(cube)
            # target position from the WCS of the aperture extension
            from astropy.wcs import WCS
            w = WCS(h[2].header)
            cx, cy = w.all_world2pix([[h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"]]], 0)[0]
        row["spoc_ap_npix"] = int(spoc_ap.sum())
        row["target_px"] = [float(cx), float(cy)]
        f = cube[:, spoc_ap].sum(1)
        row["tpf_spoc_ap_depth"], row["tpf_spoc_ap_err"], _, _ = depth_local(tt, f)
        for r in RADII:
            m = circ_mask(cube.shape[1:], cx, cy, r)
            f = cube[:, m].sum(1)
            dd, ee, _, _ = depth_local(tt, f)
            row[f"r{r}_depth"], row[f"r{r}_err"], row[f"r{r}_npix"] = dd, ee, int(m.sum())
        rows.append(row)
        print(sec, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items() if "depth" in k}, flush=True)
    Path(out).write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
