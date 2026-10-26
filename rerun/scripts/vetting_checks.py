"""Odd/even depth comparison and secondary-eclipse search for TOI-1258.01.

Uses SPOC PDCSAP 2-min data with the same per-event local linear baseline as
the other rerun scripts. Depth = mean flux decrement over the central 80% of
the transit duration; uncertainty from the out-of-event scatter. The secondary
search uses a box of the same duration at phase 0.5 (circular orbit).

Usage: python vetting_checks.py DATA_DIR TRANSIT_FIT_JSON OUT_JSON
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

data, fitp, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
fit = json.loads(fitp.read_text())
P, T0 = fit["P_days"]["median"], fit["T0_BTJD"]["median"]
HALF = fit["T14_hours"]["median"] / 48
WIN = 0.5


def events(shift):
    """Per-event (epoch, in-event residuals, out-of-event residuals) around T0 + shift."""
    res = []
    for p in sorted((data / "tess" / "lc").rglob("*lc.fits")):
        d = fits.getdata(p, 1)
        q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
        t, f = d["TIME"][q], d["PDCSAP_FLUX"][q]
        n = np.round((t - T0 - shift) / P)
        ph = t - (T0 + shift + n * P)
        for k in np.unique(n):
            m = (n == k) & (np.abs(ph) < WIN)
            oot = m & (np.abs(ph) > HALF + 0.03)
            it = m & (np.abs(ph) < 0.8 * HALF)
            if it.sum() < 30 or (oot & (ph < 0)).sum() < 20 or (oot & (ph > 0)).sum() < 20:
                continue
            c = np.polyfit(ph[oot], f[oot], 1)
            r = f / np.polyval(c, ph) - 1
            res.append((int(k), r[it], r[oot]))
    return res


def depth(evts):
    rin = np.concatenate([e[1] for e in evts])
    rout = np.concatenate([e[2] for e in evts])
    return float(-rin.mean() * 1e6), float(rout.std() / math.sqrt(len(rin)) * 1e6), len(evts)


prim = events(0.0)
odd = [e for e in prim if e[0] % 2]
even = [e for e in prim if not e[0] % 2]
do, eo, no = depth(odd)
de, ee, ne = depth(even)
ds, es, ns = depth(events(P / 2))
result = dict(
    primary_ppm=depth(prim)[:2], odd_ppm=[do, eo], n_odd=no, even_ppm=[de, ee], n_even=ne,
    odd_even_difference_sigma=abs(do - de) / math.hypot(eo, ee),
    secondary_ppm=[ds, es], n_secondary_windows=ns,
    secondary_3sigma_upper_ppm=max(ds, 0) + 3 * es,
)
out.write_text(json.dumps(result, indent=1) + "\n")
print(json.dumps(result, indent=1))
