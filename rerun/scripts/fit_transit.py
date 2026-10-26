"""MCMC transit fit of TOI-1258.01 to 12 sectors of SPOC 2-min PDCSAP photometry.

Each transit window (+/-0.5 d) is normalised by a linear out-of-transit baseline,
as in prep_inputs.py. Free parameters: T0, P, k = Rp/R*, b, stellar density
rho* (g/cm^3), Kipping q1/q2 limb darkening (uniform), and a white-noise term.
Circular orbit. The fitted rho* is compared with the TIC value as an
on-target consistency check.

Usage: python fit_transit.py DATA_DIR OUT_JSON [--steps N]
"""
import argparse
import json
from pathlib import Path

import emcee
import numpy as np
import scipy.integrate
from astropy.io import fits

# pytransit 2.2 (pinned by TRICERATOPS) uses names removed in NumPy 2 / SciPy 1.14;
# the same shims as triceratops/likelihoods.py, each a no-op where the name exists
for _name, _value in (("int", int), ("trapz", getattr(np, "trapezoid", None)), ("NaN", np.nan), ("Inf", np.inf)):
    if not hasattr(np, _name):
        setattr(np, _name, _value)
if not hasattr(scipy.integrate, "trapz"):
    scipy.integrate.trapz = scipy.integrate.trapezoid

from pytransit import QuadraticModel  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("data", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("--steps", type=int, default=6000)
ap.add_argument("--seed", type=int, default=1258)
args = ap.parse_args()

P0, T00 = 10.2613397, 3333.648786
HALF, WIN = 0.2023 / 2, 0.5
G = 6.674e-8
# TIC v8.2 host star (ExoFOP), used only for derived quantities and the density check
RSTAR, RSTAR_E = 1.49192, 0.0699799
RHO_TIC, RHO_TIC_E = 0.491264 * 1.41, 0.107459 * 1.41  # solar units -> g/cm^3

ts, fs = [], []
for p in sorted((args.data / "tess" / "lc").rglob("*lc.fits")):
    d = fits.getdata(p, 1)
    q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
    t, f = d["TIME"][q], d["PDCSAP_FLUX"][q]
    n = np.round((t - T00) / P0)
    ph = t - (T00 + n * P0)
    for k in np.unique(n):
        m = (n == k) & (np.abs(ph) < WIN)
        oot = m & (np.abs(ph) > HALF + 0.03)
        if oot.sum() < 60 or (oot & (ph < 0)).sum() < 20 or (oot & (ph > 0)).sum() < 20:
            continue
        c = np.polyfit(ph[oot], f[oot], 1)
        ts.append(t[m])
        fs.append(f[m] / np.polyval(c, ph[m]))
t = np.concatenate(ts).astype("d")
f = np.concatenate(fs).astype("d")
n_tr = len(ts)
# reference epoch near the middle of the baseline minimises the T0-P covariance
EPOCH_REF = int(np.round((np.median(t) - T00) / P0))
TREF = T00 + EPOCH_REF * P0

tm = QuadraticModel()
tm.set_data(t)
LABELS = ["T0", "P", "k", "b", "rho", "q1", "q2", "log_sigma"]


def unpack(th):
    t0, per, k, b, rho, q1, q2, lsig = th
    a = (G * rho * (per * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)  # a/R*
    inc = np.arccos(b / a)
    u1, u2 = 2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)
    return t0, per, k, b, rho, a, inc, u1, u2, np.exp(lsig)


def lnprob(th):
    t0, per, k, b, rho, q1, q2, lsig = th
    if not (0 < k < 0.2 and 0 <= b < 1 + k and 0.01 < rho < 10 and 0 < q1 < 1 and 0 < q2 < 1
            and -12 < lsig < -3 and abs(per - P0) < 0.01 and abs(t0 - TREF) < 0.1):
        return -np.inf
    t0, per, k, b, rho, a, inc, u1, u2, sig = unpack(th)
    if not np.isfinite(inc) or a <= 1:
        return -np.inf
    model = tm.evaluate(k, [u1, u2], t0, per, a, inc, 0.0, 0.0)
    r = f - model
    return -0.5 * np.sum(r ** 2 / sig ** 2 + 2 * np.log(sig))


p0 = np.array([TREF, P0, 0.028, 0.4, 0.7, 0.3, 0.35, np.log(np.std(f))])
rng = np.random.default_rng(args.seed)
nw = 40
start = p0 + rng.normal(0, 1, (nw, len(p0))) * np.array([5e-4, 2e-6, 5e-4, 0.05, 0.05, 0.05, 0.05, 0.01])
start[:, 3] = np.abs(start[:, 3])
sampler = emcee.EnsembleSampler(nw, len(p0), lnprob)
sampler.run_mcmc(start, args.steps, progress=False)
burn = args.steps // 2
chain = sampler.get_chain(discard=burn, flat=True)
try:
    tau = sampler.get_autocorr_time(discard=burn, quiet=True).tolist()
except Exception:  # noqa: BLE001  emcee raises when chains are short
    tau = None

der = np.array([unpack(th) for th in chain])
_, _, k, b, rho, a, inc, u1, u2, sig = der.T
rp = k * rng.normal(RSTAR, RSTAR_E, len(k)) * 109.076  # Earth radii
dur = (P0 / np.pi) * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / (a * np.sin(inc))) * 24


def summ(x):
    lo, mid, hi = np.percentile(x, [15.865, 50, 84.135])
    return dict(median=float(mid), minus=float(mid - lo), plus=float(hi - mid))


out = dict(
    n_points=int(len(t)), n_transits=n_tr, steps=args.steps, walkers=nw, burn=burn,
    autocorr_time=tau, acceptance=float(np.mean(sampler.acceptance_fraction)),
    T0_BTJD=summ(chain[:, 0]), P_days=summ(chain[:, 1]), k=summ(k), depth_ppm=summ(k ** 2 * 1e6),
    b=summ(b), rho_star_cgs=summ(rho), a_over_rstar=summ(a), inc_deg=summ(np.degrees(inc)),
    u1=summ(u1), u2=summ(u2), sigma_ppm=summ(sig * 1e6), T14_hours=summ(dur), Rp_earth=summ(rp),
    rho_tic_cgs=dict(value=RHO_TIC, err=RHO_TIC_E),
    rho_tension_sigma=float((np.median(rho) - RHO_TIC) / np.hypot(RHO_TIC_E, np.std(rho))),
)
args.out.write_text(json.dumps(out, indent=1) + "\n")
np.save(args.out.with_suffix(".chain.npy"), chain[:: max(1, len(chain) // 20000)])
print(json.dumps(out, indent=1))
