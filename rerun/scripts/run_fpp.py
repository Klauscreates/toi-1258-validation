"""TOI-1258.01 FPP/NFPP with current upstream TRICERATOPS.

Inputs: freshly folded SPOC PDCSAP light curve (dilution corrected), SPOC
apertures from the local TPFs, and the Gemini 562/832 nm speckle plus
Palomar Kcont AO contrast curves. Runs several seeds; optional frozen
TRILEGAL population for comparison with the Gaia default.
"""
import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ap = argparse.ArgumentParser()
ap.add_argument("--data", type=Path, required=True)
ap.add_argument("--inputs", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--N", type=int, default=1_000_000)
ap.add_argument("--seeds", default="1,2,3,4,5")
ap.add_argument("--trilegal", type=Path, default=None)
ap.add_argument("--curves", default="all", choices=["all", "none", "gemini832", "gemini"])
args = ap.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

import triceratops.triceratops as tr  # noqa: E402  (heavy import after argument parsing)

P = 10.2613397
lc = np.loadtxt(args.inputs / "folded_pdcsap_10min.csv", delimiter=",")
meta = json.loads((args.inputs / "inputs_meta.json").read_text())
depth = meta["central_depth_ppm"] * 1e-6

aps, sectors = {}, []
for p in sorted((args.data / "tess" / "tpf").rglob("*tp.fits")):
    with fits.open(p) as h:
        sec = h[0].header["SECTOR"]
        a = h["APERTURE"].data.astype(int)
        rows, cols = np.nonzero(a & 2)
        aps[sec] = np.column_stack([cols + h["APERTURE"].header["CRVAL1P"], rows + h["APERTURE"].header["CRVAL2P"]])
sectors = np.array(sorted(aps))

cache = args.out / "target.pkl"
if cache.exists():
    target = pickle.loads(cache.read_bytes())
else:
    kw = dict(trilegal_fname=str(args.trilegal)) if args.trilegal else {}
    target = tr.target(ID=288132261, sectors=sectors, search_radius=10, **kw)
    target.calc_depths(tdepth=depth, all_ap_pixels=[aps[s] for s in sectors])
    cache.write_bytes(pickle.dumps(target))
target.stars.to_csv(args.out / "stars.csv", index=False)

if args.curves == "all":
    cc = [str(args.inputs / f) for f in ("cc_gemini_562.csv", "cc_gemini_832.csv", "cc_palomar_K.csv")]
    filt = ["Vis", "TESS", "K"]
elif args.curves == "gemini":
    cc = [str(args.inputs / f) for f in ("cc_gemini_562.csv", "cc_gemini_832.csv")]
    filt = ["Vis", "TESS"]
elif args.curves == "gemini832":
    cc, filt = str(args.inputs / "cc_gemini_832.csv"), "TESS"
else:
    cc, filt = None, "TESS"

results = []
for seed in (int(s) for s in args.seeds.split(",")):
    np.random.seed(seed)
    t0 = time.monotonic()
    target.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=meta["flux_err_per_bin"], P_orb=P,
                      contrast_curve_file=cc, filt=filt, N=args.N, parallel=False, verbose=0,
                      exptime=meta["bin_minutes"] / 1440, nsamples=20)
    target.probs.to_csv(args.out / f"probs_seed{seed}.csv", index=False)
    r = dict(seed=seed, N=args.N, FPP=float(target.FPP), NFPP=float(target.NFPP), seconds=round(time.monotonic() - t0))
    results.append(r)
    print(r, flush=True)
    (args.out / "summary.json").write_text(json.dumps(dict(
        triceratops_commit="d20875c27d", curves=args.curves, filters=filt,
        background="trilegal-frozen" if args.trilegal else "gaia-dr3", depth_ppm=depth * 1e6,
        sectors=sectors.tolist(), runs=results), indent=1))
