"""Rerun historical posterior sampling from frozen stars, priors and photometry.

The saved flux-ratio/depth table is an input, not a fresh aperture reconstruction.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="1258,1259,1260,1261,1262")
    parser.add_argument("--samples", type=int, default=1_000_000)
    parser.add_argument("--binned-lightcurve", type=Path, help="Optional regenerated folded table from run_raw.py")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    if args.output.exists() or args.output.resolve().is_relative_to(repo):
        raise ValueError("Choose a new external output directory")
    import numpy as np
    import pandas as pd
    import scipy.integrate
    np.int = int
    np.float = float
    np.trapz = np.trapezoid
    scipy.integrate.trapz = scipy.integrate.trapezoid
    from triceratops import likelihoods

    def scalar(value):
        a = np.asarray(value)
        return float(a.reshape(-1)[0]) if a.size == 1 else value

    original_tp = likelihoods.simulate_TP_transit
    original_eb = likelihoods.simulate_EB_transit

    def tp(t, *values, **kwargs):
        return original_tp(t, *(scalar(v) if i < 10 else v for i, v in enumerate(values)), **kwargs)

    def eb(t, *values, **kwargs):
        return original_eb(t, *(scalar(v) if i < 11 else v for i, v in enumerate(values)), **kwargs)

    likelihoods.simulate_TP_transit = tp
    likelihoods.simulate_EB_transit = eb
    patch = repo / "research/third_party/triceratops/marginal_likelihoods_py312_compat_patched.py"
    spec = importlib.util.spec_from_file_location("triceratops.marginal_likelihoods", patch)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    from triceratops import triceratops as tri
    data = repo / "reproduce/historical-inputs"
    # Verify preserved input tables before using a frozen historical catalog.
    from raw_inputs import verify
    for row in json.loads((data / "PROVENANCE.json").read_text()):
        if row["status"] != "PRESERVED":
            raise ValueError("Missing historical model input")
        verify(data / row["name"], row)
    lightcurve = args.binned_lightcurve or data / "triceratops_folded_binned_lightcurve.csv"
    binned = pd.read_csv(lightcurve)
    from raw_inputs import digest
    lightcurve_sha = digest(lightcurve)
    oot = binned[np.abs(binned["time_from_mid_days"]) > 1.5 * 0.202292]
    flux_err = float(np.nanstd(oot["flux_rel_norm"]))
    period = json.loads((repo / "reproduce/target.json").read_text())["period_days"]
    args.output.mkdir(parents=True)
    results = []
    for seed in (int(s) for s in args.seeds.split(",")):
        np.random.seed(seed)
        target = tri.target.__new__(tri.target)
        target.ID = 288132261
        target.mission = "TESS"
        target.stars = pd.read_csv(data / "triceratops_sector14_searchradius2_stars.csv", dtype={"ID": str})
        target.trilegal_fname = str(data / "288132261_TRILEGAL.csv")
        start = time.monotonic()
        target.calc_probs(binned["time_from_mid_days"].to_numpy(float),
                          binned["flux_rel_norm"].to_numpy(float), flux_err, P_orb=period,
                          contrast_curve_file=str(data / "triceratops_contrast_curve_832nm.csv"),
                          filt="Vis", N=args.samples, parallel=False, verbose=1,
                          exptime=0.00139, nsamples=10, drop_scenario=[])
        target.probs.to_csv(args.output / f"probabilities_full_seed_{seed}.csv", index=False)
        results.append(dict(seed=seed, N=args.samples, FPP=float(target.FPP), NFPP=float(target.NFPP),
                            seconds=time.monotonic() - start))
        (args.output / "RUN.json").write_text(json.dumps(dict(
            scope="Frozen-input posterior resampling, not independent validation or fresh catalog retrieval",
            lightcurve_sha256=lightcurve_sha,
            runs=results), indent=2) + "\n")
        print(results[-1], flush=True)


if __name__ == "__main__":
    main()
