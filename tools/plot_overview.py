"""Render a measured-data overview, with an explicit claim boundary."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toi1258.analysis import analyze
from toi1258.common import number, read_csv, verify_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("figures/overview.png"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    verify_manifest(root)
    analyze(root / "data")
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    fig.suptitle("TOI-1258.01 | Transit Validation", fontsize=17, fontweight="bold", color="#18181b")
    summary, rows = analyze(root / "data")
    axes[0].plot([r["seed"] for r in rows], [r["FPP"] for r in rows], "o-", color="#087f8c", label="FPP")
    axes[0].plot([r["seed"] for r in rows], [r["NFPP"] for r in rows], "s-", color="#a63a50", label="NFPP")
    axes[0].axhline(0.001, color="#71717a", ls="--", label="Reported threshold")
    axes[0].set(
        yscale="log",
        xlabel="Random seed",
        ylabel="Saved posterior probability",
        title="Full-scenario posterior reconstruction",
    )
    axes[0].ticklabel_format(axis="x", useOffset=False)
    axes[0].set_xticks([r["seed"] for r in rows])
    axes[0].legend(fontsize=8)
    folded = read_csv(root / "data/triceratops_folded_binned_lightcurve.csv")
    axes[1].errorbar(
        [number(r["time_from_mid_days"]) for r in folded],
        [(number(r["flux_rel_norm"]) - 1) * 1e6 for r in folded],
        yerr=[number(r["flux_err_rel"]) * 1e6 for r in folded],
        fmt=".",
        ms=4,
        color="#087f8c",
        ecolor="#c1d9da",
        lw=0.7,
    )
    axes[1].set(
        xlim=(-0.4, 0.4),
        xlabel="Days from mid-transit",
        ylabel="Relative flux - 1 (ppm)",
        title="Historical folded, binned light curve",
    )
    footnote = "July 2026 evidence snapshot. Probabilities are conditional on the historical model; no RV mass measured here."

    for ax in axes:
        ax.grid(alpha=0.15)
    fig.supxlabel(footnote, fontsize=8, color="#52525b")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, facecolor="white")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
