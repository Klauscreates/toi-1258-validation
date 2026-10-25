"""Reconstruct FPP/NFPP from saved posterior scenario tables."""

from __future__ import annotations

import math
from pathlib import Path

from .common import number, read_csv

SCENARIOS = {
    "TP",
    "EB",
    "EBx2P",
    "PTP",
    "PEB",
    "PEBx2P",
    "STP",
    "SEB",
    "SEBx2P",
    "DTP",
    "DEB",
    "DEBx2P",
    "BTP",
    "BEB",
    "BEBx2P",
    "NTP",
    "NEB",
    "NEBx2P",
}


def scenario_probabilities(rows: list[dict]) -> tuple[float, float]:
    if not rows:
        raise ValueError("No scenario probabilities")
    probabilities = [(r["scenario"], number(r["prob"])) for r in rows]
    if any(s not in SCENARIOS or not 0 <= p <= 1 for s, p in probabilities):
        raise ValueError("Unknown scenario or invalid probability")
    if not math.isclose(sum(p for _, p in probabilities), 1.0, abs_tol=1e-8):
        raise ValueError("Scenario probabilities do not sum to one")
    fpp = sum(p for s, p in probabilities if s not in {"TP", "PTP", "DTP"})
    nfpp = sum(p for s, p in probabilities if s in {"NTP", "NEB", "NEBx2P"})
    return fpp, nfpp


def analyze(data: Path) -> tuple[dict, list[dict]]:
    seeds = read_csv(data / "triceratops_seed_sweep_N1000000_summary.csv")
    seen = set()
    audit = []
    for row in seeds:
        seed = int(row["seed"])
        if seed in seen or row["full_status"] != "COMPLETED":
            raise ValueError("Duplicate or incomplete seed")
        seen.add(seed)
        fpp, nfpp = scenario_probabilities(read_csv(data / f"probabilities_full_seed_{seed}.csv"))
        if not math.isclose(fpp, number(row["full_FPP"]), abs_tol=1e-10):
            raise ValueError("FPP does not match the saved seed summary")
        if not math.isclose(nfpp, number(row["full_NFPP"]), abs_tol=1e-10):
            raise ValueError("NFPP does not match the saved seed summary")
        audit.append(
            {
                "seed": seed,
                "samples": int(row["N"]),
                "FPP": fpp,
                "NFPP": nfpp,
                "below_reported_0p001_thresholds": fpp < 0.001 and nfpp < 0.001,
            }
        )
    bls = read_csv(data / "stitched_bls_summary.csv")[0]
    period = number(bls["best_bls_period_days"])
    catalog = number(bls["catalog_period_days"])
    if catalog <= 0:
        raise ValueError("Catalog period must be positive")
    return {
        "target": "TOI-1258.01",
        "tic": "288132261",
        "snapshot": "2026-07-14",
        "interpretation": "Archived posterior scenarios reproduce the reported low FPP/NFPP values under the historical model. This audit does not rerun TRICERATOPS, confirm a planet independently, measure mass, or verify current catalog status.",
        "seed_count": len(audit),
        "FPP_range": [min(r["FPP"] for r in audit), max(r["FPP"] for r in audit)],
        "NFPP_range": [min(r["NFPP"] for r in audit), max(r["NFPP"] for r in audit)],
        "all_below_reported_thresholds": all(r["below_reported_0p001_thresholds"] for r in audit),
        "bls_period_days": period,
        "catalog_period_days_snapshot": catalog,
        "relative_period_difference_percent": abs(period - catalog) / catalog * 100,
        "bls_depth_ppm": number(bls["best_bls_depth_ppm"]),
        "sectors": int(bls["n_sectors"]),
        "mass_measured_here": False,
        "current_archive_status_checked": False,
    }, audit
