from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


VERDICT_UNRESOLVED = "unresolved_but_constrained"


def build_claim_gate(
    mission2_dir: str | Path,
    mission4_dir: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Merge JWST audit artifacts and target mining into a claim-gate packet."""

    mission2_path = Path(mission2_dir)
    mission4_path = Path(mission4_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    direct_summary = _read_json(mission2_path / "jwst_archive_summary.json")
    audit_rows = _read_audit_rows(mission2_path)
    download_rows = _read_download_rows(mission2_path)
    toi_rows = _read_csv_rows(mission4_path / "toi_candidate_shortlist.csv")
    confirmed_rows = _read_csv_rows(mission4_path / "confirmed_atmosphere_targets.csv")
    tess_rows = _read_csv_rows(mission4_path / "tess_archive_observations.csv")

    feature_rows = [_feature_row(row) for row in audit_rows]
    direct_status = _direct_status(direct_summary)
    download_status_counts = Counter(row.get("status", "unknown") or "unknown" for row in download_rows)
    instrument_counts = Counter(row.get("instrume", "UNKNOWN") or "UNKNOWN" for row in audit_rows)
    target_counts = Counter(row.get("targname", "UNKNOWN") or "UNKNOWN" for row in audit_rows)

    claim_level = (
        "proxima_spectroscopy_ready"
        if direct_status == "retrievable_proxima_spectroscopy"
        else "proxy_pipeline_validated_only"
    )
    verdict = (
        "requires_retrieval"
        if direct_status == "retrievable_proxima_spectroscopy"
        else VERDICT_UNRESOLVED
    )

    packet: dict[str, Any] = {
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "verdict": verdict,
        "claim_level": claim_level,
        "go_probe": False,
        "direct_proxima": {
            "status": direct_status,
            "observation_count": _as_int(direct_summary.get("observation_count")),
            "spectroscopy_observation_count": _as_int(
                direct_summary.get("spectroscopy_observation_count")
            ),
            "science_product_count": _as_int(direct_summary.get("science_product_count")),
            "nirspec_or_niriss_science_product_count": _as_int(
                direct_summary.get("nirspec_or_niriss_science_product_count")
            ),
            "retrievable_product_count": _as_int(
                direct_summary.get("retrievable_product_count")
            ),
            "verdict_impact": direct_summary.get("verdict_impact", ""),
        },
        "analog_audit": {
            "audited_product_count": len(audit_rows),
            "instrument_counts": dict(sorted(instrument_counts.items())),
            "target_counts": dict(sorted(target_counts.items())),
            "feature_rows": feature_rows,
        },
        "download_budget": {
            "attempt_count": len(download_rows),
            "status_counts": dict(sorted(download_status_counts.items())),
            "large_products_blocked": [
                row
                for row in download_rows
                if row.get("status") in {"skipped_size_limit", "error"}
                and _as_int(row.get("product_size")) >= 100_000_000
            ],
        },
        "molecule_target_mining": {
            "top_toi_candidate": toi_rows[0] if toi_rows else {},
            "top_confirmed_atmosphere_target": confirmed_rows[0] if confirmed_rows else {},
            "toi_shortlist_count": len(toi_rows),
            "confirmed_target_count": len(confirmed_rows),
            "tess_observation_count": len(tess_rows),
        },
        "claims_allowed": _claims_allowed(claim_level, direct_status, feature_rows),
        "claims_blocked": _claims_blocked(claim_level),
        "next_actions": _next_actions(direct_status, download_status_counts),
    }

    _write_csv(output_path / "feature_screen_summary.csv", feature_rows)
    (output_path / "claim_gate.json").write_text(
        json.dumps(packet, indent=2),
        encoding="utf-8",
    )
    (output_path / "feature_screen_report.md").write_text(
        render_claim_gate_report(packet),
        encoding="utf-8",
    )
    return packet


def render_claim_gate_report(packet: dict[str, Any]) -> str:
    direct = packet["direct_proxima"]
    audit = packet["analog_audit"]
    mining = packet["molecule_target_mining"]
    download = packet["download_budget"]
    top_toi = mining.get("top_toi_candidate") or {}
    top_known = mining.get("top_confirmed_atmosphere_target") or {}

    lines = [
        "# Mission 2 Claim Gate: Proxima b Habitability Verdict",
        "",
        f"Generated: {packet['generated_at_utc']}",
        f"Verdict: {packet['verdict']}",
        f"Claim level: {packet['claim_level']}",
        f"Future probe gate: {'go' if packet['go_probe'] else 'no-go'}",
        "",
        "## Direct Proxima JWST Gate",
        "",
        f"- Status: {direct['status']}",
        f"- JWST observations: {direct['observation_count']}",
        f"- Spectroscopy observations: {direct['spectroscopy_observation_count']}",
        f"- Science products: {direct['science_product_count']}",
        f"- NIRSpec/NIRISS science products: {direct['nirspec_or_niriss_science_product_count']}",
        f"- Retrievable products: {direct['retrievable_product_count']}",
        f"- Verdict impact: {direct['verdict_impact']}",
        "",
        "## Audited Proxy Data",
        "",
        f"- Audited JWST products: {audit['audited_product_count']}",
        f"- Instruments: {_format_counts(audit['instrument_counts'])}",
        f"- Targets: {_format_counts(audit['target_counts'])}",
        "",
        "| Product | Target | Instrument | Setup | Inferred feature use | Claim use |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in audit["feature_rows"]:
        lines.append(
            "| {filename} | {target} | {instrument} | {setup} | {features} | {claim_use} |".format(
                filename=row["filename"],
                target=row["target"],
                instrument=row["instrument"],
                setup=row["setup"],
                features=row["inferred_feature_use"],
                claim_use=row["claim_use"],
            )
        )

    lines.extend(
        [
            "",
            "## Download Budget",
            "",
            f"- Download attempts: {download['attempt_count']}",
            f"- Status counts: {_format_counts(download['status_counts'])}",
            f"- Large products blocked by disk/size budget: {len(download['large_products_blocked'])}",
            "",
            "## Methane/Water Mining Feed",
            "",
            f"- TOI shortlist rows: {mining['toi_shortlist_count']}",
            f"- Confirmed atmosphere targets: {mining['confirmed_target_count']}",
            f"- TESS observations found for top TOIs: {mining['tess_observation_count']}",
            f"- Top discovery candidate: {_top_toi_label(top_toi)}",
            f"- Top known atmosphere target: {_top_known_label(top_known)}",
            "",
            "## Claims Allowed Now",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in packet["claims_allowed"])
    lines.extend(["", "## Claims Blocked Now", ""])
    lines.extend(f"- {item}" for item in packet["claims_blocked"])
    lines.extend(["", "## Next Actions", ""])
    lines.extend(f"{index}. {item}" for index, item in enumerate(packet["next_actions"], start=1))
    lines.append("")
    return "\n".join(lines)


def _read_audit_rows(mission2_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for audit_path in sorted(mission2_path.glob("spectral_audit*")):
        if not audit_path.is_dir():
            continue
        for row in _read_csv_rows(audit_path / "fits_inventory.csv"):
            row["audit_dir"] = str(audit_path)
            rows.append(row)
    return rows


def _read_download_rows(mission2_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for download_path in sorted(mission2_path.glob("calibration_analogs*/jwst_downloads.json")):
        for row in _read_json(download_path, default=[]):
            row["download_manifest"] = str(download_path)
            rows.append(row)
    return rows


def _feature_row(row: dict[str, Any]) -> dict[str, str]:
    instrument = row.get("instrume", "") or "UNKNOWN"
    exp_type = row.get("exp_type", "") or ""
    grating = row.get("grating", "") or ""
    filter_name = row.get("filter", "") or ""
    setup = " ".join(item for item in (exp_type, filter_name, grating) if item).strip()

    return {
        "filename": row.get("filename", ""),
        "target": row.get("targname", "") or row.get("targprop", "") or "UNKNOWN",
        "instrument": instrument,
        "setup": setup or "UNKNOWN",
        "date_obs": row.get("date_obs", ""),
        "mast_data_uri": row.get("mast_data_uri", ""),
        "inferred_feature_use": _infer_feature_use(instrument, exp_type, grating, filter_name),
        "claim_use": _claim_use(instrument, exp_type, grating),
    }


def _infer_feature_use(instrument: str, exp_type: str, grating: str, filter_name: str) -> str:
    setup = " ".join([instrument, exp_type, grating, filter_name]).upper()
    if "NIRSPEC" in setup and "PRISM" in setup:
        return "near-IR H2O/CH4/CO2/CO method validation"
    if "MIRI" in setup and ("P750L" in setup or "SLITLESS" in setup):
        return "mid-IR thermal/spectral pipeline validation"
    if "MIRI" in setup:
        return "mid-IR target-acquisition/support sanity check"
    return "generic detector/profile sanity check"


def _claim_use(instrument: str, exp_type: str, grating: str) -> str:
    setup = " ".join([instrument, exp_type, grating]).upper()
    if "NIRSPEC" in setup and "PRISM" in setup:
        return "proxy spectral extraction only; not a Proxima molecule claim"
    if "MIRI" in setup and "TACQ" in setup:
        return "support data only; not science spectroscopy"
    if "MIRI" in setup and ("SLITLESS" in setup or "P750L" in setup):
        return "proxy mid-IR method development only"
    return "audit provenance only"


def _direct_status(summary: dict[str, Any]) -> str:
    if _as_int(summary.get("nirspec_or_niriss_science_product_count")) > 0:
        return "retrievable_proxima_spectroscopy"
    if _as_int(summary.get("spectroscopy_observation_count")) > 0:
        return "proxima_spectroscopy_metadata_only"
    if _as_int(summary.get("observation_count")) > 0:
        return "proxima_metadata_no_retrievable_spectroscopy"
    return "no_proxima_jwst_match"


def _claims_allowed(
    claim_level: str,
    direct_status: str,
    feature_rows: list[dict[str, str]],
) -> list[str]:
    if claim_level == "proxima_spectroscopy_ready":
        return [
            "A Proxima-specific retrieval run is justified because retrievable spectroscopy is present.",
            "Habitability probabilities may be updated after wavelength-calibrated extraction and stellar-contamination controls.",
        ]
    claims = [
        "Proxima b remains scientifically unresolved with current local archive products.",
        "The archive search found public Proxima JWST metadata but no retrievable NIRSpec/NIRISS spectroscopy products.",
        "Audited analog products can validate FITS handling, provenance capture, detector profiles, and cleanup behavior.",
    ]
    if any("near-IR" in row.get("inferred_feature_use", "") for row in feature_rows):
        claims.append(
            "TRAPPIST-1 NIRSpec PRISM analog data can exercise the near-IR H2O/CH4/CO2/CO screening path."
        )
    if direct_status == "proxima_metadata_no_retrievable_spectroscopy":
        claims.append("The next Proxima gate is product release or a new observing proposal, not a rover/probe launch.")
    return claims


def _claims_blocked(claim_level: str) -> list[str]:
    if claim_level == "proxima_spectroscopy_ready":
        return [
            "No life detection claim is allowed from spectroscopy alone.",
            "No probe go-decision is allowed until retrieval posteriors cross the configured evidence threshold.",
        ]
    return [
        "No Proxima H2O, CH4, CO2, or CO detection claim.",
        "No habitable or uninhabitable verdict beyond unresolved-but-constrained.",
        "No thermal phase-curve claim without Proxima or properly downloaded analog time-series spectroscopy.",
        "No probe/rover go-decision from analog or calibration products.",
    ]


def _next_actions(
    direct_status: str,
    download_status_counts: Counter[str],
) -> list[str]:
    actions = [
        "Rerun the direct Proxima JWST archive search on a cadence and ingest any newly released products.",
        "Keep the FITS-prune audit path as the default so downloaded products become compact CSV/PNG artifacts immediately.",
    ]
    if download_status_counts.get("error", 0) or download_status_counts.get("skipped_size_limit", 0):
        actions.append(
            "Use a larger external cache or raise the download cap only for selected x1dints products needed for wavelength-calibrated retrieval tests."
        )
    if direct_status != "retrievable_proxima_spectroscopy":
        actions.append(
            "Prepare the JWST/ELT follow-up justification from the unresolved Proxima gate and the validated proxy pipeline."
        )
    return actions


def _read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path, default: Any | None = None) -> Any:
    if not path.exists():
        return {} if default is None else default
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    if not fields:
        fields = ["status"]
        rows = [{"status": "no_rows"}]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _as_int(value: Any) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _format_counts(counts: dict[str, Any]) -> str:
    if not counts:
        return "none"
    return ", ".join(f"{key}={value}" for key, value in counts.items())


def _top_toi_label(row: dict[str, Any]) -> str:
    if not row:
        return "none"
    name = row.get("working_name") or "unnamed"
    score = row.get("score") or row.get("total_score") or "unknown"
    return (
        f"{name} TOI {row.get('toi', 'unknown')} TIC {row.get('tic_id', 'unknown')} "
        f"score={score}"
    )


def _top_known_label(row: dict[str, Any]) -> str:
    if not row:
        return "none"
    score = row.get("molecule_score") or row.get("score") or row.get("total_score") or "unknown"
    return (
        f"{row.get('planet_name', 'unknown')} around {row.get('host_name', 'unknown')} "
        f"score={score}"
    )
