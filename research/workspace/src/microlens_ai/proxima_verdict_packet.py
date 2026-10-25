from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


PROVEN_HABITABLE = "proven_habitable"
PROVEN_UNINHABITABLE = "proven_uninhabitable"
UNRESOLVED = "unresolved_but_constrained"


def build_proxima_verdict_packet(
    mission_package: dict[str, Any],
    archive_summary: dict[str, Any],
) -> dict[str, Any]:
    """Build a Proxima-only habitability verdict from mission and archive state."""

    verdict = mission_package["verdict"]
    probabilities = verdict["probabilities"]
    direct_spectroscopy_products = int(
        archive_summary.get("nirspec_or_niriss_science_product_count", 0) or 0
    )
    retrievable_products = int(archive_summary.get("retrievable_product_count", 0) or 0)
    spectroscopy_observations = int(archive_summary.get("spectroscopy_observation_count", 0) or 0)
    atmosphere_data_ready = direct_spectroscopy_products > 0 and retrievable_products > 0

    proof_status = _proof_status(verdict["label"], probabilities, atmosphere_data_ready)
    packet = {
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "mission": "mission2_proxima_b_habitability",
        "target": mission_package["target"],
        "verdict": {
            "label": verdict["label"],
            "proof_status": proof_status,
            "probabilities": probabilities,
            "go_probe": bool(verdict.get("go_probe", False)),
        },
        "direct_jwst_evidence": {
            "observation_count": int(archive_summary.get("observation_count", 0) or 0),
            "public_observation_count": int(
                archive_summary.get("public_observation_count", 0) or 0
            ),
            "imaging_observation_count": int(
                archive_summary.get("imaging_observation_count", 0) or 0
            ),
            "spectroscopy_observation_count": spectroscopy_observations,
            "science_product_count": int(archive_summary.get("science_product_count", 0) or 0),
            "retrievable_product_count": retrievable_products,
            "nirspec_or_niriss_science_product_count": direct_spectroscopy_products,
            "miri_observation_count": int(archive_summary.get("miri_observation_count", 0) or 0),
            "nirspec_observation_count": int(
                archive_summary.get("nirspec_observation_count", 0) or 0
            ),
            "instruments": list(archive_summary.get("instruments", [])),
            "filters": list(archive_summary.get("filters", [])),
            "proposal_ids": list(archive_summary.get("proposal_ids", [])),
            "verdict_impact": archive_summary.get("verdict_impact", ""),
        },
        "proof_standard": {
            "habitable_requires": [
                "direct or strongly constraining atmosphere-retention evidence",
                "climate or chemistry constraints compatible with surface liquid water",
                "stellar activity and contamination controls",
                "data quality high enough to move posterior probability over the decision threshold",
            ],
            "uninhabitable_requires": [
                "strong atmosphere-loss, airless, runaway-greenhouse, desiccation, or stellar-stripping evidence",
                "data quality high enough to move posterior probability over the decision threshold",
            ],
        },
        "claims_allowed": _claims_allowed(proof_status, atmosphere_data_ready),
        "claims_blocked": _claims_blocked(proof_status),
        "next_actions": _next_actions(atmosphere_data_ready),
    }
    return packet


def write_proxima_verdict_packet(
    mission_package_path: str | Path,
    archive_summary_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Path]:
    mission_package = json.loads(Path(mission_package_path).read_text(encoding="utf-8"))
    archive_summary = json.loads(Path(archive_summary_path).read_text(encoding="utf-8"))
    packet = build_proxima_verdict_packet(mission_package, archive_summary)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    json_path = output_path / "proxima_b_verdict_packet.json"
    report_path = output_path / "proxima_b_verdict_packet.md"
    json_path.write_text(json.dumps(packet, indent=2), encoding="utf-8")
    report_path.write_text(render_proxima_verdict_packet(packet), encoding="utf-8")
    return {"json": json_path, "report": report_path}


def render_proxima_verdict_packet(packet: dict[str, Any]) -> str:
    target = packet["target"]
    verdict = packet["verdict"]
    archive = packet["direct_jwst_evidence"]
    probabilities = verdict["probabilities"]
    lines = [
        "# Proxima b Habitability Verdict Packet",
        "",
        f"Generated: {packet['generated_at_utc']}",
        f"Target: {target['planet_name']} around {target['star_name']}",
        f"Verdict: {verdict['label']}",
        f"Proof status: {verdict['proof_status']}",
        f"Future probe gate: {'go' if verdict['go_probe'] else 'no-go'}",
        "",
        "## Posterior State",
        "",
        f"- likely_habitable: {probabilities['likely_habitable']:.3f}",
        f"- likely_uninhabitable: {probabilities['likely_uninhabitable']:.3f}",
        f"- unresolved_but_constrained: {probabilities['unresolved_but_constrained']:.3f}",
        "",
        "## Direct JWST Evidence",
        "",
        f"- observations: {archive['observation_count']}",
        f"- public observations: {archive['public_observation_count']}",
        f"- imaging observations: {archive['imaging_observation_count']}",
        f"- spectroscopy observations: {archive['spectroscopy_observation_count']}",
        f"- science products: {archive['science_product_count']}",
        f"- retrievable products: {archive['retrievable_product_count']}",
        f"- NIRSpec/NIRISS science products: {archive['nirspec_or_niriss_science_product_count']}",
        f"- MIRI observations: {archive['miri_observation_count']}",
        f"- NIRSpec observations: {archive['nirspec_observation_count']}",
        f"- instruments: {_join_or_none(archive['instruments'])}",
        f"- filters: {_join_or_none(archive['filters'])}",
        f"- proposal ids: {_join_or_none(archive['proposal_ids'])}",
        f"- verdict impact: {archive['verdict_impact']}",
        "",
        "## Claims Allowed",
        "",
    ]
    lines.extend(f"- {item}" for item in packet["claims_allowed"])
    lines.extend(["", "## Claims Blocked", ""])
    lines.extend(f"- {item}" for item in packet["claims_blocked"])
    lines.extend(["", "## Next Actions", ""])
    lines.extend(f"{index}. {item}" for index, item in enumerate(packet["next_actions"], start=1))
    lines.append("")
    return "\n".join(lines)


def _proof_status(
    label: str,
    probabilities: dict[str, float],
    atmosphere_data_ready: bool,
) -> str:
    if label == PROVEN_HABITABLE and atmosphere_data_ready:
        return PROVEN_HABITABLE
    if label == PROVEN_UNINHABITABLE and atmosphere_data_ready:
        return PROVEN_UNINHABITABLE
    if (
        atmosphere_data_ready
        and probabilities.get("likely_habitable", 0.0) >= 0.75
        and label != UNRESOLVED
    ):
        return PROVEN_HABITABLE
    if (
        atmosphere_data_ready
        and probabilities.get("likely_uninhabitable", 0.0) >= 0.75
        and label != UNRESOLVED
    ):
        return PROVEN_UNINHABITABLE
    return UNRESOLVED


def _claims_allowed(proof_status: str, atmosphere_data_ready: bool) -> list[str]:
    if proof_status == PROVEN_HABITABLE:
        return ["Proxima b meets the configured likely-habitable evidence threshold."]
    if proof_status == PROVEN_UNINHABITABLE:
        return ["Proxima b meets the configured likely-uninhabitable evidence threshold."]
    claims = [
        "Proxima b habitability is not proven by the current JWST archive state.",
        "The correct mission verdict is unresolved but constrained.",
        "The successful local archive pull found public Proxima JWST metadata, not retrievable atmosphere spectra.",
    ]
    if not atmosphere_data_ready:
        claims.append(
            "Current JWST products do not support a Proxima H2O, CH4, CO2, CO, phase-curve, or atmosphere-retention claim."
        )
    return claims


def _claims_blocked(proof_status: str) -> list[str]:
    if proof_status != UNRESOLVED:
        return [
            "A life-detection claim remains blocked.",
            "A physical probe launch remains blocked until a separate engineering gate is passed.",
        ]
    return [
        "Do not claim Proxima b is habitable.",
        "Do not claim Proxima b is uninhabitable.",
        "Do not claim JWST detected water, methane, carbon dioxide, carbon monoxide, or a thermal phase curve for Proxima b.",
        "Do not open the rover/probe gate from this evidence state.",
    ]


def _next_actions(atmosphere_data_ready: bool) -> list[str]:
    if atmosphere_data_ready:
        return [
            "Run wavelength-calibrated retrievals on the Proxima-specific products.",
            "Update posterior probabilities only after stellar-contamination and flare controls.",
        ]
    return [
        "Keep rerunning a Proxima-only JWST/MAST archive monitor for newly released MIRI, NIRSpec, or NIRISS spectroscopy products.",
        "Prepare an observing proposal or next-observatory request focused on atmosphere retention and thermal constraints.",
        "Keep the verdict as unresolved until Proxima-specific atmospheric data exists.",
    ]


def _join_or_none(values: list[Any]) -> str:
    return ", ".join(str(value) for value in values) if values else "none"
