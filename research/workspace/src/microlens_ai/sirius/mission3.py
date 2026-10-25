from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .gaia import (
    SIRIUS_DEC_DEGREES,
    SIRIUS_RA_DEGREES,
    build_gaia_dr3_field_adql,
    build_gaia_dr3_nss_adql,
    build_raw_astrometry_schema_adql,
    build_source_hit_count_adql,
    build_source_keyed_candidate_tables_adql,
    gaia_interpretation_notes,
    gaia_source_from_row,
    query_gaia_tap,
)
from .orbit import astrometric_epoch_from_dict, fit_orbit_residuals
from .records import AstrometricEpoch, GaiaSourceRecord, OrbitFitResult, ResidualCandidate
from .report import audit_verdict_to_markdown, evaluate_audit_verdict

SIRIUS_C_PERIOD_YEARS = 6.3
DAGON_HYPOTHESIS_PERIOD_YEARS = 400.0


def run_mission3_package(
    output_dir: str | Path,
    live_gaia: bool = True,
    astrometric_epochs_path: str | Path | None = None,
    gaia_radius_degrees: float = 0.02,
    gaia_limit: int = 100,
) -> dict[str, Path]:
    """Execute the Mission 3 Sirius C strict audit package."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    field_query = build_gaia_dr3_field_adql(
        SIRIUS_RA_DEGREES,
        SIRIUS_DEC_DEGREES,
        gaia_radius_degrees,
        gaia_limit,
    )
    (output_path / "gaia_dr3_field_query.sql").write_text(field_query, encoding="utf-8")

    if live_gaia:
        gaia_rows, gaia_metadata = query_gaia_tap(field_query)
    else:
        gaia_rows, gaia_metadata = [], {"endpoint": None, "metadata": [], "row_count": 0}
    _write_json(output_path / "gaia_dr3_field_rows.json", gaia_rows)
    _write_csv(output_path / "gaia_dr3_field_rows.csv", gaia_rows)
    _write_json(output_path / "gaia_dr3_field_metadata.json", gaia_metadata)

    source_ids = [str(row["source_id"]) for row in gaia_rows if row.get("source_id") is not None]
    raw_astrometry_hunt = run_raw_astrometry_hunt(output_path, source_ids, live_gaia)
    if source_ids and live_gaia:
        nss_query = build_gaia_dr3_nss_adql(source_ids)
        nss_rows, nss_metadata = query_gaia_tap(nss_query)
    else:
        nss_query = "-- No Gaia DR3 source IDs were available for NSS lookup."
        nss_rows, nss_metadata = [], {"endpoint": None, "metadata": [], "row_count": 0}
    (output_path / "gaia_dr3_nss_query.sql").write_text(nss_query, encoding="utf-8")
    _write_json(output_path / "gaia_dr3_nss_rows.json", nss_rows)
    _write_json(output_path / "gaia_dr3_nss_metadata.json", nss_metadata)

    gaia_records = [gaia_source_from_row(row) for row in gaia_rows]
    gaia_summary = summarize_gaia_records(gaia_records, gaia_rows, nss_rows)
    _write_json(output_path / "gaia_astrometric_residual_proxy_summary.json", gaia_summary)

    epochs = _load_astrometric_epochs(astrometric_epochs_path)
    if not epochs:
        _write_json(
            output_path / "astrometric_epochs_required_template.json",
            [
                {
                    "epoch_year": 2018.0,
                    "observed_delta_ra_arcsec": 0.0,
                    "observed_delta_dec_arcsec": 0.0,
                    "uncertainty_mas": 1.0,
                    "instrument": "WFC3/UVIS",
                    "filter_name": "replace_with_real_measurement",
                }
            ],
        )

    fit_6p3 = fit_orbit_residuals(epochs, target_period_years=SIRIUS_C_PERIOD_YEARS)
    fit_400 = fit_orbit_residuals(epochs, target_period_years=DAGON_HYPOTHESIS_PERIOD_YEARS)
    _write_json(output_path / "orbit_fit_6p3_year_sirius_c.json", asdict(fit_6p3))
    _write_json(output_path / "orbit_fit_400_year_dagon_hypothesis.json", asdict(fit_400))

    candidates: list[ResidualCandidate] = []
    verdict = evaluate_audit_verdict(fit_6p3, candidates)
    dagon_audit = evaluate_period_hypothesis(fit_400, "400_year_dagon_hypothesis")
    _write_json(output_path / "candidate_ledger.json", [])
    _write_json(output_path / "audit_verdict.json", asdict(verdict))
    _write_json(output_path / "dagon_hypothesis_audit.json", dagon_audit)

    package = build_mission3_payload(
        gaia_summary=gaia_summary,
        gaia_metadata=gaia_metadata,
        nss_metadata=nss_metadata,
        raw_astrometry_hunt=raw_astrometry_hunt,
        orbit_6p3=fit_6p3,
        orbit_400=fit_400,
        dagon_audit=dagon_audit,
        verdict=asdict(verdict),
        epoch_count=len(epochs),
    )
    _write_json(output_path / "mission3_discovery_package.json", package)

    audit_markdown = render_mission3_audit(package, fit_6p3, fit_400)
    (output_path / "scientific_audit.md").write_text(audit_markdown, encoding="utf-8")

    signature = write_signature_files(output_path)
    _write_json(output_path / "signature.json", signature)

    return {
        "package": output_path / "mission3_discovery_package.json",
        "audit": output_path / "scientific_audit.md",
        "signature": output_path / "signature.json",
        "gaia_rows": output_path / "gaia_dr3_field_rows.json",
        "sirius_c_fit": output_path / "orbit_fit_6p3_year_sirius_c.json",
        "dagon_fit": output_path / "orbit_fit_400_year_dagon_hypothesis.json",
        "raw_astrometry_hunt": output_path / "raw_astrometry_hunt.json",
    }


def run_raw_astrometry_hunt(
    output_path: Path,
    source_ids: list[str],
    live_gaia: bool,
) -> dict[str, Any]:
    """Interrogate Gaia DR3 for raw stellar epoch astrometry surfaces."""

    schema_query = build_raw_astrometry_schema_adql()
    table_query = build_source_keyed_candidate_tables_adql()
    (output_path / "raw_astrometry_schema_hunt.sql").write_text(schema_query, encoding="utf-8")
    (output_path / "raw_astrometry_source_keyed_tables.sql").write_text(table_query, encoding="utf-8")

    if live_gaia:
        schema_rows, schema_metadata = query_gaia_tap(schema_query)
        table_rows, table_metadata = query_gaia_tap(table_query)
    else:
        schema_rows, schema_metadata = [], {"row_count": 0, "endpoint": None}
        table_rows, table_metadata = [], {"row_count": 0, "endpoint": None}

    source_keyed_tables = [row["table_name"] for row in table_rows if row.get("table_name")]
    hit_counts = []
    for table_name in source_keyed_tables:
        if not source_ids:
            hit_counts.append(
                {
                    "table_name": table_name,
                    "source_id_hit_count": 0,
                    "query_status": "skipped_no_source_ids",
                    "surface_class": classify_gaia_astrometry_surface(table_name),
                }
            )
            continue
        try:
            count_query = build_source_hit_count_adql(table_name, source_ids)
            count_rows, _ = query_gaia_tap(count_query) if live_gaia else ([{"n": 0}], {})
            count = int(count_rows[0]["n"]) if count_rows else 0
            status = "ok"
        except Exception as exc:  # pragma: no cover - live archive defensive path
            count = 0
            status = f"error: {exc}"
        hit_counts.append(
            {
                "table_name": table_name,
                "source_id_hit_count": count,
                "query_status": status,
                "surface_class": classify_gaia_astrometry_surface(table_name),
            }
        )

    schema_tables = sorted({row["table_name"] for row in schema_rows if row.get("table_name")})
    raw_candidate_tables = [
        table
        for table in schema_tables
        if classify_gaia_astrometry_surface(table) == "possible_raw_stellar_epoch_astrometry"
    ]
    public_raw_available = bool(raw_candidate_tables)
    hunt = {
        "public_raw_stellar_epoch_astrometry_available": public_raw_available,
        "schema_candidate_column_count": len(schema_rows),
        "schema_candidate_table_count": len(schema_tables),
        "source_keyed_candidate_tables": hit_counts,
        "possible_raw_epoch_astrometry_tables": raw_candidate_tables,
        "conclusion": (
            "A Gaia DR3 TAP_SCHEMA attack found no public source-keyed raw stellar epoch astrometry "
            "or along-scan residual table for Sirius-field sources."
            if not public_raw_available
            else "Potential raw stellar epoch astrometry surfaces were found and require manual validation."
        ),
        "data_surfaces_found": {
            "scan_law": [table for table in schema_tables if table.endswith("commanded_scan_law")],
            "catalog_summary_astrometry": [table for table in schema_tables if table.endswith("gaia_source") or table.endswith("gaia_source_lite")],
            "nss_fitted_solutions": [table for table in schema_tables if ".nss_" in table],
            "solar_system_observations": [table for table in schema_tables if table.endswith("sso_observation")],
            "photometric_or_spectroscopic_epoch_products": [
                table
                for table in schema_tables
                if "epoch" in table.lower() or "radial_velocity" in table.lower()
            ],
        },
        "why_this_does_not_confirm_wobble": [
            "Scan-law tables provide spacecraft pointing geometry, not per-source measured residuals.",
            "gaia_source astrometric fields are compressed summary diagnostics at reference epoch J2016.0.",
            "NSS tables are fitted model products, not raw along-scan epoch measurements.",
            "epoch_photometry and epoch radial velocity are not astrometric sky-position residuals.",
        ],
        "schema_metadata": schema_metadata,
        "table_metadata": table_metadata,
    }
    _write_json(output_path / "raw_astrometry_hunt.json", hunt)
    _write_json(output_path / "raw_astrometry_schema_candidates.json", schema_rows)
    _write_csv(output_path / "raw_astrometry_schema_candidates.csv", schema_rows)
    _write_json(output_path / "raw_astrometry_source_table_hits.json", hit_counts)
    return hunt


def classify_gaia_astrometry_surface(table_name: str) -> str:
    """Classify Gaia DR3 archive surfaces found during the raw-data hunt."""

    lower = table_name.lower()
    if lower.endswith("commanded_scan_law"):
        return "spacecraft_scan_law_not_source_astrometry"
    if lower.endswith("gaia_source") or lower.endswith("gaia_source_lite"):
        return "catalog_summary_astrometry_not_epoch_residuals"
    if ".nss_" in lower:
        return "nss_fitted_solution_not_raw_epoch_astrometry"
    if lower.endswith("sso_observation"):
        return "solar_system_observation_not_stellar_astrometry"
    if "epoch_photometry" in lower:
        return "epoch_photometry_not_astrometry"
    if "radial_velocity" in lower:
        return "epoch_radial_velocity_not_astrometry"
    if "epoch" in lower and "astrom" in lower and "source" in lower:
        return "possible_raw_stellar_epoch_astrometry"
    return "not_raw_stellar_epoch_astrometry"


def summarize_gaia_records(
    records: list[GaiaSourceRecord],
    rows: list[dict[str, Any]],
    nss_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize Gaia DR3 rows without treating quality flags as discoveries."""

    source_rows_by_id = {str(row.get("source_id")): row for row in rows}
    probable_b = probable_sirius_b_record(records, source_rows_by_id)
    bright_rows = [record for record in records if record.is_too_bright_for_standard_gaia]
    quality_anomalies = [
        {
            "source_id": record.source_id,
            "designation": record.designation,
            "astrometric_excess_noise_mas": record.astrometric_excess_noise_mas,
            "astrometric_excess_noise_sig": record.astrometric_excess_noise_sig,
            "ruwe": record.ruwe,
            "interpretation": "catalog_quality_anomaly_not_gravitational_confirmation",
        }
        for record in records
        if (record.astrometric_excess_noise_sig or 0.0) > 5.0 or (record.ruwe or 0.0) > 1.4
    ]
    return {
        "row_count": len(records),
        "sirius_a_standard_gaia_row_found": bool(bright_rows),
        "probable_sirius_b_gaia_dr3_source": asdict(probable_b) if probable_b else None,
        "gaia_nss_solution_count": len(nss_rows),
        "gaia_quality_anomaly_count": len(quality_anomalies),
        "gaia_quality_anomalies": quality_anomalies,
        "interpretation_notes": gaia_interpretation_notes(records)
        + [
            "Gaia DR3 TAP gaia_source rows do not provide the raw along-scan residual time series needed for a standalone wobble fit.",
            "Sirius A is expected to be problematic or absent in ordinary Gaia DR3 source treatment because it is far brighter than G~=3.",
        ],
    }


def probable_sirius_b_record(
    records: list[GaiaSourceRecord],
    source_rows_by_id: dict[str, dict[str, Any]],
) -> GaiaSourceRecord | None:
    """Pick the most plausible Sirius B Gaia source by parallax/proper motion."""

    candidates = []
    for record in records:
        parallax = record.parallax_mas or 0.0
        pmra = record.pmra_mas_per_year or 0.0
        pmdec = record.pmdec_mas_per_year or 0.0
        total_pm = (pmra * pmra + pmdec * pmdec) ** 0.5
        if 250.0 <= parallax <= 450.0 and total_pm >= 500.0:
            distance = float(source_rows_by_id.get(record.source_id, {}).get("angular_distance_degrees") or 99.0)
            candidates.append((distance, record))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def evaluate_period_hypothesis(
    fit: OrbitFitResult,
    hypothesis_name: str,
) -> dict[str, Any]:
    """Evaluate a named fixed-period gravitational anomaly hypothesis."""

    match = fit.wobble_match_ok()
    reasons = list(fit.reasons)
    if not match:
        reasons.append("hypothesis_not_supported_by_available_astrometric_epochs")
    return {
        "hypothesis": hypothesis_name,
        "period_years": fit.wobble_period_years,
        "matched": match,
        "status": "matched_requires_independent_review" if match else "not_matched",
        "reasons": tuple(dict.fromkeys(reasons)),
        "metrics": {
            "epoch_count": fit.epoch_count,
            "delta_bic": fit.delta_bic,
            "wobble_significance_sigma": fit.wobble_significance_sigma,
            "wobble_amplitude_mas": fit.wobble_amplitude_mas,
            "leave_one_out_stable": fit.leave_one_out_stable,
        },
    }


def build_mission3_payload(
    gaia_summary: dict[str, Any],
    gaia_metadata: dict[str, Any],
    nss_metadata: dict[str, Any],
    raw_astrometry_hunt: dict[str, Any],
    orbit_6p3: OrbitFitResult,
    orbit_400: OrbitFitResult,
    dagon_audit: dict[str, Any],
    verdict: dict[str, Any],
    epoch_count: int,
) -> dict[str, Any]:
    generated_at = datetime.now(UTC).isoformat()
    return {
        "mission": {
            "id": "mission3_sirius_c_strict_audit",
            "generated_at_utc": generated_at,
            "target": "Sirius AB / Sirius C claim",
            "scope": "Gaia DR3 catalog-frame check plus fixed-period residual audit",
            "claim_standard": "strict_null_test",
        },
        "archive_queries": {
            "gaia_tap_endpoint": gaia_metadata.get("endpoint"),
            "gaia_field_row_count": gaia_metadata.get("row_count", 0),
            "gaia_nss_row_count": nss_metadata.get("row_count", 0),
        },
        "gaia_astrometry": gaia_summary,
        "raw_astrometry_hunt": {
            "public_raw_stellar_epoch_astrometry_available": raw_astrometry_hunt[
                "public_raw_stellar_epoch_astrometry_available"
            ],
            "schema_candidate_column_count": raw_astrometry_hunt["schema_candidate_column_count"],
            "schema_candidate_table_count": raw_astrometry_hunt["schema_candidate_table_count"],
            "conclusion": raw_astrometry_hunt["conclusion"],
            "source_keyed_candidate_tables": raw_astrometry_hunt["source_keyed_candidate_tables"],
        },
        "orbital_data": {
            "epoch_count": epoch_count,
            "residual_time_series_available": epoch_count > 0,
            "limitation": "No raw Gaia DR3 A/B astrometric residual time series is available from gaia_source; real wobble confirmation requires independent epoch astrometry.",
        },
        "wobble_tests": {
            "sirius_c_6p3_year": asdict(orbit_6p3),
            "dagon_400_year": asdict(orbit_400),
            "dagon_hypothesis_audit": dagon_audit,
        },
        "verdict": verdict,
        "groundbreaking_claim_control": {
            "claim": "No groundbreaking gravitational or artificial-signal claim is made from Gaia catalog diagnostics alone.",
            "reason": "The public Gaia DR3 archive surfaces found by the schema hunt are insufficient to confirm a coherent wobble or 400-year anomaly under the configured thresholds.",
        },
    }


def render_mission3_audit(
    package: dict[str, Any],
    fit_6p3: OrbitFitResult,
    fit_400: OrbitFitResult,
) -> str:
    """Render the signed scientific audit text."""

    base_report = audit_verdict_to_markdown(
        evaluate_audit_verdict(fit_6p3, []),
        fit_6p3,
        [],
    )
    gaia = package["gaia_astrometry"]
    raw = package["raw_astrometry_hunt"]
    probable_b = gaia.get("probable_sirius_b_gaia_dr3_source") or {}
    lines = [
        "# Mission 3: Sirius C Scientific Audit",
        "",
        f"Generated UTC: {package['mission']['generated_at_utc']}",
        "",
        "## Gaia DR3 Archive Result",
        "",
        f"- Gaia field rows returned: {package['archive_queries']['gaia_field_row_count']}",
        f"- Gaia NSS rows returned: {package['archive_queries']['gaia_nss_row_count']}",
        f"- Sirius A standard Gaia row found: {gaia['sirius_a_standard_gaia_row_found']}",
        f"- Probable Sirius B Gaia DR3 source: {probable_b.get('source_id', 'not identified')}",
        f"- Gaia quality anomaly rows: {gaia['gaia_quality_anomaly_count']}",
        "",
        "Gaia quality anomalies are treated as catalog diagnostics, not gravitational-anomaly confirmations.",
        "",
        "## Raw Astrometric Time-Series Hunt",
        "",
        f"- Public raw stellar epoch astrometry found: {raw['public_raw_stellar_epoch_astrometry_available']}",
        f"- Schema candidate columns inspected: {raw['schema_candidate_column_count']}",
        f"- Schema candidate tables inspected: {raw['schema_candidate_table_count']}",
        f"- Conclusion: {raw['conclusion']}",
        "",
        "Source-keyed candidate tables checked:",
        "",
        "| Table | Sirius-field hits | Surface class |",
        "| --- | ---: | --- |",
        *[
            f"| {row['table_name']} | {row['source_id_hit_count']} | {row['surface_class']} |"
            for row in raw["source_keyed_candidate_tables"]
        ],
        "",
        "## Fixed-Period Wobble Tests",
        "",
        f"- 6.3-year Sirius C fit: sigma={fit_6p3.wobble_significance_sigma:.3f}, delta_bic={fit_6p3.delta_bic:.3f}, matched={fit_6p3.wobble_match_ok()}",
        f"- 400-year Dagon hypothesis fit: sigma={fit_400.wobble_significance_sigma:.3f}, delta_bic={fit_400.delta_bic:.3f}, matched={fit_400.wobble_match_ok()}",
        "",
        "## Core Verdict",
        "",
        base_report,
        "",
        "## Signed Audit Statement",
        "",
        "This package is signed with a SHA-256 checksum manifest for tamper evidence. It is not an institutional endorsement or a discovery announcement.",
        "",
    ]
    return "\n".join(lines)


def write_signature_files(output_dir: Path) -> dict[str, Any]:
    """Write checksums and return a package-level SHA-256 signature payload."""

    files = sorted(
        path
        for path in output_dir.iterdir()
        if path.is_file() and path.name not in {"signature.json", "checksums.sha256"}
    )
    checksum_lines = []
    digest = hashlib.sha256()
    for path in files:
        file_digest = hashlib.sha256(path.read_bytes()).hexdigest()
        checksum_lines.append(f"{file_digest}  {path.name}")
        digest.update(path.name.encode("utf-8"))
        digest.update(file_digest.encode("ascii"))
    (output_dir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    return {
        "signature_type": "sha256_checksum_manifest",
        "signed_at_utc": datetime.now(UTC).isoformat(),
        "signed_by": "Codex Mission 3 scientific audit runner",
        "package_sha256": digest.hexdigest(),
        "checksum_file": "checksums.sha256",
        "discovery_claim": "not_confirmed",
        "caveat": "Tamper-evident package digest only; not a legal identity signature.",
    }


def _load_astrometric_epochs(path: str | Path | None) -> list[AstrometricEpoch]:
    if path is None:
        return []
    epoch_path = Path(path)
    if not epoch_path.exists():
        return []
    if epoch_path.suffix.lower() == ".csv":
        with epoch_path.open("r", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    else:
        rows = json.loads(epoch_path.read_text(encoding="utf-8"))
    return [astrometric_epoch_from_dict(row) for row in rows]


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    build_raw_astrometry_schema_adql,
    build_source_hit_count_adql,
    build_source_keyed_candidate_tables_adql,
