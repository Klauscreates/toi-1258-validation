from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Any

from .records import GaiaSourceRecord

SIRIUS_RA_DEGREES = 101.28715533
SIRIUS_DEC_DEGREES = -16.71611586
GAIA_TAP_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"

GAIA_FIELD_COLUMNS = (
    "source_id",
    "designation",
    "ra",
    "dec",
    "parallax",
    "pmra",
    "pmdec",
    "ruwe",
    "astrometric_excess_noise",
    "astrometric_excess_noise_sig",
    "duplicated_source",
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
)

RAW_ASTROMETRY_SCHEMA_WHERE = (
    "table_name LIKE 'gaiadr3.%' AND ("
    "LOWER(table_name) LIKE '%epoch%' OR "
    "LOWER(table_name) LIKE '%astrom%' OR "
    "LOWER(table_name) LIKE '%transit%' OR "
    "LOWER(table_name) LIKE '%observation%' OR "
    "LOWER(table_name) LIKE '%nss%' OR "
    "LOWER(column_name) LIKE '%epoch%' OR "
    "LOWER(column_name) LIKE '%astrom%' OR "
    "LOWER(column_name) LIKE '%transit%' OR "
    "LOWER(column_name) LIKE '%scan%' OR "
    "LOWER(column_name) LIKE '%residual%' OR "
    "LOWER(column_name) LIKE '%obs%'"
    ")"
)


def build_gaia_dr3_field_adql(
    ra_degrees: float = SIRIUS_RA_DEGREES,
    dec_degrees: float = SIRIUS_DEC_DEGREES,
    radius_degrees: float = 0.0833333333,
    limit: int = 200,
) -> str:
    """Build the Gaia DR3 field query used by the Sirius audit.

    The query deliberately gathers quality indicators and photometry rather
    than treating Gaia as a direct proof of Sirius C. Sirius A is brighter than
    Gaia's ordinary bright-limit regime, so these rows are catalog/frame
    constraints and nearby-source context.
    """

    if radius_degrees <= 0:
        raise ValueError("radius_degrees must be positive")
    if limit <= 0:
        raise ValueError("limit must be positive")

    columns = ",\n       ".join(GAIA_FIELD_COLUMNS)
    point = f"POINT('ICRS', ra, dec)"
    center = f"POINT('ICRS', {ra_degrees:.8f}, {dec_degrees:.8f})"
    circle = f"CIRCLE('ICRS', {ra_degrees:.8f}, {dec_degrees:.8f}, {radius_degrees:.10f})"
    return (
        f"SELECT TOP {int(limit)}\n"
        f"       {columns},\n"
        f"       DISTANCE({point}, {center}) AS angular_distance_degrees\n"
        "FROM gaiadr3.gaia_source\n"
        f"WHERE 1 = CONTAINS({point}, {circle})\n"
        "ORDER BY angular_distance_degrees ASC"
    )


def build_gaia_dr3_nss_adql(source_ids: list[str] | tuple[str, ...]) -> str:
    """Build a Gaia DR3 non-single-star lookup query for known source IDs."""

    if not source_ids:
        raise ValueError("source_ids must not be empty")
    ids = ", ".join(str(int(source_id)) for source_id in source_ids)
    return (
        "SELECT source_id, nss_solution_type, period, period_error, eccentricity, "
        "eccentricity_error, a_thiele_innes, b_thiele_innes, f_thiele_innes, "
        "g_thiele_innes, inclination, inclination_error, astrometric_jitter, "
        "significance, flags\n"
        "FROM gaiadr3.nss_two_body_orbit\n"
        f"WHERE source_id IN ({ids})"
    )


def build_raw_astrometry_schema_adql() -> str:
    """Build a Gaia DR3 schema query for raw/epoch/astrometric surfaces."""

    return (
        "SELECT table_name, column_name, datatype, description\n"
        "FROM TAP_SCHEMA.columns\n"
        f"WHERE {RAW_ASTROMETRY_SCHEMA_WHERE}\n"
        "ORDER BY table_name, column_name"
    )


def build_source_keyed_candidate_tables_adql() -> str:
    """Find Gaia DR3 candidate time-series/NSS tables keyed by source_id."""

    return (
        "SELECT DISTINCT table_name\n"
        "FROM TAP_SCHEMA.columns\n"
        "WHERE table_name LIKE 'gaiadr3.%' AND column_name = 'source_id' AND ("
        "LOWER(table_name) LIKE '%nss%' OR "
        "LOWER(table_name) LIKE '%epoch%' OR "
        "LOWER(table_name) LIKE '%astrom%' OR "
        "LOWER(table_name) LIKE '%transit%' OR "
        "LOWER(table_name) LIKE '%observation%' OR "
        "LOWER(table_name) LIKE '%rv%'"
        ")\n"
        "ORDER BY table_name"
    )


def build_source_hit_count_adql(table_name: str, source_ids: list[str] | tuple[str, ...]) -> str:
    """Build a count query for one trusted Gaia DR3 source-keyed table."""

    if not table_name.startswith("gaiadr3.") or not table_name.replace(".", "_").replace("_", "").isalnum():
        raise ValueError(f"Unsafe Gaia table name: {table_name}")
    if not source_ids:
        raise ValueError("source_ids must not be empty")
    ids = ", ".join(str(int(source_id)) for source_id in source_ids)
    return f"SELECT COUNT(*) AS n FROM {table_name} WHERE source_id IN ({ids})"


def query_gaia_tap(
    adql: str,
    endpoint: str = GAIA_TAP_SYNC_URL,
    timeout_seconds: float = 60.0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Execute a synchronous Gaia TAP query and return row dictionaries.

    This intentionally uses the standard library so the Mission 3 audit can run
    without requiring astroquery. The returned metadata is persisted in the
    discovery package for provenance.
    """

    payload = urllib.parse.urlencode(
        {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "json",
            "QUERY": adql,
        }
    ).encode("utf-8")
    request = urllib.request.Request(endpoint, data=payload)
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            raw = response.read().decode("utf-8")
    except Exception as exc:  # pragma: no cover - exercised only in live archive failures
        raise RuntimeError(f"Gaia TAP query failed: {exc}") from exc

    document = json.loads(raw)
    return tap_json_to_rows(document), {
        "endpoint": endpoint,
        "metadata": document.get("metadata", []),
        "row_count": len(document.get("data", [])),
    }


def tap_json_to_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    """Convert a Gaia TAP JSON response into dictionaries keyed by column name."""

    columns = [entry["name"] for entry in document.get("metadata", [])]
    rows = []
    for raw_row in document.get("data", []):
        rows.append({column: raw_row[index] for index, column in enumerate(columns)})
    return rows


def gaia_source_from_row(row: dict[str, Any]) -> GaiaSourceRecord:
    """Normalize a Gaia archive row into the public audit dataclass."""

    record = GaiaSourceRecord(
        source_id=str(row["source_id"]),
        designation=_optional_str(row.get("designation")),
        ra_degrees=float(row["ra"]),
        dec_degrees=float(row["dec"]),
        parallax_mas=_optional_float(row.get("parallax")),
        pmra_mas_per_year=_optional_float(row.get("pmra")),
        pmdec_mas_per_year=_optional_float(row.get("pmdec")),
        ruwe=_optional_float(row.get("ruwe")),
        astrometric_excess_noise_mas=_optional_float(row.get("astrometric_excess_noise")),
        astrometric_excess_noise_sig=_optional_float(row.get("astrometric_excess_noise_sig")),
        duplicated_source=_optional_bool(row.get("duplicated_source")),
        phot_g_mean_mag=_optional_float(row.get("phot_g_mean_mag")),
        phot_bp_mean_mag=_optional_float(row.get("phot_bp_mean_mag")),
        phot_rp_mean_mag=_optional_float(row.get("phot_rp_mean_mag")),
        non_single_star_solution=_optional_str(row.get("non_single_star_solution")),
    )
    warnings = tuple(gaia_interpretation_notes([record]))
    return GaiaSourceRecord(**{**record.__dict__, "warnings": warnings})


def gaia_interpretation_notes(records: list[GaiaSourceRecord]) -> list[str]:
    """Return strict-audit notes that prevent overclaiming Gaia rows."""

    notes = [
        "Gaia DR3 is used as catalog/frame context, not as standalone Sirius C proof.",
    ]
    if any(record.is_too_bright_for_standard_gaia for record in records):
        notes.append(
            "At least one Gaia row is brighter than G=3; interpret ordinary Gaia astrometry with the bright-star caveat."
        )
    if any((record.astrometric_excess_noise_sig or 0.0) > 2.0 for record in records):
        notes.append(
            "Astrometric excess noise is a diagnostic flag and can absorb calibration or modeling errors."
        )
    return notes


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _optional_str(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value)
    return text if text else None


def _optional_bool(value: Any) -> bool | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"true", "t", "1", "yes", "y"}
    return bool(value)
