from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterable


NACO_S27_PIXEL_SCALE_ARCSEC = 0.02715


@dataclass(frozen=True)
class AoArchiveLead:
    """Ground adaptive-optics archive group relevant to the Sirius audit."""

    instrument: str
    program_id: str
    program_title: str
    date_start: str
    date_end: str
    row_count: int
    science_object_count: int
    filters: tuple[str, ...]
    total_bytes: int = 0


@dataclass(frozen=True)
class AoQuicklookVerdict:
    """Conservative interpretation of a ground-AO quicklook pass."""

    status: str
    anomaly_detected: bool
    wobble_test_possible: bool
    reasons: tuple[str, ...]


def summarize_ao_manifest(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Summarize normalized ESO archive rows or download-manifest rows."""

    materialized = list(rows)
    by_cat_type_filter: Counter[tuple[str, str, str]] = Counter()
    by_filter: Counter[str] = Counter()
    bytes_by_filter: Counter[str] = Counter()
    science_object_counts: Counter[str] = Counter()
    total_bytes = 0

    for row in materialized:
        dp_cat = str(row.get("dp_cat") or row.get("dataproduct_type") or "")
        dp_type = str(row.get("dp_type") or "")
        filter_name = str(row.get("filter_path") or row.get("filter") or "")
        bytes_value = int(float(row.get("bytes") or row.get("access_estsize") or 0))
        by_cat_type_filter[(dp_cat, dp_type, filter_name)] += 1
        by_filter[filter_name] += 1
        bytes_by_filter[filter_name] += bytes_value
        total_bytes += bytes_value
        if dp_cat.upper() == "SCIENCE" and dp_type.upper() == "OBJECT":
            science_object_counts[filter_name] += 1

    return {
        "row_count": len(materialized),
        "total_bytes": total_bytes,
        "by_cat_type_filter": {
            " | ".join(key): value for key, value in sorted(by_cat_type_filter.items())
        },
        "by_filter": {
            filter_name: {"n": by_filter[filter_name], "bytes": bytes_by_filter[filter_name]}
            for filter_name in sorted(by_filter)
        },
        "science_object_counts": dict(sorted(science_object_counts.items())),
    }


def build_eso_sirius_raw_adql(
    instrument: str,
    radius_degrees: float = 0.08,
    limit: int | None = None,
) -> str:
    """Build an ESO TAP query for raw VLT rows near Sirius."""

    top = f"TOP {int(limit)} " if limit is not None else ""
    return f"""
SELECT {top}
  dp_id, instrument, telescope, object, target, date_obs, mjd_obs,
  dp_cat, dp_type, dp_tech, ins_mode, filter_path, exposure,
  det_dit, det_ndit, prog_id, prog_title, ra, dec,
  access_url, datalink_url, release_date
FROM ivoa.ObsCore
WHERE instrument = '{instrument}'
  AND CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 101.28715533, -16.71611586, {radius_degrees})) = 1
  AND dataproduct_type = 'image'
ORDER BY date_obs
""".strip()


def pixel_offset_from_sky_offset(
    east_arcsec: float,
    north_arcsec: float,
    pixel_scale_arcsec: float = NACO_S27_PIXEL_SCALE_ARCSEC,
    cd1_sign: int = -1,
    cd2_sign: int = 1,
) -> tuple[float, float]:
    """Convert east/north sky offset into NACO-like x/y pixel offsets."""

    if pixel_scale_arcsec <= 0:
        raise ValueError("pixel_scale_arcsec must be positive")
    return (
        cd1_sign * east_arcsec / pixel_scale_arcsec,
        cd2_sign * north_arcsec / pixel_scale_arcsec,
    )


def evaluate_ao_quicklook(
    epoch_count: int,
    compact_candidate_count: int,
    calibration_veto_count: int = 0,
) -> AoQuicklookVerdict:
    """Return the strict audit interpretation for ground-AO quicklook output."""

    reasons: list[str] = []
    if epoch_count < 3:
        reasons.append("single_or_sparse_epoch_cannot_test_6p3_year_or_400_year_wobble")
    if compact_candidate_count <= 0:
        reasons.append("no_compact_off_axis_residual_survived_quicklook")
    if calibration_veto_count:
        reasons.append("candidate_rows_have_raw_calibration_or_artifact_vetoes")

    wobble_test_possible = epoch_count >= 3
    anomaly_detected = compact_candidate_count > 0 and calibration_veto_count == 0
    status = "unexplained_ao_anomaly" if anomaly_detected and wobble_test_possible else "not_confirmed"
    return AoQuicklookVerdict(
        status=status,
        anomaly_detected=anomaly_detected,
        wobble_test_possible=wobble_test_possible,
        reasons=tuple(reasons),
    )
