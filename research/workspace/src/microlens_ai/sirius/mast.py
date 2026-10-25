from __future__ import annotations

from typing import Any

from .gaia import SIRIUS_DEC_DEGREES, SIRIUS_RA_DEGREES
from .records import Wfc3Observation

WFC3_INSTRUMENTS = ("WFC3/UVIS", "WFC3/IR")


def build_mast_sirius_criteria(
    ra_degrees: float = SIRIUS_RA_DEGREES,
    dec_degrees: float = SIRIUS_DEC_DEGREES,
    radius_degrees: float = 0.05,
) -> dict[str, Any]:
    """Return an astroquery/MAST criteria payload for public Sirius WFC3 data."""

    if radius_degrees <= 0:
        raise ValueError("radius_degrees must be positive")
    return {
        "obs_collection": "HST",
        "coordinates": f"{ra_degrees:.8f} {dec_degrees:.8f}",
        "radius": f"{radius_degrees:.8f} deg",
        "instrument_name": list(WFC3_INSTRUMENTS),
        "dataproduct_type": "image",
        "intentType": "science",
    }


def normalize_wfc3_observation(row: dict[str, Any]) -> Wfc3Observation:
    """Normalize a MAST observation row and annotate Sirius-specific risks."""

    instrument = str(row.get("instrument_name") or row.get("instrument") or "")
    filter_name = str(row.get("filters") or row.get("filter_name") or "")
    target_name = str(row.get("target_name") or row.get("target") or "")
    obs_id = str(row.get("obs_id") or row.get("obsid") or "")
    proposal_id = str(row.get("proposal_id") or "")
    calibration_level = int(row.get("calib_level") or row.get("calibration_level") or 0)
    intent_type = str(row.get("intentType") or row.get("intent_type") or "science")
    data_product_type = str(row.get("dataproduct_type") or row.get("data_product_type") or "image")
    public = _public_from_row(row)
    return Wfc3Observation(
        obs_id=obs_id,
        proposal_id=proposal_id,
        instrument=instrument,
        filter_name=filter_name,
        target_name=target_name,
        t_min_mjd=float(row.get("t_min") or row.get("t_min_mjd") or 0.0),
        ra_degrees=float(row.get("s_ra") or row.get("ra") or SIRIUS_RA_DEGREES),
        dec_degrees=float(row.get("s_dec") or row.get("dec") or SIRIUS_DEC_DEGREES),
        calibration_level=calibration_level,
        data_product_type=data_product_type,
        intent_type=intent_type,
        product_filenames=tuple(row.get("product_filenames") or ()),
        crds_context=_optional_str(row.get("crds_context")),
        artifact_risks=artifact_risks_for_observation(instrument, target_name),
        public=public,
    )


def filter_preferred_wfc3_observations(
    observations: list[Wfc3Observation],
    require_public: bool = True,
) -> list[Wfc3Observation]:
    """Keep science WFC3 image rows and sort best calibrated products first."""

    filtered = [
        observation
        for observation in observations
        if observation.instrument in WFC3_INSTRUMENTS
        and observation.data_product_type.lower() == "image"
        and observation.intent_type.lower() == "science"
        and observation.calibration_level >= 2
        and (observation.public or not require_public)
    ]
    return sorted(
        filtered,
        key=lambda item: (
            item.calibration_level,
            1 if item.instrument == "WFC3/UVIS" else 0,
            -item.t_min_mjd,
        ),
        reverse=True,
    )


def artifact_risks_for_observation(instrument: str, target_name: str) -> tuple[str, ...]:
    """List detector/target risks that must be vetoed before a claim."""

    risks = ["saturation", "diffraction_spikes", "scattered_light", "geometric_distortion"]
    if instrument == "WFC3/UVIS":
        risks.extend(["cte_loss", "under_sampled_psf"])
    if instrument == "WFC3/IR":
        risks.extend(["ir_persistence", "detector_nonlinearity"])
    if "SIRIUS" not in target_name.upper():
        risks.append("target_name_not_sirius")
    return tuple(dict.fromkeys(risks))


def _public_from_row(row: dict[str, Any]) -> bool:
    value = row.get("dataRights") or row.get("data_rights") or row.get("public")
    if value is None:
        return True
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"public", "true", "t", "1", "yes"}


def _optional_str(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value)
    return text if text else None
