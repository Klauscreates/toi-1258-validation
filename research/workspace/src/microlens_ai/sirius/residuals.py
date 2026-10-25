from __future__ import annotations

from .records import ResidualCandidate


def score_residual_candidate(candidate: ResidualCandidate) -> float:
    """Score a WFC3 residual candidate for human audit, not discovery."""

    score = candidate.peak_sigma
    score += 0.8 * candidate.sky_repeat_count
    score -= 1.2 * candidate.detector_repeat_count
    score -= 5.0 * candidate.detector_correlation_abs
    score -= 3.0 * len(candidate.artifact_flags)
    score -= 4.0 * len(candidate.calibration_flags)
    if candidate.separation_from_sirius_a_arcsec < 0.7:
        score -= 5.0
    return score


def classify_residual_candidate(candidate: ResidualCandidate) -> str:
    """Classify a candidate with strict calibration-first labels."""

    if candidate.is_detector_fixed:
        return "rejected_detector_fixed"
    if candidate.has_calibration_veto:
        return "rejected_calibration_or_artifact_veto"
    if candidate.is_persistent_sky_residual:
        return "persistent_sky_residual"
    if candidate.epoch_count < 2:
        return "single_epoch_residual"
    return "weak_or_incoherent_residual"


def candidate_survives_strict(candidate: ResidualCandidate) -> bool:
    """Return True only for residuals allowed to survive to the final audit."""

    return (
        classify_residual_candidate(candidate) == "persistent_sky_residual"
        and score_residual_candidate(candidate) >= 5.0
    )


def wfc3_artifact_flags(
    x: float,
    y: float,
    width: int,
    height: int,
    star_x: float,
    star_y: float,
    edge_margin: float = 10.0,
    saturated_radius: float = 16.0,
    diffraction_spike_width: float = 3.0,
) -> tuple[str, ...]:
    """Return simple detector-position flags for bright-star WFC3 residuals."""

    flags = []
    if x < edge_margin or y < edge_margin or x > width - edge_margin or y > height - edge_margin:
        flags.append("detector_edge")
    dx = x - star_x
    dy = y - star_y
    radius = (dx * dx + dy * dy) ** 0.5
    if radius <= saturated_radius:
        flags.append("saturated_core")
    if abs(dx) <= diffraction_spike_width or abs(dy) <= diffraction_spike_width:
        flags.append("diffraction_spike_zone")
    return tuple(flags)
