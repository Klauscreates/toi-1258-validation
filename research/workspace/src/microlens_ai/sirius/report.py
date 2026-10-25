from __future__ import annotations

from .records import AuditVerdict, OrbitFitResult, ResidualCandidate
from .residuals import candidate_survives_strict, classify_residual_candidate


def evaluate_audit_verdict(
    orbit_fit: OrbitFitResult | None,
    candidates: list[ResidualCandidate],
    threshold_sigma: float = 5.0,
    threshold_delta_bic: float = 10.0,
) -> AuditVerdict:
    """Apply the strict Sirius C and artificial-signal decision rules."""

    surviving = [candidate for candidate in candidates if candidate_survives_strict(candidate)]
    astrometric_pass = bool(
        orbit_fit and orbit_fit.wobble_match_ok(threshold_sigma, threshold_delta_bic)
    )

    reasons = []
    if orbit_fit is None:
        reasons.append("no_astrometric_fit_available")
    elif not astrometric_pass:
        reasons.extend(orbit_fit.reasons or ("coherent_wobble_thresholds_not_met",))

    if not surviving:
        reasons.append("no_persistent_wfc3_sky_residual_survived_artifact_vetoes")

    reasons.append(
        "artificial_outpost_signal_scope_is_astrometry_image_only_no_external_technosignature_data"
    )

    metrics = {}
    if orbit_fit:
        metrics = {
            "epoch_count": orbit_fit.epoch_count,
            "delta_bic": orbit_fit.delta_bic,
            "wobble_significance_sigma": orbit_fit.wobble_significance_sigma,
            "wobble_amplitude_mas": orbit_fit.wobble_amplitude_mas,
            "leave_one_out_stable": orbit_fit.leave_one_out_stable,
        }

    if astrometric_pass:
        status = "sirius_c_confirmed_candidate_requires_independent_review"
        summary = (
            "A coherent fixed-period astrometric residual passed the configured strict thresholds; "
            "this remains a candidate requiring independent review, not an artificial-signal claim."
        )
        sirius_confirmed = True
    elif surviving:
        status = "unexplained_astrometric_image_anomaly"
        summary = (
            "One or more WFC3 residuals survived the image-artifact veto, but the Sirius C wobble "
            "confirmation thresholds were not met."
        )
        sirius_confirmed = False
    else:
        status = "not_confirmed"
        summary = (
            "Sirius C is not confirmed under the strict audit rules, and the artificial outpost "
            "claim has no evidence in the astrometry/image-only scope."
        )
        sirius_confirmed = False

    return AuditVerdict(
        status=status,
        sirius_c_confirmed=sirius_confirmed,
        artificial_signal_evidence=False,
        summary=summary,
        reasons=tuple(dict.fromkeys(reasons)),
        surviving_candidate_ids=tuple(candidate.candidate_id for candidate in surviving),
        confirmation_metrics=metrics,
    )


def audit_verdict_to_markdown(
    verdict: AuditVerdict,
    orbit_fit: OrbitFitResult | None,
    candidates: list[ResidualCandidate],
) -> str:
    """Render a compact reproducible audit report."""

    lines = [
        "# Sirius C / WFC3 Residual Audit",
        "",
        f"## Verdict: {verdict.status}",
        "",
        verdict.summary,
        "",
        "## Astrometric Wobble Test",
    ]
    if orbit_fit is None:
        lines.append("- No astrometric fit result was available.")
    else:
        lines.extend(
            [
                f"- Epochs: {orbit_fit.epoch_count}",
                f"- Candidate residual period: {orbit_fit.wobble_period_years:.3f} yr",
                f"- Wobble amplitude: {orbit_fit.wobble_amplitude_mas:.3f} mas",
                f"- Significance proxy: {orbit_fit.wobble_significance_sigma:.3f} sigma",
                f"- Delta BIC: {orbit_fit.delta_bic:.3f}",
                f"- Leave-one-epoch-out stable: {orbit_fit.leave_one_out_stable}",
            ]
        )
    lines.extend(["", "## WFC3 Residual Candidates"])
    if not candidates:
        lines.append("- No residual candidates were supplied or generated.")
    else:
        for candidate in candidates:
            lines.append(
                "- "
                f"{candidate.candidate_id}: {classify_residual_candidate(candidate)}, "
                f"peak={candidate.peak_sigma:.2f} sigma, "
                f"sky_repeats={candidate.sky_repeat_count}, "
                f"detector_repeats={candidate.detector_repeat_count}"
            )
    lines.extend(
        [
            "",
            "## Decision Notes",
            *[f"- {reason}" for reason in verdict.reasons],
            "",
            "Artificial outpost signal: no evidence in this astrometry/image-only audit scope.",
            "",
        ]
    )
    return "\n".join(lines)
