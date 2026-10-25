from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class GaiaSourceRecord:
    """Gaia DR3 source row used for Sirius field/frame context."""

    source_id: str
    ra_degrees: float
    dec_degrees: float
    parallax_mas: float | None = None
    pmra_mas_per_year: float | None = None
    pmdec_mas_per_year: float | None = None
    ruwe: float | None = None
    astrometric_excess_noise_mas: float | None = None
    astrometric_excess_noise_sig: float | None = None
    duplicated_source: bool | None = None
    phot_g_mean_mag: float | None = None
    phot_bp_mean_mag: float | None = None
    phot_rp_mean_mag: float | None = None
    designation: str | None = None
    non_single_star_solution: str | None = None
    warnings: tuple[str, ...] = ()

    @property
    def is_too_bright_for_standard_gaia(self) -> bool:
        return self.phot_g_mean_mag is not None and self.phot_g_mean_mag < 3.0


@dataclass(frozen=True)
class Wfc3Observation:
    """Normalized MAST/HST observation row for the Sirius WFC3 audit."""

    obs_id: str
    proposal_id: str
    instrument: str
    filter_name: str
    target_name: str
    t_min_mjd: float
    ra_degrees: float
    dec_degrees: float
    calibration_level: int
    data_product_type: str = "image"
    intent_type: str = "science"
    product_filenames: tuple[str, ...] = ()
    crds_context: str | None = None
    artifact_risks: tuple[str, ...] = ()
    public: bool = True


@dataclass(frozen=True)
class AstrometricEpoch:
    """Measured Sirius B relative position at one epoch."""

    epoch_year: float
    observed_delta_ra_arcsec: float
    observed_delta_dec_arcsec: float
    uncertainty_mas: float
    instrument: str
    filter_name: str = ""
    roll_angle_deg: float | None = None
    detector_x: float | None = None
    detector_y: float | None = None
    flags: tuple[str, ...] = ()


@dataclass(frozen=True)
class OrbitFitResult:
    """Baseline Sirius AB residual fit plus fixed-period wobble search."""

    model_name: str
    epoch_count: int
    baseline_chi2: float
    wobble_chi2: float
    delta_bic: float
    wobble_period_years: float
    wobble_amplitude_mas: float
    wobble_phase_degrees: float
    wobble_significance_sigma: float
    leave_one_out_stable: bool
    residual_rms_mas: float
    reasons: tuple[str, ...] = ()
    residuals: tuple[dict[str, float], ...] = ()

    def wobble_match_ok(
        self,
        threshold_sigma: float = 5.0,
        threshold_delta_bic: float = 10.0,
    ) -> bool:
        return (
            self.epoch_count >= 6
            and self.wobble_significance_sigma >= threshold_sigma
            and self.delta_bic >= threshold_delta_bic
            and self.leave_one_out_stable
        )


@dataclass(frozen=True)
class ResidualCandidate:
    """WFC3 image residual or astrometric-image anomaly candidate."""

    candidate_id: str
    ra_degrees: float
    dec_degrees: float
    peak_sigma: float
    epoch_count: int
    sky_repeat_count: int
    detector_repeat_count: int
    detector_correlation_abs: float
    separation_from_sirius_a_arcsec: float
    separation_from_sirius_b_arcsec: float
    artifact_flags: tuple[str, ...] = ()
    calibration_flags: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def is_detector_fixed(self) -> bool:
        return (
            self.detector_repeat_count > self.sky_repeat_count
            and self.detector_correlation_abs >= 0.7
        )

    @property
    def has_calibration_veto(self) -> bool:
        return bool(self.artifact_flags or self.calibration_flags)

    @property
    def is_persistent_sky_residual(self) -> bool:
        return self.peak_sigma >= 5.0 and self.sky_repeat_count >= 2


@dataclass(frozen=True)
class AuditVerdict:
    """Final strict-audit verdict for Sirius C and artificial-signal claims."""

    status: str
    sirius_c_confirmed: bool
    artificial_signal_evidence: bool
    summary: str
    reasons: tuple[str, ...]
    surviving_candidate_ids: tuple[str, ...] = ()
    confirmation_metrics: dict[str, float | bool | str] = field(default_factory=dict)


def dataclass_to_dict(value: Any) -> dict[str, Any]:
    """Convert one of the audit dataclasses to JSON-friendly dictionaries."""

    return asdict(value)
