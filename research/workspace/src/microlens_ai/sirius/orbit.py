from __future__ import annotations

import math
from collections.abc import Iterable
from typing import Any

from .records import AstrometricEpoch, OrbitFitResult

DEFAULT_SIRIUS_AB_ORBIT: dict[str, float] = {
    "period_years": 50.1284,
    "semimajor_axis_arcsec": 7.4957,
    "eccentricity": 0.59142,
    "inclination_degrees": 136.336,
    "longitude_node_degrees": 44.57,
    "argument_periastron_degrees": 149.161,
    "periastron_epoch_year": 1994.5715,
}


def solve_kepler(mean_anomaly_rad: float, eccentricity: float, iterations: int = 30) -> float:
    """Solve Kepler's equation for eccentric anomaly in radians."""

    if not 0.0 <= eccentricity < 1.0:
        raise ValueError("eccentricity must be in [0, 1)")
    anomaly = mean_anomaly_rad % (2.0 * math.pi)
    estimate = anomaly if eccentricity < 0.8 else math.pi
    for _ in range(iterations):
        delta = (estimate - eccentricity * math.sin(estimate) - anomaly) / (
            1.0 - eccentricity * math.cos(estimate)
        )
        estimate -= delta
        if abs(delta) < 1e-12:
            break
    return estimate


def projected_binary_offset_arcsec(
    epoch_year: float,
    orbit_elements: dict[str, float] | None = None,
) -> tuple[float, float]:
    """Return projected Sirius B minus Sirius A offset as east/north arcsec."""

    orbit = {**DEFAULT_SIRIUS_AB_ORBIT, **(orbit_elements or {})}
    period = orbit["period_years"]
    eccentricity = orbit["eccentricity"]
    semimajor_axis = orbit["semimajor_axis_arcsec"]
    mean_anomaly = 2.0 * math.pi * ((epoch_year - orbit["periastron_epoch_year"]) / period)
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    orbital_x = semimajor_axis * (math.cos(eccentric_anomaly) - eccentricity)
    orbital_y = semimajor_axis * math.sqrt(1.0 - eccentricity * eccentricity) * math.sin(eccentric_anomaly)

    node = math.radians(orbit["longitude_node_degrees"])
    arg = math.radians(orbit["argument_periastron_degrees"])
    inc = math.radians(orbit["inclination_degrees"])
    cos_node = math.cos(node)
    sin_node = math.sin(node)
    cos_arg = math.cos(arg)
    sin_arg = math.sin(arg)
    cos_inc = math.cos(inc)

    east = orbital_x * (cos_node * cos_arg - sin_node * sin_arg * cos_inc)
    east -= orbital_y * (cos_node * sin_arg + sin_node * cos_arg * cos_inc)
    north = orbital_x * (sin_node * cos_arg + cos_node * sin_arg * cos_inc)
    north += orbital_y * (-sin_node * sin_arg + cos_node * cos_arg * cos_inc)
    return east, north


def compute_astrometric_residuals(
    epochs: Iterable[AstrometricEpoch],
    orbit_elements: dict[str, float] | None = None,
) -> list[dict[str, float]]:
    """Compute observed-minus-calculated residuals against the baseline AB orbit."""

    residuals = []
    for epoch in epochs:
        model_east, model_north = projected_binary_offset_arcsec(epoch.epoch_year, orbit_elements)
        residuals.append(
            {
                "epoch_year": float(epoch.epoch_year),
                "observed_delta_ra_arcsec": float(epoch.observed_delta_ra_arcsec),
                "observed_delta_dec_arcsec": float(epoch.observed_delta_dec_arcsec),
                "model_delta_ra_arcsec": float(model_east),
                "model_delta_dec_arcsec": float(model_north),
                "residual_ra_mas": 1000.0 * (epoch.observed_delta_ra_arcsec - model_east),
                "residual_dec_mas": 1000.0 * (epoch.observed_delta_dec_arcsec - model_north),
                "uncertainty_mas": max(float(epoch.uncertainty_mas), 1e-9),
            }
        )
    return residuals


def fit_orbit_residuals(
    epochs: list[AstrometricEpoch],
    target_period_years: float = 6.3,
    orbit_elements: dict[str, float] | None = None,
) -> OrbitFitResult:
    """Fit fixed-period residual wobble after subtracting the Sirius AB orbit."""

    residuals = compute_astrometric_residuals(epochs, orbit_elements)
    if not residuals:
        return OrbitFitResult(
            model_name="sirius_ab_fixed_orbit_plus_fixed_period_wobble",
            epoch_count=0,
            baseline_chi2=0.0,
            wobble_chi2=0.0,
            delta_bic=0.0,
            wobble_period_years=target_period_years,
            wobble_amplitude_mas=0.0,
            wobble_phase_degrees=0.0,
            wobble_significance_sigma=0.0,
            leave_one_out_stable=False,
            residual_rms_mas=0.0,
            reasons=("no_astrometric_epochs",),
            residuals=(),
        )

    baseline_chi2 = _chi2_zero(residuals)
    wobble = _fit_wobble_to_residuals(residuals, target_period_years)
    n_values = max(1, 2 * len(residuals))
    wobble_parameter_count = 4
    baseline_bic = baseline_chi2
    wobble_bic = wobble["chi2"] + wobble_parameter_count * math.log(n_values)
    delta_bic = baseline_bic - wobble_bic
    improvement = max(0.0, baseline_chi2 - wobble["chi2"])
    significance = math.sqrt(improvement)
    residual_rms = _residual_rms_mas(residuals)
    leave_one_out_stable = _leave_one_out_stability(
        residuals,
        target_period_years,
        wobble["amplitude_mas"],
        wobble["phase_degrees"],
    )

    reasons = []
    if significance < 5.0:
        reasons.append("wobble_significance_below_5_sigma")
    if delta_bic < 10.0:
        reasons.append("delta_bic_below_10")
    if not leave_one_out_stable:
        reasons.append("leave_one_epoch_out_not_stable")
    if len(residuals) < 6:
        reasons.append("fewer_than_6_independent_epochs")

    return OrbitFitResult(
        model_name="sirius_ab_fixed_orbit_plus_fixed_period_wobble",
        epoch_count=len(residuals),
        baseline_chi2=float(baseline_chi2),
        wobble_chi2=float(wobble["chi2"]),
        delta_bic=float(delta_bic),
        wobble_period_years=float(target_period_years),
        wobble_amplitude_mas=float(wobble["amplitude_mas"]),
        wobble_phase_degrees=float(wobble["phase_degrees"]),
        wobble_significance_sigma=float(significance),
        leave_one_out_stable=leave_one_out_stable,
        residual_rms_mas=float(residual_rms),
        reasons=tuple(reasons),
        residuals=tuple(residuals),
    )


def inject_synthetic_wobble(
    epochs: list[AstrometricEpoch],
    amplitude_mas: float,
    period_years: float = 6.3,
    phase_degrees: float = 0.0,
    orbit_elements: dict[str, float] | None = None,
) -> list[AstrometricEpoch]:
    """Return epochs with a deterministic circular sky wobble injected."""

    phase = math.radians(phase_degrees)
    injected = []
    for epoch in epochs:
        baseline_east, baseline_north = projected_binary_offset_arcsec(epoch.epoch_year, orbit_elements)
        angle = 2.0 * math.pi * epoch.epoch_year / period_years + phase
        wobble_east = amplitude_mas * math.sin(angle) / 1000.0
        wobble_north = amplitude_mas * math.cos(angle) / 1000.0
        injected.append(
            AstrometricEpoch(
                epoch_year=epoch.epoch_year,
                observed_delta_ra_arcsec=baseline_east + wobble_east,
                observed_delta_dec_arcsec=baseline_north + wobble_north,
                uncertainty_mas=epoch.uncertainty_mas,
                instrument=epoch.instrument,
                filter_name=epoch.filter_name,
                roll_angle_deg=epoch.roll_angle_deg,
                detector_x=epoch.detector_x,
                detector_y=epoch.detector_y,
                flags=epoch.flags,
            )
        )
    return injected


def astrometric_epoch_from_dict(row: dict[str, Any]) -> AstrometricEpoch:
    """Parse a JSON/CSV row into an astrometric epoch."""

    return AstrometricEpoch(
        epoch_year=float(row["epoch_year"]),
        observed_delta_ra_arcsec=float(row["observed_delta_ra_arcsec"]),
        observed_delta_dec_arcsec=float(row["observed_delta_dec_arcsec"]),
        uncertainty_mas=float(row["uncertainty_mas"]),
        instrument=str(row.get("instrument", "")),
        filter_name=str(row.get("filter_name", "")),
        roll_angle_deg=_optional_float(row.get("roll_angle_deg")),
        detector_x=_optional_float(row.get("detector_x")),
        detector_y=_optional_float(row.get("detector_y")),
        flags=tuple(row.get("flags", ())),
    )


def _fit_wobble_to_residuals(
    residuals: list[dict[str, float]],
    period_years: float,
) -> dict[str, float]:
    ra_sin, ra_cos = _fit_axis(residuals, "residual_ra_mas", period_years)
    dec_sin, dec_cos = _fit_axis(residuals, "residual_dec_mas", period_years)
    chi2 = 0.0
    for row in residuals:
        theta = 2.0 * math.pi * row["epoch_year"] / period_years
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        model_ra = ra_sin * sin_theta + ra_cos * cos_theta
        model_dec = dec_sin * sin_theta + dec_cos * cos_theta
        sigma = row["uncertainty_mas"]
        chi2 += ((row["residual_ra_mas"] - model_ra) / sigma) ** 2
        chi2 += ((row["residual_dec_mas"] - model_dec) / sigma) ** 2

    amplitude = math.sqrt(
        max(0.0, ra_sin * ra_sin + ra_cos * ra_cos + dec_sin * dec_sin + dec_cos * dec_cos) / 2.0
    )
    phase = math.degrees(math.atan2(ra_cos, ra_sin)) % 360.0
    return {
        "chi2": float(chi2),
        "amplitude_mas": float(amplitude),
        "phase_degrees": float(phase),
    }


def _fit_axis(
    residuals: list[dict[str, float]],
    key: str,
    period_years: float,
) -> tuple[float, float]:
    ss = sc = cc = sy = cy = 0.0
    for row in residuals:
        theta = 2.0 * math.pi * row["epoch_year"] / period_years
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        weight = 1.0 / (row["uncertainty_mas"] ** 2)
        value = row[key]
        ss += weight * sin_theta * sin_theta
        sc += weight * sin_theta * cos_theta
        cc += weight * cos_theta * cos_theta
        sy += weight * sin_theta * value
        cy += weight * cos_theta * value
    determinant = ss * cc - sc * sc
    if abs(determinant) < 1e-18:
        return 0.0, 0.0
    sin_coeff = (sy * cc - cy * sc) / determinant
    cos_coeff = (ss * cy - sc * sy) / determinant
    return sin_coeff, cos_coeff


def _chi2_zero(residuals: list[dict[str, float]]) -> float:
    chi2 = 0.0
    for row in residuals:
        sigma = row["uncertainty_mas"]
        chi2 += (row["residual_ra_mas"] / sigma) ** 2
        chi2 += (row["residual_dec_mas"] / sigma) ** 2
    return chi2


def _residual_rms_mas(residuals: list[dict[str, float]]) -> float:
    values = []
    for row in residuals:
        values.extend([row["residual_ra_mas"], row["residual_dec_mas"]])
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / len(values))


def _leave_one_out_stability(
    residuals: list[dict[str, float]],
    period_years: float,
    full_amplitude_mas: float,
    full_phase_degrees: float,
) -> bool:
    if len(residuals) < 6 or full_amplitude_mas <= 0.0:
        return False
    for index in range(len(residuals)):
        subset = residuals[:index] + residuals[index + 1 :]
        fit = _fit_wobble_to_residuals(subset, period_years)
        if fit["amplitude_mas"] < 0.5 * full_amplitude_mas:
            return False
        if _phase_distance_degrees(fit["phase_degrees"], full_phase_degrees) > 45.0:
            return False
    return True


def _phase_distance_degrees(a: float, b: float) -> float:
    raw = abs((a - b + 180.0) % 360.0 - 180.0)
    return min(raw, 360.0 - raw)


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)
