from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .config import _load_yaml


VERDICT_HABITABLE = "likely_habitable"
VERDICT_UNINHABITABLE = "likely_uninhabitable"
VERDICT_UNRESOLVED = "unresolved_but_constrained"


@dataclass(frozen=True)
class ProximaTarget:
    """Known target facts needed to keep the mission honest."""

    star_name: str
    planet_name: str
    discovery_method: str
    transit_status: str
    mass_earth: float
    radius_earth_estimate: float
    orbital_period_days: float
    semi_major_axis_au: float
    eccentricity: float
    distance_light_years: float
    host_star_type: str
    approximate_j_mag: float
    non_transit_constraints: tuple[str, ...]


@dataclass(frozen=True)
class JwstMode:
    """One candidate JWST observing mode and its immediate caveats."""

    key: str
    instrument: str
    mode: str
    wavelength_min_um: float
    wavelength_max_um: float
    resolving_power: int
    observation_role: str
    optimized_for: str
    bright_limit_j_mag: float | None
    bright_limit_note: str


@dataclass(frozen=True)
class ScenarioModel:
    """Atmospheric or surface state to keep in the retrieval ensemble."""

    key: str
    prior_weight: float
    expected_verdict_bias: str


@dataclass(frozen=True)
class VerdictThresholds:
    decision_probability: float
    unresolved_probability: float
    probe_gate_probability: float
    min_data_quality_for_probe: float


@dataclass(frozen=True)
class HabitabilityEvidence:
    """Evidence summary used by the mission verdict classifier."""

    retained_atmosphere_probability: float
    liquid_water_compatibility_probability: float
    atmospheric_escape_risk: float
    stellar_contamination_risk: float
    data_quality: float


@dataclass(frozen=True)
class FeasibilityResult:
    mode_key: str
    instrument: str
    mode: str
    status: str
    saturation_margin_mag: float | None
    science_role: str
    blockers: tuple[str, ...]
    next_checks: tuple[str, ...]


@dataclass(frozen=True)
class VerdictResult:
    label: str
    probabilities: dict[str, float]
    go_probe: bool
    rationale: tuple[str, ...]


@dataclass(frozen=True)
class Mission2Config:
    """Typed Mission 2 configuration."""

    project_name: str
    science_goal: str
    generated_product: str
    target: ProximaTarget
    sources: dict[str, str]
    jwst_modes: tuple[JwstMode, ...]
    retrieval_molecules: tuple[str, ...]
    retrieval_thermal_or_orbital: tuple[str, ...]
    retrieval_stellar_activity: tuple[str, ...]
    scenario_models: tuple[ScenarioModel, ...]
    verdict_thresholds: VerdictThresholds
    baseline_evidence: HabitabilityEvidence
    raw: dict[str, Any]

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Mission2Config":
        raw = _load_yaml(Path(path))
        target = raw["target"]
        thresholds = raw["verdict_thresholds"]
        evidence = raw["baseline_evidence"]
        retrieval = raw["retrieval_targets"]
        project = raw["project"]

        return cls(
            project_name=project["name"],
            science_goal=project["science_goal"],
            generated_product=project["generated_product"],
            target=ProximaTarget(
                star_name=target["star_name"],
                planet_name=target["planet_name"],
                discovery_method=target["discovery_method"],
                transit_status=target["transit_status"],
                mass_earth=float(target["mass_earth"]),
                radius_earth_estimate=float(target["radius_earth_estimate"]),
                orbital_period_days=float(target["orbital_period_days"]),
                semi_major_axis_au=float(target["semi_major_axis_au"]),
                eccentricity=float(target["eccentricity"]),
                distance_light_years=float(target["distance_light_years"]),
                host_star_type=target["host_star_type"],
                approximate_j_mag=float(target["approximate_j_mag"]),
                non_transit_constraints=tuple(target.get("non_transit_constraints", ())),
            ),
            sources=dict(raw["sources"]),
            jwst_modes=tuple(
                _mode_from_dict(key, value)
                for key, value in raw["jwst_modes"].items()
            ),
            retrieval_molecules=tuple(retrieval["molecules"]),
            retrieval_thermal_or_orbital=tuple(retrieval["thermal_or_orbital"]),
            retrieval_stellar_activity=tuple(retrieval["stellar_activity"]),
            scenario_models=tuple(
                ScenarioModel(
                    key=key,
                    prior_weight=float(value["prior_weight"]),
                    expected_verdict_bias=value["expected_verdict_bias"],
                )
                for key, value in raw["scenario_models"].items()
            ),
            verdict_thresholds=VerdictThresholds(
                decision_probability=float(thresholds["decision_probability"]),
                unresolved_probability=float(thresholds["unresolved_probability"]),
                probe_gate_probability=float(thresholds["probe_gate_probability"]),
                min_data_quality_for_probe=float(thresholds["min_data_quality_for_probe"]),
            ),
            baseline_evidence=HabitabilityEvidence(
                retained_atmosphere_probability=float(evidence["retained_atmosphere_probability"]),
                liquid_water_compatibility_probability=float(
                    evidence["liquid_water_compatibility_probability"]
                ),
                atmospheric_escape_risk=float(evidence["atmospheric_escape_risk"]),
                stellar_contamination_risk=float(evidence["stellar_contamination_risk"]),
                data_quality=float(evidence["data_quality"]),
            ),
            raw=raw,
        )


def _mode_from_dict(key: str, raw: dict[str, Any]) -> JwstMode:
    return JwstMode(
        key=key,
        instrument=raw["instrument"],
        mode=raw["mode"],
        wavelength_min_um=float(raw["wavelength_min_um"]),
        wavelength_max_um=float(raw["wavelength_max_um"]),
        resolving_power=int(raw["resolving_power"]),
        observation_role=raw["observation_role"],
        optimized_for=raw["optimized_for"],
        bright_limit_j_mag=(
            None if raw.get("bright_limit_j_mag") is None else float(raw["bright_limit_j_mag"])
        ),
        bright_limit_note=raw["bright_limit_note"],
    )


def assess_jwst_feasibility(config: Mission2Config) -> list[FeasibilityResult]:
    """Produce an offline JWST feasibility matrix.

    The matrix is intentionally conservative. It does not replace an official
    ETC run; it identifies the proposal blockers that must be cleared there.
    """

    results: list[FeasibilityResult] = []
    for mode in config.jwst_modes:
        blockers: list[str] = []
        next_checks = [
            "Run the current JWST ETC with Proxima's spectrum, exact subarray, groups, and readout.",
            "Model stellar activity and contamination jointly with any planet signal.",
        ]
        saturation_margin = None

        if mode.bright_limit_j_mag is None:
            status = "requires_etc"
            blockers.append("No offline bright-limit benchmark is encoded for this setup.")
        else:
            saturation_margin = config.target.approximate_j_mag - mode.bright_limit_j_mag
            if saturation_margin < 0:
                status = "saturation_risk"
                blockers.append(
                    "Target is brighter than the encoded benchmark bright limit."
                )
            elif saturation_margin < 0.5:
                status = "tight_brightness_margin"
                blockers.append("Brightness margin is small; ETC and acquisition checks are mandatory.")
            else:
                status = "brightness_plausible"

        if (
            config.target.transit_status == "likely_non_transiting"
            and mode.optimized_for == "transit_spectroscopy"
        ):
            blockers.append("Likely non-transiting geometry weakens classic transmission spectroscopy.")
            if status == "brightness_plausible":
                status = "geometry_limited"
            elif status == "tight_brightness_margin":
                status = "geometry_and_brightness_limited"

        if mode.observation_role.startswith("stellar_context"):
            next_checks.append("Use as support data for stellar and contamination priors, not as a direct habitability verdict.")

        results.append(
            FeasibilityResult(
                mode_key=mode.key,
                instrument=mode.instrument,
                mode=mode.mode,
                status=status,
                saturation_margin_mag=saturation_margin,
                science_role=mode.observation_role,
                blockers=tuple(blockers),
                next_checks=tuple(next_checks),
            )
        )
    return results


def classify_habitability(
    evidence: HabitabilityEvidence,
    thresholds: VerdictThresholds,
) -> VerdictResult:
    """Classify the current evidence into the three Mission 2 verdict classes."""

    retained = _clamp01(evidence.retained_atmosphere_probability)
    water = _clamp01(evidence.liquid_water_compatibility_probability)
    escape = _clamp01(evidence.atmospheric_escape_risk)
    contamination = _clamp01(evidence.stellar_contamination_risk)
    quality = _clamp01(evidence.data_quality)

    habitable_raw = retained * water * (1.0 - escape) * max(0.25, quality)
    uninhabitable_raw = (0.45 * (1.0 - retained)) + (0.35 * escape) + (0.20 * (1.0 - water))
    unresolved_raw = max(0.05, 1.0 - quality) + (0.5 * contamination)

    total = habitable_raw + uninhabitable_raw + unresolved_raw
    probabilities = {
        VERDICT_HABITABLE: habitable_raw / total,
        VERDICT_UNINHABITABLE: uninhabitable_raw / total,
        VERDICT_UNRESOLVED: unresolved_raw / total,
    }

    top_label = max(probabilities, key=probabilities.__getitem__)
    top_probability = probabilities[top_label]
    if (
        probabilities[VERDICT_UNRESOLVED] >= thresholds.unresolved_probability
        or top_probability < thresholds.decision_probability
    ):
        label = VERDICT_UNRESOLVED
    else:
        label = top_label

    go_probe = (
        probabilities[VERDICT_HABITABLE] >= thresholds.probe_gate_probability
        and quality >= thresholds.min_data_quality_for_probe
    )
    rationale = _verdict_rationale(label, probabilities, evidence, go_probe)
    return VerdictResult(
        label=label,
        probabilities={key: round(value, 6) for key, value in probabilities.items()},
        go_probe=go_probe,
        rationale=tuple(rationale),
    )


def build_mission_package(config: Mission2Config) -> dict[str, Any]:
    """Build the full mission package payload used by docs and scripts."""

    feasibility = assess_jwst_feasibility(config)
    verdict = classify_habitability(config.baseline_evidence, config.verdict_thresholds)
    return {
        "mission": {
            "project_name": config.project_name,
            "science_goal": config.science_goal,
            "generated_product": config.generated_product,
        },
        "target": asdict(config.target),
        "sources": config.sources,
        "jwst_feasibility": [asdict(row) for row in feasibility],
        "retrieval_targets": {
            "molecules": list(config.retrieval_molecules),
            "thermal_or_orbital": list(config.retrieval_thermal_or_orbital),
            "stellar_activity": list(config.retrieval_stellar_activity),
        },
        "scenario_models": [asdict(row) for row in config.scenario_models],
        "baseline_evidence": asdict(config.baseline_evidence),
        "verdict": asdict(verdict),
        "rover_probe_gate": {
            "status": "go" if verdict.go_probe else "no_go",
            "rule": "Open the interstellar rover/probe track only if the habitable posterior and data quality clear the configured thresholds.",
        },
    }


def build_mission_package_with_archive(
    config: Mission2Config,
    archive_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the mission package and attach live JWST archive state."""

    package = build_mission_package(config)
    if archive_summary is not None:
        package["jwst_archive"] = archive_summary
        package["verdict"]["archive_limited"] = (
            archive_summary.get("nirspec_or_niriss_science_product_count", 0) == 0
        )
        if package["verdict"]["archive_limited"]:
            package["verdict"]["rationale"] = list(package["verdict"]["rationale"]) + [
                "Live MAST pull did not find retrievable NIRSpec/NIRISS science products for a stronger atmospheric verdict."
            ]
    return package


def write_mission_outputs(
    config: Mission2Config,
    output_dir: str | Path,
    archive_summary: dict[str, Any] | None = None,
) -> dict[str, Path]:
    """Write JSON, CSV, and Markdown products for Mission 2."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    package = build_mission_package_with_archive(config, archive_summary=archive_summary)

    package_path = output_path / "mission_package.json"
    package_path.write_text(json.dumps(package, indent=2), encoding="utf-8")

    feasibility_path = output_path / "jwst_feasibility_matrix.csv"
    _write_csv(feasibility_path, package["jwst_feasibility"])

    scenarios_path = output_path / "retrieval_scenario_priors.csv"
    _write_csv(scenarios_path, package["scenario_models"])

    report_path = output_path / "verdict_report.md"
    report_path.write_text(render_verdict_report(package), encoding="utf-8")

    return {
        "package": package_path,
        "feasibility": feasibility_path,
        "scenarios": scenarios_path,
        "report": report_path,
    }


def render_verdict_report(package: dict[str, Any]) -> str:
    verdict = package["verdict"]
    target = package["target"]
    probabilities = verdict["probabilities"]
    lines = [
        "# Mission 2: Proxima b Habitability Verdict",
        "",
        f"Target: {target['planet_name']} around {target['star_name']}.",
        f"Current verdict: **{verdict['label']}**.",
        f"Future rover/probe gate: **{package['rover_probe_gate']['status']}**.",
        "",
        "## Posterior Probabilities",
        "",
    ]
    for key in (VERDICT_HABITABLE, VERDICT_UNINHABITABLE, VERDICT_UNRESOLVED):
        lines.append(f"- {key}: {probabilities[key]:.3f}")

    lines.extend(
        [
            "",
            "## JWST Feasibility",
            "",
            "| Mode | Status | Science role | Primary blockers |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in package["jwst_feasibility"]:
        blockers = "; ".join(row["blockers"]) if row["blockers"] else "None encoded"
        lines.append(
            f"| {row['instrument']} {row['mode']} | {row['status']} | {row['science_role']} | {blockers} |"
        )

    lines.extend(
        [
            "",
            "## Retrieval Ensemble",
            "",
            "Use the scenario priors as an initial bookkeeping scaffold, not as fixed astrophysical truth.",
        ]
    )
    for scenario in package["scenario_models"]:
        lines.append(
            f"- {scenario['key']}: prior {scenario['prior_weight']:.2f}, bias {scenario['expected_verdict_bias']}"
        )

    lines.extend(
        [
            "",
            "## Rationale",
            "",
        ]
    )
    for item in verdict["rationale"]:
        lines.append(f"- {item}")

    if "jwst_archive" in package:
        archive = package["jwst_archive"]
        lines.extend(
            [
                "",
                "## Live JWST Archive Pull",
                "",
                f"- observations_found: {archive['observation_count']}",
                f"- public_observations: {archive['public_observation_count']}",
                f"- science_products_found: {archive['science_product_count']}",
                f"- retrievable_products: {archive.get('retrievable_product_count', 0)}",
                f"- nirspec_or_niriss_science_products: {archive['nirspec_or_niriss_science_product_count']}",
                f"- miri_observations: {archive.get('miri_observation_count', 0)}",
                f"- nirspec_observations: {archive.get('nirspec_observation_count', 0)}",
                f"- phase_curve_candidates: {archive.get('phase_curve_candidate_count', 0)}",
                f"- water_spectroscopy_candidates: {archive.get('water_spectroscopy_candidate_count', 0)}",
                f"- proposal_types: {', '.join(archive.get('proposal_types', ())) if archive.get('proposal_types') else 'none'}",
                f"- instruments: {', '.join(archive['instruments']) if archive['instruments'] else 'none'}",
                f"- verdict_impact: {archive['verdict_impact']}",
                f"- next_action: {archive['next_action']}",
            ]
        )

    lines.extend(
        [
            "",
            "## Sources",
            "",
        ]
    )
    for name, url in package["sources"].items():
        lines.append(f"- {name}: {url}")
    lines.append("")
    return "\n".join(lines)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key)) for key in fields})


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict)):
        return json.dumps(value)
    return value


def _verdict_rationale(
    label: str,
    probabilities: dict[str, float],
    evidence: HabitabilityEvidence,
    go_probe: bool,
) -> list[str]:
    rationale = [
        "The classifier separates habitability from biosignature detection; life is not claimed.",
        "Likely non-transiting geometry keeps JWST transmission spectroscopy high risk.",
    ]
    if label == VERDICT_UNRESOLVED:
        rationale.append(
            "Current evidence quality is not strong enough to choose habitable or uninhabitable."
        )
    elif label == VERDICT_HABITABLE:
        rationale.append("Atmosphere retention and liquid-water compatibility dominate the posterior.")
    else:
        rationale.append("Atmosphere-loss and climate-incompatibility terms dominate the posterior.")

    if evidence.stellar_contamination_risk >= 0.4:
        rationale.append("Stellar activity remains a major false-positive and false-negative control.")
    if not go_probe:
        rationale.append("The future rover/probe track remains gated off until the habitable posterior clears the threshold.")
    else:
        rationale.append("The future probe track is open for concept study, not immediate launch.")
    rationale.append(
        "Posterior probabilities: "
        + ", ".join(f"{key}={value:.3f}" for key, value in probabilities.items())
    )
    return rationale


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))
