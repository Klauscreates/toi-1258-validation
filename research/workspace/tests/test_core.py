from pathlib import Path

from microlens_ai.config import ArchitectureConfig
from microlens_ai.mast_archive import (
    MastObservation,
    MastProduct,
    build_retrieval_manifest,
    summarize_jwst_archive,
)
from microlens_ai.model import describe_model
from microlens_ai.proxima_mission import (
    VERDICT_HABITABLE,
    VERDICT_UNINHABITABLE,
    VERDICT_UNRESOLVED,
    HabitabilityEvidence,
    Mission2Config,
    build_mission_package,
    classify_habitability,
)
from microlens_ai.scoring import calibrated_score
from microlens_ai.simulation import paczynski_magnification


ROOT = Path(__file__).resolve().parents[1]


def test_config_loads() -> None:
    config = ArchitectureConfig.from_yaml(ROOT / "configs" / "roman_hubble_m31.yaml")

    assert config.project_name == "roman_hubble_m31_microlensing"
    assert config.roman.archive == "MAST"
    assert config.hubble.starting_collections == ("PHAT", "Hubble Legacy Archive M31 programs with repeated visits")
    assert config.candidate_generation.min_epochs == 3
    assert "image_artifact" in config.model.classes


def test_paczynski_peak_exceeds_baseline() -> None:
    peak = paczynski_magnification(time=10.0, t0=10.0, t_e=30.0, u0=0.1)
    wing = paczynski_magnification(time=100.0, t0=10.0, t_e=30.0, u0=0.1)

    assert peak > wing
    assert wing > 1.0


def test_scoring_prefers_lensing_logit() -> None:
    score = calibrated_score(
        microlensing_logit=4.0,
        variable_logit=1.0,
        artifact_logit=0.0,
        uncertainty=0.2,
    )

    assert score.microlensing_probability > score.variable_probability
    assert score.vetting_priority() > 0


def test_model_description_is_serializable() -> None:
    description = describe_model()

    assert description["name"] == "hybrid_microlensing_transformer"
    assert "image" in description["branches"]


def test_mission2_config_loads() -> None:
    config = Mission2Config.from_yaml(ROOT / "configs" / "proxima_b_mission2.yaml")

    assert config.project_name == "mission2_proxima_b_habitability"
    assert config.target.discovery_method == "radial_velocity"
    assert config.target.transit_status == "likely_non_transiting"
    assert {mode.instrument for mode in config.jwst_modes} == {"NIRISS", "NIRSpec"}
    assert len(config.scenario_models) == 6


def test_mission2_feasibility_tracks_non_transit_constraint() -> None:
    config = Mission2Config.from_yaml(ROOT / "configs" / "proxima_b_mission2.yaml")
    package = build_mission_package(config)

    soss = next(row for row in package["jwst_feasibility"] if row["mode_key"] == "niriss_soss_multistripe")
    bots = next(row for row in package["jwst_feasibility"] if row["mode_key"] == "nirspec_bots_prism")

    assert soss["status"] == "geometry_limited"
    assert any("non-transiting" in blocker for blocker in soss["blockers"])
    assert bots["status"] == "saturation_risk"
    assert bots["saturation_margin_mag"] < 0


def test_mission2_baseline_verdict_is_unresolved() -> None:
    config = Mission2Config.from_yaml(ROOT / "configs" / "proxima_b_mission2.yaml")
    package = build_mission_package(config)

    assert package["verdict"]["label"] == VERDICT_UNRESOLVED
    assert package["rover_probe_gate"]["status"] == "no_go"


def test_mission2_verdict_can_promote_strong_evidence() -> None:
    config = Mission2Config.from_yaml(ROOT / "configs" / "proxima_b_mission2.yaml")

    habitable = classify_habitability(
        HabitabilityEvidence(
            retained_atmosphere_probability=0.92,
            liquid_water_compatibility_probability=0.88,
            atmospheric_escape_risk=0.08,
            stellar_contamination_risk=0.05,
            data_quality=0.92,
        ),
        config.verdict_thresholds,
    )
    uninhabitable = classify_habitability(
        HabitabilityEvidence(
            retained_atmosphere_probability=0.08,
            liquid_water_compatibility_probability=0.12,
            atmospheric_escape_risk=0.9,
            stellar_contamination_risk=0.05,
            data_quality=0.92,
        ),
        config.verdict_thresholds,
    )

    assert habitable.label == VERDICT_HABITABLE
    assert habitable.go_probe
    assert uninhabitable.label == VERDICT_UNINHABITABLE
    assert not uninhabitable.go_probe


def test_jwst_archive_summary_flags_missing_spectroscopy_products() -> None:
    observations = [
        MastObservation(
            obsid=1,
            obs_id="jw_test_miri",
            obs_collection="JWST",
            provenance_name="APT",
            proposal_id="8581",
            proposal_type="SURVEY",
            proposal_pi="Bowens-Rubin",
            target_name="NAME-Proxima-Centauri",
            target_classification="STAR",
            instrument_name="MIRI/IMAGE",
            dataproduct_type="image",
            filters="F2100W",
            wave_region="Infrared",
            wavelength_region="Infrared",
            em_min=18.0,
            em_max=24.0,
            wave_min=18.0,
            wave_max=24.0,
            calib_level=-1,
            t_min=None,
            t_max=None,
            t_exptime=190.0,
            t_obs_release=None,
            data_rights="PUBLIC",
            obs_title="The HOTH Survey",
            data_url=None,
            s_ra=217.42894,
            s_dec=-62.67949,
        )
    ]

    summary = summarize_jwst_archive(observations, [], query_patterns=("%PROXIMA%",))

    assert summary.observation_count == 1
    assert summary.public_observation_count == 1
    assert summary.science_product_count == 0
    assert summary.nirspec_or_niriss_science_product_count == 0
    assert "cannot move the habitability verdict" in summary.verdict_impact


def test_jwst_archive_summary_recognizes_nirspec_products() -> None:
    observations = [
        MastObservation(
            obsid=2,
            obs_id="jw_test_nirspec",
            obs_collection="JWST",
            provenance_name="CALJWST",
            proposal_id="9999",
            proposal_type="ERS",
            proposal_pi="Example",
            target_name="NAME-Proxima-Centauri",
            target_classification="EXOPLANET",
            instrument_name="NIRSpec/BOTS",
            dataproduct_type="spectrum",
            filters="PRISM/CLEAR",
            wave_region="Infrared",
            wavelength_region="Infrared",
            em_min=0.6,
            em_max=5.3,
            wave_min=0.6,
            wave_max=5.3,
            calib_level=2,
            t_min=60000.0,
            t_max=60000.1,
            t_exptime=1000.0,
            t_obs_release=60030.0,
            data_rights="PUBLIC",
            obs_title="Proxima spectrum",
            data_url="mast:JWST/product/jw_test_nirspec_x1d.fits",
            s_ra=217.42894,
            s_dec=-62.67949,
        )
    ]
    products = [
        MastProduct(
            obsid=2,
            obs_id="jw_test_nirspec",
            product_filename="jw_test_nirspec_x1d.fits",
            product_type="SCIENCE",
            product_group_description="Minimum Recommended Products",
            calib_level=2,
            size=123,
            data_uri="mast:JWST/product/jw_test_nirspec_x1d.fits",
            data_rights="PUBLIC",
        )
    ]

    summary = summarize_jwst_archive(observations, products, query_patterns=("%PROXIMA%",))

    assert summary.science_product_count == 1
    assert summary.nirspec_or_niriss_science_product_count == 1
    assert "available for retrieval" in summary.verdict_impact


def test_retrieval_manifest_marks_nirspec_timeseries_as_water_screen() -> None:
    observations = [
        MastObservation(
            obsid=3,
            obs_id="jw_trappist_nirspec",
            obs_collection="JWST",
            provenance_name="CALJWST",
            proposal_id="9256",
            proposal_type="GO",
            proposal_pi="Example",
            target_name="TRAPPIST-1",
            target_classification="STAR",
            instrument_name="NIRSPEC/SLIT",
            dataproduct_type="timeseries",
            filters="CLEAR;PRISM",
            wave_region="Infrared",
            wavelength_region="Infrared",
            em_min=0.6,
            em_max=5.3,
            wave_min=0.6,
            wave_max=5.3,
            calib_level=2,
            t_min=61000.0,
            t_max=61000.2,
            t_exptime=25000.0,
            t_obs_release=61001.0,
            data_rights="PUBLIC",
            obs_title="Using stellar contamination proxy",
            data_url=None,
            s_ra=None,
            s_dec=None,
        )
    ]
    products = [
        MastProduct(
            obsid=3,
            obs_id="jw_trappist_nirspec",
            product_filename="jw_trappist_nirspec_rate.fits",
            product_type="SCIENCE",
            product_group_description="Minimum Recommended Products",
            calib_level=2,
            size=403200,
            data_uri="mast:JWST/product/jw_trappist_nirspec_rate.fits",
            data_rights="PUBLIC",
        )
    ]

    manifest = build_retrieval_manifest(observations, products)

    assert manifest[0].target_name == "TRAPPIST-1"
    assert manifest[0].analysis_role == "spectroscopy_water_screen"
    assert "H2O_1.4um" in manifest[0].molecular_targets
    assert manifest[0].download_recommended
