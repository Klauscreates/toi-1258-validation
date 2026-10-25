from microlens_ai.sirius import (
    AstrometricEpoch,
    GaiaSourceRecord,
    NACO_S27_PIXEL_SCALE_ARCSEC,
    ResidualCandidate,
    build_eso_sirius_raw_adql,
    build_gaia_dr3_field_adql,
    build_raw_astrometry_schema_adql,
    build_mast_sirius_criteria,
    build_source_hit_count_adql,
    classify_gaia_astrometry_surface,
    classify_residual_candidate,
    evaluate_audit_verdict,
    evaluate_ao_quicklook,
    evaluate_period_hypothesis,
    filter_preferred_wfc3_observations,
    fit_orbit_residuals,
    gaia_interpretation_notes,
    inject_synthetic_wobble,
    normalize_wfc3_observation,
    pixel_offset_from_sky_offset,
    projected_binary_offset_arcsec,
    run_mission3_package,
    summarize_ao_manifest,
    tap_json_to_rows,
)


def test_gaia_adql_includes_quality_columns_and_bright_caveat() -> None:
    query = build_gaia_dr3_field_adql(limit=25)

    assert "gaiadr3.gaia_source" in query
    assert "ruwe" in query
    assert "astrometric_excess_noise_sig" in query
    assert "duplicated_source" in query

    record = GaiaSourceRecord(
        source_id="123",
        ra_degrees=101.0,
        dec_degrees=-16.0,
        phot_g_mean_mag=-1.4,
    )

    assert record.is_too_bright_for_standard_gaia
    assert any("bright" in note for note in gaia_interpretation_notes([record]))


def test_mast_criteria_and_wfc3_filtering() -> None:
    criteria = build_mast_sirius_criteria()

    assert criteria["obs_collection"] == "HST"
    assert "WFC3/UVIS" in criteria["instrument_name"]
    assert criteria["dataproduct_type"] == "image"

    rows = [
        {
            "obs_id": "wfc3-good",
            "proposal_id": "12345",
            "instrument_name": "WFC3/UVIS",
            "filters": "F814W",
            "target_name": "SIRIUS",
            "t_min": 56000.0,
            "calib_level": 3,
            "dataproduct_type": "image",
            "intentType": "science",
            "dataRights": "PUBLIC",
        },
        {
            "obs_id": "acs-skip",
            "proposal_id": "12345",
            "instrument_name": "ACS/WFC",
            "filters": "F814W",
            "target_name": "SIRIUS",
            "t_min": 56000.0,
            "calib_level": 3,
            "dataproduct_type": "image",
            "intentType": "science",
            "dataRights": "PUBLIC",
        },
    ]

    observations = [normalize_wfc3_observation(row) for row in rows]
    preferred = filter_preferred_wfc3_observations(observations)

    assert [observation.obs_id for observation in preferred] == ["wfc3-good"]
    assert "geometric_distortion" in preferred[0].artifact_risks
    assert "diffraction_spikes" in preferred[0].artifact_risks


def test_synthetic_6_3_year_wobble_is_recovered() -> None:
    base_epochs = []
    for index in range(16):
        epoch_year = 2010.0 + 0.55 * index
        east, north = projected_binary_offset_arcsec(epoch_year)
        base_epochs.append(
            AstrometricEpoch(
                epoch_year=epoch_year,
                observed_delta_ra_arcsec=east,
                observed_delta_dec_arcsec=north,
                uncertainty_mas=0.4,
                instrument="WFC3/UVIS",
            )
        )
    injected = inject_synthetic_wobble(base_epochs, amplitude_mas=8.0, period_years=6.3)
    result = fit_orbit_residuals(injected, target_period_years=6.3)

    assert result.wobble_amplitude_mas > 7.0
    assert result.wobble_significance_sigma > 5.0
    assert result.delta_bic > 10.0
    assert result.leave_one_out_stable


def test_detector_fixed_residual_is_rejected() -> None:
    candidate = ResidualCandidate(
        candidate_id="detector-fixed",
        ra_degrees=101.0,
        dec_degrees=-16.0,
        peak_sigma=12.0,
        epoch_count=4,
        sky_repeat_count=1,
        detector_repeat_count=4,
        detector_correlation_abs=0.95,
        separation_from_sirius_a_arcsec=2.0,
        separation_from_sirius_b_arcsec=5.0,
    )

    assert classify_residual_candidate(candidate) == "rejected_detector_fixed"


def test_null_epochs_produce_not_confirmed_verdict() -> None:
    epochs = []
    for index in range(8):
        epoch_year = 2011.0 + index
        east, north = projected_binary_offset_arcsec(epoch_year)
        epochs.append(
            AstrometricEpoch(
                epoch_year=epoch_year,
                observed_delta_ra_arcsec=east,
                observed_delta_dec_arcsec=north,
                uncertainty_mas=1.0,
                instrument="WFC3/UVIS",
            )
        )
    result = fit_orbit_residuals(epochs, target_period_years=6.3)
    verdict = evaluate_audit_verdict(result, [])

    assert verdict.status == "not_confirmed"
    assert not verdict.sirius_c_confirmed
    assert not verdict.artificial_signal_evidence


def test_tap_json_to_rows_maps_gaia_metadata() -> None:
    document = {
        "metadata": [{"name": "source_id"}, {"name": "ra"}],
        "data": [[123, 101.2]],
    }

    assert tap_json_to_rows(document) == [{"source_id": 123, "ra": 101.2}]


def test_mission3_offline_package_is_signed_and_not_confirmed(tmp_path) -> None:
    outputs = run_mission3_package(tmp_path, live_gaia=False)

    assert outputs["package"].exists()
    assert outputs["audit"].exists()
    assert outputs["signature"].exists()
    assert (tmp_path / "checksums.sha256").exists()
    assert "not_confirmed" in (tmp_path / "audit_verdict.json").read_text(encoding="utf-8")


def test_400_year_hypothesis_rejected_without_epochs() -> None:
    result = fit_orbit_residuals([], target_period_years=400.0)
    audit = evaluate_period_hypothesis(result, "400_year_dagon_hypothesis")

    assert audit["status"] == "not_matched"
    assert "hypothesis_not_supported_by_available_astrometric_epochs" in audit["reasons"]


def test_raw_astrometry_schema_queries_are_guarded() -> None:
    query = build_raw_astrometry_schema_adql()
    count_query = build_source_hit_count_adql(
        "gaiadr3.nss_two_body_orbit",
        ["2947050466531873024"],
    )

    assert "TAP_SCHEMA.columns" in query
    assert "residual" in query.lower()
    assert "COUNT(*)" in count_query


def test_raw_astrometry_surface_classification_rejects_summaries() -> None:
    assert (
        classify_gaia_astrometry_surface("gaiadr3.gaia_source")
        == "catalog_summary_astrometry_not_epoch_residuals"
    )
    assert (
        classify_gaia_astrometry_surface("gaiadr3.nss_two_body_orbit")
        == "nss_fitted_solution_not_raw_epoch_astrometry"
    )
    assert (
        classify_gaia_astrometry_surface("gaiadr3.commanded_scan_law")
        == "spacecraft_scan_law_not_source_astrometry"
    )


def test_eso_ao_query_targets_sirius_and_instrument() -> None:
    query = build_eso_sirius_raw_adql("NAOS+CONICA", limit=25)

    assert "ivoa.ObsCore" in query
    assert "NAOS+CONICA" in query
    assert "101.28715533" in query
    assert "TOP 25" in query


def test_ao_manifest_summary_counts_science_objects() -> None:
    rows = [
        {"dp_cat": "SCIENCE", "dp_type": "OBJECT", "filter_path": "ND_SHORT,NB_2.17", "bytes": 10},
        {"dp_cat": "SCIENCE", "dp_type": "SKY", "filter_path": "ND_SHORT,NB_2.17", "bytes": 5},
        {"dp_cat": "ACQUISITION", "dp_type": "SKY", "filter_path": "ND_SHORT,NB_1.64", "bytes": 2},
    ]

    summary = summarize_ao_manifest(rows)

    assert summary["row_count"] == 3
    assert summary["total_bytes"] == 17
    assert summary["science_object_counts"] == {"ND_SHORT,NB_2.17": 1}


def test_ao_quicklook_stays_not_confirmed_for_single_epoch() -> None:
    dx, dy = pixel_offset_from_sky_offset(0.2715, 0.543, NACO_S27_PIXEL_SCALE_ARCSEC)

    assert round(dx, 3) == -10.0
    assert round(dy, 3) == 20.0

    verdict = evaluate_ao_quicklook(epoch_count=1, compact_candidate_count=1)

    assert verdict.status == "not_confirmed"
    assert not verdict.wobble_test_possible
    assert "single_or_sparse_epoch_cannot_test_6p3_year_or_400_year_wobble" in verdict.reasons
