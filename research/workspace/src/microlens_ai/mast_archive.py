from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


MAST_INVOKE = "https://mast.stsci.edu/api/v0/invoke"
PROXIMA_RA_DEG = 217.4289420833333
PROXIMA_DEC_DEG = -62.6794888888889
PROXIMA_B_PERIOD_DAYS = 11.2
DIRECTOR_OR_ERS_TYPES = {"DD", "DDT", "GO/DD", "ERS"}


@dataclass(frozen=True)
class MastObservation:
    """A compact public-archive observation row from MAST CAOM."""

    obsid: int
    obs_id: str
    obs_collection: str
    provenance_name: str | None
    proposal_id: str
    proposal_type: str | None
    proposal_pi: str | None
    target_name: str
    target_classification: str | None
    instrument_name: str
    dataproduct_type: str
    filters: str | None
    wave_region: str | None
    wavelength_region: str | None
    em_min: float | None
    em_max: float | None
    wave_min: float | None
    wave_max: float | None
    calib_level: int | None
    t_min: float | None
    t_max: float | None
    t_exptime: float | None
    t_obs_release: float | None
    data_rights: str | None
    obs_title: str | None
    data_url: str | None
    s_ra: float | None
    s_dec: float | None


@dataclass(frozen=True)
class MastProduct:
    """A compact related-product row for a MAST observation."""

    obsid: int | None
    obs_id: str | None
    product_filename: str | None
    product_type: str | None
    product_group_description: str | None
    calib_level: int | None
    size: int | None
    data_uri: str | None
    data_rights: str | None


@dataclass(frozen=True)
class RetrievalManifestRow:
    """A mission-facing row that says whether a record can feed analysis."""

    source_kind: str
    obsid: int | None
    obs_id: str | None
    target_name: str | None
    proposal_id: str | None
    proposal_type: str | None
    instrument_name: str | None
    filter_name: str | None
    dataproduct_type: str | None
    product_filename: str | None
    product_type: str | None
    product_size: int | None
    data_uri: str | None
    data_rights: str | None
    access_class: str
    analysis_role: str
    molecular_targets: str
    orbital_phase: float | None
    download_recommended: bool
    blocker: str | None


@dataclass(frozen=True)
class JwstArchiveSummary:
    """Mission-facing summary of the live JWST archive pull."""

    query_patterns: tuple[str, ...]
    instrument_patterns: tuple[str, ...]
    program_type_filters: tuple[str, ...]
    search_radius_degrees: float | None
    observation_count: int
    public_observation_count: int
    restricted_observation_count: int
    science_product_count: int
    public_product_count: int
    restricted_product_count: int
    retrievable_product_count: int
    miri_observation_count: int
    nirspec_observation_count: int
    miri_science_product_count: int
    nirspec_science_product_count: int
    nirspec_or_niriss_observation_count: int
    nirspec_or_niriss_science_product_count: int
    spectroscopy_observation_count: int
    imaging_observation_count: int
    phase_curve_candidate_count: int
    water_spectroscopy_candidate_count: int
    director_or_ers_observation_count: int
    planned_or_unproductized_count: int
    instruments: tuple[str, ...]
    filters: tuple[str, ...]
    proposal_ids: tuple[str, ...]
    proposal_types: tuple[str, ...]
    verdict_impact: str
    next_action: str


def mast_request(
    service: str,
    params: dict[str, Any],
    pagesize: int = 500,
    page: int = 1,
    timeout_seconds: int = 90,
    max_polls: int = 6,
) -> dict[str, Any]:
    """Run a MAST API request using only the Python standard library."""

    request = {
        "service": service,
        "format": "json",
        "params": params,
        "pagesize": pagesize,
        "page": page,
        "removenullcolumns": False,
        "timeout": min(timeout_seconds, 60),
    }
    for poll_index in range(max_polls):
        payload = _post_mast_request(request, timeout_seconds=timeout_seconds)
        status = payload.get("status")
        if status == "COMPLETE":
            return payload
        if status != "EXECUTING":
            raise RuntimeError(f"MAST request failed: {payload.get('msg', payload)}")
        time.sleep(min(2 + poll_index, 8))
    raise TimeoutError(f"MAST request did not complete after {max_polls} polls")


def query_jwst_proxima_observations(
    target_patterns: tuple[str, ...] = ("%PROXIMA%",),
    instrument_patterns: tuple[str, ...] = ("%MIRI%", "%NIRSPEC%"),
    program_types: tuple[str, ...] = (),
    include_coordinate_search: bool = True,
    ra_degrees: float = PROXIMA_RA_DEG,
    dec_degrees: float = PROXIMA_DEC_DEG,
    radius_degrees: float = 0.05,
    pagesize: int = 200,
    timeout_seconds: int = 90,
) -> list[MastObservation]:
    """Pull JWST observation metadata matching Proxima target-name aliases.

    Data-rights filters are deliberately omitted. The mission wants to screen
    immediately public, restricted, ERS, DDT, and planned rows separately rather
    than hiding any class of archive evidence.
    """

    rows_by_obsid: dict[int, MastObservation] = {}
    columns = (
        "obsid,obs_id,obs_collection,provenance_name,proposal_id,proposal_type,"
        "proposal_pi,target_name,target_classification,instrument_name,filters,"
        "wave_region,wavelength_region,em_min,em_max,wave_min,wave_max,"
        "dataproduct_type,calib_level,t_min,t_max,t_exptime,t_obs_release,"
        "dataRights,obs_title,dataURL,s_ra,s_dec"
    )
    base_filters = [{"paramName": "obs_collection", "values": ["JWST"]}]
    if program_types:
        base_filters.append({"paramName": "proposal_type", "values": list(program_types)})

    instrument_filter_sets = _instrument_filter_sets(instrument_patterns)
    for pattern in target_patterns:
        for instrument_filters in instrument_filter_sets:
            filters = [
                *base_filters,
                *instrument_filters,
                {"paramName": "target_name", "values": [], "freeText": pattern},
            ]
            for row in _paged_mast_rows(
                "Mast.Caom.Filtered",
                {"columns": columns, "filters": filters},
                pagesize=pagesize,
                timeout_seconds=timeout_seconds,
            ):
                observation = _observation_from_row(row)
                rows_by_obsid[observation.obsid] = observation

    if include_coordinate_search:
        position = f"{ra_degrees}, {dec_degrees}, {radius_degrees}"
        for instrument_filters in instrument_filter_sets:
            filters = [*base_filters, *instrument_filters]
            for row in _paged_mast_rows(
                "Mast.Caom.Filtered.Position",
                {"columns": columns, "filters": filters, "position": position},
                pagesize=pagesize,
                timeout_seconds=timeout_seconds,
            ):
                observation = _observation_from_row(row)
                rows_by_obsid[observation.obsid] = observation

    return sorted(rows_by_obsid.values(), key=lambda row: (row.proposal_id, row.obs_id))


def query_jwst_target_observations(
    target_patterns: tuple[str, ...],
    instrument_patterns: tuple[str, ...] = ("%MIRI%", "%NIRSPEC%"),
    program_types: tuple[str, ...] = (),
    pagesize: int = 200,
    timeout_seconds: int = 90,
) -> list[MastObservation]:
    """Pull JWST observations for target-name patterns without a coordinate cone."""

    return query_jwst_proxima_observations(
        target_patterns=target_patterns,
        instrument_patterns=instrument_patterns,
        program_types=program_types,
        include_coordinate_search=False,
        pagesize=pagesize,
        timeout_seconds=timeout_seconds,
    )


def query_jwst_keyword_screen(
    instrument: str,
    target_patterns: tuple[str, ...] = ("%PROXIMA%",),
    title_patterns: tuple[str, ...] = (),
    pagesize: int = 200,
    timeout_seconds: int = 90,
) -> list[dict[str, Any]]:
    """Screen JWST instrument keyword tables for processed exposure products."""

    service = {
        "MIRI": "Mast.Jwst.Filtered.Miri",
        "NIRSPEC": "Mast.Jwst.Filtered.Nirspec",
        "NIRISS": "Mast.Jwst.Filtered.Niriss",
    }.get(instrument.upper())
    if service is None:
        raise ValueError(f"Unsupported JWST keyword instrument: {instrument}")

    rows_by_filename: dict[str, dict[str, Any]] = {}
    for pattern in target_patterns:
        for row in _paged_mast_rows(
            service,
            {
                "columns": "*",
                "filters": [{"paramName": "targname", "values": [], "freeText": pattern}],
            },
            pagesize=pagesize,
            timeout_seconds=timeout_seconds,
        ):
            rows_by_filename[str(row.get("filename") or row.get("obs_id") or id(row))] = row

    for pattern in title_patterns:
        for row in _paged_mast_rows(
            service,
            {
                "columns": "*",
                "filters": [{"paramName": "title", "values": [], "freeText": pattern}],
            },
            pagesize=pagesize,
            timeout_seconds=timeout_seconds,
        ):
            rows_by_filename[str(row.get("filename") or row.get("obs_id") or id(row))] = row
    return list(rows_by_filename.values())


def query_jwst_keyword_screens(
    instruments: tuple[str, ...],
    target_patterns: tuple[str, ...] = ("%PROXIMA%",),
    title_patterns: tuple[str, ...] = (),
    timeout_seconds: int = 90,
) -> dict[str, list[dict[str, Any]]]:
    return {
        instrument.upper(): query_jwst_keyword_screen(
            instrument,
            target_patterns=target_patterns,
            title_patterns=title_patterns,
            timeout_seconds=timeout_seconds,
        )
        for instrument in instruments
    }


def _paged_mast_rows(
    service: str,
    params: dict[str, Any],
    pagesize: int,
    timeout_seconds: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 1
    while True:
        payload = mast_request(
            service,
            params,
            pagesize=pagesize,
            page=page,
            timeout_seconds=timeout_seconds,
        )
        page_rows = payload.get("data", [])
        rows.extend(page_rows)
        paging = payload.get("paging", {})
        total = int(paging.get("rowsTotal") or len(page_rows))
        if page * pagesize >= total or not page_rows:
            break
        page += 1
    return rows


def _instrument_filter_sets(instrument_patterns: tuple[str, ...]) -> list[list[dict[str, Any]]]:
    if not instrument_patterns:
        return [[]]
    return [
        [{"paramName": "instrument_name", "values": [], "freeText": pattern}]
        for pattern in instrument_patterns
    ]


def query_products_for_observations(
    observations: list[MastObservation],
    timeout_seconds: int = 90,
) -> list[MastProduct]:
    """Fetch related products for discovered observations without data-rights filtering."""

    products_by_key: dict[tuple[int | None, str | None, str | None], MastProduct] = {}
    for observation in observations:
        payload = mast_request(
            "Mast.Caom.Products",
            {"obsid": observation.obsid},
            pagesize=1000,
            page=1,
            timeout_seconds=timeout_seconds,
        )
        for row in payload.get("data", []):
            product = _product_from_row(row)
            key = (product.obsid, product.obs_id, product.product_filename)
            products_by_key[key] = product
    return list(products_by_key.values())


def build_retrieval_manifest(
    observations: list[MastObservation],
    products: list[MastProduct],
    orbital_period_days: float = PROXIMA_B_PERIOD_DAYS,
    reference_epoch_mjd: float = 0.0,
) -> list[RetrievalManifestRow]:
    products_by_obsid: dict[int | None, list[MastProduct]] = {}
    for product in products:
        products_by_obsid.setdefault(product.obsid, []).append(product)

    rows: list[RetrievalManifestRow] = []
    for observation in observations:
        observation_products = products_by_obsid.get(observation.obsid, [])
        if not observation_products:
            role, molecules = _analysis_role_for_observation(observation, None)
            rows.append(
                RetrievalManifestRow(
                    source_kind="observation_metadata",
                    obsid=observation.obsid,
                    obs_id=observation.obs_id,
                    target_name=observation.target_name,
                    proposal_id=observation.proposal_id,
                    proposal_type=observation.proposal_type,
                    instrument_name=observation.instrument_name,
                    filter_name=observation.filters,
                    dataproduct_type=observation.dataproduct_type,
                    product_filename=None,
                    product_type=None,
                    product_size=None,
                    data_uri=observation.data_url,
                    data_rights=observation.data_rights,
                    access_class=_access_class(observation.data_rights),
                    analysis_role=role,
                    molecular_targets=molecules,
                    orbital_phase=_orbital_phase(
                        observation.t_min,
                        orbital_period_days=orbital_period_days,
                        reference_epoch_mjd=reference_epoch_mjd,
                    ),
                    download_recommended=False,
                    blocker="No related MAST product rows were returned for this observation.",
                )
            )
            continue

        for product in observation_products:
            role, molecules = _analysis_role_for_observation(observation, product)
            blocker = None
            if not product.data_uri:
                blocker = "Product has no dataURI."
            elif product.data_rights and product.data_rights != "PUBLIC":
                blocker = "Product may require authenticated archive access."
            elif product.size is not None and product.size > 250_000_000:
                blocker = "Large product; increase download size limit if raw integrations are required."
            rows.append(
                RetrievalManifestRow(
                    source_kind="product",
                    obsid=observation.obsid,
                    obs_id=observation.obs_id,
                    target_name=observation.target_name,
                    proposal_id=observation.proposal_id,
                    proposal_type=observation.proposal_type,
                    instrument_name=observation.instrument_name,
                    filter_name=observation.filters,
                    dataproduct_type=observation.dataproduct_type,
                    product_filename=product.product_filename,
                    product_type=product.product_type,
                    product_size=product.size,
                    data_uri=product.data_uri,
                    data_rights=product.data_rights or observation.data_rights,
                    access_class=_access_class(product.data_rights or observation.data_rights),
                    analysis_role=role,
                    molecular_targets=molecules,
                    orbital_phase=_orbital_phase(
                        observation.t_min,
                        orbital_period_days=orbital_period_days,
                        reference_epoch_mjd=reference_epoch_mjd,
                    ),
                    download_recommended=bool(
                        product.data_uri
                        and product.product_type == "SCIENCE"
                        and _is_download_preferred_product(product)
                    ),
                    blocker=blocker,
                )
            )
    return rows


def download_recommended_products(
    manifest_rows: list[RetrievalManifestRow],
    output_dir: str | Path,
    limit: int | None = None,
    max_size_bytes: int = 250_000_000,
    timeout_seconds: int = 180,
) -> list[dict[str, Any]]:
    """Download recommended MAST products in the retrieval manifest."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    selected = [
        row
        for row in manifest_rows
        if row.download_recommended and row.data_uri and row.product_filename
    ]
    if limit is not None:
        selected = selected[:limit]

    downloads: list[dict[str, Any]] = []
    for row in selected:
        destination = output_path / _safe_filename(row.product_filename or row.data_uri or "product.fits")
        uri = row.data_uri or ""
        url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urllib.parse.urlencode({"uri": uri})
        status = "downloaded"
        error = None
        if row.product_size is not None and row.product_size > max_size_bytes:
            status = "skipped_size_limit"
            error = f"size={row.product_size} exceeds max_size_bytes={max_size_bytes}"
        else:
            try:
                with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
                    destination.write_bytes(response.read())
            except Exception as exc:  # pragma: no cover - network failure path
                status = "error"
                error = str(exc)
        if status != "downloaded" and destination.exists():
            destination.unlink()
        downloads.append(
            {
                "product_filename": row.product_filename,
                "data_uri": row.data_uri,
                "product_size": row.product_size,
                "destination": str(destination),
                "status": status,
                "error": error,
            }
        )
    return downloads


def summarize_jwst_archive(
    observations: list[MastObservation],
    products: list[MastProduct],
    retrieval_manifest: list[RetrievalManifestRow] | None = None,
    query_patterns: tuple[str, ...] = ("%PROXIMA%",),
    instrument_patterns: tuple[str, ...] = ("%MIRI%", "%NIRSPEC%"),
    program_type_filters: tuple[str, ...] = (),
    search_radius_degrees: float | None = 0.05,
) -> JwstArchiveSummary:
    instruments = tuple(sorted({row.instrument_name for row in observations if row.instrument_name}))
    filters = tuple(sorted({row.filters for row in observations if row.filters}))
    proposal_ids = tuple(sorted({row.proposal_id for row in observations if row.proposal_id}))
    proposal_types = tuple(sorted({row.proposal_type for row in observations if row.proposal_type}))
    science_products = [row for row in products if row.product_type == "SCIENCE"]
    public_products = [row for row in products if row.data_rights == "PUBLIC"]
    restricted_products = [
        row for row in products if row.data_rights not in {None, "", "PUBLIC"}
    ]
    retrievable_products = [row for row in products if row.data_uri]
    miri_observations = [
        row for row in observations if "miri" in row.instrument_name.lower()
    ]
    nirspec_observations = [
        row for row in observations if "nirspec" in row.instrument_name.lower()
    ]
    miri_products = [
        row
        for row in science_products
        if (row.obs_id and "miri" in row.obs_id.lower())
        or (row.product_filename and "miri" in row.product_filename.lower())
    ]
    nirspec_products = [
        row
        for row in science_products
        if (row.obs_id and "nirspec" in row.obs_id.lower())
        or (row.product_filename and "nirspec" in row.product_filename.lower())
    ]
    nirspec_or_niriss_observations = [
        row for row in observations if _is_nirspec_or_niriss(row.instrument_name)
    ]
    nirspec_or_niriss_products = [
        row
        for row in science_products
        if row.obs_id and ("nirspec" in row.obs_id.lower() or "niriss" in row.obs_id.lower())
    ]
    spectroscopy_observations = [
        row for row in observations if _looks_like_spectroscopy(row.instrument_name, row.dataproduct_type)
    ]
    imaging_observations = [
        row for row in observations if row.dataproduct_type == "image" or "IMAGE" in row.instrument_name
    ]
    retrieval_manifest = retrieval_manifest or build_retrieval_manifest(observations, products)
    phase_curve_candidates = [
        row for row in retrieval_manifest if row.analysis_role == "thermal_phase_curve_candidate"
    ]
    water_spectroscopy_candidates = [
        row for row in retrieval_manifest if row.analysis_role == "spectroscopy_water_screen"
    ]
    director_or_ers_observations = [
        row
        for row in observations
        if (row.proposal_type or "").upper() in DIRECTOR_OR_ERS_TYPES
        or "early release" in (row.obs_title or "").lower()
        or "director" in (row.obs_title or "").lower()
    ]
    planned_or_unproductized = [
        row
        for row in observations
        if (row.calib_level is not None and row.calib_level < 0)
        or not any(product.obsid == row.obsid for product in products)
    ]

    if not observations:
        impact = "No JWST Proxima observations were found by the configured target-name patterns."
        next_action = "Broaden aliases and coordinate searches, then check exo.MAST and proposal text."
    elif not science_products:
        impact = (
            "MAST returned JWST Proxima observation metadata, but no retrievable science "
            "products; this cannot move the habitability verdict."
        )
        next_action = "Monitor proposal IDs for product release, keep screening instrument-keyword tables, and rerun retrieval when products appear."
    elif not nirspec_or_niriss_products and not water_spectroscopy_candidates:
        impact = (
            "JWST products exist, but not NIRSpec/NIRISS spectroscopy products "
            "that can constrain Proxima b atmosphere in this mission."
        )
        next_action = "Use MIRI products for stellar/thermal context only and keep seeking spectroscopy or phase-curve data."
    elif water_spectroscopy_candidates and not nirspec_or_niriss_products:
        impact = "MIRI spectroscopy products are available for thermal and mid-infrared molecular-method development."
        next_action = "Download extracted spectral time-series products and test thermal/H2O/CO2/CH4 screening on analog data."
    else:
        impact = "NIRSpec/NIRISS science products are available for retrieval and atmospheric-screening analysis."
        next_action = "Download minimum recommended products and run spectrum extraction plus stellar contamination checks."

    return JwstArchiveSummary(
        query_patterns=query_patterns,
        instrument_patterns=instrument_patterns,
        program_type_filters=program_type_filters,
        search_radius_degrees=search_radius_degrees,
        observation_count=len(observations),
        public_observation_count=sum(1 for row in observations if row.data_rights == "PUBLIC"),
        restricted_observation_count=sum(
            1 for row in observations if row.data_rights not in {None, "", "PUBLIC"}
        ),
        science_product_count=len(science_products),
        public_product_count=len(public_products),
        restricted_product_count=len(restricted_products),
        retrievable_product_count=len(retrievable_products),
        miri_observation_count=len(miri_observations),
        nirspec_observation_count=len(nirspec_observations),
        miri_science_product_count=len(miri_products),
        nirspec_science_product_count=len(nirspec_products),
        nirspec_or_niriss_observation_count=len(nirspec_or_niriss_observations),
        nirspec_or_niriss_science_product_count=len(nirspec_or_niriss_products),
        spectroscopy_observation_count=len(spectroscopy_observations),
        imaging_observation_count=len(imaging_observations),
        phase_curve_candidate_count=len(phase_curve_candidates),
        water_spectroscopy_candidate_count=len(water_spectroscopy_candidates),
        director_or_ers_observation_count=len(director_or_ers_observations),
        planned_or_unproductized_count=len(planned_or_unproductized),
        instruments=instruments,
        filters=filters,
        proposal_ids=proposal_ids,
        proposal_types=proposal_types,
        verdict_impact=impact,
        next_action=next_action,
    )


def write_archive_outputs(
    observations: list[MastObservation],
    products: list[MastProduct],
    summary: JwstArchiveSummary,
    output_dir: str | Path,
    retrieval_manifest: list[RetrievalManifestRow] | None = None,
    keyword_screens: dict[str, list[dict[str, Any]]] | None = None,
    downloads: list[dict[str, Any]] | None = None,
) -> dict[str, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    observations_json = output_path / "jwst_archive_observations.json"
    observations_csv = output_path / "jwst_archive_observations.csv"
    products_json = output_path / "jwst_archive_products.json"
    products_csv = output_path / "jwst_archive_products.csv"
    summary_json = output_path / "jwst_archive_summary.json"
    retrieval_json = output_path / "jwst_retrieval_manifest.json"
    retrieval_csv = output_path / "jwst_retrieval_manifest.csv"
    download_queue = output_path / "jwst_download_queue.txt"
    keyword_json = output_path / "jwst_keyword_screens.json"
    downloads_json = output_path / "jwst_downloads.json"

    observation_rows = [asdict(row) for row in observations]
    product_rows = [asdict(row) for row in products]
    retrieval_manifest = retrieval_manifest or build_retrieval_manifest(observations, products)
    retrieval_rows = [asdict(row) for row in retrieval_manifest]
    observations_json.write_text(json.dumps(observation_rows, indent=2), encoding="utf-8")
    products_json.write_text(json.dumps(product_rows, indent=2), encoding="utf-8")
    summary_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    retrieval_json.write_text(json.dumps(retrieval_rows, indent=2), encoding="utf-8")
    download_queue.write_text(
        "\n".join(
            row.data_uri or ""
            for row in retrieval_manifest
            if row.download_recommended and row.data_uri
        )
        + ("\n" if any(row.download_recommended and row.data_uri for row in retrieval_manifest) else ""),
        encoding="utf-8",
    )
    _write_csv(observations_csv, observation_rows)
    _write_csv(products_csv, product_rows)
    _write_csv(retrieval_csv, retrieval_rows)

    outputs = {
        "archive_observations_json": observations_json,
        "archive_observations_csv": observations_csv,
        "archive_products_json": products_json,
        "archive_products_csv": products_csv,
        "archive_summary": summary_json,
        "retrieval_manifest_json": retrieval_json,
        "retrieval_manifest_csv": retrieval_csv,
        "download_queue": download_queue,
    }
    if keyword_screens is not None:
        keyword_json.write_text(json.dumps(keyword_screens, indent=2), encoding="utf-8")
        outputs["keyword_screens"] = keyword_json
    if downloads is not None:
        downloads_json.write_text(json.dumps(downloads, indent=2), encoding="utf-8")
        outputs["downloads"] = downloads_json
    return outputs


def _post_mast_request(request: dict[str, Any], timeout_seconds: int) -> dict[str, Any]:
    data = urllib.parse.urlencode({"request": json.dumps(request)}).encode("utf-8")
    http_request = urllib.request.Request(MAST_INVOKE, data=data, method="POST")
    with urllib.request.urlopen(http_request, timeout=timeout_seconds) as response:
        return json.loads(response.read().decode("utf-8"))


def _observation_from_row(row: dict[str, Any]) -> MastObservation:
    return MastObservation(
        obsid=int(row["obsid"]),
        obs_id=str(row.get("obs_id") or ""),
        obs_collection=str(row.get("obs_collection") or ""),
        provenance_name=row.get("provenance_name"),
        proposal_id=str(row.get("proposal_id") or ""),
        proposal_type=row.get("proposal_type"),
        proposal_pi=row.get("proposal_pi"),
        target_name=str(row.get("target_name") or ""),
        target_classification=row.get("target_classification"),
        instrument_name=str(row.get("instrument_name") or ""),
        dataproduct_type=str(row.get("dataproduct_type") or ""),
        filters=row.get("filters"),
        wave_region=row.get("wave_region"),
        wavelength_region=row.get("wavelength_region"),
        em_min=_optional_float(row.get("em_min")),
        em_max=_optional_float(row.get("em_max")),
        wave_min=_optional_float(row.get("wave_min")),
        wave_max=_optional_float(row.get("wave_max")),
        calib_level=_optional_int(row.get("calib_level")),
        t_min=_optional_float(row.get("t_min")),
        t_max=_optional_float(row.get("t_max")),
        t_exptime=_optional_float(row.get("t_exptime")),
        t_obs_release=_optional_float(row.get("t_obs_release")),
        data_rights=row.get("dataRights"),
        obs_title=row.get("obs_title"),
        data_url=row.get("dataURL"),
        s_ra=_optional_float(row.get("s_ra")),
        s_dec=_optional_float(row.get("s_dec")),
    )


def _product_from_row(row: dict[str, Any]) -> MastProduct:
    return MastProduct(
        obsid=_optional_int(row.get("obsID") or row.get("obsid")),
        obs_id=row.get("obs_id"),
        product_filename=row.get("productFilename"),
        product_type=row.get("productType"),
        product_group_description=row.get("productGroupDescription"),
        calib_level=_optional_int(row.get("calib_level")),
        size=_optional_int(row.get("size")),
        data_uri=row.get("dataURI"),
        data_rights=row.get("dataRights"),
    )


def _looks_like_spectroscopy(instrument_name: str, dataproduct_type: str) -> bool:
    text = f"{instrument_name} {dataproduct_type}".lower()
    return "spec" in text or "spectrum" in text or "x1d" in text


def _is_nirspec_or_niriss(instrument_name: str) -> bool:
    text = instrument_name.lower()
    return "nirspec" in text or "niriss" in text


def _analysis_role_for_observation(
    observation: MastObservation,
    product: MastProduct | None,
) -> tuple[str, str]:
    instrument = observation.instrument_name.lower()
    product_name = (product.product_filename if product else "") or ""
    filters = (observation.filters or "").lower()
    data_kind = f"{observation.dataproduct_type} {filters} {product_name}".lower()
    is_science_product = product is not None and product.product_type == "SCIENCE"

    if product is None:
        if "miri" in instrument and (observation.dataproduct_type == "image" or "image" in instrument):
            return "metadata_only_unretrievable", "thermal_emission;day_night_contrast;phase_curve_shape"
        if "nirspec" in instrument or "niriss" in instrument:
            return "metadata_only_unretrievable", "H2O_1.4um;H2O_1.9um;H2O_2.7um;CO2;CH4;CO"
        return "metadata_only_unretrievable", "none"

    if ("nirspec" in instrument or "niriss" in instrument) and is_science_product and (
        "spectrum" in data_kind
        or "timeseries" in data_kind
        or "prism" in data_kind
        or "g395" in data_kind
        or "x1d" in data_kind
        or "_rate.fits" in data_kind
    ):
        return "spectroscopy_water_screen", "H2O_1.4um;H2O_1.9um;H2O_2.7um;CO2;CH4;CO"
    if "miri" in instrument and is_science_product and (
        "spectrum" in data_kind or "mrs" in data_kind or "lrs" in data_kind or "x1d" in data_kind
    ):
        return "spectroscopy_water_screen", "H2O_6.3um;CO2_15um;CH4_7.7um;thermal_emission"
    if "miri" in instrument and (observation.dataproduct_type == "image" or "image" in instrument):
        role = "thermal_phase_curve_candidate" if is_science_product else "metadata_only_unretrievable"
        return role, "thermal_emission;day_night_contrast;phase_curve_shape"
    return "support_product", "stellar_context;calibration_context"


def _is_download_preferred_product(product: MastProduct) -> bool:
    filename = (product.product_filename or "").lower()
    if not filename.endswith(".fits"):
        return False
    if filename.endswith("_x1d.fits"):
        return True
    if filename.endswith("_x1dints.fits"):
        return True
    if filename.endswith("_rate.fits"):
        return True
    if filename.endswith("_cal.fits") and (product.size or 0) <= 250_000_000:
        return True
    if filename.endswith("_i2d.fits") and (product.size or 0) <= 250_000_000:
        return True
    return False


def _access_class(data_rights: str | None) -> str:
    if data_rights == "PUBLIC":
        return "public"
    if data_rights in {None, ""}:
        return "unknown"
    return "restricted_or_auth_required"


def _orbital_phase(
    mjd: float | None,
    orbital_period_days: float,
    reference_epoch_mjd: float,
) -> float | None:
    if mjd is None or orbital_period_days <= 0:
        return None
    return round(((mjd - reference_epoch_mjd) % orbital_period_days) / orbital_period_days, 6)


def _safe_filename(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)


def _optional_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def _optional_int(value: Any) -> int | None:
    if value in {None, ""}:
        return None
    return int(value)


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
