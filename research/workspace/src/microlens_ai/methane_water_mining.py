from __future__ import annotations

import csv
import json
import math
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .config import _load_yaml


EXOPLANET_TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

DEFAULT_TOI_LIMIT = 1500
DEFAULT_CONFIRMED_LIMIT = 800

TOI_COLUMNS = (
    "toi",
    "tid",
    "ctoi_alias",
    "tfopwg_disp",
    "ra",
    "dec",
    "pl_orbper",
    "pl_tranmid",
    "pl_trandurh",
    "pl_trandep",
    "pl_rade",
    "pl_insol",
    "pl_eqt",
    "st_tmag",
    "st_dist",
    "st_teff",
    "st_logg",
    "st_rad",
    "toi_created",
    "rowupdate",
)

CONFIRMED_COLUMNS = (
    "pl_name",
    "hostname",
    "discoverymethod",
    "disc_year",
    "disc_facility",
    "ra",
    "dec",
    "pl_orbper",
    "pl_rade",
    "pl_bmasse",
    "pl_eqt",
    "pl_trandep",
    "st_teff",
    "st_rad",
    "sy_tmag",
    "sy_jmag",
    "sy_dist",
    "pl_ntranspec",
    "pl_nespec",
    "pl_ndispec",
    "st_nspec",
)


@dataclass(frozen=True)
class MethaneWaterConfig:
    """Configuration for the Milky Way methane/water mining mission."""

    project_name: str
    science_goal: str
    working_candidate_name: str
    output_tag: str
    archive_sources: dict[str, str]
    max_toi_rows: int = DEFAULT_TOI_LIMIT
    max_confirmed_rows: int = DEFAULT_CONFIRMED_LIMIT
    shortlist_size: int = 25
    top_mast_targets: int = 5
    raw: dict[str, Any] | None = None

    @classmethod
    def from_yaml(cls, path: str | Path) -> "MethaneWaterConfig":
        raw = _load_yaml(Path(path))
        project = raw["project"]
        ranking = raw.get("ranking", {})
        return cls(
            project_name=project["name"],
            science_goal=project["science_goal"],
            working_candidate_name=project.get("working_candidate_name", "Haven-1"),
            output_tag=project.get("output_tag", "mission4_methane_water"),
            archive_sources=dict(raw.get("archive_sources", {})),
            max_toi_rows=int(ranking.get("max_toi_rows", DEFAULT_TOI_LIMIT)),
            max_confirmed_rows=int(
                ranking.get("max_confirmed_rows", DEFAULT_CONFIRMED_LIMIT)
            ),
            shortlist_size=int(ranking.get("shortlist_size", 25)),
            top_mast_targets=int(ranking.get("top_mast_targets", 5)),
            raw=raw,
        )


@dataclass(frozen=True)
class ToiCandidate:
    """Ranked TESS Object of Interest candidate from the public archive."""

    rank: int
    working_name: str | None
    toi: str
    tic_id: str
    ctoi_alias: str | None
    disposition: str
    ra: float | None
    dec: float | None
    orbital_period_days: float | None
    transit_midpoint_bjd: float | None
    transit_duration_hours: float | None
    transit_depth_ppm: float | None
    planet_radius_earth: float | None
    equilibrium_temperature_k: float | None
    insolation_earth: float | None
    stellar_tmag: float | None
    stellar_distance_pc: float | None
    stellar_teff_k: float | None
    stellar_logg: float | None
    stellar_radius_solar: float | None
    toi_created: str | None
    rowupdate: str | None
    methane_priority: float
    water_priority: float
    discovery_priority: float
    followup_priority: float
    total_score: float
    signal_proxy_ppm: float | None
    flags: tuple[str, ...]
    blockers: tuple[str, ...]
    exofop_url: str
    mast_portal_url: str


@dataclass(frozen=True)
class ConfirmedAtmosphereTarget:
    """Known transiting planet with archive spectroscopy mining value."""

    rank: int
    planet_name: str
    host_name: str
    discovery_method: str | None
    discovery_year: int | None
    discovery_facility: str | None
    ra: float | None
    dec: float | None
    orbital_period_days: float | None
    planet_radius_earth: float | None
    planet_mass_earth: float | None
    equilibrium_temperature_k: float | None
    transit_depth_percent: float | None
    stellar_teff_k: float | None
    stellar_radius_solar: float | None
    tess_mag: float | None
    j_mag: float | None
    distance_pc: float | None
    transmission_spectra: int
    eclipse_spectra: int
    direct_imaging_spectra: int
    stellar_spectra: int
    methane_priority: float
    water_priority: float
    archive_priority: float
    total_score: float
    molecular_targets: str
    archive_url: str


@dataclass(frozen=True)
class TessArchiveObservation:
    """Compact MAST observation row for a top TIC target."""

    tic_id: str
    toi: str
    obsid: int | None
    obs_id: str | None
    target_name: str | None
    dataproduct_type: str | None
    instrument_name: str | None
    calib_level: int | None
    data_rights: str | None
    t_min: float | None
    t_max: float | None
    exptime_seconds: float | None
    product_count: int
    lightcurve_product_count: int
    targetpixel_product_count: int
    recommended_downloads: int


@dataclass(frozen=True)
class MiningPackage:
    """Complete mission output payload."""

    generated_at_utc: str
    project_name: str
    science_goal: str
    working_candidate_name: str
    query_status: dict[str, str]
    candidate_count: int
    confirmed_atmosphere_count: int
    tess_observation_count: int
    top_candidates: list[ToiCandidate]
    confirmed_atmosphere_targets: list[ConfirmedAtmosphereTarget]
    tess_observations: list[TessArchiveObservation]
    archive_sources: dict[str, str]


def tap_query_csv(sql: str, timeout_seconds: int = 90) -> list[dict[str, str]]:
    """Run a NASA Exoplanet Archive TAP query and return CSV rows."""

    params = urllib.parse.urlencode({"query": sql, "format": "csv"})
    url = f"{EXOPLANET_TAP_SYNC}?{params}"
    with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
        text = response.read().decode("utf-8")
    if text.lstrip().lower().startswith("<"):
        raise RuntimeError(f"TAP returned non-CSV response for query: {sql[:120]}")
    return list(csv.DictReader(text.splitlines()))


def build_toi_sql(limit: int = DEFAULT_TOI_LIMIT) -> str:
    columns = ",".join(TOI_COLUMNS)
    return (
        f"select top {int(limit)} {columns} from toi "
        "where tfopwg_disp in ('PC','APC') "
        "and pl_rade is not null "
        "and st_tmag is not null "
        "and st_rad is not null "
        "order by rowupdate desc"
    )


def build_confirmed_sql(limit: int = DEFAULT_CONFIRMED_LIMIT) -> str:
    columns = ",".join(CONFIRMED_COLUMNS)
    return (
        f"select top {int(limit)} {columns} from pscomppars "
        "where tran_flag=1 "
        "and pl_rade is not null "
        "and st_rad is not null "
        "and (pl_ntranspec > 0 or pl_nespec > 0 or pl_ndispec > 0 or st_nspec > 0) "
        "order by pl_ntranspec desc"
    )


def rank_toi_rows(
    rows: list[dict[str, Any]],
    working_candidate_name: str = "Haven-1",
    shortlist_size: int = 25,
) -> list[ToiCandidate]:
    ranked = [_candidate_from_toi_row(row, working_candidate_name=None) for row in rows]
    ranked = sorted(ranked, key=lambda row: row.total_score, reverse=True)
    output: list[ToiCandidate] = []
    for index, row in enumerate(ranked[:shortlist_size], start=1):
        output.append(
            ToiCandidate(
                **{
                    **asdict(row),
                    "rank": index,
                    "working_name": working_candidate_name if index == 1 else None,
                }
            )
        )
    return output


def rank_confirmed_rows(
    rows: list[dict[str, Any]],
    shortlist_size: int = 25,
) -> list[ConfirmedAtmosphereTarget]:
    ranked = [_confirmed_from_row(row) for row in rows]
    ranked = sorted(ranked, key=lambda row: row.total_score, reverse=True)
    output: list[ConfirmedAtmosphereTarget] = []
    for index, row in enumerate(ranked[:shortlist_size], start=1):
        output.append(ConfirmedAtmosphereTarget(**{**asdict(row), "rank": index}))
    return output


def build_mining_package(
    config: MethaneWaterConfig,
    live_archive: bool = False,
    include_mast_tess: bool = False,
    timeout_seconds: int = 90,
) -> MiningPackage:
    """Build a mission package from live archive rows or an offline skeleton."""

    query_status: dict[str, str] = {}
    toi_rows: list[dict[str, str]] = []
    confirmed_rows: list[dict[str, str]] = []
    tess_observations: list[TessArchiveObservation] = []

    if live_archive:
        toi_sql = build_toi_sql(config.max_toi_rows)
        confirmed_sql = build_confirmed_sql(config.max_confirmed_rows)
        try:
            toi_rows = tap_query_csv(toi_sql, timeout_seconds=timeout_seconds)
            query_status["toi"] = f"ok rows={len(toi_rows)}"
        except Exception as exc:  # pragma: no cover - live archive failure path
            query_status["toi"] = f"error {exc}"
        try:
            confirmed_rows = tap_query_csv(confirmed_sql, timeout_seconds=timeout_seconds)
            query_status["confirmed_atmospheres"] = f"ok rows={len(confirmed_rows)}"
        except Exception as exc:  # pragma: no cover - live archive failure path
            query_status["confirmed_atmospheres"] = f"error {exc}"
    else:
        query_status["toi"] = "offline mode: no TAP query executed"
        query_status["confirmed_atmospheres"] = "offline mode: no TAP query executed"

    top_candidates = rank_toi_rows(
        toi_rows,
        working_candidate_name=config.working_candidate_name,
        shortlist_size=config.shortlist_size,
    )
    confirmed_targets = rank_confirmed_rows(
        confirmed_rows,
        shortlist_size=config.shortlist_size,
    )

    if live_archive and include_mast_tess and top_candidates:
        try:
            tess_observations = query_tess_observations_for_candidates(
                top_candidates[: config.top_mast_targets],
                timeout_seconds=timeout_seconds,
            )
            query_status["mast_tess"] = f"ok observations={len(tess_observations)}"
        except Exception as exc:  # pragma: no cover - live MAST failure path
            query_status["mast_tess"] = f"error {exc}"
    elif include_mast_tess:
        query_status["mast_tess"] = "skipped: no top candidates available"
    else:
        query_status["mast_tess"] = "skipped: use --mast-tess to query MAST products"

    return MiningPackage(
        generated_at_utc=datetime.now(UTC).replace(microsecond=0).isoformat(),
        project_name=config.project_name,
        science_goal=config.science_goal,
        working_candidate_name=config.working_candidate_name,
        query_status=query_status,
        candidate_count=len(top_candidates),
        confirmed_atmosphere_count=len(confirmed_targets),
        tess_observation_count=len(tess_observations),
        top_candidates=top_candidates,
        confirmed_atmosphere_targets=confirmed_targets,
        tess_observations=tess_observations,
        archive_sources=config.archive_sources,
    )


def query_tess_observations_for_candidates(
    candidates: list[ToiCandidate],
    timeout_seconds: int = 90,
) -> list[TessArchiveObservation]:
    """Query MAST for TESS observation/product metadata for the top TICs."""

    from .mast_archive import mast_request

    observations: list[TessArchiveObservation] = []
    columns = (
        "obsid,obs_id,target_name,dataproduct_type,instrument_name,calib_level,"
        "dataRights,t_min,t_max,t_exptime"
    )
    for candidate in candidates:
        payload = mast_request(
            "Mast.Caom.Filtered",
            {
                "columns": columns,
                "filters": [
                    {"paramName": "obs_collection", "values": ["TESS"]},
                    {
                        "paramName": "target_name",
                        "values": [],
                        "freeText": f"%{candidate.tic_id}%",
                    },
                ],
            },
            pagesize=100,
            timeout_seconds=timeout_seconds,
        )
        for row in payload.get("data", []):
            obsid = _optional_int(row.get("obsid"))
            product_counts = _mast_product_counts(obsid, timeout_seconds=timeout_seconds)
            observations.append(
                TessArchiveObservation(
                    tic_id=candidate.tic_id,
                    toi=candidate.toi,
                    obsid=obsid,
                    obs_id=_optional_str(row.get("obs_id")),
                    target_name=_optional_str(row.get("target_name")),
                    dataproduct_type=_optional_str(row.get("dataproduct_type")),
                    instrument_name=_optional_str(row.get("instrument_name")),
                    calib_level=_optional_int(row.get("calib_level")),
                    data_rights=_optional_str(row.get("dataRights")),
                    t_min=_optional_float(row.get("t_min")),
                    t_max=_optional_float(row.get("t_max")),
                    exptime_seconds=_optional_float(row.get("t_exptime")),
                    product_count=product_counts["products"],
                    lightcurve_product_count=product_counts["lightcurves"],
                    targetpixel_product_count=product_counts["target_pixels"],
                    recommended_downloads=product_counts["recommended"],
                )
            )
    return observations


def write_mining_outputs(package: MiningPackage, output_dir: str | Path) -> dict[str, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    package_path = output_path / "mission_package.json"
    candidates_csv = output_path / "toi_candidate_shortlist.csv"
    confirmed_csv = output_path / "confirmed_atmosphere_targets.csv"
    tess_csv = output_path / "tess_archive_observations.csv"
    report_path = output_path / "discovery_shortlist.md"
    working_name_path = output_path / "working_name_registry.json"

    payload = _package_to_jsonable(package)
    package_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    _write_csv(candidates_csv, [asdict(row) for row in package.top_candidates])
    _write_csv(
        confirmed_csv,
        [asdict(row) for row in package.confirmed_atmosphere_targets],
    )
    _write_csv(tess_csv, [asdict(row) for row in package.tess_observations])
    report_path.write_text(render_discovery_report(package), encoding="utf-8")
    working_name_path.write_text(
        json.dumps(
            {
                "working_candidate_name": package.working_candidate_name,
                "status": "provisional codename, not an official IAU exoplanet name",
                "assigned_to": (
                    package.top_candidates[0].toi if package.top_candidates else None
                ),
                "generated_at_utc": package.generated_at_utc,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return {
        "package": package_path,
        "candidates": candidates_csv,
        "confirmed_atmospheres": confirmed_csv,
        "tess_observations": tess_csv,
        "report": report_path,
        "working_name_registry": working_name_path,
    }


def render_discovery_report(package: MiningPackage) -> str:
    lines = [
        "# Mission 4: Milky Way Methane and Water Planet Hunt",
        "",
        f"Generated: {package.generated_at_utc}",
        f"Science goal: {package.science_goal}",
        f"Working codename for the current top candidate: {package.working_candidate_name}",
        "",
        "## Reality Gate",
        "",
        "- This run does not claim a new confirmed planet.",
        "- The discovery track starts with unconfirmed TESS candidates and public TESS light-curve products.",
        "- The molecule track starts with known transiting planets that already have public spectroscopy counts.",
        "- Methane and water priorities are triage scores; they identify targets worth retrieval analysis.",
        "",
        "## Query Status",
        "",
    ]
    for name, status in package.query_status.items():
        lines.append(f"- {name}: {status}")

    lines.extend(
        [
            "",
            "## Top Discovery Candidates",
            "",
            "| Rank | Work name | TOI | TIC | Disp. | R_Earth | Teq K | Tmag | Score | Flags | Blockers |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in package.top_candidates:
        flags = "; ".join(row.flags) if row.flags else "none"
        blockers = "; ".join(row.blockers) if row.blockers else "none"
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.rank),
                    row.working_name or "",
                    row.toi,
                    row.tic_id,
                    row.disposition,
                    _fmt(row.planet_radius_earth),
                    _fmt(row.equilibrium_temperature_k),
                    _fmt(row.stellar_tmag),
                    f"{row.total_score:.3f}",
                    flags,
                    blockers,
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Known Atmosphere Mining Targets",
            "",
            "| Rank | Planet | Host | Teq K | Spectra | Score | Molecular targets |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in package.confirmed_atmosphere_targets:
        spectra = (
            f"trans={row.transmission_spectra}, eclipse={row.eclipse_spectra}, "
            f"direct={row.direct_imaging_spectra}, stellar={row.stellar_spectra}"
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.rank),
                    row.planet_name,
                    row.host_name,
                    _fmt(row.equilibrium_temperature_k),
                    spectra,
                    f"{row.total_score:.3f}",
                    row.molecular_targets,
                ]
            )
            + " |"
        )

    if package.tess_observations:
        lines.extend(
            [
                "",
                "## TESS Product Recon",
                "",
                "| TIC | TOI | ObsID | Type | Products | Light curves | Target pixels | Recommended |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in package.tess_observations:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row.tic_id,
                        row.toi,
                        str(row.obsid or ""),
                        row.dataproduct_type or "",
                        str(row.product_count),
                        str(row.lightcurve_product_count),
                        str(row.targetpixel_product_count),
                        str(row.recommended_downloads),
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Next Execution Steps",
            "",
            "1. Download the top TIC light curves and target-pixel files from MAST.",
            "2. Run a transit refit and odd/even, secondary-eclipse, centroid, and contamination vetting.",
            "3. Cross-check the top signal against TOI, ExoFOP, Gaia, nearby-source, and known-variable records.",
            "4. For any survivor, build a validation packet before treating it as a planet candidate.",
            "5. For molecule mining, retrieve spectra and run H2O/CH4 model comparison against stellar contamination controls.",
            "",
            "## Source Routes",
            "",
        ]
    )
    for name, url in package.archive_sources.items():
        lines.append(f"- {name}: {url}")
    lines.append("")
    return "\n".join(lines)


def _candidate_from_toi_row(
    row: dict[str, Any],
    working_candidate_name: str | None,
) -> ToiCandidate:
    disposition = str(row.get("tfopwg_disp") or "").strip()
    radius = _optional_float(row.get("pl_rade"))
    teq = _optional_float(row.get("pl_eqt"))
    tmag = _optional_float(row.get("st_tmag"))
    stellar_radius = _optional_float(row.get("st_rad"))
    stellar_teff = _optional_float(row.get("st_teff"))
    depth_ppm = _optional_float(row.get("pl_trandep"))

    methane = _methane_priority(teq, radius)
    water = _water_priority(teq, radius)
    discovery = _discovery_priority(disposition, radius, tmag, stellar_radius, depth_ppm)
    followup = _followup_priority(tmag, stellar_radius, depth_ppm, stellar_teff)
    signal_proxy = _signal_proxy_ppm(radius, stellar_radius, teq)
    total = (
        0.34 * discovery
        + 0.28 * followup
        + 0.20 * max(methane, water)
        + 0.10 * methane
        + 0.08 * water
    )

    flags, blockers = _candidate_flags_and_blockers(
        disposition=disposition,
        radius=radius,
        teq=teq,
        tmag=tmag,
        stellar_radius=stellar_radius,
        stellar_teff=stellar_teff,
        depth_ppm=depth_ppm,
    )
    tic_id = str(row.get("tid") or "").strip()
    toi = str(row.get("toi") or "").strip()
    return ToiCandidate(
        rank=0,
        working_name=working_candidate_name,
        toi=toi,
        tic_id=tic_id,
        ctoi_alias=_optional_str(row.get("ctoi_alias")),
        disposition=disposition,
        ra=_optional_float(row.get("ra")),
        dec=_optional_float(row.get("dec")),
        orbital_period_days=_optional_float(row.get("pl_orbper")),
        transit_midpoint_bjd=_optional_float(row.get("pl_tranmid")),
        transit_duration_hours=_optional_float(row.get("pl_trandurh")),
        transit_depth_ppm=depth_ppm,
        planet_radius_earth=radius,
        equilibrium_temperature_k=teq,
        insolation_earth=_optional_float(row.get("pl_insol")),
        stellar_tmag=tmag,
        stellar_distance_pc=_optional_float(row.get("st_dist")),
        stellar_teff_k=stellar_teff,
        stellar_logg=_optional_float(row.get("st_logg")),
        stellar_radius_solar=stellar_radius,
        toi_created=_optional_str(row.get("toi_created")),
        rowupdate=_optional_str(row.get("rowupdate")),
        methane_priority=round(methane, 6),
        water_priority=round(water, 6),
        discovery_priority=round(discovery, 6),
        followup_priority=round(followup, 6),
        total_score=round(total, 6),
        signal_proxy_ppm=None if signal_proxy is None else round(signal_proxy, 3),
        flags=flags,
        blockers=blockers,
        exofop_url=f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic_id}",
        mast_portal_url=f"https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery={urllib.parse.quote('TIC ' + tic_id)}",
    )


def _confirmed_from_row(row: dict[str, Any]) -> ConfirmedAtmosphereTarget:
    radius = _optional_float(row.get("pl_rade"))
    teq = _optional_float(row.get("pl_eqt"))
    tmag = _optional_float(row.get("sy_tmag"))
    jmag = _optional_float(row.get("sy_jmag"))
    stellar_radius = _optional_float(row.get("st_rad"))
    depth_percent = _optional_float(row.get("pl_trandep"))
    ntrans = _optional_int(row.get("pl_ntranspec")) or 0
    nespec = _optional_int(row.get("pl_nespec")) or 0
    ndispec = _optional_int(row.get("pl_ndispec")) or 0
    nspec = _optional_int(row.get("st_nspec")) or 0

    methane = _methane_priority(teq, radius)
    water = _water_priority(teq, radius)
    archive = _archive_priority(ntrans, nespec, ndispec, nspec)
    brightness = _brightness_score(jmag if jmag is not None else tmag, bright=7.0, faint=14.0)
    star = _small_star_score(stellar_radius)
    total = 0.36 * archive + 0.22 * max(methane, water) + 0.18 * brightness + 0.14 * star + 0.10 * methane

    molecular_targets = _molecular_targets(teq, methane, water)
    planet_name = str(row.get("pl_name") or "").strip()
    return ConfirmedAtmosphereTarget(
        rank=0,
        planet_name=planet_name,
        host_name=str(row.get("hostname") or "").strip(),
        discovery_method=_optional_str(row.get("discoverymethod")),
        discovery_year=_optional_int(row.get("disc_year")),
        discovery_facility=_optional_str(row.get("disc_facility")),
        ra=_optional_float(row.get("ra")),
        dec=_optional_float(row.get("dec")),
        orbital_period_days=_optional_float(row.get("pl_orbper")),
        planet_radius_earth=radius,
        planet_mass_earth=_optional_float(row.get("pl_bmasse")),
        equilibrium_temperature_k=teq,
        transit_depth_percent=depth_percent,
        stellar_teff_k=_optional_float(row.get("st_teff")),
        stellar_radius_solar=stellar_radius,
        tess_mag=tmag,
        j_mag=jmag,
        distance_pc=_optional_float(row.get("sy_dist")),
        transmission_spectra=ntrans,
        eclipse_spectra=nespec,
        direct_imaging_spectra=ndispec,
        stellar_spectra=nspec,
        methane_priority=round(methane, 6),
        water_priority=round(water, 6),
        archive_priority=round(archive, 6),
        total_score=round(total, 6),
        molecular_targets=molecular_targets,
        archive_url=f"https://exoplanetarchive.ipac.caltech.edu/overview/{urllib.parse.quote(planet_name)}",
    )


def _mast_product_counts(obsid: int | None, timeout_seconds: int) -> dict[str, int]:
    if obsid is None:
        return {"products": 0, "lightcurves": 0, "target_pixels": 0, "recommended": 0}
    from .mast_archive import mast_request

    payload = mast_request(
        "Mast.Caom.Products",
        {"obsid": obsid},
        pagesize=1000,
        timeout_seconds=timeout_seconds,
    )
    rows = payload.get("data", [])
    lightcurves = 0
    target_pixels = 0
    recommended = 0
    for row in rows:
        filename = str(row.get("productFilename") or "").lower()
        group = str(row.get("productSubGroupDescription") or "").lower()
        product_type = str(row.get("productType") or "").lower()
        if "lightcurve" in product_type or filename.endswith("_lc.fits") or group == "lc":
            lightcurves += 1
        if "targetpixel" in product_type or filename.endswith("_tp.fits") or group == "tp":
            target_pixels += 1
        if row.get("dataURI") and row.get("dataRights") == "PUBLIC" and (
            filename.endswith("_lc.fits") or filename.endswith("_tp.fits")
        ):
            recommended += 1
    return {
        "products": len(rows),
        "lightcurves": lightcurves,
        "target_pixels": target_pixels,
        "recommended": recommended,
    }


def _methane_priority(teq: float | None, radius_earth: float | None) -> float:
    if teq is None:
        temp_score = 0.35
    elif 250 <= teq <= 700:
        temp_score = 1.0
    elif 700 < teq <= 1000:
        temp_score = 0.72
    elif 150 <= teq < 250:
        temp_score = 0.62
    elif 1000 < teq <= 1400:
        temp_score = 0.32
    else:
        temp_score = 0.12

    if radius_earth is None:
        radius_score = 0.55
    elif 1.0 <= radius_earth <= 4.0:
        radius_score = 1.0
    elif 4.0 < radius_earth <= 8.0:
        radius_score = 0.58
    elif radius_earth < 1.0:
        radius_score = 0.62
    else:
        radius_score = 0.28
    return _clamp01(0.72 * temp_score + 0.28 * radius_score)


def _water_priority(teq: float | None, radius_earth: float | None) -> float:
    if teq is None:
        temp_score = 0.45
    elif 180 <= teq <= 450:
        temp_score = 1.0
    elif 450 < teq <= 900:
        temp_score = 0.78
    elif 900 < teq <= 1400:
        temp_score = 0.52
    elif 120 <= teq < 180:
        temp_score = 0.50
    else:
        temp_score = 0.22

    if radius_earth is None:
        radius_score = 0.55
    elif 0.7 <= radius_earth <= 3.5:
        radius_score = 1.0
    elif 3.5 < radius_earth <= 8.0:
        radius_score = 0.67
    elif radius_earth < 0.7:
        radius_score = 0.45
    else:
        radius_score = 0.35
    return _clamp01(0.68 * temp_score + 0.32 * radius_score)


def _discovery_priority(
    disposition: str,
    radius_earth: float | None,
    tmag: float | None,
    stellar_radius: float | None,
    depth_ppm: float | None,
) -> float:
    disposition_score = {"PC": 1.0, "APC": 0.72}.get(disposition.upper(), 0.35)
    radius_score = _candidate_radius_score(radius_earth)
    brightness = _brightness_score(tmag, bright=8.0, faint=14.5)
    star = _small_star_score(stellar_radius)
    depth = _depth_score(depth_ppm)
    return _clamp01(
        0.36 * disposition_score
        + 0.18 * radius_score
        + 0.18 * brightness
        + 0.16 * star
        + 0.12 * depth
    )


def _followup_priority(
    tmag: float | None,
    stellar_radius: float | None,
    depth_ppm: float | None,
    stellar_teff: float | None,
) -> float:
    brightness = _brightness_score(tmag, bright=7.5, faint=14.0)
    star = _small_star_score(stellar_radius)
    depth = _depth_score(depth_ppm)
    cool_star = _cool_star_score(stellar_teff)
    return _clamp01(0.34 * brightness + 0.27 * star + 0.22 * depth + 0.17 * cool_star)


def _candidate_flags_and_blockers(
    disposition: str,
    radius: float | None,
    teq: float | None,
    tmag: float | None,
    stellar_radius: float | None,
    stellar_teff: float | None,
    depth_ppm: float | None,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    flags: list[str] = []
    blockers: list[str] = []
    if disposition.upper() == "PC":
        flags.append("planetary_candidate")
    elif disposition.upper() == "APC":
        flags.append("ambiguous_candidate")
        blockers.append("ambiguous TFOPWG disposition")
    else:
        blockers.append("not a PC/APC candidate")

    if radius is not None and radius <= 4.0:
        flags.append("rocky_or_sub_neptune_size")
    elif radius is not None and radius > 8.0:
        blockers.append("giant radius weakens first-small-planet objective")

    if teq is not None and 250 <= teq <= 1000:
        flags.append("methane_temperature_window")
    if teq is not None and 180 <= teq <= 900:
        flags.append("water_temperature_window")
    if stellar_teff is not None and stellar_teff <= 4200:
        flags.append("cool_host_star")
    if stellar_radius is not None and stellar_radius <= 0.75:
        flags.append("small_star_deeper_transit")
    if tmag is not None and tmag <= 12.5:
        flags.append("bright_tess_host")
    if depth_ppm is not None and depth_ppm >= 600:
        flags.append("deep_tess_transit")

    for label, value in (
        ("planet radius", radius),
        ("equilibrium temperature", teq),
        ("TESS magnitude", tmag),
        ("stellar radius", stellar_radius),
        ("transit depth", depth_ppm),
    ):
        if value is None:
            blockers.append(f"missing {label}")
    return tuple(flags), tuple(blockers)


def _molecular_targets(teq: float | None, methane: float, water: float) -> str:
    targets: list[str] = []
    if water >= 0.55:
        targets.extend(["H2O_1.4um", "H2O_1.9um", "H2O_2.7um"])
    if methane >= 0.55:
        targets.extend(["CH4_1.6um", "CH4_2.3um", "CH4_3.3um"])
    if teq is not None and teq >= 700:
        targets.extend(["CO2", "CO"])
    if not targets:
        targets.append("stellar_contamination_control")
    return ";".join(targets)


def _archive_priority(ntrans: int, nespec: int, ndispec: int, nspec: int) -> float:
    value = (1.4 * ntrans) + (1.0 * nespec) + (0.8 * ndispec) + (0.25 * nspec)
    return _clamp01(math.log1p(value) / math.log(18.0))


def _candidate_radius_score(radius_earth: float | None) -> float:
    if radius_earth is None:
        return 0.45
    if 0.8 <= radius_earth <= 2.2:
        return 1.0
    if 2.2 < radius_earth <= 4.0:
        return 0.86
    if 4.0 < radius_earth <= 8.0:
        return 0.48
    if radius_earth < 0.8:
        return 0.52
    return 0.18


def _brightness_score(magnitude: float | None, bright: float, faint: float) -> float:
    if magnitude is None:
        return 0.35
    if magnitude <= bright:
        return 1.0
    if magnitude >= faint:
        return 0.0
    return _clamp01((faint - magnitude) / (faint - bright))


def _small_star_score(stellar_radius: float | None) -> float:
    if stellar_radius is None:
        return 0.4
    if stellar_radius <= 0.35:
        return 1.0
    if stellar_radius <= 0.75:
        return 0.82
    if stellar_radius <= 1.1:
        return 0.45
    if stellar_radius <= 1.6:
        return 0.20
    return 0.05


def _cool_star_score(stellar_teff: float | None) -> float:
    if stellar_teff is None:
        return 0.4
    if stellar_teff <= 3900:
        return 1.0
    if stellar_teff <= 4700:
        return 0.78
    if stellar_teff <= 5600:
        return 0.38
    return 0.14


def _depth_score(depth_ppm: float | None) -> float:
    if depth_ppm is None:
        return 0.35
    return _clamp01(math.log1p(max(depth_ppm, 0.0)) / math.log(6000.0))


def _signal_proxy_ppm(
    radius_earth: float | None,
    stellar_radius_solar: float | None,
    teq: float | None,
) -> float | None:
    if radius_earth is None or stellar_radius_solar is None or stellar_radius_solar <= 0:
        return None
    temp_factor = 0.8 if teq is None else _clamp01((teq - 120.0) / 900.0)
    radius_factor = min(max(radius_earth, 0.5), 8.0) ** 0.7
    star_factor = 1.0 / max(stellar_radius_solar, 0.12) ** 2
    return 12.0 * temp_factor * radius_factor * star_factor


def _package_to_jsonable(package: MiningPackage) -> dict[str, Any]:
    return {
        **asdict(package),
        "top_candidates": [asdict(row) for row in package.top_candidates],
        "confirmed_atmosphere_targets": [
            asdict(row) for row in package.confirmed_atmosphere_targets
        ],
        "tess_observations": [asdict(row) for row in package.tess_observations],
    }


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


def _fmt(value: float | int | None) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def _optional_str(value: Any) -> str | None:
    if value in {None, ""}:
        return None
    return str(value)


def _optional_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def _optional_int(value: Any) -> int | None:
    if value in {None, ""}:
        return None
    return int(float(value))


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))
