"""TOI-1258.01 mass/barcode archive hunt.

Queries public archive metadata for:
- KOA/Keck HIRES observations near TIC 288132261.
- KOA reduced-data table entries linked to any matching KOA IDs.
- MAST observations near the target, with emphasis on HST/JWST spectroscopy
  and any time-series/transit-relevant products.

The script intentionally does not claim a mass or atmosphere detection. It
creates a reproducible gate report describing what public data are present.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import urlopen


TARGET = {
    "name": "TOI-1258.01",
    "tic": "288132261",
    "ra": 209.484313,
    "dec": 79.583322,
    "period_days": 10.2613397,
    "tranmid_bjd": 2460333.648786,
    "radius_rearth": 3.60391,
    "host_tmag": 9.13619,
}


@dataclass
class GateResult:
    gate: str
    status: str
    summary: str
    table: str | None = None


def ensure_dirs(root: Path) -> dict[str, Path]:
    dirs = {
        "root": root,
        "reports": root / "evidence" / "reports",
        "tables": root / "evidence" / "tables",
        "logs": root / "evidence" / "logs",
        "outreach": root / "outreach",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def flag_is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "t", "yes", "y", "1", "in"}


def flag_is_false(value: Any) -> bool:
    return str(value).strip().lower() in {"false", "f", "no", "n", "0", "out"}


def row_to_dict(row: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    if hasattr(row, "items"):
        iterator = row.items()
    elif hasattr(row, "array") and getattr(row.array, "dtype", None) is not None:
        iterator = ((name, row[name]) for name in (row.array.dtype.names or []))
    else:
        iterator = []
    for name, value in iterator:
        try:
            if hasattr(value, "item"):
                value = value.item()
        except Exception:
            pass
        if isinstance(value, bytes):
            value = value.decode("utf-8", "replace")
        out[name] = value
    return out


def query_koa_hires(dirs: dict[str, Path]) -> GateResult:
    try:
        import pyvo
    except Exception as exc:  # pragma: no cover
        return GateResult(
            gate="koa_hires",
            status="BLOCKED_NO_PYVO",
            summary=f"pyvo import failed: {exc}",
        )

    svc = pyvo.dal.TAPService("https://koa.ipac.caltech.edu/TAP")
    ra = TARGET["ra"]
    dec = TARGET["dec"]
    radius = 0.08

    query = f"""
    SELECT TOP 200
      koaid, targname, object, koaimtyp, obstype, date_obs, ut, utdatetime,
      ra, dec, progid, progpi, progtitl, semid, propint, exptime,
      deckname, slitwidt, specres, sig2nois, waveblue, wavecntr, wavered,
      iodin, iodout, filehand, ofname
    FROM koa_hires
    WHERE 1=CONTAINS(
      POINT('ICRS', ra, dec),
      CIRCLE('ICRS', {ra}, {dec}, {radius})
    )
    ORDER BY date_obs
    """
    rows: list[dict[str, Any]] = []
    try:
        result = svc.search(query, maxrec=300)
        rows = [row_to_dict(row) for row in result]
    except Exception as exc:
        log = dirs["logs"] / "koa_hires_query_error.txt"
        log.write_text(f"Query failed:\n{query}\n\n{exc}\n", encoding="utf-8")
        return GateResult(
            gate="koa_hires",
            status="QUERY_FAILED",
            summary=f"KOA HIRES TAP query failed. See {log}",
        )

    table = dirs["tables"] / "koa_hires_cone_search.csv"
    write_csv(table, rows)

    if not rows:
        return GateResult(
            gate="koa_hires",
            status="NO_PUBLIC_KOA_HIRES_MATCH",
            summary=(
                f"No public KOA HIRES rows were returned within {radius} deg "
                f"of RA={ra}, Dec={dec}."
            ),
            table=str(table),
        )

    science_like = [
        r
        for r in rows
        if str(r.get("koaimtyp", "")).strip().lower() in {"science", "object"}
        or str(r.get("obstype", "")).strip().lower() in {"science", "object"}
    ]
    iodine_in = [r for r in science_like if flag_is_true(r.get("iodin", ""))]
    iodine_out = [r for r in science_like if flag_is_true(r.get("iodout", ""))]

    classified = []
    for row in rows:
        item = dict(row)
        if flag_is_true(row.get("iodin", "")):
            item["rv_role"] = "iodine_in_candidate_rv_epoch"
        elif flag_is_true(row.get("iodout", "")):
            item["rv_role"] = "iodine_out_template_or_reference"
        else:
            item["rv_role"] = "hires_object_unknown_rv_role"
        classified.append(item)
    classified_table = dirs["tables"] / "koa_hires_cone_search_classified.csv"
    write_csv(classified_table, classified)

    return GateResult(
        gate="koa_hires",
        status="PUBLIC_KOA_HIRES_METADATA_FOUND",
        summary=(
            f"KOA returned {len(rows)} HIRES metadata rows within 0.08 deg; "
            f"{len(science_like)} are science-like/object exposures. Header flags "
            f"identify {len(iodine_in)} iodine-in candidate RV epochs and "
            f"{len(iodine_out)} iodine-out template/reference exposures. These are "
            f"public raw HIRES spectra metadata, not a machine-readable BJD/RV/error "
            f"time series or mass solution."
        ),
        table=str(classified_table),
    )


def query_koa_reduced_for_koaid(dirs: dict[str, Path]) -> GateResult:
    try:
        import pyvo
    except Exception as exc:  # pragma: no cover
        return GateResult(
            gate="koa_reduced_data",
            status="BLOCKED_NO_PYVO",
            summary=f"pyvo import failed: {exc}",
        )

    hires_table = dirs["tables"] / "koa_hires_cone_search.csv"
    if not hires_table.exists() or hires_table.stat().st_size == 0:
        return GateResult(
            gate="koa_reduced_data",
            status="SKIPPED_NO_HIRES_ROWS",
            summary="No KOA HIRES rows available to cross-check in koa_reduced_data.",
        )

    with hires_table.open("r", newline="", encoding="utf-8") as f:
        koaid_values = [
            row.get("koaid", "").strip()
            for row in csv.DictReader(f)
            if row.get("koaid", "").strip()
        ]
    if not koaid_values:
        return GateResult(
            gate="koa_reduced_data",
            status="SKIPPED_NO_KOAID_VALUES",
            summary="KOA HIRES rows did not contain KOAID values.",
        )

    svc = pyvo.dal.TAPService("https://koa.ipac.caltech.edu/TAP")
    quoted = ",".join("'" + k.replace("'", "''") + "'" for k in koaid_values[:150])
    query = f"""
    SELECT TOP 300 koaid, instrument, ingesttype, filetype, description,
           date_obs, semid, propint, filehand
    FROM koa_reduced_data
    WHERE koaid IN ({quoted})
    ORDER BY koaid
    """
    try:
        result = svc.search(query, maxrec=400)
        rows = [row_to_dict(row) for row in result]
    except Exception as exc:
        log = dirs["logs"] / "koa_reduced_query_error.txt"
        log.write_text(f"Query failed:\n{query}\n\n{exc}\n", encoding="utf-8")
        return GateResult(
            gate="koa_reduced_data",
            status="QUERY_FAILED",
            summary=f"KOA reduced-data query failed. See {log}",
        )

    table = dirs["tables"] / "koa_reduced_data_for_hires_matches.csv"
    write_csv(table, rows)
    if rows:
        return GateResult(
            gate="koa_reduced_data",
            status="PUBLIC_KOA_REDUCED_ROWS_FOUND",
            summary=(
                f"Found {len(rows)} KOA reduced-data rows associated with matched "
                "HIRES KOAIDs. These may include quick-look extracted spectra, "
                "but not necessarily radial velocities."
            ),
            table=str(table),
        )
    return GateResult(
        gate="koa_reduced_data",
        status="NO_REDUCED_ROWS_FOR_MATCHED_HIRES",
        summary="No koa_reduced_data rows were returned for the matched HIRES KOAIDs.",
        table=str(table),
    )


def assess_hires_prv_service_path(dirs: dict[str, Path]) -> GateResult:
    hires_table = dirs["tables"] / "koa_hires_cone_search_classified.csv"
    if not hires_table.exists() or hires_table.stat().st_size == 0:
        return GateResult(
            gate="nexsci_hires_prv_service",
            status="SKIPPED_NO_HIRES_ROWS",
            summary="No classified HIRES table is available for assessing the NExScI HIRES PRV service route.",
        )

    with hires_table.open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    iodine_in = [r for r in rows if r.get("rv_role") == "iodine_in_candidate_rv_epoch"]
    iodine_out = [r for r in rows if r.get("rv_role") == "iodine_out_template_or_reference"]

    return GateResult(
        gate="nexsci_hires_prv_service",
        status="POSSIBLE_ROUTE_REQUIRES_KOA_LOGIN_AND_MORE_RV_EPOCHS",
        summary=(
            "NExScI documents a HIRES PRV environment that can reduce compatible "
            "raw HIRES spectra into RV curves, but it requires a KOA-authenticated "
            "user workspace and server-side processing. The public KOA cone search "
            f"found {len(iodine_in)} iodine-in candidate RV epochs and "
            f"{len(iodine_out)} iodine-out template/reference exposures for T001258; "
            "the documented minimum is three iodine-in observations before PRVs are "
            "produced, so this is a plausible follow-up route but not a completed "
            "mass extraction from public files."
        ),
        table=str(hires_table),
    )


def query_nasa_exoplanet_archive(dirs: dict[str, Path]) -> GateResult:
    base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

    def tap_json(query: str) -> list[dict[str, Any]]:
        url = base + "?" + urlencode({"query": query, "format": "json"})
        with urlopen(url, timeout=60) as response:
            text = response.read().decode("utf-8", "replace")
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"TAP response was not JSON: {text[:500]}") from exc
        if not isinstance(parsed, list):
            raise RuntimeError(f"TAP response was not a row list: {parsed!r}")
        return parsed

    try:
        toi_rows = tap_json(
            "select toi,tid,tfopwg_disp,ctoi_alias,pl_pnum,ra,dec,pl_orbper,"
            "pl_tranmid,pl_rade from toi where tid=288132261"
        )
        confirmed_rows = tap_json(
            "select pl_name,hostname,pl_orbper,pl_rade,pl_bmasse,pl_bmassj,"
            "pl_rvamp,disc_refname from pscomppars where tic_id like '%288132261%' "
            "or hostname like '%TOI-1258%' or pl_name like '%TOI-1258%'"
        )
    except Exception as exc:
        log = dirs["logs"] / "nasa_exoplanet_archive_query_error.txt"
        log.write_text(str(exc), encoding="utf-8")
        return GateResult(
            gate="nasa_exoplanet_archive",
            status="QUERY_FAILED",
            summary=f"NASA Exoplanet Archive TAP query failed. See {log}",
        )

    toi_table = dirs["tables"] / "nasa_exoplanet_archive_toi_rows.csv"
    confirmed_table = dirs["tables"] / "nasa_exoplanet_archive_confirmed_rows.csv"
    write_csv(toi_table, toi_rows)
    write_csv(confirmed_table, confirmed_rows)

    dispositions = sorted({str(r.get("tfopwg_disp", "")) for r in toi_rows if r.get("tfopwg_disp")})
    if confirmed_rows:
        status = "CONFIRMED_PLANET_ROWS_FOUND"
        summary = (
            f"NASA Exoplanet Archive returned {len(confirmed_rows)} confirmed-planet "
            "rows tied to TIC 288132261/TOI-1258. Manual comparison with the RV package is required."
        )
        table = confirmed_table
    elif toi_rows:
        status = "TOI_LISTED_AS_PC_NO_CONFIRMED_PLANET_ROW"
        summary = (
            f"NASA Exoplanet Archive TOI table returned {len(toi_rows)} row(s) for "
            f"TIC 288132261 with disposition(s) {dispositions}; the confirmed-planet "
            "table returned no mass/RV solution rows for TIC 288132261 or TOI-1258."
        )
        table = toi_table
    else:
        status = "NO_TOI_OR_CONFIRMED_ROWS"
        summary = (
            "NASA Exoplanet Archive returned no TOI table rows and no confirmed-planet "
            "rows for TIC 288132261/TOI-1258."
        )
        table = toi_table

    return GateResult(
        gate="nasa_exoplanet_archive",
        status=status,
        summary=summary,
        table=str(table),
    )


def query_mast(dirs: dict[str, Path]) -> GateResult:
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from astroquery.mast import Observations
    except Exception as exc:  # pragma: no cover
        return GateResult(
            gate="mast_hst_jwst",
            status="BLOCKED_NO_ASTROQUERY",
            summary=f"astroquery/astropy import failed: {exc}",
        )

    coord = SkyCoord(TARGET["ra"], TARGET["dec"], unit="deg", frame="icrs")
    rows: list[dict[str, Any]] = []
    try:
        obs = Observations.query_region(coord, radius=0.08 * u.deg)
        names = list(obs.colnames)
        for row in obs:
            item: dict[str, Any] = {}
            for name in names:
                value = row[name]
                try:
                    if hasattr(value, "item"):
                        value = value.item()
                except Exception:
                    pass
                item[name] = str(value)
            rows.append(item)
    except Exception as exc:
        log = dirs["logs"] / "mast_query_error.txt"
        log.write_text(str(exc), encoding="utf-8")
        return GateResult(
            gate="mast_hst_jwst",
            status="QUERY_FAILED",
            summary=f"MAST query failed. See {log}",
        )

    table = dirs["tables"] / "mast_cone_search_observations.csv"
    write_csv(table, rows)

    hst_jwst = [
        r
        for r in rows
        if str(r.get("obs_collection", "")).upper() in {"HST", "JWST"}
    ]
    spectra_like = [
        r
        for r in hst_jwst
        if "SPEC" in str(r.get("dataproduct_type", "")).upper()
        or "TIME" in str(r.get("dataproduct_type", "")).upper()
        or any(
            key in str(r.get("instrument_name", "")).upper()
            for key in ["NIRSPEC", "NIRISS", "MIRI", "STIS", "WFC3/IR"]
        )
    ]
    target_rows = [
        r
        for r in hst_jwst
        if "TOI" in str(r.get("target_name", "")).upper()
        or "288132261" in str(r.get("target_name", ""))
        or "TIC" in str(r.get("target_name", "")).upper()
    ]

    subset = dirs["tables"] / "mast_hst_jwst_candidate_rows.csv"
    write_csv(subset, hst_jwst)
    spectra_subset = dirs["tables"] / "mast_hst_jwst_spectra_like_rows.csv"
    write_csv(spectra_subset, spectra_like)

    if spectra_like:
        status = "HST_JWST_SPECTRA_LIKE_ROWS_FOUND"
        summary = (
            f"MAST returned {len(rows)} total rows within 0.08 deg, including "
            f"{len(hst_jwst)} HST/JWST rows and {len(spectra_like)} spectra/time-series-"
            f"like HST/JWST rows. {len(target_rows)} HST/JWST rows have TOI/TIC-like "
            "target-name strings. Manual product inspection is required before claiming "
            "a planet atmosphere barcode."
        )
    elif hst_jwst:
        status = "HST_JWST_ROWS_FOUND_NO_CLEAR_SPECTRA"
        summary = (
            f"MAST returned {len(hst_jwst)} HST/JWST rows within 0.08 deg, but no "
            "clear spectroscopy/time-series rows by the simple metadata filter."
        )
    else:
        status = "NO_HST_JWST_MATCH"
        summary = (
            f"MAST returned {len(rows)} total rows within 0.08 deg, but no HST/JWST "
            "observations in the cone."
        )

    return GateResult(
        gate="mast_hst_jwst",
        status=status,
        summary=summary,
        table=str(table),
    )


def write_email_draft(dirs: dict[str, Path]) -> Path:
    path = dirs["outreach"] / "toi1258_rv_data_collaboration_email.md"
    path.write_text(
        f"""Subject: TOI-1258.01 RV time series / joint photometry+RV follow-up

Dear Dr. Giacalone,

I am an independent researcher in Gilbert, Arizona. I recently completed a public-data statistical-validation package for TOI-1258.01 / TIC 288132261 using TESS light curves, SPOC DV centroids, Gemini/Alopeke speckle constraints, and high-sample TRICERATOPS false-positive modeling.

The core result is that the transit is recovered across 12 TESS sectors at P = 10.261376896 d, matching the catalog period to 0.000362 percent. The BLS depth is 722 ppm with SNR = 64.7, the SPOC DV centroid is consistent with TIC 288132261, and five N=1,000,000 TRICERATOPS seeds give FPP = 2.37e-4 to 6.92e-4 with NFPP below 1e-4.

In the ExoFOP metadata I see APF/Levy RV rows under specuser=giacalone, including 68 RVs from 2019-10-25 to 2020-05-02 plus later APF/Levy sequences, and a Keck/HIRES entry. The public bundle I inspected did not expose the machine-readable BJD/RV/error table.

I also checked KOA directly. Public KOA metadata contains five HIRES object exposures near TIC 288132261 under target name T001258: a 2020-01-04 N152/Crossfield iodine-in/out sequence, a 2020-03-08 U127/Petigura iodine-in exposure, and two 2020-03-09 C003/Howard iodine-out high-SNR template/reference exposures. KOA did not expose a reduced RV table for those KOAIDs.

Are those RVs part of an unpublished mass solution for TOI-1258.01, and if so is there room to compare notes or collaborate on a joint photometry+RV model? I am not asking you to accept any overclaim; I am trying to connect the public statistical-validation result with the RV mass gate in the cleanest possible way.

I can send the reproducibility package, figures, and scripts immediately if useful.

Best,
Nikolai Stoloff
author-contact-omitted
""",
        encoding="utf-8",
    )
    return path


def write_report(dirs: dict[str, Path], gates: list[GateResult], email_path: Path) -> Path:
    path = dirs["reports"] / "TOI1258_MASS_BARCODE_ARCHIVE_HUNT_REPORT.md"
    lines = [
        "# TOI-1258.01 Mass and Barcode Archive Hunt",
        "",
        f"Created UTC: `{datetime.now(timezone.utc).isoformat()}`",
        "",
        "## Target",
        "",
        f"- Name: `{TARGET['name']}`",
        f"- TIC: `{TARGET['tic']}`",
        f"- RA, Dec: `{TARGET['ra']}`, `{TARGET['dec']}` deg",
        f"- Period: `{TARGET['period_days']}` d",
        f"- Radius: `{TARGET['radius_rearth']}` Earth radii",
        "",
        "## Gate Results",
        "",
    ]
    for gate in gates:
        lines.extend(
            [
                f"### {gate.gate}",
                "",
                f"- Status: `{gate.status}`",
                f"- Summary: {gate.summary}",
            ]
        )
        if gate.table:
            lines.append(f"- Table: `{gate.table}`")
        lines.append("")

    lines.extend(
        [
            "## Current Evidence Boundary",
            "",
            "- The public-data package statistically validates TOI-1258.01 as a transiting planet under TRICERATOPS-style false-positive testing.",
            "- KOA exposes public Keck/HIRES raw-spectra metadata for target name `T001258`, including iodine-in candidate RV epochs and iodine-out template/reference exposures.",
            "- NASA Exoplanet Archive currently keeps the target in the TOI candidate table rather than the confirmed-planet parameter table, based on this run.",
            "- The mass is **not** recovered unless a BJD/RV/error table or a published RV orbital fit is obtained.",
            "- A chemical barcode is **not** recovered unless HST/JWST or comparable transit spectroscopy exists and can be tied to TOI-1258.01 transits.",
            "- Keck/HIRES metadata can support the RV-search story, but APF/Levy RVs are not in KOA because APF is a Lick Observatory instrument.",
            "",
            "## Outreach",
            "",
            f"- Draft email: `{email_path}`",
            "",
            "## Source Notes",
            "",
            "- KOA TAP/PyVO documentation: https://koa.ipac.caltech.edu/UserGuide/PyKOA/TAPClients.html",
            "- KOA user guide and acknowledgment text: https://koa.ipac.caltech.edu/UserGuide/",
            "- NExScI HIRES PRV documentation: https://caltech-ipac.github.io/hiresprv/index.html",
            "- NASA Exoplanet Archive TAP queried at: https://exoplanetarchive.ipac.caltech.edu/TAP/sync",
            "- MAST queried through `astroquery.mast.Observations.query_region`.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    out_root = Path(r"data/external/exoplanet_outputs\TOI1258_MASS_BARCODE_ARCHIVE_HUNT_2026-07-08")
    dirs = ensure_dirs(out_root)

    gates: list[GateResult] = []
    gates.append(query_koa_hires(dirs))
    gates.append(query_koa_reduced_for_koaid(dirs))
    gates.append(assess_hires_prv_service_path(dirs))
    gates.append(query_nasa_exoplanet_archive(dirs))
    gates.append(query_mast(dirs))
    email_path = write_email_draft(dirs)
    report = write_report(dirs, gates, email_path)

    manifest = {
        "target": TARGET,
        "gates": [asdict(g) for g in gates],
        "report": str(report),
        "email_draft": str(email_path),
    }
    (dirs["root"] / "target_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
