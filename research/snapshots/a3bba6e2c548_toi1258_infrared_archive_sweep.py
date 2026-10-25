"""TOI-1258.01 infrared / spectroscopy archive sweep.

Searches public metadata around TIC 288132261 / TOI-1258.01 for:
- MAST HST/JWST observations,
- IRSA infrared catalog detections (2MASS/WISE/NEOWISE/CatWISE),
- Keck/KOA near-IR instrument metadata,
- Gemini public archive region metadata.

This does not claim an atmospheric barcode; it identifies whether suitable
wavelength-resolved or infrared observations already exist publicly.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

TARGET = {
    "name": "TOI-1258.01",
    "tic": "288132261",
    "ra": 209.484313,
    "dec": 79.583322,
    "radius_deg": 0.08,
}

OUT = Path(r"data/external/exoplanet_outputs\TOI1258_INFRARED_ARCHIVE_SWEEP_2026-07-08")
REPORTS = OUT / "evidence" / "reports"
TABLES = OUT / "evidence" / "tables"
LOGS = OUT / "evidence" / "logs"
SCRIPTS = OUT / "scripts"
for d in [REPORTS, TABLES, LOGS, SCRIPTS]:
    d.mkdir(parents=True, exist_ok=True)

@dataclass
class Gate:
    gate: str
    status: str
    summary: str
    table: str | None = None


def row_to_dict(row: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    if hasattr(row, "items"):
        iterator = row.items()
    elif hasattr(row, "keys"):
        iterator = ((k, row[k]) for k in row.keys())
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
        out[str(name)] = value
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields=[]
    for r in rows:
        for k in r.keys():
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def mast_gate() -> Gate:
    try:
        from astroquery.mast import Observations
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except Exception as e:
        return Gate("mast_hst_jwst", "BLOCKED_IMPORT", f"astroquery/astropy import failed: {e}")
    try:
        coord = SkyCoord(TARGET["ra"], TARGET["dec"], unit="deg", frame="icrs")
        obs = Observations.query_region(coord, radius=TARGET["radius_deg"] * u.deg)
        rows=[]
        for row in obs:
            rows.append({c: str(row[c]) for c in obs.colnames})
    except Exception as e:
        p=LOGS/"mast_error.txt"; p.write_text(str(e), encoding="utf-8")
        return Gate("mast_hst_jwst", "QUERY_FAILED", f"MAST query failed: {p}")
    table=TABLES/"mast_region_observations.csv"
    write_csv(table, rows)
    hst_jwst=[r for r in rows if str(r.get("obs_collection","")).upper() in {"HST","JWST"}]
    spectra=[r for r in hst_jwst if "SPEC" in str(r.get("dataproduct_type","")).upper() or any(x in str(r.get("instrument_name","")).upper() for x in ["NIRSPEC","NIRISS","MIRI","STIS","WFC3"])]
    write_csv(TABLES/"mast_hst_jwst_rows.csv", hst_jwst)
    write_csv(TABLES/"mast_hst_jwst_spectra_like_rows.csv", spectra)
    if spectra:
        return Gate("mast_hst_jwst", "HST_JWST_SPECTRA_LIKE_FOUND", f"MAST returned {len(rows)} total rows; {len(hst_jwst)} HST/JWST and {len(spectra)} spectra-like rows. Product inspection required.", str(table))
    if hst_jwst:
        return Gate("mast_hst_jwst", "HST_JWST_IMAGES_ONLY_OR_UNCLEAR", f"MAST returned {len(rows)} total rows and {len(hst_jwst)} HST/JWST rows, but no clear spectroscopy rows.", str(table))
    collections={}
    for r in rows:
        k=r.get("obs_collection","")
        collections[k]=collections.get(k,0)+1
    return Gate("mast_hst_jwst", "NO_HST_JWST_MATCH", f"MAST returned {len(rows)} total rows but no HST/JWST rows. Collections: {collections}", str(table))


def irsa_gate() -> Gate:
    try:
        from astroquery.ipac.irsa import Irsa
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except Exception as e:
        return Gate("irsa_infrared_catalogs", "BLOCKED_IMPORT", f"astroquery/astropy import failed: {e}")
    coord=SkyCoord(TARGET["ra"], TARGET["dec"], unit="deg", frame="icrs")
    catalogs = {
        "2MASS_PSC_fp_psc": "fp_psc",
        "AllWISE_allwise_p3as_psd": "allwise_p3as_psd",
        "NEOWISE_neowiser_p1bs_psd": "neowiser_p1bs_psd",
        "CatWISE2020_catwise_2020": "catwise_2020",
    }
    summaries=[]
    errors=[]
    total=0
    for label, cat in catalogs.items():
        try:
            tab=Irsa.query_region(coord, catalog=cat, spatial="Cone", radius=5*u.arcsec)
            rows=[{c: (None if getattr(row[c], 'mask', False) else str(row[c])) for c in tab.colnames} for row in tab]
            out=TABLES/(label+".csv")
            write_csv(out, rows)
            total += len(rows)
            summaries.append(f"{label}: {len(rows)} rows")
        except Exception as e:
            errors.append(f"{label}: {e}")
            (LOGS/(label+"_error.txt")).write_text(str(e), encoding="utf-8")
    status="IR_CATALOG_DETECTIONS_FOUND" if total else "NO_IR_CATALOG_ROWS_OR_QUERIES_FAILED"
    summary="; ".join(summaries) if summaries else "No catalog summaries produced."
    if errors:
        summary += "; errors: " + " | ".join(errors)
    return Gate("irsa_infrared_catalogs", status, summary, str(TABLES))


def koa_ir_gate() -> Gate:
    try:
        import pyvo
    except Exception as e:
        return Gate("koa_ground_ir", "BLOCKED_NO_PYVO", f"pyvo import failed: {e}")
    svc=pyvo.dal.TAPService("https://koa.ipac.caltech.edu/TAP")
    ra=TARGET["ra"]; dec=TARGET["dec"]; radius=TARGET["radius_deg"]
    tables={
        "koa_nirc":"koaid,targname,koaimtyp,date_obs,ut,progid,progpi,progtitl,ra,dec,filehand",
        "koa_nirc2":"koaid,targname,koaimtyp,date_obs,ut,filter,slitname,specres,waveblue,wavecntr,wavered,progid,progpi,progtitl,ra,dec,filehand",
        "koa_nirspec":"koaid,targname,object,koaimtyp,date_obs,ut,filter,slitname,specres,waveblue,wavecntr,wavered,progid,progpi,progtitl,ra,dec,filehand",
        "koa_mosfire":"koaid,targname,object,koaimtyp,date_obs,ut,filter,obsmode,waveblue,wavecntr,wavered,progid,progpi,progtitl,ra,dec,filehand",
        "koa_osiris":"koaid,targname,object,koaimtyp,date_obs,ut,filter,waveblue,wavecntr,wavered,progid,progpi,progtitl,ra,dec,filehand",
    }
    counts={}
    all_rows=[]
    for table, cols in tables.items():
        q=f"""
        SELECT TOP 200 {cols}
        FROM {table}
        WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))
        ORDER BY date_obs
        """
        try:
            res=svc.search(q, maxrec=250)
            rows=[]
            for row in res:
                d=row_to_dict(row); d["koa_table"]=table; rows.append(d)
            counts[table]=len(rows)
            write_csv(TABLES/(table+"_region_rows.csv"), rows)
            all_rows.extend(rows)
        except Exception as e:
            counts[table]=f"ERROR: {e}"
            (LOGS/(table+"_error.txt")).write_text(q+"\n\n"+str(e), encoding="utf-8")
    combined=TABLES/"koa_ground_ir_combined_rows.csv"
    write_csv(combined, all_rows)
    spectroscopy=[r for r in all_rows if r.get("koa_table") in {"koa_nirspec","koa_mosfire","koa_osiris"}]
    status="KOA_GROUND_IR_ROWS_FOUND" if all_rows else "NO_KOA_GROUND_IR_ROWS"
    return Gate("koa_ground_ir", status, f"KOA IR search counts: {counts}. Spectroscopy-capable rows: {len(spectroscopy)}. These are host-coordinate archive rows, not planet-isolated barcodes.", str(combined))


def gemini_gate() -> Gate:
    try:
        from astroquery.gemini import Observations
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except Exception as e:
        return Gate("gemini_archive", "BLOCKED_IMPORT", f"astroquery/astropy import failed: {e}")
    coord=SkyCoord(TARGET["ra"], TARGET["dec"], unit="deg", frame="icrs")
    try:
        tab=Observations.query_region(coord, radius=TARGET["radius_deg"]*u.deg)
        rows=[]
        for row in tab:
            rows.append({c: str(row[c]) for c in tab.colnames})
    except Exception as e:
        p=LOGS/"gemini_error.txt"; p.write_text(str(e), encoding="utf-8")
        return Gate("gemini_archive", "QUERY_FAILED", f"Gemini query failed: {p}")
    table=TABLES/"gemini_region_rows.csv"
    write_csv(table, rows)
    instruments={}
    for r in rows:
        k=r.get("instrument", r.get("Instrument", ""))
        instruments[k]=instruments.get(k,0)+1
    status="GEMINI_ROWS_FOUND" if rows else "NO_GEMINI_ROWS"
    return Gate("gemini_archive", status, f"Gemini returned {len(rows)} rows. Instrument counts: {instruments}", str(table))


def main():
    gates=[mast_gate(), irsa_gate(), koa_ir_gate(), gemini_gate()]
    report_lines=[
        "# TOI-1258.01 Infrared and Spectroscopy Archive Sweep",
        "",
        f"Created UTC: `{datetime.now(timezone.utc).isoformat()}`",
        "",
        "## Target",
        "",
        f"- Target: `{TARGET['name']} / TIC {TARGET['tic']}`",
        f"- RA, Dec: `{TARGET['ra']}`, `{TARGET['dec']}` deg",
        f"- Search radius: `{TARGET['radius_deg']}` deg",
        "",
        "## Gate Results",
        "",
    ]
    for g in gates:
        report_lines += [f"### {g.gate}", "", f"- Status: `{g.status}`", f"- Summary: {g.summary}"]
        if g.table: report_lines.append(f"- Table: `{g.table}`")
        report_lines.append("")
    report_lines += [
        "## Interpretation Boundary",
        "",
        "- Infrared catalog detections of the host star do not equal a planet atmosphere barcode.",
        "- A barcode requires wavelength-resolved in-transit/out-of-transit spectroscopy, eclipse spectroscopy, phase spectroscopy, or high-resolution Doppler cross-correlation data that can isolate the planet signal.",
        "- The planet is expected to be separated from the host star by less than 1 milliarcsecond on the sky, so normal imaging cannot resolve it as a separate object.",
        "- If no JWST/HST/Gemini/Keck spectroscopy rows exist at the target coordinate, the barcode is not in public archive data yet and must be obtained through new observations or hidden/non-public data release.",
    ]
    report=REPORTS/"TOI1258_INFRARED_ARCHIVE_SWEEP_REPORT.md"
    report.write_text("\n".join(report_lines)+"\n", encoding="utf-8")
    manifest={"target":TARGET,"gates":[asdict(g) for g in gates],"report":str(report)}
    (OUT/"target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    try:
        src=Path(__file__)
        (SCRIPTS/src.name).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    except Exception:
        pass
    print(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()
