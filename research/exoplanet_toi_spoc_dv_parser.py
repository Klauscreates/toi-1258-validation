#!/usr/bin/env python
"""Extract SPOC Data Validation centroid metrics from TESS DVR XML files."""

from __future__ import annotations

import argparse
import csv
import xml.etree.ElementTree as ET
from pathlib import Path


NS = {"dv": "http://www.nasa.gov/2018/TESS/DV"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml-root", required=True, help="Directory containing downloaded *_dvr.xml files.")
    parser.add_argument("--output-dir", required=True, help="Directory where CSV tables should be written.")
    parser.add_argument("--period-days", type=float, required=True)
    parser.add_argument("--period-tolerance-days", type=float, default=0.02)
    args = parser.parse_args()

    xml_root = Path(args.xml_root)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary, sector_rows = parse_xml_tree(xml_root, args.period_days, args.period_tolerance_days)
    write_csv(out_dir / "spoc_dv_centroid_summary.csv", summary)
    write_csv(out_dir / "spoc_dv_sector_difference_image_offsets.csv", sector_rows)
    print(
        {
            "xml_files": len(list(xml_root.glob("**/*_dvr.xml"))),
            "summary_rows": len(summary),
            "sector_rows": len(sector_rows),
            "summary_csv": str(out_dir / "spoc_dv_centroid_summary.csv"),
            "sector_csv": str(out_dir / "spoc_dv_sector_difference_image_offsets.csv"),
        }
    )
    return 0


def parse_xml_tree(xml_root: Path, period_days: float, period_tolerance_days: float):
    summary = []
    sector_rows = []
    for xml_path in sorted(xml_root.glob("**/*_dvr.xml")):
        root = ET.parse(xml_path).getroot()
        root_info = {
            "xml_file": str(xml_path),
            "filename": xml_path.name,
            "startCadence": root.attrib.get("startCadence", ""),
            "endCadence": root.attrib.get("endCadence", ""),
            "ticId": root.attrib.get("ticId", ""),
            "toiId_root": root.attrib.get("toiId", ""),
            "matchedToiId": root.findtext("dv:matchedToiId", default="", namespaces=NS) or "",
        }
        for planet_result in root.findall(".//dv:planetResults", NS):
            fit = planet_result.find("dv:allTransitsFit", NS)
            period, period_unc = model_param(fit, "orbitalPeriodDays")
            try:
                measured_period = float(period)
            except Exception:
                continue
            if abs(measured_period - period_days) > period_tolerance_days:
                continue

            epoch, epoch_unc = model_param(fit, "transitEpochBtjd")
            duration, duration_unc = model_param(fit, "transitDurationHours")
            depth, depth_unc = model_param(fit, "transitDepthPpm")
            radius, radius_unc = model_param(fit, "planetRadiusEarthRadii")
            motion = planet_result.find(".//dv:differenceImageMotionResults", NS)
            ms_tic = motion.find("dv:msTicCentroidOffsets", NS) if motion is not None else None
            ms_control = motion.find("dv:msControlCentroidOffsets", NS) if motion is not None else None
            quality = motion.find("dv:summaryQualityMetric", NS) if motion is not None else None
            overlap = motion.find("dv:summaryOverlapMetric", NS) if motion is not None else None

            row = dict(root_info)
            row.update(
                {
                    "planetNumber": planet_result.attrib.get("planetNumber", ""),
                    "toiId_planet": planet_result.attrib.get("toiId", ""),
                    "toiCorrelation": planet_result.attrib.get("toiCorrelation", ""),
                    "modelFitSnr": fit.attrib.get("modelFitSnr", "") if fit is not None else "",
                    "period_days": period,
                    "period_uncertainty_days": period_unc,
                    "epoch_btjd": epoch,
                    "epoch_uncertainty_days": epoch_unc,
                    "duration_hours": duration,
                    "duration_uncertainty_hours": duration_unc,
                    "depth_ppm": depth,
                    "depth_uncertainty_ppm": depth_unc,
                    "planet_radius_rearth": radius,
                    "planet_radius_uncertainty_rearth": radius_unc,
                    "ms_tic_mean_ra_offset_arcsec": value(child(ms_tic, "meanRaOffset")),
                    "ms_tic_mean_ra_offset_unc_arcsec": value(child(ms_tic, "meanRaOffset"), "uncertainty"),
                    "ms_tic_mean_dec_offset_arcsec": value(child(ms_tic, "meanDecOffset")),
                    "ms_tic_mean_dec_offset_unc_arcsec": value(child(ms_tic, "meanDecOffset"), "uncertainty"),
                    "ms_tic_mean_sky_offset_arcsec": value(child(ms_tic, "meanSkyOffset")),
                    "ms_tic_mean_sky_offset_unc_arcsec": value(child(ms_tic, "meanSkyOffset"), "uncertainty"),
                    "ms_control_mean_sky_offset_arcsec": value(child(ms_control, "meanSkyOffset")),
                    "ms_control_mean_sky_offset_unc_arcsec": value(child(ms_control, "meanSkyOffset"), "uncertainty"),
                    "quality_good_metrics": attr(quality, "numberOfGoodMetrics"),
                    "quality_total_metrics": attr(quality, "numberOfMetrics"),
                    "quality_fraction_good_metrics": attr(quality, "fractionOfGoodMetrics"),
                    "overlap_image_count": attr(overlap, "imageCount"),
                    "overlap_no_overlap_count": attr(overlap, "imageCountNoOverlap"),
                    "overlap_fraction_no_overlap": attr(overlap, "imageCountFractionNoOverlap"),
                }
            )
            summary.append(row)

            for difference_image in planet_result.findall("dv:differenceImageResults", NS):
                control_offsets = difference_image.find("dv:controlCentroidOffsets", NS)
                diff_centroid = difference_image.find("dv:differenceImageCentroid", NS)
                control_centroid = difference_image.find("dv:controlImageCentroid", NS)
                sector_rows.append(
                    {
                        "xml_file": str(xml_path),
                        "filename": xml_path.name,
                        "sector": difference_image.attrib.get("sector", ""),
                        "numberOfTransits": difference_image.attrib.get("numberOfTransits", ""),
                        "numberOfCadencesInTransit": difference_image.attrib.get("numberOfCadencesInTransit", ""),
                        "numberOfCadencesOutOfTransit": difference_image.attrib.get("numberOfCadencesOutOfTransit", ""),
                        "overlappedTransits": difference_image.attrib.get("overlappedTransits", ""),
                        "control_row_offset_pix": value(child(control_offsets, "rowOffset")),
                        "control_col_offset_pix": value(child(control_offsets, "columnOffset")),
                        "control_focal_plane_offset_pix": value(child(control_offsets, "focalPlaneOffset")),
                        "control_sky_offset_arcsec": value(child(control_offsets, "skyOffset")),
                        "control_sky_offset_unc_arcsec": value(child(control_offsets, "skyOffset"), "uncertainty"),
                        "difference_centroid_row": value(child(diff_centroid, "row")),
                        "difference_centroid_col": value(child(diff_centroid, "column")),
                        "control_centroid_row": value(child(control_centroid, "row")),
                        "control_centroid_col": value(child(control_centroid, "column")),
                    }
                )
    return summary, sector_rows


def model_param(fit, name: str) -> tuple[str, str]:
    if fit is None:
        return "", ""
    for param in fit.findall(".//dv:modelParameter", NS):
        if param.attrib.get("name") == name:
            return param.attrib.get("value", ""), param.attrib.get("uncertainty", "")
    return "", ""


def child(parent, tag: str):
    return parent.find(f"dv:{tag}", NS) if parent is not None else None


def value(element, name: str = "value") -> str:
    return element.attrib.get(name, "") if element is not None else ""


def attr(element, name: str) -> str:
    return element.attrib.get(name, "") if element is not None else ""


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
