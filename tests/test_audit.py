"""Deterministic tests on real historical evidence and explicit failure modes."""

import json
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from toi1258.analysis import analyze
from toi1258.common import number, ratio, read_csv, verify_manifest

ROOT = Path(__file__).resolve().parents[1]


class AuditTests(unittest.TestCase):
    def test_source_inventory_is_complete(self):
        inventory = json.loads((ROOT / "research/SOURCE_INVENTORY.json").read_text())
        self.assertEqual(
            inventory["workspace_project_scripts_discovered"], inventory["workspace_project_scripts_exported"]
        )
        self.assertGreater(inventory["workspace_project_scripts_exported"], 10)
        seen = set()
        for item in inventory["files"]:
            self.assertNotIn(item["path"], seen)
            seen.add(item["path"])
            path = (ROOT / item["path"]).resolve()
            self.assertTrue(path.is_relative_to(ROOT.resolve()))
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), item["sha256"], item["path"])

    def test_export_does_not_corrupt_urls(self):
        for path in (ROOT / "research").rglob("*.py"):
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("httpdata/external", text)
            self.assertNotIn("httpdata/", text)

    def test_preserved_original_empty_files_are_explicit(self):
        inventory = json.loads((ROOT / "research/SOURCE_INVENTORY.json").read_text())
        flagged = {r["path"] for r in inventory["files"] if r["empty_original"]}
        self.assertEqual(flagged, set(inventory["empty_original_files"]))
        for path in flagged:
            self.assertFalse((ROOT / path).read_bytes().strip())

    def test_example_matches_regenerated_summary(self):
        summary, _ = analyze(ROOT / "data")
        summary["verified_provenance_files"] = verify_manifest(ROOT)
        expected = json.loads((ROOT / "examples/summary.json").read_text())
        self.assertEqual(summary, expected)

    def test_manifest(self):
        self.assertGreater(verify_manifest(ROOT), 5)

    def test_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            (root / "data").mkdir()
            (root / "data/example.csv").write_text("x\n1\n")
            (root / "data/manifest.json").write_text(
                json.dumps({"files": [{"path": "data/example.csv", "sha256": "0" * 64}]})
            )
            with self.assertRaisesRegex(ValueError, "checksum"):
                verify_manifest(root)

    def test_manifest_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            (root / "data").mkdir()
            (root / "data/manifest.json").write_text(
                json.dumps({"files": [{"path": "../outside.txt", "sha256": "0" * 64}]})
            )
            with self.assertRaisesRegex(ValueError, "escapes"):
                verify_manifest(root)

    def test_invalid_numbers(self):
        for value in ("nan", "inf", "-inf"):
            with self.assertRaises(ValueError):
                number(value)
        for value in (0, -1):
            with self.assertRaises(ValueError):
                ratio(1, value)

    def test_malformed_csv(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "bad.csv"
            path.write_text("a,b\n1\n")
            with self.assertRaises(ValueError):
                read_csv(path)

    def test_cli_from_other_directory(self):
        with tempfile.TemporaryDirectory() as name:
            result = subprocess.run(
                [sys.executable, "-m", "toi1258", "--root", str(ROOT)],
                cwd=name,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            output = Path(name) / "results/summary.json"
            self.assertTrue(output.exists())
            self.assertIn("interpretation", json.loads(output.read_text()))

    def test_cli_missing_data_fails(self):
        with tempfile.TemporaryDirectory() as name:
            result = subprocess.run(
                [sys.executable, "-m", "toi1258", "--root", name], cwd=name, capture_output=True, text=True
            )
            self.assertNotEqual(result.returncode, 0)

    def test_cli_protects_inputs(self):
        result = subprocess.run(
            [sys.executable, "-m", "toi1258", "--root", str(ROOT), "--output", str(ROOT / "data")],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_posterior_reconstruction(self):
        summary, rows = analyze(ROOT / "data")
        self.assertEqual(len(rows), 5)
        self.assertTrue(summary["all_below_reported_thresholds"])
        self.assertAlmostEqual(summary["FPP_range"][0], 0.000236863454368796, places=12)
        self.assertEqual(summary["sectors"], 12)
        self.assertFalse(summary["mass_measured_here"])
        self.assertFalse(summary["current_archive_status_checked"])

    def test_probability_partition(self):
        from toi1258.analysis import scenario_probabilities

        fpp, nfpp = scenario_probabilities(
            [
                {"scenario": "TP", "prob": "0.8"},
                {"scenario": "PTP", "prob": "0.1"},
                {"scenario": "EB", "prob": "0.06"},
                {"scenario": "NTP", "prob": "0.04"},
            ]
        )
        self.assertAlmostEqual(fpp, 0.1)
        self.assertAlmostEqual(nfpp, 0.04)

    def test_invalid_probability_sets(self):
        from toi1258.analysis import scenario_probabilities

        for rows in (
            [],
            [{"scenario": "TP", "prob": "0.5"}],
            [{"scenario": "UNKNOWN", "prob": "1"}],
            [{"scenario": "TP", "prob": "-1"}],
        ):
            with self.assertRaises(ValueError):
                scenario_probabilities(rows)


if __name__ == "__main__":
    unittest.main()
