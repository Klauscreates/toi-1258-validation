"""Historical model inputs must retain their exact exported bytes."""
import hashlib
import json
from pathlib import Path
import unittest


class HistoricalInputTests(unittest.TestCase):
    def test_inputs_match_export_provenance(self):
        root = Path(__file__).resolve().parents[1] / "reproduce/historical-inputs"
        records = json.loads((root / "PROVENANCE.json").read_text())
        self.assertTrue(records)
        for row in records:
            self.assertEqual(row["status"], "PRESERVED")
            data = (root / row["name"]).read_bytes()
            self.assertEqual(len(data), row["bytes"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), row["sha256"])


if __name__ == "__main__":
    unittest.main()
