import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("raw_inputs", ROOT / "tools/raw_inputs.py")
raw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(raw)


class RawInputTests(unittest.TestCase):
    def test_paths_cannot_escape(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            for bad in ("../outside", "/outside", "C:/outside", "a\\b", "."):
                with self.assertRaises(ValueError):
                    raw.safe_path(root, bad)

    def test_hash_is_exact(self):
        with tempfile.TemporaryDirectory() as name:
            p = Path(name) / "input.fits"
            p.write_bytes(b"evidence")
            record = dict(path=p.name, bytes=8, sha256=hashlib.sha256(b"evidence").hexdigest())
            raw.verify(p, record)
            p.write_bytes(b"Evidence")
            with self.assertRaises(ValueError):
                raw.verify(p, record)

    def test_zip_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            z = root / "bad.zip"
            with zipfile.ZipFile(z, "w") as handle:
                handle.writestr("../escape", "bad")
            with self.assertRaises(ValueError):
                raw.checked_extract(z, root / "data")

    def test_zip_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            z = root / "data.zip"
            with zipfile.ZipFile(z, "w") as handle:
                handle.writestr("input.txt", "original")
            raw.checked_extract(z, root / "data")
            with self.assertRaises(FileExistsError):
                raw.checked_extract(z, root / "data")

    def test_manifest_paths_unique(self):
        import json
        manifest = json.loads((ROOT / "reproduce/inputs.json").read_text())
        names = [r["path"] for r in manifest["inputs"]]
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(names)
        for row in manifest["inputs"]:
            self.assertEqual(len(row["sha256"]), 64)
            self.assertGreater(row["bytes"], 0)
            self.assertTrue(row["url"] or row["release_asset"])


if __name__ == "__main__":
    unittest.main()
