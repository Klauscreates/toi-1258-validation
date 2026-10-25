"""Acquire or verify exact research inputs. Downloads are explicit and hash checked."""
from pathlib import Path, PurePosixPath
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import urllib.request
import zipfile


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def safe_path(root, name):
    rel = PurePosixPath(name)
    if rel.is_absolute() or ".." in rel.parts or "\\" in name or ":" in name:
        raise ValueError(f"Unsafe input path: {name}")
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()) or path == root.resolve():
        raise ValueError(f"Input escapes root: {name}")
    return path


def verify(path, record):
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
        raise ValueError(f"Exact input checksum mismatch: {record.get('path', record.get('name'))}")


def checked_extract(archive, root):
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            target = safe_path(root, info.filename.rstrip("/"))
            if (info.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Symlink ZIP member rejected")
            if not info.is_dir() and target.exists():
                raise FileExistsError(f"Refusing to overwrite existing input: {target}")
        z.extractall(root)


def fetch(record, target):
    url = record["url"]
    if not url or not url.startswith("https://"):
        raise ValueError("No HTTPS archive retrieval URL")
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".download-part")
    if part.exists():
        raise FileExistsError(f"Previous partial download exists: {part}")
    req = urllib.request.Request(url, headers={"User-Agent": "Research-input-reproduction/1.0"})
    with urllib.request.urlopen(req, timeout=180) as response, part.open("xb") as output:
        if not response.geturl().startswith("https://"):
            raise ValueError("Insecure redirect")
        shutil.copyfileobj(response, output, length=1024 * 1024)
    # ESO may return .Z or gzip despite a .fits filename. Hash the decoded file.
    if target.suffix == ".fits" or ".eso.org/" in url:
        with part.open("rb") as f:
            magic = f.read(2)
        decoded = part.with_name(part.name + ".decoded")
        if magic == b"\x1f\x8b":
            with gzip.open(part, "rb") as source, decoded.open("xb") as output:
                shutil.copyfileobj(source, output)
        elif magic == b"\x1f\x9d":
            from unlzw3 import unlzw
            decoded.write_bytes(unlzw(part.read_bytes()))
        else:
            decoded = part
        verify(decoded, record)
        decoded.replace(target)
        if part.exists():
            part.unlink()
    else:
        verify(part, record)
        part.replace(target)


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=repo / "reproduce/inputs.json")
    parser.add_argument("--root", type=Path, required=True, help="External data directory, not Git-tracked")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--release-assets", action="store_true", help="Download supplemental ZIPs using authenticated gh")
    parser.add_argument("--select", default="", help="Substring of input path")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--local-map", type=Path, help="Optional private path map for verifying an existing cache")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if root.is_relative_to(repo):
        raise ValueError("Keep external data outside this Git repository")
    manifest = json.loads(args.manifest.read_text())
    mapping = json.loads(args.local_map.read_text())[manifest["project"]] if args.local_map else {}
    if args.release_assets:
        downloads = root / ".downloads"
        downloads.mkdir(parents=True, exist_ok=True)
        for asset in manifest["assets"]:
            archive = downloads / asset["name"]
            if not archive.exists():
                subprocess.run(["gh", "release", "download", asset["tag"], "--repo", asset["repository"],
                                "--pattern", asset["name"], "--dir", str(downloads)], check=True)
            verify(archive, asset)
            members = [r for r in manifest["inputs"] if r["release_asset"] == asset["name"]]
            if all(safe_path(root, r["path"]).exists() for r in members):
                for r in members:
                    verify(safe_path(root, r["path"]), r)
            else:
                checked_extract(archive, root)
    records = [r for r in manifest["inputs"] if args.select in r["path"]]
    if args.limit:
        records = records[:args.limit]
    if not records:
        raise ValueError("No matching inputs")
    results = []
    for r in records:
        path = Path(mapping[r["path"]]) if r["path"] in mapping else safe_path(root, r["path"])
        try:
            if not path.exists() and args.download:
                fetch(r, path)
            verify(path, r)
            results.append(dict(path=r["path"], status="PASS"))
        except (OSError, ValueError) as exc:
            # Do not leak local machine paths into the portable verification report.
            results.append(dict(path=r["path"], status="FAIL", error_type=type(exc).__name__))
            print(r["path"], type(exc).__name__, flush=True)
    report = dict(project=manifest["project"], tested=len(results),
                  passed=sum(r["status"] == "PASS" for r in results), results=results)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "results"}))
    if report["passed"] != report["tested"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
