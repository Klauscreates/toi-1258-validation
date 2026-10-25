"""Repository hygiene checks; no credentials, raw archives, or broken doc links."""

from __future__ import annotations

import ast
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SKIP = {".git", ".venv", "__pycache__", "build", "dist", "results"}
BLOCKED = {".fits", ".fit", ".fz", ".exe", ".dll", ".pfx", ".pem"}
TEXT = {
    ".py",
    ".md",
    ".json",
    ".csv",
    ".toml",
    ".cff",
    ".yml",
    ".yaml",
    ".txt",
    ".ps1",
    ".sh",
    ".tex",
    ".bib",
}


def main() -> int:
    errors = []
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT)
        if any(p in SKIP or p.endswith(".egg-info") for p in rel.parts) or not path.is_file():
            continue
        if path.suffix.lower() in BLOCKED or path.stat().st_size > 10_000_000:
            errors.append(f"{rel}: forbidden file type or >10 MB")
        if path.suffix.lower() not in TEXT:
            continue
        text = path.read_text(encoding="utf-8")
        if re.search(r"[A-Z]:[\\/](?:Users|microlens|exoplanet|sirius)", text, re.I):
            errors.append(f"{rel}: local machine path")
        private_links = ("ms_id_" + "key=", "form_type=" + "object_log")
        if any(term in text for term in private_links) or re.search(
            r"(?:github_pat_|ghp_)[A-Za-z0-9_]{20,}|-{5}BEGIN [A-Z ]*PRIVATE KEY-{5}", text
        ):
            errors.append(f"{rel}: possible credential or private submission link")
        if re.search(r"[\w.+-]+@(?:gmail|outlook|hotmail)\.com", text, re.I):
            errors.append(f"{rel}: personal email")
        if path.suffix == ".py":
            ast.parse(text, filename=str(rel))
        if path.suffix == ".md":
            for link in re.findall(r"\]\(([^)]+)\)", text):
                link = link.split("#")[0].strip("<>")
                if not link or "://" in link or link.startswith("mailto:"):
                    continue
                if not (path.parent / link).exists():
                    errors.append(f"{rel}: broken relative link {link}")
    if errors:
        print("\n".join(errors))
        return 1
    print("Repository hygiene: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
