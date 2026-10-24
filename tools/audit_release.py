"""Check the intended release file boundary, local document links and identifiers."""

import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = {".gitignore", "README.md", "LICENSE.md", "MANIFEST.in", "pyproject.toml"}
# These two are standard setuptools additions in a source-distribution round trip.
SDIST_FILES = {"PKG-INFO", "setup.cfg"}
DIRECTORIES = {"short_read_aligner", "docs", "tests", "tools", "validation", "benchmarks", "evidence", "examples"}
IGNORED = {".git", "output", "benchmark-output", "build", "dist", ".venv", "__pycache__"}


def main():
    files, errors = [], []
    for path in sorted(ROOT.rglob("*")):
        rel = path.relative_to(ROOT)
        if any(part in IGNORED or part.endswith(".egg-info") for part in rel.parts):
            continue
        if not path.is_file():
            continue
        if path.is_symlink():
            errors.append(f"unexpected symbolic link: {rel}")
        if (len(rel.parts) == 1 and rel.name not in ROOT_FILES | SDIST_FILES) or (len(rel.parts) > 1 and rel.parts[0] not in DIRECTORIES):
            errors.append(f"unexpected release path: {rel}")
        if path.suffix in (".pdf", ".tex", ".pkl", ".bam", ".bt2", ".npz", ".log", ".pyc"):
            errors.append(f"excluded artifact type: {rel}")
        text = path.read_text(encoding="utf-8")
        for pattern in (r"mail\.huji\.ac\.il", r"/cs[/]cbio/", r"\b76553\b", r"[A-Za-z]:[\\/](?:Users|CS|Genomics)[\\/]", r"\b(?:student[_ -]?id|teudat[_ -]?zehut)\s*[:=]\s*\d+"):
            if re.search(pattern, text, re.IGNORECASE):
                errors.append(f"identifier or institutional/local path pattern in {rel}")
        if path.suffix == ".md":
            for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", text):
                target = target.split("#", 1)[0]
                if not target or "://" in target:
                    continue
                if not (path.parent / target).resolve().is_file():
                    errors.append(f"broken link in {rel}: {target}")
        files.append({"path": rel.as_posix(), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    result = {"success": not errors, "files": len(files), "bytes": sum(row["bytes"] for row in files), "errors": errors,
              "scope": "intended source release only; generated packaging/git/temporary files excluded; pattern audit is not a proof of authorship or rights"}
    print(json.dumps(result, indent=2))
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
