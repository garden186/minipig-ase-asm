"""Fingerprint source files used by an analysis run."""

import hashlib
import json
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parent


def code_identity(*directories):
    files = set(ANALYSIS.glob("*.py"))
    for directory in directories:
        for path in Path(directory).rglob("*"):
            if path.is_file() and path.suffix in {".py", ".R", ".json"}:
                files.add(path)
    hashes = {}
    for path in sorted(files):
        with path.open("rb") as handle:
            hashes[path.relative_to(ANALYSIS).as_posix()] = hashlib.file_digest(
                handle, "sha256"
            ).hexdigest()
    return hashlib.sha256(
        json.dumps(hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
