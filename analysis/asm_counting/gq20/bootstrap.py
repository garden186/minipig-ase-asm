"""Load the fragment-counting implementation."""

from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
FROZEN = ROOT.parent / "engine"
sys.path.insert(0, str(FROZEN))
from common import (
    read_json,
    write_json as _write_json,
    rows,
    write_table,
    sha,
    signature,
    stamp,
    CALLS,
    STATES,
)

ANALYSIS_ID = "gq20"
POLICY = {
    "min_gq": 20,
    "minimum_haplotype_depth": None,
    "mask": "Retain the upstream genotype-disruption exclusion unchanged",
    "phase": "Exact within-animal phase set; never pool disconnected phase sets",
    "depth": "Distinct paired fragments with valid M/U calls at one CpG",
    "inference": "Descriptive only; no ASM tests or GQ30 comparison",
}


def write_json(path, value):
    """Retry transient Windows sync/indexer locks; preserve other IO failures."""
    for attempt in range(6):
        try:
            return _write_json(path, value)
        except PermissionError as exc:
            if getattr(exc, "winerror", None) not in (5, 32, 33) or attempt == 5:
                raise
            time.sleep(0.1 * (attempt + 1))


def verify_release():
    import sys

    sys.path.insert(0, str(ROOT.parents[1]))
    from source_identity import code_identity

    return code_identity(ROOT.parent)


def contained(root, relative):
    """Resolve a relative source/output member without allowing path escape."""
    from pathlib import PurePosixPath

    p = PurePosixPath(relative)
    if p.is_absolute() or ".." in p.parts or "\\" in relative or ":" in relative:
        raise ValueError("Unsafe member path: " + relative)
    result = root.joinpath(*p.parts).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError("Member escaped its root")
    return result
