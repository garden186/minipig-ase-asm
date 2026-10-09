"""I/O, provenance, and fixed scientific contracts for the single-unit recount."""
from __future__ import annotations
import contextlib
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import time

PACKAGE = Path(__file__).resolve().parent
VERSION = "0.1.0"
CHROMS = [str(i) for i in range(1, 19)]
PATTERNS = [f"{h}_nM{i}" for h in ("H1", "H2") for i in range(6)]
FIELDS = ["sample", "tissue", "chrom", "k", "window_id", "first_cpg_index",
          "last_cpg_index", "cpg_positions_1based", "start0", "end0", "span_bp", "ps",
          *PATTERNS, "H1_all", "H2_all", "H1_U", "H1_X", "H1_M", "H2_U", "H2_X", "H2_M",
          "H1_UM", "H2_UM", "T_UM", "minor_UM", "H1_M_fraction", "H2_M_fraction",
          "delta_M", "abs_delta_M", "within_span_200", "depth_T10m4",
          "testable_T10m4_span200"]
CONTRACT = {
    "assembly": "Sscrofa11.1",
    "chromosomes": CHROMS, "k": 4, "reference_consecutive_CpGs": True,
    "states": "q35: nM0/1=U; nM2=X; nM3/4=M; X is excluded from UM denominators",
    "delta_M": "H1_M/(H1_U+H1_M) - H2_M/(H2_U+H2_M); NA if either denominator is zero",
    "retention": "All observed window-by-PS rows, including low depth and spans over 200 bp",
    "testability_annotation": "span<=200 AND total_UM>=10 AND minor_UM>=4; not a row filter",
    "effect_filter": None, "p_values": False, "multiple_testing": False,
    "assignment": "Unique allele support allowing C->T on XG=CT and G->A on XG=GA before mate merge",
    "ambiguous_link_SNPs": "Global C/T and A/G exclusion retained",
    "min_mapq": 0, "base_quality_filter": None, "min_link_snps": 1, "pair_wait_bp": 2000,
    "phase_scope": "One exact WGS PS and one non-conflicting haplotype per fragment",
    "genotype_CpG_disruption": "Unmodified v0.3.3 helper applied to variant-only DeepVariant VCF",
    "unobserved_rows": "Not materialized; absence is not a zero or a negative ASM observation",
    "inference_scope": "Counts and effect estimates only; no final ASM calls or tissue inference",
}


def timestamp():
    return dt.datetime.now(dt.timezone(dt.timedelta(hours=9))).strftime("%Y%m%d_%H%M%S_%f")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def fingerprint(path):
    p = Path(path)
    s = p.stat()
    if not p.is_file() or not s.st_size or not os.access(p, os.R_OK):
        raise ValueError("Missing, empty, or unreadable input: " + str(p))
    return {"path": str(p.resolve()), "size_bytes": s.st_size, "mtime_ns": s.st_mtime_ns,
            "inode": s.st_ino}


def check_inputs(request):
    current = {k: fingerprint(v["path"]) for k, v in request["inputs"].items()}
    if current != request["inputs"]:
        raise ValueError("Input fingerprint changed. Do not reuse results after modifying inputs.")


def verify_package():
    manifest = read_json(PACKAGE / "package_manifest.json")
    for rel, expected in manifest["sha256"].items():
        if digest(PACKAGE / rel) != expected:
            raise ValueError("Package checksum mismatch: " + rel)
    return signature(manifest)


@contextlib.contextmanager
def table_writer(path, fields):
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as txt:
                w = csv.DictWriter(txt, fieldnames=fields, delimiter="\t", lineterminator="\n")
                w.writeheader()
                yield w
    os.replace(tmp, path)


def table_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def run_capture(command, stderr_path, timeout=120):
    with Path(stderr_path).open("w", encoding="utf-8") as log:
        p = subprocess.run(command, stdout=subprocess.PIPE, stderr=log,
                           text=True, timeout=timeout, check=False)
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}); see {stderr_path}")
    return p.stdout


def stream(command, stderr_path):
    """Drain stdout with stderr on disk; terminate a producer on early consumer exit."""
    with Path(stderr_path).open("w", encoding="utf-8") as err:
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err,
                             text=True, bufsize=1024 * 1024)
        exhausted = False
        try:
            assert p.stdout is not None
            yield from p.stdout
            exhausted = True
        finally:
            if p.stdout:
                p.stdout.close()
            if not exhausted and p.poll() is None:
                p.terminate()
            try:
                rc = p.wait(timeout=15)
            except subprocess.TimeoutExpired:
                p.kill()
                rc = p.wait()
            if exhausted and rc:
                raise RuntimeError(f"Command failed ({rc}); see {stderr_path}")


def rss_kib():
    try:
        import resource
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except ImportError:
        return 0


def derived(vector, span):
    if len(vector) != 12 or vector[5] or vector[11] or any(n < 0 for n in vector):
        raise ValueError("Invalid CpG4 pattern vector")
    row = {}
    for h, offset in (("H1", 0), ("H2", 6)):
        v = vector[offset:offset + 5]
        u, x, m = v[0] + v[1], v[2], v[3] + v[4]
        row.update({h + "_all": sum(v), h + "_U": u, h + "_X": x, h + "_M": m,
                    h + "_UM": u + m, h + "_M_fraction": m / (u + m) if u + m else "NA"})
    a, b = row["H1_UM"], row["H2_UM"]
    delta = row["H1_M_fraction"] - row["H2_M_fraction"] if a and b else "NA"
    row.update(T_UM=a + b, minor_UM=min(a, b), delta_M=delta,
               abs_delta_M=abs(delta) if delta != "NA" else "NA",
               within_span_200=int(span <= 200), depth_T10m4=int(a + b >= 10 and min(a, b) >= 4),
               testable_T10m4_span200=int(span <= 200 and a + b >= 10 and min(a, b) >= 4))
    return row
