#!/usr/bin/env python3
"""Independently validate cohort-aware GATK SNP-level ASE outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, Sequence

from call_gatk_snp_level_ase import (
    AnalysisError,
    bh_adjust,
    exact_binom_two_sided_half_counts,
)


class ValidationError(RuntimeError):
    """Raised when an output invariant is violated."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    if path.suffix == ".gz":
        context = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        context = path.open("r", encoding="utf-8", newline="")
    with context as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def same_float(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=5e-10, abs_tol=1e-300)


def row_key(row: Dict[str, str]) -> tuple[str, str, str, str, str, str]:
    return (
        row["sample"],
        row["tissue"],
        row["contig"],
        row["position"],
        row["ref_allele"],
        row["alt_allele"],
    )


def validate(out_dir: Path, output: Path) -> dict[str, object]:
    required = [
        "cohort_audit.tsv",
        "effective_cohort.tsv",
        "site_tests.tsv.gz",
        "snp_level_ase_calls.tsv.gz",
        "sample_tissue_summary.tsv",
        "run_metadata.json",
    ]
    for name in required:
        path = out_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValidationError(f"missing or empty output: {path}")

    metadata = json.loads((out_dir / "run_metadata.json").read_text(encoding="utf-8"))
    if metadata.get("status") != "PASS":
        raise ValidationError("run metadata status is not PASS")
    definition = metadata["call_definition"]
    minimum_total = int(definition["minimum_total_count"])
    minimum_fraction = float(definition["minimum_major_allele_fraction"])
    fdr = float(definition["maximum_q_bh_within_animal_tissue"])

    cohort_audit = list(iter_tsv(out_dir / "cohort_audit.tsv"))
    effective = list(iter_tsv(out_dir / "effective_cohort.tsv"))
    effective_jobs = {row["job_id"] for row in effective}
    if not effective_jobs or len(effective_jobs) != len(effective):
        raise ValidationError("effective cohort is empty or has duplicate jobs")
    if any(row["inclusion_status"] != "INCLUDED" for row in effective):
        raise ValidationError("effective cohort contains an excluded row")
    audit_included = {
        row["job_id"] for row in cohort_audit if row["inclusion_status"] == "INCLUDED"
    }
    if effective_jobs != audit_included:
        raise ValidationError("effective cohort does not match cohort-audit inclusion flags")

    summary = {row["job_id"]: row for row in iter_tsv(out_dir / "sample_tissue_summary.tsv")}
    if set(summary) != effective_jobs:
        raise ValidationError("summary jobs do not match the effective cohort")

    calls = iter(iter_tsv(out_dir / "snp_level_ase_calls.tsv.gz"))
    next_call = next(calls, None)
    total_tests = 0
    total_calls = 0
    observed_summary: Dict[str, Dict[str, int]] = {}
    call_direction: Counter[str] = Counter()
    monoallelic_calls: Counter[str] = Counter()

    def validate_group(job_id: str, rows: list[Dict[str, str]]) -> None:
        nonlocal next_call, total_calls
        if job_id in observed_summary:
            raise ValidationError(f"non-contiguous repeated BH family: {job_id}")
        pvalues = [float(row["p_exact_two_sided"]) for row in rows]
        expected_qvalues = bh_adjust(pvalues)
        n_fdr = 0
        n_effect = 0
        n_calls = 0
        for row, expected_q in zip(rows, expected_qvalues):
            observed_q = float(row["q_bh_within_animal_tissue"])
            if not same_float(expected_q, observed_q):
                raise ValidationError(f"BH q mismatch: {row_key(row)}")
            expected_fdr = int(observed_q <= fdr)
            expected_effect = int(float(row["major_allele_fraction"]) >= minimum_fraction)
            expected_call = expected_fdr * expected_effect
            if int(row["is_fdr_significant"]) != expected_fdr:
                raise ValidationError(f"FDR flag mismatch: {row_key(row)}")
            if int(row["is_strong_effect"]) != expected_effect:
                raise ValidationError(f"effect flag mismatch: {row_key(row)}")
            if int(row["is_snp_level_ase"]) != expected_call:
                raise ValidationError(f"call flag mismatch: {row_key(row)}")
            if int(row["bh_family_size"]) != len(rows):
                raise ValidationError(f"BH family-size mismatch: {row_key(row)}")
            n_fdr += expected_fdr
            n_effect += expected_effect
            n_calls += expected_call
            if expected_call:
                if next_call is None or row_key(next_call) != row_key(row):
                    raise ValidationError(f"significant-only output mismatch: {row_key(row)}")
                if next_call != row:
                    raise ValidationError(
                        f"significant-only row fields differ from all-test row: {row_key(row)}"
                    )
                call_direction[row["imbalance_direction"]] += 1
                monoallelic_calls["ref_zero"] += int(int(row["ref_count"]) == 0)
                monoallelic_calls["alt_zero"] += int(int(row["alt_count"]) == 0)
                next_call = next(calls, None)
                total_calls += 1
        observed_summary[job_id] = {
            "n_testable_sites": len(rows),
            "n_fdr_significant_sites": n_fdr,
            "n_strong_effect_sites": n_effect,
            "n_snp_level_ase_calls": n_calls,
        }

    current_job = ""
    current_rows: list[Dict[str, str]] = []
    for row in iter_tsv(out_dir / "site_tests.tsv.gz"):
        job_id = row["job_id"]
        if job_id not in effective_jobs:
            raise ValidationError(f"site test outside effective cohort: {job_id}")
        if current_job and job_id != current_job:
            validate_group(current_job, current_rows)
            current_rows = []
        current_job = job_id
        if int(row["total_count"]) < minimum_total:
            raise ValidationError(f"site below depth threshold: {row_key(row)}")
        ref_count = int(row["ref_count"])
        alt_count = int(row["alt_count"])
        total_count = int(row["total_count"])
        if ref_count + alt_count != total_count:
            raise ValidationError(f"count invariant failed: {row_key(row)}")
        expected_p = exact_binom_two_sided_half_counts(
            total_count, min(ref_count, alt_count)
        )
        if not same_float(expected_p, float(row["p_exact_two_sided"])):
            raise ValidationError(f"exact p mismatch: {row_key(row)}")
        current_rows.append(row)
        total_tests += 1
    if current_job:
        validate_group(current_job, current_rows)

    if set(observed_summary) != effective_jobs:
        raise ValidationError("at least one retained unit has no testable SNPs")
    if next_call is not None:
        raise ValidationError("significant-only output contains an extra row")

    for job_id, observed in observed_summary.items():
        expected = summary[job_id]
        for field, value in observed.items():
            if int(expected[field]) != value:
                raise ValidationError(
                    f"summary mismatch for {job_id} {field}: {expected[field]} != {value}"
                )

    expected_counts = metadata["counts"]
    if int(expected_counts["testable_sites"]) != total_tests:
        raise ValidationError("metadata testable-site count mismatch")
    if int(expected_counts["snp_level_ase_calls"]) != total_calls:
        raise ValidationError("metadata SNP-level ASE call count mismatch")
    if int(metadata["cohort"]["retained_units"]) != len(effective_jobs):
        raise ValidationError("metadata retained-unit count mismatch")

    caller_path = Path(__file__).with_name("call_gatk_snp_level_ase.py")
    validator_path = Path(__file__)
    fingerprints = {}
    for name in required:
        path = out_dir / name
        fingerprints[name] = {
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }

    directional_total = call_direction["REF"] + call_direction["ALT"]
    result: dict[str, object] = {
        "status": "PASS",
        "checks": {
            "required_outputs_present": True,
            "effective_cohort_matches_audit": True,
            "all_test_rows_belong_to_effective_cohort": True,
            "exact_binomial_p_values_reproduced": True,
            "within_animal_tissue_bh_reproduced": True,
            "effect_and_call_flags_reproduced": True,
            "significant_only_output_matches_call_flags": True,
            "summary_and_metadata_counts_reproduced": True,
        },
        "observed": {
            "manifest_units": len(cohort_audit),
            "retained_units": len(effective_jobs),
            "testable_sites": total_tests,
            "snp_level_ase_calls": total_calls,
            "call_direction": dict(sorted(call_direction.items())),
            "ref_fraction_among_directional_calls": (
                call_direction["REF"] / directional_total if directional_total else None
            ),
            "monoallelic_calls": dict(sorted(monoallelic_calls.items())),
        },
        "source_code": {
            "caller": {"path": str(caller_path.resolve()), "sha256": sha256(caller_path)},
            "validator": {
                "path": str(validator_path.resolve()),
                "sha256": sha256(validator_path),
            },
        },
        "output_fingerprints": fingerprints,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output = args.output or args.out_dir / "validation.json"
    try:
        result = validate(args.out_dir, output)
    except (OSError, ValueError, KeyError, AnalysisError, ValidationError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] {result['observed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
