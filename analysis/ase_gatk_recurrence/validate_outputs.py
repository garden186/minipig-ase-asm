#!/usr/bin/env python3
"""Independently validate GATK SNP-level recurrence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Sequence, Set, Tuple


REQUIRED_FILES = {
    "gatk_snp_recurrence_all.tsv.gz",
    "statistically_recurrent_gatk_snp_ase.tsv",
    "sample_tissue_snp_burden.tsv",
    "physical_snp_summary.tsv.gz",
    "run_metadata.json",
}

REQUIRED_COLUMNS = {
    "tissue",
    "physical_snp_id",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "n_animals_testable",
    "animals_testable",
    "n_animals_with_snp_level_ase",
    "animals_with_snp_level_ase",
    "recurrence_fraction",
    "expected_animals_leave_one_site_out",
    "null_probabilities_leave_one_site_out",
    "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
}


class ValidationError(RuntimeError):
    """Raised when validation cannot proceed."""


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle = opener("r", encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise ValidationError(f"{path}: missing columns: {', '.join(missing)}")
        for row in reader:
            yield {key: value or "" for key, value in row.items()}


def split_values(value: str) -> List[str]:
    return [item for item in value.split(",") if item]


def poisson_binomial_tail(probabilities: Sequence[float], observed: int) -> float:
    if observed <= 0:
        return 1.0
    distribution = [1.0] + [0.0] * len(probabilities)
    populated = 0
    for probability in probabilities:
        for count in range(populated + 1, 0, -1):
            distribution[count] = (
                distribution[count] * (1.0 - probability)
                + distribution[count - 1] * probability
            )
        distribution[0] *= 1.0 - probability
        populated += 1
    return min(1.0, max(0.0, math.fsum(distribution[observed:])))


def close(a: float, b: float, tolerance: float = 2e-10) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(output: Path) -> Dict[str, object]:
    missing_files = sorted(name for name in REQUIRED_FILES if not (output / name).is_file())
    if missing_files:
        raise ValidationError(f"missing required outputs: {', '.join(missing_files)}")
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    fdr = float(metadata["parameters"]["fdr"])
    minimum_observed = int(metadata["parameters"]["minimum_observed_animals"])

    p_values: List[float] = []
    q_values: List[float] = []
    success_counts: List[int] = []
    flags: List[int] = []
    recurrent_keys: Set[Tuple[str, str]] = set()
    errors: List[str] = []

    for index, row in enumerate(iter_tsv(output / "gatk_snp_recurrence_all.tsv.gz")):
        try:
            n_testable = int(row["n_animals_testable"])
            n_success = int(row["n_animals_with_snp_level_ase"])
            testable = split_values(row["animals_testable"])
            positive = split_values(row["animals_with_snp_level_ase"])
            probabilities = [
                float(value)
                for value in split_values(row["null_probabilities_leave_one_site_out"])
            ]
            p_value = float(row["recurrence_p_poisson_binomial"])
            q_value = float(row["recurrence_q_bh_global"])
            flag = int(row["is_statistically_recurrent_snp_ase"])
            expected_id = (
                f"{row['contig']}_{int(row['position'])}_"
                f"{row['ref_allele']}_{row['alt_allele']}"
            )
            if row["physical_snp_id"] != expected_id:
                raise ValueError("physical_snp_id mismatch")
            if len(testable) != n_testable or len(probabilities) != n_testable:
                raise ValueError("testable/probability count mismatch")
            if len(positive) != n_success or not set(positive).issubset(testable):
                raise ValueError("positive animal count/subset mismatch")
            if not close(float(row["recurrence_fraction"]), n_success / n_testable):
                raise ValueError("recurrence fraction mismatch")
            if not close(
                float(row["expected_animals_leave_one_site_out"]), math.fsum(probabilities)
            ):
                raise ValueError("expected count mismatch")
            if not close(p_value, poisson_binomial_tail(probabilities, n_success)):
                raise ValueError("Poisson-binomial P mismatch")
            p_values.append(p_value)
            q_values.append(q_value)
            success_counts.append(n_success)
            flags.append(flag)
            if flag:
                recurrent_keys.add((row["tissue"], row["physical_snp_id"]))
        except (ValueError, KeyError) as error:
            if len(errors) < 100:
                errors.append(f"row {index + 2}: {error}")

    n = len(p_values)
    expected_q = [1.0] * n
    running = 1.0
    for rank_from_end, index in enumerate(
        sorted(range(n), key=lambda i: (p_values[i], i), reverse=True), start=1
    ):
        rank = n - rank_from_end + 1
        running = min(running, min(1.0, p_values[index] * n / rank))
        expected_q[index] = running
    n_q_mismatch = sum(not close(observed, expected) for observed, expected in zip(q_values, expected_q))
    n_flag_mismatch = sum(
        flag != int(success >= minimum_observed and q_value <= fdr)
        for flag, success, q_value in zip(flags, success_counts, q_values)
    )

    significant_keys = {
        (row["tissue"], row["physical_snp_id"])
        for row in iter_tsv(output / "statistically_recurrent_gatk_snp_ase.tsv")
    }
    metadata_counts = metadata["counts"]
    observed = {
        "n_evaluable_hypotheses": n,
        "n_statistically_recurrent_hypotheses": len(recurrent_keys),
        "n_q_mismatches": n_q_mismatch,
        "n_flag_mismatches": n_flag_mismatch,
        "n_row_errors": len(errors),
    }
    checks = {
        "required_outputs_present": True,
        "row_invariants_and_poisson_binomial_p_reproduced": not errors,
        "global_bh_reproduced": n_q_mismatch == 0,
        "recurrent_definition_reproduced": n_flag_mismatch == 0,
        "significant_subset_exact": significant_keys == recurrent_keys,
        "metadata_hypothesis_count_matches": int(metadata_counts["n_evaluable_hypotheses"]) == n,
        "metadata_recurrent_count_matches": int(
            metadata_counts["n_statistically_recurrent_hypotheses"]
        ) == len(recurrent_keys),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed": observed,
        "errors_first_100": errors,
        "output_fingerprints": {
            name: {
                "size_bytes": (output / name).stat().st_size,
                "sha256": sha256(output / name),
            }
            for name in sorted(REQUIRED_FILES)
        },
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.out_dir.resolve())
    except (ValidationError, OSError, ValueError, KeyError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    path = args.out_dir.resolve() / "validation.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[{result['status']}] {path}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
