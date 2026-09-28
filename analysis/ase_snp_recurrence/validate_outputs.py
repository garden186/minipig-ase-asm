#!/usr/bin/env python3
"""Validate SNP-level recurrence outputs without loading them into memory."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterator, Sequence


class ValidationError(RuntimeError):
    pass


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle = opener("r", encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {key: value if value is not None else "" for key, value in row.items()}


def sample_count(value: str) -> int:
    return len([item for item in value.split(",") if item])


def validate(out_dir: Path) -> Dict[str, object]:
    required = [
        "eligible_gene_variant_measurements.tsv.gz",
        "sample_tissue_snp_burden.tsv",
        "sample_tissue_gene_snp_burden.tsv.gz",
        "snp_recurrence_all.tsv.gz",
        "statistically_recurrent_block_concordant_snp_ase.tsv",
        "analysis_qc.json",
        "run_metadata.json",
    ]
    for name in required:
        path = out_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValidationError(f"missing or empty output: {path}")

    with (out_dir / "run_metadata.json").open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    with (out_dir / "analysis_qc.json").open(encoding="utf-8") as handle:
        qc = json.load(handle)
    if metadata.get("status") != "PASS" or qc.get("status") != "PASS":
        raise ValidationError("metadata or analysis QC status is not PASS")
    fdr = float(metadata["statistics"]["fdr"])

    n_measurements = 0
    n_success_measurements = 0
    for row in iter_tsv(out_dir / "eligible_gene_variant_measurements.tsv.gz"):
        n_measurements += 1
        success = int(row["block_concordant_snp_level_ase"])
        site_signal = int(row["site_level_ase_signal"])
        if success and not site_signal:
            raise ValidationError("a block-concordant success is not a site-level ASE signal")
        n_success_measurements += success

    n_hypotheses = 0
    n_recurrent = 0
    n_ge2 = 0
    n_ge3 = 0
    recurrent_keys = set()
    for row in iter_tsv(out_dir / "snp_recurrence_all.tsv.gz"):
        n_hypotheses += 1
        n_testable = int(row["n_animals_testable"])
        n_success = int(row["n_animals_with_block_concordant_snp_ase"])
        q_value = float(row["recurrence_q_bh_global"])
        recurrent = int(row["is_statistically_recurrent_snp_ase"])
        if n_testable < 2 or not 0 <= n_success <= n_testable:
            raise ValidationError(f"invalid testable/success counts: {row}")
        if sample_count(row["animals_testable"]) != n_testable:
            raise ValidationError(f"testable sample-list mismatch: {row}")
        if sample_count(row["animals_with_block_concordant_snp_ase"]) != n_success:
            raise ValidationError(f"success sample-list mismatch: {row}")
        fraction = float(row["recurrence_fraction"])
        if not math.isclose(fraction, n_success / n_testable, rel_tol=1e-10, abs_tol=1e-12):
            raise ValidationError(f"recurrence fraction mismatch: {row}")
        expected_flag = int(n_success >= 2 and q_value <= fdr)
        if recurrent != expected_flag:
            raise ValidationError(f"recurrent flag mismatch: {row}")
        n_ge2 += int(n_success >= 2)
        n_ge3 += int(n_success >= 3)
        n_recurrent += recurrent
        if recurrent:
            recurrent_keys.add((row["gene_id"], row["tissue"], row["variant_id"]))

    significant_keys = set()
    for row in iter_tsv(out_dir / "statistically_recurrent_block_concordant_snp_ase.tsv"):
        if int(row["is_statistically_recurrent_snp_ase"]) != 1:
            raise ValidationError("non-recurrent row found in significant-only output")
        significant_keys.add((row["gene_id"], row["tissue"], row["variant_id"]))
    if recurrent_keys != significant_keys:
        raise ValidationError("significant-only keys do not equal flagged keys in all hypotheses")

    expected = metadata["counts"]
    observed = {
        "n_testable_gene_snp_measurements": n_measurements,
        "n_gene_linked_block_concordant_successes": n_success_measurements,
        "n_evaluable_hypotheses": n_hypotheses,
        "n_hypotheses_with_at_least_two_observed": n_ge2,
        "n_hypotheses_with_at_least_three_observed": n_ge3,
        "n_statistically_recurrent_gene_tissue_snp_hypotheses": n_recurrent,
    }
    for key, value in observed.items():
        if int(expected[key]) != value:
            raise ValidationError(
                f"metadata count mismatch for {key}: expected={expected[key]}, observed={value}"
            )

    result = {
        "status": "PASS",
        "checks": {
            "all_required_outputs_present": True,
            "success_subset_of_site_level_signal": True,
            "exact_testable_and_observed_sample_counts": True,
            "recurrent_definition_is_n_observed_ge2_and_global_q_le_fdr": True,
            "significant_output_matches_all_hypothesis_flags": True,
            "metadata_counts_reproduced": True,
        },
        "observed_counts": observed,
    }
    with (out_dir / "validation.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = validate(args.out_dir)
    except (ValidationError, KeyError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] validation: {result['observed_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
