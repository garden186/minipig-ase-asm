#!/usr/bin/env python3
"""Rebuild and exactly validate ASE v0.2.6 cohort recurrence outputs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Mapping, Sequence

from build_ase_cohort import (
    DEFAULT_SAMPLES,
    DEFAULT_TISSUE_ALIASES,
    OUTPUT_SPECS,
    CohortError,
    build_expected,
    normalize_tissue,
    parse_csv_list,
    parse_tissue_aliases,
)


VERSION = "ase_candidates_v0.2.6_cohort_validator_v0.3.0"


def compare_tsv(
    path: Path,
    header: Sequence[str],
    expected: Sequence[Mapping[str, object]],
) -> int:
    if not path.is_file():
        raise CohortError(f"missing cohort output: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(header):
            raise CohortError(
                f"unexpected output header in {path.name}; "
                f"expected={list(header)}, observed={reader.fieldnames}"
            )
        count = 0
        for count, (observed, expected_row) in enumerate(
            zip(reader, expected), start=1
        ):
            canonical_expected = {
                field: str(expected_row.get(field, "")) for field in header
            }
            if observed != canonical_expected:
                raise CohortError(
                    f"{path.name}: row {count} differs; "
                    f"observed={observed}, expected={canonical_expected}"
                )
        else:
            observed_extra = next(reader, None)
            if observed_extra is not None:
                raise CohortError(
                    f"{path.name}: observed output has extra rows after {count}"
                )
        if count != len(expected):
            raise CohortError(
                f"{path.name}: row count differs; observed={count}, "
                f"expected={len(expected)}"
            )
    return count


def validate_recurrence_invariants(
    recurrence: Sequence[Mapping[str, object]],
    recurrent: Sequence[Mapping[str, object]],
    q_threshold: float,
) -> None:
    expected_recurrent = []
    previous_p = None
    previous_q = None
    ordered_by_p = sorted(
        recurrence,
        key=lambda row: (
            float(row["recurrence_p"]),
            str(row["gene_id"]),
            str(row["tissue"]),
        ),
    )
    for row in recurrence:
        n_testable = int(row["n_samples_testable"])
        n_core = int(row["n_samples_core"])
        p_value = float(row["recurrence_p"])
        q_value = float(row["recurrence_q_bh"])
        flag = int(row["is_statistically_recurrent_core"])
        if n_testable < 2 or not 0 <= n_core <= n_testable:
            raise CohortError(f"invalid recurrence counts: {row}")
        if not 0.0 <= p_value <= 1.0 or not 0.0 <= q_value <= 1.0:
            raise CohortError(f"invalid recurrence probability: {row}")
        expected_flag = int(n_core >= 2 and q_value <= q_threshold)
        if flag != expected_flag:
            raise CohortError(f"invalid recurrence flag: {row}")
        if flag:
            expected_recurrent.append(dict(row))
    if list(recurrent) != expected_recurrent:
        raise CohortError(
            "statistically recurrent table is not the exact recurrence=1 subset"
        )
    for row in ordered_by_p:
        p_value = float(row["recurrence_p"])
        q_value = float(row["recurrence_q_bh"])
        if previous_p is not None and p_value >= previous_p and q_value + 1e-12 < previous_q:
            raise CohortError("BH q-values are not monotone in p-value order")
        previous_p = p_value
        previous_q = q_value


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate ASE v0.2.6 cohort")
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--samples", default=",".join(DEFAULT_SAMPLES))
    parser.add_argument("--exclude-tissues", default="Cranial")
    parser.add_argument(
        "--tissue-aliases",
        default=",".join(
            f"{source}={target}" for source, target in DEFAULT_TISSUE_ALIASES.items()
        ),
    )
    parser.add_argument("--expect-sample-tissue-combinations", type=int, default=0)
    parser.add_argument("--expect-tissues", type=int, default=0)
    parser.add_argument("--recurrence-q-threshold", type=float, default=0.05)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    samples = parse_csv_list(args.samples)
    aliases = parse_tissue_aliases(args.tissue_aliases)
    excluded = {
        normalize_tissue(tissue, aliases)
        for tissue in parse_csv_list(args.exclude_tissues)
    }
    report = {
        "version": VERSION,
        "status": "FAIL",
        "input_root": str(args.input_root.resolve()),
        "out_dir": str(args.out_dir.resolve()),
        "samples": samples,
        "errors": [],
        "warnings": [],
    }
    try:
        expected, metadata = build_expected(
            args.input_root,
            samples,
            excluded,
            tissue_aliases=aliases,
            expect_sample_tissue_combinations=args.expect_sample_tissue_combinations,
            expect_tissues=args.expect_tissues,
            recurrence_q_threshold=args.recurrence_q_threshold,
        )
        observed_counts = {}
        for name, header in OUTPUT_SPECS:
            observed_counts[name] = compare_tsv(
                args.out_dir / name, header, expected[name]
            )
        qc_path = args.out_dir / "cohort_qc.json"
        if not qc_path.is_file():
            raise CohortError(f"missing cohort QC: {qc_path}")
        with qc_path.open("r", encoding="utf-8") as handle:
            observed_qc = json.load(handle)
        if observed_qc != metadata:
            raise CohortError("cohort_qc.json differs from exact recomputation")
        validate_recurrence_invariants(
            expected["gene_tissue_core_recurrence_test_results.tsv"],
            expected["statistically_recurrent_core_gene_tissue.tsv"],
            args.recurrence_q_threshold,
        )
        report.update(metadata)
        report["status"] = "PASS"
        report["observed_row_counts"] = observed_counts
        warning_samples = [
            sample
            for sample, entry in metadata["source_validation"].items()
            if int(entry.get("n_warnings") or 0) > 0
        ]
        if warning_samples:
            report["warnings"] = [
                "source validation warnings retained for: "
                + ",".join(warning_samples)
            ]
    except Exception as exc:
        report["errors"] = [str(exc)]
    report["n_errors"] = len(report["errors"])
    report["n_warnings"] = len(report["warnings"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "cohort_validation.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
