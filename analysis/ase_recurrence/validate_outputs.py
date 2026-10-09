#!/usr/bin/env python3
"""Independently validate recurrence output counts and subset relations."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Mapping, Set, Tuple


VERSION = "0.1.0"


def normalize_tissue(tissue: str) -> str:
    return {
        "blood": "Blood",
        "L-blood": "Blood",
        "Blood": "Blood",
        "Tenderlo": "Tenderloin",
        "tenderloin": "Tenderloin",
        "Tenderloin": "Tenderloin",
    }.get(tissue, tissue)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def count_and_keys(
    path: Path,
    fields: Tuple[str, ...],
    normalize_tissue_field: bool = False,
) -> Tuple[int, Set[Tuple[str, ...]]]:
    count = 0
    keys: Set[Tuple[str, ...]] = set()
    for row in rows(path):
        count += 1
        key = tuple(
            normalize_tissue(row[field])
            if normalize_tissue_field and field == "tissue"
            else row[field]
            for field in fields
        )
        if key in keys:
            raise ValueError(f"duplicate key in {path.name}: {key}")
        keys.add(key)
    return count, keys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--canonical-core", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report_path = args.report or args.out_dir / "validation.json"
    errors = []
    checks: Dict[str, object] = {}

    try:
        qc = json.loads((args.out_dir / "analysis_qc.json").read_text(encoding="utf-8"))
        if qc.get("status") != "PASS":
            errors.append("analysis_qc status is not PASS")
        if qc.get("version") != VERSION:
            errors.append(f"unexpected analysis version: {qc.get('version')}")
    except Exception as exc:
        errors.append(f"cannot read analysis_qc.json: {exc}")
        qc = {}

    summary_rows = {
        int(row["fragment_threshold"]): row
        for row in rows(args.out_dir / "threshold_summary.tsv")
    }
    for threshold in (15, 20, 30):
        prefix = f"threshold_{threshold}"
        summary = summary_rows.get(threshold)
        if summary is None:
            errors.append(f"missing threshold {threshold} summary")
            continue
        measurement_path = args.out_dir / f"{prefix}_eligible_gene_measurements.tsv.gz"
        recurrence_path = args.out_dir / f"{prefix}_recurrence_all.tsv.gz"
        recurrent_path = args.out_dir / f"{prefix}_statistically_recurrent.tsv"
        burden_path = args.out_dir / f"{prefix}_sample_tissue_burden.tsv"
        try:
            n_measurements, measurement_keys = count_and_keys(
                measurement_path, ("sample", "tissue", "gene_id")
            )
            n_core = 0
            eligible_sum = 0
            for row in rows(measurement_path):
                n_core += int(row["is_core_ase"])
                eligible_sum += 1
                if int(row["n_eligible_blocks"]) < 1:
                    errors.append(f"{prefix}: output contains a non-eligible measurement")
                    break
                if int(row["is_core_ase"]) and int(row["n_core_blocks"]) < 1:
                    errors.append(f"{prefix}: core flag without a core block")
                    break
            n_hypotheses, all_keys = count_and_keys(
                recurrence_path, ("gene_id", "tissue")
            )
            n_recurrent, recurrent_keys = count_and_keys(
                recurrent_path, ("gene_id", "tissue")
            )
            if not recurrent_keys.issubset(all_keys):
                errors.append(f"{prefix}: recurrent table is not a subset of all hypotheses")
            flagged_keys = set()
            for row in rows(recurrence_path):
                flag = int(row["is_statistically_recurrent_core"])
                expected = int(
                    int(row["n_samples_core"]) >= 2
                    and float(row["recurrence_q_bh_leave_one_gene_out"]) <= 0.05
                )
                if flag != expected:
                    errors.append(f"{prefix}: recurrence flag disagrees with primary rule")
                    break
                if flag:
                    flagged_keys.add((row["gene_id"], row["tissue"]))
            if flagged_keys != recurrent_keys:
                errors.append(f"{prefix}: recurrent subset does not match flagged hypotheses")
            burden_eligible = 0
            burden_core = 0
            burden_units = set()
            for row in rows(burden_path):
                unit = (row["sample"], row["tissue"])
                if unit in burden_units:
                    errors.append(f"{prefix}: duplicate burden unit {unit}")
                    break
                burden_units.add(unit)
                burden_eligible += int(row["n_eligible_genes"])
                burden_core += int(row["n_core_genes"])
            comparisons = {
                "n_eligible_gene_measurements": n_measurements,
                "n_core_gene_measurements": n_core,
                "n_recurrence_hypotheses": n_hypotheses,
                "n_recurrent_gene_tissue_pairs": n_recurrent,
            }
            for field, observed in comparisons.items():
                if int(summary[field]) != observed:
                    errors.append(
                        f"{prefix}: {field} summary={summary[field]} observed={observed}"
                    )
            if burden_eligible != n_measurements or burden_core != n_core:
                errors.append(
                    f"{prefix}: burden totals ({burden_eligible},{burden_core}) "
                    f"!= measurement totals ({n_measurements},{n_core})"
                )
            checks[prefix] = {**comparisons, "n_sample_tissue_units": len(burden_units)}
        except Exception as exc:
            errors.append(f"{prefix}: validation exception: {exc}")

    try:
        canonical_count, canonical_keys = count_and_keys(
            args.canonical_core,
            ("sample", "tissue", "gene_id"),
            normalize_tissue_field=True,
        )
        core15_keys = {
            (row["sample"], normalize_tissue(row["tissue"]), row["gene_id"])
            for row in rows(args.out_dir / "threshold_15_eligible_gene_measurements.tsv.gz")
            if row["is_core_ase"] == "1"
        }
        if canonical_keys != core15_keys:
            errors.append(
                f"canonical core mismatch: missing={len(canonical_keys-core15_keys)} "
                f"extra={len(core15_keys-canonical_keys)}"
            )
        checks["canonical_core"] = {
            "canonical": canonical_count,
            "reconstructed": len(core15_keys),
        }
    except Exception as exc:
        errors.append(f"canonical core validation exception: {exc}")

    try:
        matched_all = args.out_dir / "opportunity_matched_recurrence_all.tsv.gz"
        matched_recurrent = args.out_dir / "opportunity_matched_statistically_recurrent.tsv"
        n_matched, matched_keys = count_and_keys(matched_all, ("gene_id", "tissue"))
        n_matched_recurrent, matched_recurrent_keys = count_and_keys(
            matched_recurrent, ("gene_id", "tissue")
        )
        matched_flagged = {
            (row["gene_id"], row["tissue"])
            for row in rows(matched_all)
            if row["is_statistically_recurrent_core"] == "1"
        }
        if matched_flagged != matched_recurrent_keys:
            errors.append("matched recurrent subset does not match flags")
        if n_matched != int(summary_rows[15]["n_recurrence_hypotheses"]):
            errors.append("matched and primary hypothesis universes differ")
        checks["opportunity_matched"] = {
            "n_recurrence_hypotheses": n_matched,
            "n_recurrent_gene_tissue_pairs": n_matched_recurrent,
        }
    except Exception as exc:
        errors.append(f"matched validation exception: {exc}")

    tissue_counts = qc.get("tissue_animal_counts", {}) if isinstance(qc, Mapping) else {}
    expected_tissue_counts = {
        "Backfat": 9,
        "Blood": 5,
        "Brain": 10,
        "Heart": 10,
        "Kidney": 10,
        "Liver": 9,
        "Loin": 10,
        "Lung": 10,
        "Lymph": 10,
        "Spleen": 10,
        "Tenderloin": 10,
    }
    if tissue_counts != expected_tissue_counts:
        errors.append(f"unexpected tissue-animal counts: {tissue_counts}")

    report = {
        "status": "PASS" if not errors else "FAIL",
        "version": VERSION,
        "checks": checks,
        "errors": errors,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
