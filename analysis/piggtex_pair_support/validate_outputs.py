#!/usr/bin/env python3
"""Validate 39-pair independent SNP support outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, Sequence


REQUIRED = [
    "tissue_enriched_pair_independent_snp_support.tsv",
    "supporting_independent_recurrent_snps.tsv",
    "support_stage_summary.tsv",
    "best_support_class_summary.tsv",
    "run_metadata.json",
]


def rows(path: Path) -> Iterator[Dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    errors = []
    checks = {}
    for name in REQUIRED:
        if not (args.output_dir / name).is_file():
            errors.append(f"missing output: {name}")
    checks["required_outputs_present"] = not errors
    if errors:
        return 1
    metadata = json.loads((args.output_dir / "run_metadata.json").read_text(encoding="utf-8"))
    pair_rows = list(rows(args.output_dir / "tissue_enriched_pair_independent_snp_support.tsv"))
    site_rows = list(rows(args.output_dir / "supporting_independent_recurrent_snps.tsv"))
    stage_rows = list(rows(args.output_dir / "support_stage_summary.tsv"))
    pair_keys = {(row["gene_id"], row["target_tissue"]) for row in pair_rows}
    checks["exactly_39_unique_pairs"] = len(pair_rows) == 39 and len(pair_keys) == 39
    if not checks["exactly_39_unique_pairs"]:
        errors.append("pair universe is not exactly 39 unique pairs")
    observed = {
        "primary_tissue_enriched_pairs": len(pair_rows),
        "pairs_with_any_coordinate_assigned_independent_recurrent_snp": sum(
            int(row["n_coordinate_assigned_independent_recurrent_snps"]) > 0 for row in pair_rows
        ),
        "pairs_with_unique_protein_coding_exonic_independent_recurrent_snp": sum(
            int(row["n_unique_protein_coding_exonic_independent_recurrent_snps"]) > 0
            for row in pair_rows
        ),
        "pairs_with_piggtex_same_tissue_exact_site_support": sum(
            int(row["n_unique_exonic_snps_with_piggtex_same_tissue_exact_support"]) > 0
            for row in pair_rows
        ),
        "pairs_with_direction_comparable_piggtex_exact_site": sum(
            int(row["n_unique_exonic_snps_with_direction_comparable_piggtex_support"]) > 0
            for row in pair_rows
        ),
        "pairs_with_direction_concordant_piggtex_exact_site": sum(
            int(row["n_unique_exonic_snps_with_direction_concordant_piggtex_support"]) > 0
            for row in pair_rows
        ),
        "coordinate_assigned_independent_recurrent_snp_rows": len(site_rows),
    }
    checks["metadata_counts_match"] = all(
        metadata["counts"][key] == value for key, value in observed.items()
    )
    if not checks["metadata_counts_match"]:
        errors.append("metadata counts mismatch")
    stage_counts = [int(row["n_of_39_pairs"]) for row in stage_rows]
    checks["stage_hierarchy_nonincreasing"] = all(
        left >= right for left, right in zip(stage_counts, stage_counts[1:])
    )
    if not checks["stage_hierarchy_nonincreasing"]:
        errors.append("stage hierarchy is not nonincreasing")
    site_pair_keys = {(row["target_gene_id"], row["target_tissue"]) for row in site_rows}
    checks["supporting_sites_belong_to_pair_universe"] = site_pair_keys <= pair_keys
    if not checks["supporting_sites_belong_to_pair_universe"]:
        errors.append("supporting site outside pair universe")
    classes = Counter(row["best_available_snp_support_class"] for row in pair_rows)
    checks["exclusive_support_classes_sum_to_39"] = sum(classes.values()) == 39
    output_fingerprints = {}
    for name in REQUIRED:
        path = args.output_dir / name
        output_fingerprints[name] = {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
    result = {
        "status": "PASS" if not errors else "FAIL",
        "checks": checks,
        "errors_first_100": errors[:100],
        "observed": observed,
        "best_support_class_counts": dict(sorted(classes.items())),
        "output_fingerprints": output_fingerprints,
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": result["status"], "observed": observed}, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
