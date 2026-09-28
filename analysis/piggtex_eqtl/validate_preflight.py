#!/usr/bin/env python3
"""Validate full-universe PigGTEx gene-eQTL preflight outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterator, Sequence


REQUIRED = [
    "full_independent_recurrent_snp_gene_eqtl_input.tsv.gz",
    "all_gene_assignment_eqtl_mapping_status.tsv.gz",
    "local_permutation_lead_eqtl_lookup.tsv.gz",
    "gene_tissue_lead_eqtl_summary.tsv",
    "local_eqtl_file_audit.tsv",
    "run_metadata.json",
]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def rows(path: Path) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
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
    candidate_rows = list(rows(args.output_dir / "full_independent_recurrent_snp_gene_eqtl_input.tsv.gz"))
    mapping_rows = list(rows(args.output_dir / "all_gene_assignment_eqtl_mapping_status.tsv.gz"))
    lead_rows = list(rows(args.output_dir / "local_permutation_lead_eqtl_lookup.tsv.gz"))
    pair_rows = list(rows(args.output_dir / "gene_tissue_lead_eqtl_summary.tsv"))
    audit_rows = list(rows(args.output_dir / "local_eqtl_file_audit.tsv"))
    candidate_keys = {
        (
            row["internal_gene_id"],
            row["internal_tissue"],
            row["internal_variant_unordered_key"],
        )
        for row in candidate_rows
    }
    checks["candidate_keys_unique"] = len(candidate_keys) == len(candidate_rows)
    if not checks["candidate_keys_unique"]:
        errors.append("candidate keys are not unique")
    lead_keys = {
        (
            row["internal_gene_id"],
            row["internal_tissue"],
            row["internal_variant_unordered_key"],
        )
        for row in lead_rows
    }
    checks["lead_lookup_exact_candidate_universe"] = lead_keys == candidate_keys
    if not checks["lead_lookup_exact_candidate_universe"]:
        errors.append("lead lookup universe differs from candidate input")
    checks["all_local_eqtl_files_present"] = all(row["file_present"] == "True" for row in audit_rows)
    if not checks["all_local_eqtl_files_present"]:
        errors.append("one or more mapped-tissue permutation eQTL files are missing")
    row_errors = 0
    for row in candidate_rows:
        if not row["piggtex_gene_id"] or not row["piggtex_tissue"]:
            row_errors += 1
        if row["internal_variant_unordered_key"] != "_".join(
            [row["chrom"], row["position"], row["allele_1"], row["allele_2"]]
        ):
            row_errors += 1
    checks["candidate_row_invariants"] = row_errors == 0
    if row_errors:
        errors.append(f"candidate row errors: {row_errors}")
    observed = {
        "gene_tissue_snp_assignment_rows": len(mapping_rows),
        "v100_mapped_gene_tissue_snp_candidate_rows": len(candidate_rows),
        "candidate_gene_tissue_pairs": len(pair_rows),
        "candidate_pairs_with_piggtex_gene_eqtl_record": sum(
            int(row["piggtex_gene_eqtl_record_present"]) for row in pair_rows
        ),
        "candidate_pairs_that_are_piggtex_eGenes": sum(
            row["piggtex_is_eGene"].upper() == "TRUE" for row in pair_rows
        ),
        "candidate_gene_tissue_snp_rows_matching_piggtex_lead_variant": sum(
            int(row["internal_snp_is_exact_piggtex_lead_variant"]) for row in lead_rows
        ),
        "candidate_gene_tissue_pairs_with_exact_lead_variant_match": sum(
            int(row["n_internal_snps_matching_piggtex_lead_variant"]) > 0 for row in pair_rows
        ),
    }
    checks["metadata_counts_match"] = all(
        metadata["counts"][key] == value for key, value in observed.items()
    )
    if not checks["metadata_counts_match"]:
        errors.append("metadata counts mismatch")
    outputs = {}
    for name in REQUIRED:
        path = args.output_dir / name
        outputs[name] = {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
    result = {
        "status": "PASS" if not errors else "FAIL",
        "checks": checks,
        "errors_first_100": errors[:100],
        "observed": observed,
        "output_fingerprints": outputs,
    }
    (args.output_dir / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": result["status"], "observed": observed}, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
