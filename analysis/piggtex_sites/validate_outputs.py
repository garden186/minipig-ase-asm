#!/usr/bin/env python3
"""Independently validate GATK–PigGTEx exact-site integration outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Set, Tuple


REQUIRED_FILES = {
    "gatk_snp_piggtex_site_support.tsv.gz",
    "recurrent_gatk_snp_piggtex_site_support.tsv",
    "primary_block_concordant_gatk_snp_piggtex_site_support.tsv",
    "piggtex_site_support_stage_summary.tsv",
    "tissue_enriched_pair_gatk_piggtex_site_summary.tsv",
    "run_metadata.json",
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
        for row in reader:
            yield {key: value or "" for key, value in row.items()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_values(value: str, separator: str = ",") -> List[str]:
    return [item for item in value.split(separator) if item]


def row_key(row: Mapping[str, str]) -> Tuple[str, str]:
    return row["internal_tissue"], row["physical_snp_id"]


def expected_support_status(row: Mapping[str, str]) -> str:
    if int(row["piggtex_significant_rows_same_mapped_tissue"]) > 0:
        return "SIGNIFICANT_ASE_SITE_SAME_MAPPED_TISSUE"
    if int(row["piggtex_significant_rows_any_tissue"]) > 0:
        return "SIGNIFICANT_ASE_SITE_OTHER_TISSUE_ONLY"
    return "NOT_OBSERVED_IN_SIGNIFICANT_ONLY_SITE_TABLE"


def validate(output: Path) -> Dict[str, object]:
    missing = sorted(name for name in REQUIRED_FILES if not (output / name).is_file())
    if missing:
        raise ValidationError(f"missing outputs: {', '.join(missing)}")
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    rows: Dict[Tuple[str, str], Dict[str, str]] = {}
    errors: List[str] = []
    for index, row in enumerate(iter_tsv(output / "gatk_snp_piggtex_site_support.tsv.gz"), start=2):
        try:
            key = row_key(row)
            if key in rows:
                raise ValueError("duplicate internal tissue-SNP row")
            expected_id = (
                f"{row['contig']}_{int(row['position'])}_"
                f"{row['ref_allele']}_{row['alt_allele']}"
            )
            if row["physical_snp_id"] != expected_id:
                raise ValueError("physical_snp_id mismatch")
            if sorted((row["ref_allele"], row["alt_allele"])) != [row["allele_1"], row["allele_2"]]:
                raise ValueError("unordered allele pair mismatch")
            if len(split_values(row["internal_gatk_ase_samples"])) != int(row["n_internal_gatk_ase_measurements"]):
                raise ValueError("internal measurement/sample count mismatch")
            if row["piggtex_site_support_status"] != expected_support_status(row):
                raise ValueError("PigGTEx support status mismatch")
            if int(row["is_statistically_recurrent_gatk_snp_ase"]) and not row["gatk_recurrence_q_bh_global"]:
                raise ValueError("recurrent row lacks recurrence q")
            if int(row["n_primary_block_direction_concordant_contexts"]) > int(row["n_primary_block_gene_contexts"]):
                raise ValueError("concordant context count exceeds total")
            if int(row["n_primary_block_direction_discordant_contexts"]) > int(row["n_primary_block_gene_contexts"]):
                raise ValueError("discordant context count exceeds total")
            rows[key] = row
        except (ValueError, KeyError) as error:
            if len(errors) < 100:
                errors.append(f"row {index}: {error}")

    recurrent_keys = {key for key, row in rows.items() if int(row["is_statistically_recurrent_gatk_snp_ase"]) == 1}
    block_keys = {key for key, row in rows.items() if int(row["n_primary_block_direction_concordant_contexts"]) > 0}
    recurrent_file_keys = {row_key(row) for row in iter_tsv(output / "recurrent_gatk_snp_piggtex_site_support.tsv")}
    block_file_keys = {row_key(row) for row in iter_tsv(output / "primary_block_concordant_gatk_snp_piggtex_site_support.tsv")}

    stage_definitions = [
        ("all_independent_GATK_SNP_ASE_tissue_sites", lambda row: True),
        ("statistically_recurrent_GATK_SNP_ASE", lambda row: int(row["is_statistically_recurrent_gatk_snp_ase"]) == 1),
        ("primary_block_direction_concordant_GATK_SNP_ASE", lambda row: int(row["n_primary_block_direction_concordant_contexts"]) > 0),
        ("recurrent_and_primary_block_direction_concordant", lambda row: int(row["is_statistically_recurrent_gatk_snp_ase"]) == 1 and int(row["n_primary_block_direction_concordant_contexts"]) > 0),
        ("tissue_enriched_pair_primary_block_direction_concordant", lambda row: bool(row["tissue_enriched_gene_tissue_pairs"]) and int(row["n_primary_block_direction_concordant_contexts"]) > 0),
    ]
    expected_stages = {}
    for name, predicate in stage_definitions:
        selected = [row for row in rows.values() if predicate(row)]
        expected_stages[name] = (
            len(selected),
            len({(row["contig"], row["position"], row["allele_1"], row["allele_2"]) for row in selected}),
            sum(int(row["piggtex_significant_rows_any_tissue"]) > 0 for row in selected),
            sum(int(row["piggtex_significant_rows_same_mapped_tissue"]) > 0 for row in selected),
            sum(int(row["piggtex_significant_rows_any_tissue"]) == 0 for row in selected),
        )
    observed_stages = {}
    for row in iter_tsv(output / "piggtex_site_support_stage_summary.tsv"):
        observed_stages[row["stage"]] = tuple(
            int(row[field]) for field in (
                "n_internal_tissue_snp_sites", "n_distinct_physical_snps",
                "n_with_PigGTEx_significant_site_ASE_any_tissue",
                "n_with_PigGTEx_significant_site_ASE_same_mapped_tissue",
                "n_not_observed_in_PigGTEx_significant_only_site_table",
            )
        )

    pair_summary_rows = list(iter_tsv(output / "tissue_enriched_pair_gatk_piggtex_site_summary.tsv"))
    pair_summary_errors = 0
    for row in pair_summary_rows:
        pair = f"{row['gene_id']}|{row['target_tissue']}"
        selected = [
            site for site in rows.values()
            if pair in split_values(site["tissue_enriched_gene_tissue_pairs"], ";")
        ]
        if len(selected) != int(row["n_gatk_ase_sites_in_primary_blocks"]):
            pair_summary_errors += 1
        concordant_ids = set(split_values(row["direction_concordant_site_ids"]))
        same_ids = set(split_values(row["same_tissue_PigGTEx_supported_site_ids"]))
        if len(concordant_ids) != int(row["n_direction_concordant_gatk_ase_sites"]):
            pair_summary_errors += 1
        if len(same_ids) != int(row["n_direction_concordant_sites_with_PigGTEx_significant_ASE_same_mapped_tissue"]):
            pair_summary_errors += 1
        for site_id in same_ids:
            site = rows.get((row["target_tissue"], site_id))
            if site is None or int(site["piggtex_significant_rows_same_mapped_tissue"]) == 0:
                pair_summary_errors += 1

    counts = metadata["counts"]
    checks = {
        "required_outputs_present": True,
        "row_identity_and_support_status_invariants": not errors,
        "recurrent_subset_exact": recurrent_keys == recurrent_file_keys,
        "primary_block_concordant_subset_exact": block_keys == block_file_keys,
        "stage_summaries_reproduced": expected_stages == observed_stages,
        "pair_summary_basic_counts_reproduced": pair_summary_errors == 0 and len(pair_summary_rows) == 39,
        "metadata_GATK_measurement_count_matches": int(counts["n_gatk_ase_measurements"]) == sum(int(row["n_internal_gatk_ase_measurements"]) for row in rows.values()),
        "metadata_tissue_site_count_matches": int(counts["n_gatk_ase_tissue_snp_sites"]) == len(rows),
        "metadata_recurrent_count_matches": int(counts["n_recurrent_gatk_tissue_snp_sites"]) == len(recurrent_keys),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed": {
            "n_GATK_ASE_tissue_SNP_rows": len(rows),
            "n_recurrent_rows": len(recurrent_keys),
            "n_primary_block_concordant_rows": len(block_keys),
            "n_pair_summary_rows": len(pair_summary_rows),
            "n_pair_summary_errors": pair_summary_errors,
            "n_row_errors": len(errors),
        },
        "errors_first_100": errors,
        "output_fingerprints": {
            name: {"size_bytes": (output / name).stat().st_size, "sha256": sha256(output / name)}
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
