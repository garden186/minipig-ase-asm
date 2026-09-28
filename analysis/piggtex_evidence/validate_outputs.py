#!/usr/bin/env python3
"""Independently validate the canonical SNP/GATK/PigGTEx evidence audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Set, Tuple


REQUIRED_FILES = {
    "current_646_pair_site_evidence_audit.tsv.gz",
    "current_39_pair_site_evidence_audit.tsv",
    "current_39_pair_site_evidence_summary.tsv",
    "canonical_recurrent_snp_gatk_piggtex_audit.tsv",
    "high_pip_full_site_resolved_candidates.tsv",
    "thyn1_site_resolved_audit.tsv",
    "evidence_stage_summary.tsv",
    "run_metadata.json",
}


class ValidationError(RuntimeError):
    """Raised when output validation cannot proceed."""


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {key: value or "" for key, value in row.items()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_csv(value: str) -> List[str]:
    return [item for item in value.split(",") if item]


def site_key(row: Mapping[str, str]) -> Tuple[str, str, str, int, str, str]:
    return (
        row["gene_id"], row["internal_tissue"], row["chrom"], int(row["position"]),
        row["allele_1"], row["allele_2"],
    )


def is_full_high_pip(row: Mapping[str, str]) -> bool:
    return (
        int(row["canonical_recurrent_gene_tissue_snp"]) == 1
        and int(row["n_primary_block_direction_concordant_samples"]) > 0
        and int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0
        and int(row["piggtex_exact_significant_eqtl"]) == 1
        and int(row["piggtex_exact_finemapped_eqtl_pip_ge_0_5"]) == 1
    )


def expected_support_status(row: Mapping[str, str]) -> str:
    if int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0:
        return "SIGNIFICANT_ASE_SITE_SAME_MAPPED_TISSUE"
    if int(row["piggtex_significant_site_rows_any_tissue"]) > 0:
        return "SIGNIFICANT_ASE_SITE_OTHER_TISSUE_ONLY"
    return "NOT_OBSERVED_IN_SIGNIFICANT_ONLY_SITE_TABLE"


def scope_counts(name: str, rows: List[Mapping[str, str]]) -> Dict[str, int | str]:
    return {
        "scope": name,
        "n_sites": len(rows),
        "n_canonical_recurrent_snp": sum(int(row["canonical_recurrent_gene_tissue_snp"]) for row in rows),
        "n_with_primary_block_gatk_ase": sum(int(row["n_primary_block_gatk_ase_samples"]) > 0 for row in rows),
        "n_with_primary_block_direction_concordant_gatk_ase": sum(int(row["n_primary_block_direction_concordant_samples"]) > 0 for row in rows),
        "n_independent_gatk_recurrent": sum(int(row["independent_gatk_statistically_recurrent"]) for row in rows),
        "n_with_piggtex_significant_site_ase_same_mapped_tissue": sum(int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0 for row in rows),
        "n_with_exact_significant_eqtl": sum(int(row["piggtex_exact_significant_eqtl"]) for row in rows),
        "n_with_exact_finemapped_eqtl": sum(int(row["piggtex_exact_finemapped_eqtl"]) for row in rows),
        "n_with_exact_finemapped_eqtl_pip_ge_0_5": sum(int(row["piggtex_exact_finemapped_eqtl_pip_ge_0_5"]) for row in rows),
        "n_with_exact_finemapped_eqtl_pip_ge_0_9": sum(int(row["piggtex_exact_finemapped_eqtl_pip_ge_0_9"]) for row in rows),
        "n_full_site_resolved_candidates_pip_ge_0_5": sum(is_full_high_pip(row) for row in rows),
    }


def validate(output: Path) -> Dict[str, object]:
    missing = sorted(name for name in REQUIRED_FILES if not (output / name).is_file())
    if missing:
        raise ValidationError(f"missing outputs: {', '.join(missing)}")
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    site_rows = list(iter_tsv(output / "current_646_pair_site_evidence_audit.tsv.gz"))
    site_by_key = {}
    row_errors = []
    for index, row in enumerate(site_rows, start=2):
        try:
            key = site_key(row)
            if key in site_by_key:
                raise ValueError("duplicate site key")
            if row["variant_unordered_key"] != f"{row['chrom']}_{int(row['position'])}_{row['allele_1']}_{row['allele_2']}":
                raise ValueError("variant_unordered_key mismatch")
            if sorted((row["allele_1"], row["allele_2"])) != [row["allele_1"], row["allele_2"]]:
                raise ValueError("alleles are not in unordered-key sort order")
            if int(row["current_recurrent_gene_tissue_pair"]) != 1:
                raise ValueError("646-pair audit contains an out-of-scope pair")
            if row["piggtex_site_ase_support_status"] != expected_support_status(row):
                raise ValueError("PigGTEx significant-only support status mismatch")
            if len(split_csv(row["primary_block_direction_concordant_samples"])) != int(row["n_primary_block_direction_concordant_samples"]):
                raise ValueError("direction-concordant sample count mismatch")
            if len(split_csv(row["independent_gatk_ase_samples"])) != int(row["n_independent_gatk_ase_samples"]):
                raise ValueError("independent GATK sample count mismatch")
            site_by_key[key] = row
        except (ValueError, KeyError) as error:
            if len(row_errors) < 100:
                row_errors.append(f"row {index}: {error}")

    rows39 = list(iter_tsv(output / "current_39_pair_site_evidence_audit.tsv"))
    expected39 = {key for key, row in site_by_key.items() if int(row["current_tissue_enriched_pair"]) == 1}
    observed39 = {site_key(row) for row in rows39}

    high_rows = list(iter_tsv(output / "high_pip_full_site_resolved_candidates.tsv"))
    expected_high = {key for key, row in site_by_key.items() if is_full_high_pip(row)}
    observed_high = {site_key(row) for row in high_rows}

    thyn1_rows = list(iter_tsv(output / "thyn1_site_resolved_audit.tsv"))
    expected_thyn1 = {key for key, row in site_by_key.items() if row["gene_name"] == "THYN1"}
    observed_thyn1 = {site_key(row) for row in thyn1_rows}

    pair_summaries = list(iter_tsv(output / "current_39_pair_site_evidence_summary.tsv"))
    rows39_by_pair = defaultdict(list)
    for row in rows39:
        rows39_by_pair[(row["gene_id"], row["internal_tissue"])].append(row)
    pair_errors = []
    for row in pair_summaries:
        pair = (row["gene_id"], row["target_tissue"])
        selected = rows39_by_pair.get(pair, [])
        checks = {
            "n_server_audit_candidate_sites": len(selected),
            "n_canonical_recurrent_snp_sites": sum(int(site["canonical_recurrent_gene_tissue_snp"]) for site in selected),
            "n_sites_with_primary_block_direction_concordant_gatk_ase": sum(int(site["n_primary_block_direction_concordant_samples"]) > 0 for site in selected),
            "n_independent_gatk_recurrent_sites": sum(int(site["independent_gatk_statistically_recurrent"]) for site in selected),
            "n_sites_with_piggtex_significant_ase_same_mapped_tissue": sum(int(site["piggtex_significant_site_rows_same_mapped_tissue"]) > 0 for site in selected),
            "n_sites_with_exact_significant_eqtl": sum(int(site["piggtex_exact_significant_eqtl"]) for site in selected),
            "n_full_site_resolved_candidates_pip_ge_0_5": sum(is_full_high_pip(site) for site in selected),
        }
        for field, expected in checks.items():
            if int(row[field]) != expected:
                pair_errors.append(f"{pair}: {field}={row[field]} expected {expected}")

    canonical_rows = list(iter_tsv(output / "canonical_recurrent_snp_gatk_piggtex_audit.tsv"))
    canonical_keys = {
        (row["gene_id"], row["tissue"], row["chrom"], int(row["position"]), *sorted((row["ref"], row["alt"])))
        for row in canonical_rows
    }

    stage_rows = {row["scope"]: row for row in iter_tsv(output / "evidence_stage_summary.tsv")}
    expected_scopes = {
        "current_646_pair_server_audit_site_universe": scope_counts("current_646_pair_server_audit_site_universe", site_rows),
        "current_39_tissue_enriched_pair_site_universe": scope_counts("current_39_tissue_enriched_pair_site_universe", rows39),
        "current_646_pairs_canonical_recurrent_SNP_subset": scope_counts(
            "current_646_pairs_canonical_recurrent_SNP_subset",
            [row for row in site_rows if int(row["canonical_recurrent_gene_tissue_snp"]) == 1],
        ),
    }
    stage_errors = []
    for scope, expected in expected_scopes.items():
        observed = stage_rows.get(scope)
        if observed is None:
            stage_errors.append(f"missing stage: {scope}")
            continue
        for field, value in expected.items():
            if field != "scope" and int(observed[field]) != int(value):
                stage_errors.append(f"{scope}: {field}={observed[field]} expected {value}")

    counts = metadata["counts"]
    boundaries = metadata["canonical_boundaries"]
    checks = {
        "required_outputs_present": True,
        "site_row_identity_and_status_invariants": not row_errors,
        "current_646_pair_count_reproduced": int(counts["n_current_recurrent_pairs"]) == 646,
        "current_39_pair_count_reproduced": int(counts["n_current_tissue_enriched_pairs"]) == 39 and len(pair_summaries) == 39,
        "canonical_1619_SNP_rows_unique": len(canonical_rows) == 1619 and len(canonical_keys) == 1619,
        "current_39_site_subset_exact": expected39 == observed39,
        "high_PIP_candidate_definition_exact": expected_high == observed_high,
        "THYN1_subset_exact": expected_thyn1 == observed_thyn1,
        "pair_summary_counts_reproduced": not pair_errors and set(rows39_by_pair) == {(row["gene_id"], row["target_tissue"]) for row in pair_summaries},
        "site_stage_summaries_reproduced": not stage_errors,
        "server_candidate_external_key_count_reproduced": int(counts["n_server_candidate_unique_external_gene_tissue_site_keys"]) == 10345,
        "independent_GATK_is_not_declared_canonical_replacement": boundaries["independent_GATK_recurrence_is_validation_not_replacement"] is True,
        "retired_eight_pair_statistics_not_used": boundaries["retired_eight_pair_internal_statistics_used"] is False,
        "external_provenance_limitation_recorded": "not independently rescanned" in metadata["external_provenance_boundary"]["PigGTEx_gene_eQTL"],
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed": {
            "n_current_646_pair_site_rows": len(site_rows),
            "n_current_39_pair_site_rows": len(rows39),
            "n_current_39_pair_summaries": len(pair_summaries),
            "n_canonical_recurrent_SNP_rows": len(canonical_rows),
            "n_high_PIP_full_site_resolved_candidates": len(high_rows),
            "n_THYN1_site_rows": len(thyn1_rows),
            "n_row_errors": len(row_errors),
            "n_pair_summary_errors": len(pair_errors),
            "n_stage_errors": len(stage_errors),
        },
        "errors_first_100": row_errors + pair_errors[:100] + stage_errors[:100],
        "output_fingerprints": {
            name: {"size_bytes": (output / name).stat().st_size, "sha256": sha256(output / name)}
            for name in sorted(REQUIRED_FILES)
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.out_dir.resolve())
    except (ValidationError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    path = args.out_dir.resolve() / "validation.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[{result['status']}] {path}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
