#!/usr/bin/env python3
"""Independently validate GATK-to-primary-block linkage outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Sequence, Tuple


REQUIRED_FILES = {
    "primary_block_gatk_snp_validation.tsv.gz",
    "primary_block_validation_summary.tsv",
    "primary_gene_measurement_validation_summary.tsv",
    "validation_summary_by_tissue.tsv",
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


def expected_status(testable: int, ase: int, concordant: int, discordant: int) -> str:
    if testable == 0:
        return "NOT_GATK_TESTABLE"
    if ase == 0:
        return "GATK_TESTABLE_NO_ASE"
    if concordant > 0 and discordant == 0:
        return "DIRECTION_CONCORDANT_ONLY"
    if discordant > 0 and concordant == 0:
        return "DIRECTION_DISCORDANT_ONLY"
    if concordant > 0 and discordant > 0:
        return "MIXED_DIRECTION"
    raise ValidationError("invalid summary counts")


def aggregate(rows: Sequence[Mapping[str, str]], level: str) -> Dict[Tuple[str, ...], Counter]:
    result: Dict[Tuple[str, ...], Counter] = defaultdict(Counter)
    for row in rows:
        if level == "block":
            key = (row["sample"], row["tissue"], row["gene_id"], row["haplotype_block_id"])
        elif level == "measurement":
            key = (row["sample"], row["tissue"], row["gene_id"])
        else:
            raise ValueError(level)
        counter = result[key]
        counter["n_contexts"] += 1
        testable = int(row["gatk_site_testable"])
        ase = int(row["gatk_snp_level_ase"])
        status = row["gatk_block_direction_status"]
        counter["n_testable"] += testable
        counter["n_ase"] += ase
        counter["n_concordant"] += int(status == "DIRECTION_CONCORDANT")
        counter["n_discordant"] += int(status == "DIRECTION_DISCORDANT")
        counter["n_recurrent"] += int(row["gatk_site_is_statistically_recurrent"])
    return dict(result)


def validate(output: Path) -> Dict[str, object]:
    missing = sorted(name for name in REQUIRED_FILES if not (output / name).is_file())
    if missing:
        raise ValidationError(f"missing outputs: {', '.join(missing)}")
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    site_rows = list(iter_tsv(output / "primary_block_gatk_snp_validation.tsv.gz"))
    errors: List[str] = []
    site_counts = Counter()
    seen_contexts = set()
    for index, row in enumerate(site_rows, start=2):
        try:
            key = (
                row["sample"], row["tissue"], row["gene_id"],
                row["haplotype_block_id"], row["variant_id"],
            )
            if key in seen_contexts:
                raise ValueError("duplicate block gene-variant context")
            seen_contexts.add(key)
            ref, alt = row["ref_allele"], row["alt_allele"]
            hap_a, hap_b = row["haplotype_a_allele"], row["haplotype_b_allele"]
            if {hap_a, hap_b} != {ref, alt}:
                raise ValueError("haplotype alleles do not match REF/ALT")
            direction = row["v026_ase_direction"]
            predicted = hap_a if direction == "HAP_A" else hap_b if direction == "HAP_B" else None
            minor = hap_b if direction == "HAP_A" else hap_a if direction == "HAP_B" else None
            if predicted is None or row["predicted_major_allele"] != predicted:
                raise ValueError("predicted major allele mismatch")
            if row["predicted_minor_allele"] != minor:
                raise ValueError("predicted minor allele mismatch")
            testable = int(row["gatk_site_testable"])
            ase = int(row["gatk_snp_level_ase"])
            status = row["gatk_block_direction_status"]
            if testable == 0:
                expected = "NOT_GATK_TESTABLE"
            elif ase == 0:
                expected = "GATK_TESTABLE_NOT_ASE"
            elif row["gatk_major_allele"] == predicted:
                expected = "DIRECTION_CONCORDANT"
            elif row["gatk_major_allele"] == minor:
                expected = "DIRECTION_DISCORDANT"
            else:
                raise ValueError("GATK major allele outside block alleles")
            if status != expected:
                raise ValueError("direction status mismatch")
            recurrent = int(row["gatk_site_is_statistically_recurrent"])
            if recurrent and (ase != 1 or not row["gatk_recurrence_q_bh_global"]):
                raise ValueError("recurrent context lacks GATK ASE or recurrence q")
            site_counts["n_testable"] += testable
            site_counts["n_ase"] += ase
            site_counts["n_concordant"] += int(status == "DIRECTION_CONCORDANT")
            site_counts["n_discordant"] += int(status == "DIRECTION_DISCORDANT")
            site_counts["n_recurrent"] += recurrent
        except (ValueError, KeyError) as error:
            if len(errors) < 100:
                errors.append(f"row {index}: {error}")

    block_agg = aggregate(site_rows, "block")
    measurement_agg = aggregate(site_rows, "measurement")
    block_mismatches = 0
    for row in iter_tsv(output / "primary_block_validation_summary.tsv"):
        key = (row["sample"], row["tissue"], row["gene_id"], row["haplotype_block_id"])
        observed = block_agg.pop(key, None)
        if observed is None:
            block_mismatches += 1
            continue
        expected_values = (
            int(row["n_gene_variants"]), int(row["n_gatk_testable_variants"]),
            int(row["n_gatk_snp_level_ase_variants"]),
            int(row["n_direction_concordant_variants"]),
            int(row["n_direction_discordant_variants"]),
            int(row["n_recurrent_gatk_ase_variants"]),
        )
        calculated = tuple(observed[name] for name in (
            "n_contexts", "n_testable", "n_ase", "n_concordant", "n_discordant", "n_recurrent"
        ))
        if expected_values != calculated or row["gatk_block_validation_status"] != expected_status(
            observed["n_testable"], observed["n_ase"], observed["n_concordant"], observed["n_discordant"]
        ):
            block_mismatches += 1
    block_mismatches += len(block_agg)

    measurement_mismatches = 0
    status_counts = Counter()
    n_measurement_rows = 0
    for row in iter_tsv(output / "primary_gene_measurement_validation_summary.tsv"):
        n_measurement_rows += 1
        key = (row["sample"], row["tissue"], row["gene_id"])
        observed = measurement_agg.pop(key, None)
        if observed is None:
            measurement_mismatches += 1
            continue
        expected_values = (
            int(row["n_gene_variant_contexts"]), int(row["n_gatk_testable_variants"]),
            int(row["n_gatk_snp_level_ase_variants"]),
            int(row["n_direction_concordant_variants"]),
            int(row["n_direction_discordant_variants"]),
            int(row["n_recurrent_gatk_ase_variants"]),
        )
        calculated = tuple(observed[name] for name in (
            "n_contexts", "n_testable", "n_ase", "n_concordant", "n_discordant", "n_recurrent"
        ))
        expected_label = expected_status(
            observed["n_testable"], observed["n_ase"], observed["n_concordant"], observed["n_discordant"]
        )
        status_counts[expected_label] += 1
        if expected_values != calculated or row["gatk_measurement_validation_status"] != expected_label:
            measurement_mismatches += 1
    measurement_mismatches += len(measurement_agg)

    counts = metadata["counts"]
    checks = {
        "required_outputs_present": True,
        "site_level_allele_and_direction_invariants": not errors,
        "block_summaries_reproduced": block_mismatches == 0,
        "measurement_summaries_reproduced": measurement_mismatches == 0,
        "metadata_site_context_count_matches": int(counts["n_core_gene_variant_contexts"]) == len(site_rows),
        "metadata_block_count_matches": int(counts["n_core_block_summary_rows"]) == len(seen_contexts) - len(site_rows) + len({(r['sample'],r['tissue'],r['gene_id'],r['haplotype_block_id']) for r in site_rows}),
        "metadata_measurement_count_matches": int(counts["n_primary_measurement_summary_rows"]) == n_measurement_rows,
        "metadata_testable_count_matches": int(counts["n_gatk_testable_core_contexts"]) == site_counts["n_testable"],
        "metadata_ase_count_matches": int(counts["n_gatk_ase_core_contexts"]) == site_counts["n_ase"],
        "metadata_direction_counts_match": (
            int(counts["n_direction_concordant_core_contexts"]) == site_counts["n_concordant"]
            and int(counts["n_direction_discordant_core_contexts"]) == site_counts["n_discordant"]
        ),
        "metadata_recurrent_context_count_matches": int(counts["n_recurrent_gatk_ase_core_contexts"]) == site_counts["n_recurrent"],
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed": {
            "n_site_context_rows": len(site_rows),
            "n_primary_measurement_rows": n_measurement_rows,
            "n_block_summary_mismatches": block_mismatches,
            "n_measurement_summary_mismatches": measurement_mismatches,
            "n_site_errors": len(errors),
            "site_counts": dict(site_counts),
            "measurement_validation_status_counts": dict(sorted(status_counts.items())),
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
