#!/usr/bin/env python3
"""Validate ASE v0.2.5 output structure, counts and candidate invariants."""

from __future__ import print_function

import argparse
import csv
import gzip
import json
import os
from collections import Counter, defaultdict


VERSION = "0.2.5"
PRIMARY_CLASSES = {"PRIMARY_ROBUST", "PRIMARY_SUPPORTED"}
ALL_CLASSES = PRIMARY_CLASSES | {
    "SECONDARY_SINGLE_SNP", "NOT_PRIMARY"
}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def read_rows(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def truthy(value):
    return str(value).lower() in {"1", "true", "yes"}


def block_key(row):
    return (
        row["sample"],
        row["tissue"],
        row["bam"],
        row["haplotype_block_id"],
        row["gene_id"],
    )


def gene_key(row):
    return (
        row["sample"], row["tissue"], row["gene_id"]
    )


def parser():
    p = argparse.ArgumentParser(
        description="Validate ASE v0.2.5 outputs"
    )
    p.add_argument("--out-dir", required=True)
    p.add_argument("--expect-sample")
    p.add_argument("--require-exact-phaser-match", type=int, default=1)
    p.add_argument("--out-json")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    audit = os.path.join(args.out_dir, "audit")
    required = {
        "blocks": os.path.join(
            args.out_dir, "ase_haplotype_block_results.tsv"
        ),
        "primary_blocks": os.path.join(
            args.out_dir, "ase_primary_haplotype_blocks.tsv"
        ),
        "genes": os.path.join(args.out_dir, "ase_gene_summary.tsv"),
        "primary_genes": os.path.join(
            args.out_dir, "ase_primary_genes.tsv"
        ),
        "counts": os.path.join(audit, "ase_strand_assignment_qc.tsv"),
        "measurements": os.path.join(
            audit, "ase_measurement_unique.tsv"
        ),
        "count_validation": os.path.join(
            audit, "ase_count_validation.tsv"
        ),
        "wgs_snp_qc": os.path.join(audit, "ase_snp_qc.tsv"),
        "phase_qc": os.path.join(
            audit, "ase_local_phase_graph_qc.tsv"
        ),
        "qc_json": os.path.join(args.out_dir, "ase_v025_qc.json"),
    }
    missing = [
        path for path in required.values() if not os.path.isfile(path)
    ]
    if missing:
        raise SystemExit(
            "[ERROR] missing required outputs: {}".format(
                ",".join(missing)
            )
        )

    errors = []
    warnings = []
    tables = {}
    headers = {}
    for key, path in required.items():
        if key == "qc_json":
            continue
        tables[key], headers[key] = read_rows(path)

    blocks = tables["blocks"]
    block_by_key = {}
    measurement_ids = Counter()
    for row in blocks:
        key = block_key(row)
        if key in block_by_key:
            errors.append("duplicate_block_gene_row:{}".format(key))
        block_by_key[key] = row
        measurement_ids[row["measurement_id"]] += 1
        candidate_class = row["candidate_class"]
        if candidate_class not in ALL_CLASSES:
            errors.append(
                "unknown_candidate_class:{}:{}".format(
                    key, candidate_class
                )
            )
        is_primary = truthy(row["is_primary_candidate"])
        if is_primary != (candidate_class in PRIMARY_CLASSES):
            errors.append("primary_flag_mismatch:{}".format(key))
        pvalue = as_float(row.get("ase_p_exact"))
        qvalue = as_float(row.get("ase_q_bh"))
        if pvalue is not None and not 0 < pvalue <= 1:
            errors.append("invalid_pvalue:{}".format(key))
        if qvalue is not None and not 0 < qvalue <= 1:
            errors.append("invalid_qvalue:{}".format(key))
        if candidate_class in PRIMARY_CLASSES:
            required_primary = {
                "minor_support_class": "ADEQUATE_MINOR",
                "local_phase_status": "PASS",
            }
            for field, expected in required_primary.items():
                if row.get(field) != expected:
                    errors.append(
                        "primary_{}_mismatch:{}:{}!={}".format(
                            field, key, row.get(field), expected
                        )
                    )
            if as_int(row.get("n_wgs_balanced_gene_snps")) < 2:
                errors.append("primary_has_fewer_than_two_wgs_snps:{}".format(
                    key
                ))
            expected_support = {
                "PRIMARY_ROBUST": "MULTI_SNP_LOO_ROBUST",
                "PRIMARY_SUPPORTED": "MULTI_SNP_SUPPORTED",
            }[candidate_class]
            if row.get("snp_support_status") != expected_support:
                errors.append("primary_snp_support_mismatch:{}".format(key))
        if args.expect_sample and row["sample"] != args.expect_sample:
            errors.append(
                "unexpected_sample:{}:{}".format(key, row["sample"])
            )

    duplicate_measurements = [
        item for item, count in measurement_ids.items() if count != 1
    ]
    if duplicate_measurements:
        errors.append(
            "nonunique_measurement_ids:{}".format(
                ",".join(duplicate_measurements[:20])
            )
        )

    primary_keys = {block_key(row) for row in tables["primary_blocks"]}
    expected_primary_keys = {
        key for key, row in block_by_key.items()
        if row["candidate_class"] in PRIMARY_CLASSES
    }
    if primary_keys != expected_primary_keys:
        errors.append("primary_block_subset_mismatch")

    count_by_key = {}
    for row in tables["counts"]:
        key = block_key(row)
        if key in count_by_key:
            errors.append("duplicate_count_row:{}".format(key))
        count_by_key[key] = row
        for hap in ("a", "b"):
            strict = as_int(row["strict_{}_count".format(hap)])
            partition = (
                as_int(row["direct_{}_count".format(hap)])
                + as_int(row["strand_resolved_{}".format(hap)])
            )
            if strict != partition:
                errors.append(
                    "{}_count_conservation:{}:{}!={}".format(
                        hap, key, strict, partition
                    )
                )
        if as_int(row["strict_total_count"]) != (
            as_int(row["direct_total_count"])
            + as_int(row["strand_resolved_total"])
        ):
            errors.append("total_count_conservation:{}".format(key))

    for key, row in block_by_key.items():
        count = count_by_key.get(key)
        if count is None:
            errors.append("block_without_count_row:{}".format(key))
            continue
        for hap in ("a", "b"):
            if as_int(row["ase_{}_count".format(hap)]) != as_int(
                count["strict_{}_count".format(hap)]
            ):
                errors.append(
                    "block_count_mismatch:{}:{}".format(key, hap)
                )
        if as_int(row["ase_total_count"]) != as_int(
            count["strict_total_count"]
        ):
            errors.append("block_total_mismatch:{}".format(key))

    measurement_rows = tables["measurements"]
    if len(measurement_rows) != len(blocks):
        errors.append(
            "measurement_row_count:{}!={}".format(
                len(measurement_rows), len(blocks)
            )
        )

    validation_rows = tables["count_validation"]
    exact_matches = sum(
        truthy(row.get("count_match_exact")) for row in validation_rows
    )
    if (
        args.require_exact_phaser_match
        and exact_matches != len(validation_rows)
    ):
        errors.append(
            "phaser_reconstruction_mismatch:{}of{}".format(
                len(validation_rows) - exact_matches,
                len(validation_rows),
            )
        )

    wgs_status_counts = Counter(
        row["wgs_snp_status"] for row in tables["wgs_snp_qc"]
    )
    unknown_wgs = set(wgs_status_counts).difference(
        {"WGS_BALANCED", "WGS_UNBALANCED", "WGS_INSUFFICIENT"}
    )
    if unknown_wgs:
        errors.append(
            "unknown_wgs_snp_status:{}".format(
                ",".join(sorted(unknown_wgs))
            )
        )

    genes_by_key = {gene_key(row): row for row in tables["genes"]}
    if len(genes_by_key) != len(tables["genes"]):
        errors.append("duplicate_gene_summary_rows")
    primary_gene_keys = {
        gene_key(row) for row in tables["primary_genes"]
    }
    expected_primary_gene_keys = {
        key for key, row in genes_by_key.items()
        if row["gene_ase_status"] in {
            "ASE_GENE_ROBUST", "ASE_GENE_SUPPORTED"
        }
    }
    if primary_gene_keys != expected_primary_gene_keys:
        errors.append("primary_gene_subset_mismatch")

    fragment_path = os.path.join(
        audit, "ase_fragment_assignment.tsv.gz"
    )
    fragment_rows = 0
    fragment_counts = defaultdict(
        lambda: {
            "strict_a": set(), "strict_b": set(),
            "direct_a": set(), "direct_b": set(),
            "strand_a": set(), "strand_b": set(),
        }
    )
    if os.path.isfile(fragment_path):
        fragments, _ = read_rows(fragment_path)
        fragment_rows = len(fragments)
        assigned_seen = set()
        for row in fragments:
            if row["final_status"] != "ASSIGNED":
                continue
            key = block_key(row)
            identity = key + (row["qname"],)
            if identity in assigned_seen:
                errors.append(
                    "duplicate_assigned_fragment:{}".format(identity)
                )
            assigned_seen.add(identity)
            hap = row["haplotype"].lower()
            fragment_counts[key]["strict_" + hap].add(row["qname"])
            method = row["gene_assignment_method"]
            if method == "DIRECT":
                fragment_counts[key]["direct_" + hap].add(row["qname"])
                if truthy(row["strand_resolution_required"]):
                    errors.append(
                        "direct_marked_strand_required:{}".format(identity)
                    )
            elif method == "STRAND_RESOLVED":
                fragment_counts[key]["strand_" + hap].add(row["qname"])
                if not truthy(row["strand_resolution_required"]):
                    errors.append(
                        "strand_resolved_not_marked:{}".format(identity)
                    )
            else:
                errors.append(
                    "invalid_assigned_gene_method:{}:{}".format(
                        identity, method
                    )
                )
        for key, count in count_by_key.items():
            observed = fragment_counts.get(key, {})
            expected = {
                "strict_a": as_int(count["strict_a_count"]),
                "strict_b": as_int(count["strict_b_count"]),
                "direct_a": as_int(count["direct_a_count"]),
                "direct_b": as_int(count["direct_b_count"]),
                "strand_a": as_int(count["strand_resolved_a"]),
                "strand_b": as_int(count["strand_resolved_b"]),
            }
            for label, expected_count in expected.items():
                if len(observed.get(label, set())) != expected_count:
                    errors.append(
                        "fragment_count_mismatch:{}:{}:{}!={}".format(
                            key, label,
                            len(observed.get(label, set())),
                            expected_count,
                        )
                    )
    else:
        warnings.append(
            "fragment audit skipped: audit/ase_fragment_assignment.tsv.gz "
            "absent"
        )

    with open(required["qc_json"]) as handle:
        qc = json.load(handle)
    if qc.get("version") != VERSION:
        errors.append(
            "qc_version_mismatch:{}!={}".format(
                qc.get("version"), VERSION
            )
        )
    if qc.get("n_haplotype_block_result_rows") != len(blocks):
        errors.append("qc_block_row_count_mismatch")
    if qc.get("n_primary_haplotype_block_rows") != len(
        tables["primary_blocks"]
    ):
        errors.append("qc_primary_block_row_count_mismatch")

    report = {
        "version": VERSION,
        "out_dir": os.path.abspath(args.out_dir),
        "status": "PASS" if not errors else "FAIL",
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "errors_first_100": errors[:100],
        "warnings": warnings,
        "row_counts": {
            key: len(value) for key, value in tables.items()
        },
        "n_exact_phaser_matches": exact_matches,
        "n_fragment_rows": fragment_rows,
        "wgs_snp_status_counts": dict(sorted(wgs_status_counts.items())),
        "candidate_class_counts": dict(sorted(Counter(
            row["candidate_class"] for row in blocks
        ).items())),
    }
    out_json = args.out_json or os.path.join(
        args.out_dir, "ase_v025_validation.json"
    )
    with open(out_json, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
