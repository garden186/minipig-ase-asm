#!/usr/bin/env python3
"""Independent structural and numerical validator for ASE v0.2.6 outputs."""

import argparse
import csv
import json
import math
import os
import sqlite3
import sys
import tempfile
from collections import Counter
from pathlib import Path


CORE_CLASSES = {"ASE_CORE_ROBUST", "ASE_CORE_SUPPORTED"}
SUPPORTED_CLASSES = {
    "ASE_SUPPORTED_SINGLE_SNP",
    "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED",
    "ASE_SUPPORT_UNRESOLVED",
}
SENSITIVITY_CLASSES = {
    "ASE_FDR_MODERATE_EFFECT",
    "ASE_NEAR_FDR",
    "ASE_NOMINAL_DIRECTIONAL",
    "ASE_STATISTICAL_SMALL_EFFECT",
}

REQUIRED_OUTPUTS = [
    "ase_haplotype_block_reclassified.tsv",
    "ase_core_haplotype_blocks.tsv",
    "ase_supported_haplotype_blocks.tsv",
    "ase_sensitivity_haplotype_blocks.tsv",
    "ase_gene_summary.tsv",
    "ase_core_genes.tsv",
    "ase_supported_genes.tsv",
    "ase_sensitivity_genes.tsv",
    "ase_v026_qc.json",
    os.path.join("audit", "ase_v025_to_v026_transition.tsv"),
    os.path.join("audit", "ase_reclassification_reason_counts.tsv"),
]


def number(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "."}:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def integer(value):
    result = number(value)
    if result is None or not result.is_integer():
        return None
    return int(result)


def expected_signal(row, thresholds):
    status = row.get("ase_test_status", "")
    if status == "LOW_COVERAGE":
        return "LOW_COVERAGE"
    if status != "TESTED":
        return "NOT_TESTABLE"
    qvalue = number(row.get("ase_q_bh"))
    pvalue = number(row.get("ase_p_exact"))
    effect = number(row.get("major_haplotype_fraction"))
    if qvalue is None or effect is None:
        return "NOT_TESTABLE"
    if qvalue <= thresholds["fdr"]:
        if effect >= thresholds["strong_major_fraction"]:
            return "FDR_STRONG_EFFECT"
        if effect >= thresholds["moderate_major_fraction"]:
            return "FDR_MODERATE_EFFECT"
        return "FDR_SMALL_EFFECT"
    if qvalue <= thresholds["near_fdr"] and effect >= thresholds["strong_major_fraction"]:
        return "NEAR_FDR_STRONG_EFFECT"
    if (
        pvalue is not None
        and pvalue <= thresholds["nominal_p"]
        and effect >= thresholds["strong_major_fraction"]
    ):
        return "NOMINAL_STRONG_EFFECT"
    return "NO_SELECTED_SIGNAL"


def expected_shape(row):
    a_count = integer(row.get("ase_a_count"))
    b_count = integer(row.get("ase_b_count"))
    if a_count is None or b_count is None:
        return "NOT_ASSESSED"
    minor = min(a_count, b_count)
    if minor <= 1:
        return "MONOALLELIC_LIKE"
    if minor == 2:
        return "LOW_MINOR"
    return "TWO_HAPLOTYPE_SUPPORTED"


def expected_technical(row, thresholds):
    if row.get("analysis_set") == "EXCLUDED_TISSUE":
        return "EXCLUDED"
    reasons = []
    if row.get("analysis_set") != "MAIN":
        reasons.append("UNKNOWN_ANALYSIS_SET")
    if row.get("count_validation") != "MATCH":
        reasons.append("COUNT_MISMATCH")
    balanced = integer(row.get("n_wgs_balanced_gene_snps"))
    if balanced is None or balanced <= 0:
        reasons.append("NO_WGS_BALANCED_SNP")
    conflict = number(row.get("fragment_hap_conflict_fraction"))
    if conflict is None:
        reasons.append("FRAGMENT_CONFLICT_NOT_ASSESSED")
    elif conflict > thresholds["max_fragment_hap_conflict_fraction"]:
        reasons.append("FRAGMENT_CONFLICT")
    if not row.get("gene_id") or not row.get("haplotype_block_id"):
        reasons.append("MISSING_BLOCK_OR_GENE_ID")
    return "FAIL" if reasons else "PASS"


def expected_evidence(technical, signal, support, phase):
    if technical == "EXCLUDED":
        return "EXCLUDED_TISSUE"
    if technical != "PASS":
        return "TECHNICAL_FAIL"
    if signal == "LOW_COVERAGE":
        return "LOW_COVERAGE"
    if signal == "NOT_TESTABLE":
        return "NOT_TESTABLE"
    if signal == "FDR_STRONG_EFFECT":
        if support == "MULTI_SNP_LOO_ROBUST" and phase == "PASS":
            return "ASE_CORE_ROBUST"
        if support == "MULTI_SNP_SUPPORTED" and phase == "PASS":
            return "ASE_CORE_SUPPORTED"
        if support in {"SINGLE_SNP_ONLY", "SINGLE_SNP_DOMINANT"}:
            return "ASE_SUPPORTED_SINGLE_SNP"
        if support.startswith("MULTI_SNP") and phase in {"REVIEW", "FAIL"}:
            return "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED"
        return "ASE_SUPPORT_UNRESOLVED"
    if signal == "FDR_MODERATE_EFFECT":
        return "ASE_FDR_MODERATE_EFFECT"
    if signal == "FDR_SMALL_EFFECT":
        return "ASE_STATISTICAL_SMALL_EFFECT"
    if signal == "NEAR_FDR_STRONG_EFFECT":
        return "ASE_NEAR_FDR"
    if signal == "NOMINAL_STRONG_EFFECT":
        return "ASE_NOMINAL_DIRECTIONAL"
    return "NO_ASE_SIGNAL"


def block_key(row):
    return "\x1f".join(
        [
            row.get("sample", ""),
            row.get("tissue", ""),
            row.get("gene_id", ""),
            row.get("haplotype_block_id", ""),
        ]
    )


def gene_key(row):
    return "\x1f".join(
        [
            row.get("sample", ""),
            row.get("tissue", ""),
            row.get("gene_id", ""),
            row.get("gene_name", ""),
            row.get("gene_biotype", ""),
        ]
    )


def create_db(path):
    connection = sqlite3.connect(str(path))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute(
        "CREATE TABLE block_keys (key TEXT PRIMARY KEY, evidence TEXT NOT NULL)"
    )
    connection.execute(
        """
        CREATE TABLE measurement_obs (
            measurement_id TEXT,
            sample TEXT,
            tissue TEXT,
            pvalue REAL,
            qvalue REAL,
            UNIQUE(measurement_id,sample,tissue,pvalue,qvalue)
        )
        """
    )
    connection.execute(
        "CREATE TABLE gene_rows (gene_key TEXT, evidence TEXT, test_status TEXT, technical_status TEXT)"
    )
    return connection


def add_error(errors, message, limit=1000):
    if len(errors) < limit:
        errors.append(message)


def read_subset_keys(path, expected_classes, errors):
    keys = set()
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = block_key(row)
            if key in keys:
                add_error(errors, "duplicate_subset_block:{}".format(key))
            keys.add(key)
            if row.get("v026_evidence_class") not in expected_classes:
                add_error(
                    errors,
                    "invalid_subset_class:{}:{}".format(
                        key, row.get("v026_evidence_class")
                    ),
                )
    return keys


def validate_bh(connection, errors, tolerance):
    conflicts = connection.execute(
        """
        SELECT measurement_id,COUNT(*)
        FROM measurement_obs
        GROUP BY measurement_id
        HAVING COUNT(*) > 1
        LIMIT 100
        """
    ).fetchall()
    for measurement_id, count in conflicts:
        add_error(
            errors,
            "measurement_stat_conflict:{}:{}".format(measurement_id, count),
        )
    groups = connection.execute(
        "SELECT DISTINCT sample,tissue FROM measurement_obs ORDER BY sample,tissue"
    ).fetchall()
    n_measurements = 0
    for sample, tissue in groups:
        observations = connection.execute(
            """
            SELECT measurement_id,pvalue,qvalue
            FROM measurement_obs
            WHERE sample=? AND tissue=?
            ORDER BY pvalue,measurement_id
            """,
            (sample, tissue),
        ).fetchall()
        n = len(observations)
        n_measurements += n
        adjusted = [1.0] * n
        running = 1.0
        for index in range(n - 1, -1, -1):
            rank = index + 1
            candidate = min(1.0, observations[index][1] * n / float(rank))
            running = min(running, candidate)
            adjusted[index] = running
        for observation, expected in zip(observations, adjusted):
            measurement_id, _pvalue, observed = observation
            if observed is None or abs(observed - expected) > tolerance:
                add_error(
                    errors,
                    "bh_qvalue_mismatch:{}:{}:{}".format(
                        measurement_id, observed, expected
                    ),
                )
    return n_measurements


def expected_gene_status(classes, has_clean_tested):
    if classes.intersection(CORE_CLASSES):
        return "ASE_GENE_CORE"
    if "ASE_SUPPORTED_SINGLE_SNP" in classes:
        return "ASE_GENE_SUPPORTED_SINGLE_SNP"
    if classes.intersection(
        {"ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED", "ASE_SUPPORT_UNRESOLVED"}
    ):
        return "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED"
    if "ASE_FDR_MODERATE_EFFECT" in classes:
        return "ASE_GENE_FDR_MODERATE"
    if "ASE_NEAR_FDR" in classes:
        return "ASE_GENE_NEAR_FDR"
    if "ASE_NOMINAL_DIRECTIONAL" in classes:
        return "ASE_GENE_NOMINAL"
    if "ASE_STATISTICAL_SMALL_EFFECT" in classes:
        return "ASE_GENE_STATISTICAL_SMALL_EFFECT"
    return "NOT_ASE_GENE" if has_clean_tested else "UNTESTABLE"


def calculate_expected_gene_statuses(connection):
    expected = {}
    current = None
    classes = set()
    has_clean_tested = False
    query = """
        SELECT gene_key,evidence,test_status,technical_status
        FROM gene_rows ORDER BY gene_key
    """
    for key, evidence, test_status, technical_status in connection.execute(query):
        if current is not None and key != current:
            expected[current] = expected_gene_status(classes, has_clean_tested)
            classes = set()
            has_clean_tested = False
        current = key
        classes.add(evidence)
        if test_status == "TESTED" and technical_status == "PASS":
            has_clean_tested = True
    if current is not None:
        expected[current] = expected_gene_status(classes, has_clean_tested)
    return expected


def validate_gene_summary(path, expected, errors):
    observed = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample",
            "tissue",
            "gene_id",
            "gene_name",
            "gene_biotype",
            "gene_ase_status_v026",
            "top_evidence_class_v026",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            add_error(errors, "gene_summary_missing_fields:{}".format(",".join(sorted(missing))))
            return 0, Counter()
        for row in reader:
            key = gene_key(row)
            if key in observed:
                add_error(errors, "duplicate_gene_summary:{}".format(key))
            status = row.get("gene_ase_status_v026", "")
            observed[key] = status
            if expected.get(key) != status:
                add_error(
                    errors,
                    "gene_status_mismatch:{}:{}:{}".format(
                        key, status, expected.get(key)
                    ),
                )
    for key in set(expected).difference(observed):
        add_error(errors, "missing_gene_summary:{}".format(key))
    for key in set(observed).difference(expected):
        add_error(errors, "extra_gene_summary:{}".format(key))
    return len(observed), Counter(observed.values())


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate ASE v0.2.6 outputs")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--expect-sample")
    parser.add_argument("--skip-bh-recompute", action="store_true")
    parser.add_argument("--qvalue-tolerance", type=float, default=1e-9)
    args = parser.parse_args(argv)

    output_dir = Path(args.out_dir).resolve()
    errors = []
    warnings = []
    missing_outputs = [name for name in REQUIRED_OUTPUTS if not (output_dir / name).is_file()]
    if missing_outputs:
        report = {
            "version": "0.2.6",
            "status": "FAIL",
            "n_errors": len(missing_outputs),
            "errors_first_100": ["missing_output:{}".format(name) for name in missing_outputs],
        }
        print(json.dumps(report, indent=2, sort_keys=True))
        return 2

    with open(output_dir / "ase_v026_qc.json", "r", encoding="utf-8") as handle:
        qc = json.load(handle)
    thresholds = qc.get("thresholds", {})
    required_thresholds = {
        "fdr",
        "near_fdr",
        "nominal_p",
        "strong_major_fraction",
        "moderate_major_fraction",
        "max_fragment_hap_conflict_fraction",
    }
    missing_thresholds = required_thresholds.difference(thresholds)
    if missing_thresholds:
        add_error(errors, "qc_missing_thresholds:{}".format(",".join(sorted(missing_thresholds))))

    db_fd, db_name = tempfile.mkstemp(prefix="ase_v026_validate_", suffix=".sqlite")
    os.close(db_fd)
    connection = create_db(db_name)
    cursor = connection.cursor()
    class_counts = Counter()
    signal_counts = Counter()
    shape_counts = Counter()
    expected_core = set()
    expected_supported = set()
    expected_sensitivity = set()
    n_rows = 0
    try:
        block_path = output_dir / "ase_haplotype_block_reclassified.tsv"
        with open(block_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "sample",
                "tissue",
                "gene_id",
                "haplotype_block_id",
                "measurement_id",
                "ase_a_count",
                "ase_b_count",
                "ase_total_count",
                "major_haplotype_fraction",
                "ase_test_status",
                "ase_p_exact",
                "ase_q_bh",
                "analysis_set",
                "count_validation",
                "n_wgs_balanced_gene_snps",
                "fragment_hap_conflict_fraction",
                "snp_support_status",
                "local_phase_status",
                "v026_technical_qc_status",
                "v026_statistical_signal",
                "v026_imbalance_shape",
                "v026_evidence_class",
                "v026_is_core_ase",
            }
            missing = required.difference(reader.fieldnames or [])
            if missing:
                add_error(errors, "block_output_missing_fields:{}".format(",".join(sorted(missing))))
            for row_number, row in enumerate(reader, start=2):
                n_rows += 1
                key = block_key(row)
                if args.expect_sample and row.get("sample") != args.expect_sample:
                    add_error(
                        errors,
                        "unexpected_sample:{}:{}".format(key, row.get("sample")),
                    )
                a_count = integer(row.get("ase_a_count"))
                b_count = integer(row.get("ase_b_count"))
                total = integer(row.get("ase_total_count"))
                if (
                    a_count is None
                    or b_count is None
                    or total is None
                    or a_count + b_count != total
                ):
                    add_error(errors, "count_conservation:{}".format(key))
                elif total:
                    expected_fraction = max(a_count, b_count) / float(total)
                    observed_fraction = number(row.get("major_haplotype_fraction"))
                    if observed_fraction is None or abs(observed_fraction - expected_fraction) > 1e-9:
                        add_error(errors, "major_fraction_mismatch:{}".format(key))

                signal = expected_signal(row, thresholds) if not missing_thresholds else ""
                shape = expected_shape(row)
                technical = expected_technical(row, thresholds) if not missing_thresholds else ""
                evidence = expected_evidence(
                    technical,
                    signal,
                    row.get("snp_support_status", ""),
                    row.get("local_phase_status", ""),
                )
                if row.get("v026_statistical_signal") != signal:
                    add_error(errors, "statistical_signal_mismatch:{}".format(key))
                if row.get("v026_imbalance_shape") != shape:
                    add_error(errors, "imbalance_shape_mismatch:{}".format(key))
                if row.get("v026_technical_qc_status") != technical:
                    add_error(errors, "technical_qc_mismatch:{}".format(key))
                if row.get("v026_evidence_class") != evidence:
                    add_error(errors, "evidence_class_mismatch:{}".format(key))
                expected_core_flag = "1" if evidence in CORE_CLASSES else "0"
                if row.get("v026_is_core_ase") != expected_core_flag:
                    add_error(errors, "core_flag_mismatch:{}".format(key))

                try:
                    cursor.execute("INSERT INTO block_keys VALUES (?,?)", (key, evidence))
                except sqlite3.IntegrityError:
                    add_error(errors, "duplicate_block_gene_row:{}".format(key))
                if row.get("analysis_set") == "MAIN":
                    cursor.execute(
                        "INSERT INTO gene_rows VALUES (?,?,?,?)",
                        (
                            gene_key(row),
                            evidence,
                            row.get("ase_test_status", ""),
                            technical,
                        ),
                    )
                if (
                    row.get("analysis_set") == "MAIN"
                    and row.get("ase_test_status") == "TESTED"
                ):
                    cursor.execute(
                        "INSERT OR IGNORE INTO measurement_obs VALUES (?,?,?,?,?)",
                        (
                            row.get("measurement_id", ""),
                            row.get("sample", ""),
                            row.get("tissue", ""),
                            number(row.get("ase_p_exact")),
                            number(row.get("ase_q_bh")),
                        ),
                    )
                if evidence in CORE_CLASSES:
                    expected_core.add(key)
                elif evidence in SUPPORTED_CLASSES:
                    expected_supported.add(key)
                elif evidence in SENSITIVITY_CLASSES:
                    expected_sensitivity.add(key)
                class_counts[evidence] += 1
                signal_counts[row.get("v026_statistical_signal", "")] += 1
                shape_counts[row.get("v026_imbalance_shape", "")] += 1
                if n_rows % 10000 == 0:
                    connection.commit()
        connection.commit()

        observed_core = read_subset_keys(
            output_dir / "ase_core_haplotype_blocks.tsv", CORE_CLASSES, errors
        )
        observed_supported = read_subset_keys(
            output_dir / "ase_supported_haplotype_blocks.tsv", SUPPORTED_CLASSES, errors
        )
        observed_sensitivity = read_subset_keys(
            output_dir / "ase_sensitivity_haplotype_blocks.tsv", SENSITIVITY_CLASSES, errors
        )
        if observed_core != expected_core:
            add_error(errors, "core_block_subset_mismatch")
        if observed_supported != expected_supported:
            add_error(errors, "supported_block_subset_mismatch")
        if observed_sensitivity != expected_sensitivity:
            add_error(errors, "sensitivity_block_subset_mismatch")

        if args.skip_bh_recompute:
            warnings.append("BH_RECOMPUTATION_SKIPPED")
            n_unique_measurements = connection.execute(
                "SELECT COUNT(*) FROM measurement_obs"
            ).fetchone()[0]
        else:
            n_unique_measurements = validate_bh(
                connection, errors, args.qvalue_tolerance
            )

        expected_genes = calculate_expected_gene_statuses(connection)
        n_genes, gene_status_counts = validate_gene_summary(
            output_dir / "ase_gene_summary.tsv", expected_genes, errors
        )
    finally:
        connection.close()
        try:
            os.remove(db_name)
        except FileNotFoundError:
            pass

    if n_rows != qc.get("n_block_gene_rows"):
        add_error(
            errors,
            "qc_block_row_count_mismatch:{}:{}".format(
                n_rows, qc.get("n_block_gene_rows")
            ),
        )
    if n_genes != qc.get("n_gene_summary_rows"):
        add_error(
            errors,
            "qc_gene_row_count_mismatch:{}:{}".format(
                n_genes, qc.get("n_gene_summary_rows")
            ),
        )
    if dict(sorted(class_counts.items())) != qc.get("block_evidence_class_counts"):
        add_error(errors, "qc_evidence_class_counts_mismatch")
    if dict(sorted(gene_status_counts.items())) != qc.get("gene_status_counts"):
        add_error(errors, "qc_gene_status_counts_mismatch")

    report = {
        "version": "0.2.6",
        "status": "PASS" if not errors else "FAIL",
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "errors_first_100": errors[:100],
        "warnings": warnings[:100],
        "n_block_gene_rows": n_rows,
        "n_gene_summary_rows": n_genes,
        "n_unique_tested_measurements": n_unique_measurements,
        "block_evidence_class_counts": dict(sorted(class_counts.items())),
        "statistical_signal_counts": dict(sorted(signal_counts.items())),
        "imbalance_shape_counts": dict(sorted(shape_counts.items())),
        "gene_status_counts": dict(sorted(gene_status_counts.items())),
    }
    with open(output_dir / "ase_v026_validation.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
