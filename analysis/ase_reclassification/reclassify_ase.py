#!/usr/bin/env python3
"""Reclassify validated ASE v0.2.5 block results without recounting reads.

ASE v0.2.6 deliberately keeps the v0.2.5 measurement and multiple-testing
universe fixed.  It separates statistical signal, imbalance shape, SNP
support, local-phase confidence, assignment dependency, and technical QC,
then derives transparent reporting classes from those independent axes.
"""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import shutil
import sqlite3
import sys
import tempfile
from collections import Counter
from pathlib import Path


VERSION = "0.2.6"

DEFAULT_THRESHOLDS = {
    "fdr": 0.05,
    "near_fdr": 0.10,
    "nominal_p": 0.05,
    "strong_major_fraction": 0.65,
    "moderate_major_fraction": 0.60,
    "max_fragment_hap_conflict_fraction": 0.10,
}

REQUIRED_BLOCK_FIELDS = {
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "haplotype_block_id",
    "measurement_id",
    "analysis_set",
    "count_validation",
    "n_wgs_balanced_gene_snps",
    "fragment_hap_conflict_fraction",
    "ase_a_count",
    "ase_b_count",
    "ase_total_count",
    "major_haplotype_fraction",
    "ase_test_status",
    "ase_p_exact",
    "ase_q_bh",
    "snp_support_status",
    "local_phase_status",
    "wgs_phase_status",
    "wgs_phase_orientation",
    "wgs_phase_set",
    "candidate_class",
    "is_primary_candidate",
    "exclusion_reasons",
}

V026_FIELDS = [
    "candidate_class_v025",
    "is_primary_candidate_v025",
    "v026_technical_qc_status",
    "v026_technical_qc_reasons",
    "v026_statistical_signal",
    "v026_imbalance_shape",
    "v026_snp_support_class",
    "v026_phase_qc_status",
    "v026_assignment_dependency",
    "v026_asm_integration_readiness",
    "v026_ase_direction",
    "v026_evidence_class",
    "v026_is_fdr_strong_signal",
    "v026_is_core_ase",
    "v026_reclassification_reasons",
]

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

CLASS_PRIORITY = {
    "ASE_CORE_ROBUST": 0,
    "ASE_CORE_SUPPORTED": 1,
    "ASE_SUPPORTED_SINGLE_SNP": 2,
    "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED": 3,
    "ASE_SUPPORT_UNRESOLVED": 4,
    "ASE_FDR_MODERATE_EFFECT": 5,
    "ASE_NEAR_FDR": 6,
    "ASE_NOMINAL_DIRECTIONAL": 7,
    "ASE_STATISTICAL_SMALL_EFFECT": 8,
    "NO_ASE_SIGNAL": 9,
    "LOW_COVERAGE": 10,
    "NOT_TESTABLE": 11,
    "TECHNICAL_FAIL": 12,
    "EXCLUDED_TISSUE": 13,
}

GENE_FIELDS = [
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "n_informative_haplotype_blocks",
    "n_testable_haplotype_blocks",
    "n_fdr_strong_blocks",
    "n_core_robust_blocks",
    "n_core_supported_blocks",
    "n_supported_single_snp_blocks",
    "n_phase_unresolved_blocks",
    "n_support_unresolved_blocks",
    "n_fdr_moderate_effect_blocks",
    "n_near_fdr_blocks",
    "n_nominal_directional_blocks",
    "n_statistical_small_effect_blocks",
    "n_monoallelic_like_blocks",
    "n_low_minor_blocks",
    "n_two_haplotype_supported_blocks",
    "gene_phase_scope",
    "gene_ase_status_v026",
    "top_haplotype_block_id",
    "top_measurement_id",
    "top_evidence_class_v026",
    "top_statistical_signal_v026",
    "top_imbalance_shape_v026",
    "top_snp_support_class_v026",
    "top_phase_qc_status_v026",
    "top_assignment_dependency_v026",
    "top_asm_integration_readiness_v026",
    "top_ase_direction_v026",
    "top_ase_q_bh",
    "top_major_haplotype_fraction",
    "top_total_count",
    "core_haplotype_block_ids",
    "supported_haplotype_block_ids",
    "sensitivity_haplotype_block_ids",
]

TRANSITION_FIELDS = [
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "haplotype_block_id",
    "measurement_id",
    "candidate_class_v025",
    "exclusion_reasons_v025",
    "v026_evidence_class",
    "v026_statistical_signal",
    "v026_imbalance_shape",
    "v026_snp_support_class",
    "v026_phase_qc_status",
    "v026_technical_qc_status",
    "v026_reclassification_reasons",
]


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def parse_float(value):
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


def parse_int(value):
    number = parse_float(value)
    if number is None:
        return None
    if not float(number).is_integer():
        return None
    return int(number)


def bool_text(value):
    return "1" if value else "0"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_input_file(input_dir, explicit, names, required=True):
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise FileNotFoundError("input file not found: {}".format(path))
        return path
    if input_dir:
        root = Path(input_dir)
        for name in names:
            path = root / name
            if path.is_file():
                return path
    if required:
        raise FileNotFoundError(
            "none of the required inputs were found: {}".format(
                ", ".join(names)
            )
        )
    return None


def prepare_output_dir(path, force=False):
    path = Path(path).resolve()
    if path.exists() and any(path.iterdir()):
        if not force:
            raise RuntimeError(
                "output directory is not empty; use --force to replace it: {}".format(
                    path
                )
            )
        shutil.rmtree(str(path))
    path.mkdir(parents=True, exist_ok=True)
    (path / "audit").mkdir(parents=True, exist_ok=True)
    return path


def validate_source(validation_path, qc_path, allow_unvalidated=False):
    warnings = []
    validation = None
    qc = None
    if validation_path is None:
        if not allow_unvalidated:
            raise RuntimeError(
                "ase_v025_validation.json is required in production mode"
            )
        warnings.append("SOURCE_VALIDATION_NOT_PROVIDED")
    else:
        validation = load_json(validation_path)
        if validation.get("status") != "PASS" or validation.get("n_errors", 0) != 0:
            raise RuntimeError(
                "source ASE v0.2.5 validation did not PASS: {}".format(
                    validation_path
                )
            )
        if str(validation.get("version")) != "0.2.5":
            warnings.append("SOURCE_VALIDATION_VERSION_NOT_0.2.5")
    if qc_path is None:
        if not allow_unvalidated:
            raise RuntimeError("ase_v025_qc.json is required in production mode")
        warnings.append("SOURCE_QC_NOT_PROVIDED")
    else:
        qc = load_json(qc_path)
    return validation, qc, warnings


def load_assignment_audit(path):
    if path is None:
        return {}, []
    required = {
        "measurement_id",
        "direct_total_count",
        "strand_resolved_total",
        "rescue_fraction",
        "rescue_dependency",
    }
    result = {}
    warnings = []
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "assignment audit missing fields: {}".format(
                    ",".join(sorted(missing))
                )
            )
        for row in reader:
            measurement_id = row["measurement_id"]
            value = (
                row.get("rescue_dependency", ""),
                parse_float(row.get("rescue_fraction")),
                parse_int(row.get("direct_total_count")),
                parse_int(row.get("strand_resolved_total")),
            )
            previous = result.get(measurement_id)
            if previous is not None and previous != value:
                warnings.append(
                    "ASSIGNMENT_AUDIT_CONFLICT:{}".format(measurement_id)
                )
            result[measurement_id] = value
    return result, warnings


def assignment_dependency(measurement_id, assignment_index):
    value = assignment_index.get(measurement_id)
    if value is None:
        return "NOT_ASSESSED"
    dependency, fraction, direct_total, strand_total = value
    if dependency == "NO_STRAND_RESOLVED_FRAGMENTS" or strand_total == 0:
        return "DIRECT_ONLY"
    if dependency == "NOT_DEPENDENT":
        return "MIXED_NOT_DEPENDENT"
    if dependency == "RESCUE_DEPENDENT_REVIEW":
        return "STRAND_RESOLVED_DEPENDENT"
    if fraction is not None and fraction == 0:
        return "DIRECT_ONLY"
    if direct_total == 0 and (strand_total or 0) > 0:
        return "STRAND_RESOLVED_DEPENDENT"
    return "MIXED_DEPENDENCY_UNKNOWN"


def statistical_signal(row, thresholds):
    status = row.get("ase_test_status", "")
    if status == "LOW_COVERAGE":
        return "LOW_COVERAGE"
    if status != "TESTED":
        return "NOT_TESTABLE"
    qvalue = parse_float(row.get("ase_q_bh"))
    pvalue = parse_float(row.get("ase_p_exact"))
    effect = parse_float(row.get("major_haplotype_fraction"))
    if qvalue is None or effect is None:
        return "NOT_TESTABLE"
    if qvalue <= thresholds["fdr"]:
        if effect >= thresholds["strong_major_fraction"]:
            return "FDR_STRONG_EFFECT"
        if effect >= thresholds["moderate_major_fraction"]:
            return "FDR_MODERATE_EFFECT"
        return "FDR_SMALL_EFFECT"
    if (
        qvalue <= thresholds["near_fdr"]
        and effect >= thresholds["strong_major_fraction"]
    ):
        return "NEAR_FDR_STRONG_EFFECT"
    if (
        pvalue is not None
        and pvalue <= thresholds["nominal_p"]
        and effect >= thresholds["strong_major_fraction"]
    ):
        return "NOMINAL_STRONG_EFFECT"
    return "NO_SELECTED_SIGNAL"


def imbalance_shape(row):
    a_count = parse_int(row.get("ase_a_count"))
    b_count = parse_int(row.get("ase_b_count"))
    if a_count is None or b_count is None:
        return "NOT_ASSESSED"
    minor = min(a_count, b_count)
    if minor <= 1:
        return "MONOALLELIC_LIKE"
    if minor == 2:
        return "LOW_MINOR"
    return "TWO_HAPLOTYPE_SUPPORTED"


def technical_qc(row, thresholds):
    if row.get("analysis_set") == "EXCLUDED_TISSUE":
        return "EXCLUDED", ["EXCLUDED_TISSUE"]
    reasons = []
    if row.get("analysis_set") != "MAIN":
        reasons.append("UNKNOWN_ANALYSIS_SET")
    if row.get("count_validation") != "MATCH":
        reasons.append("COUNT_MISMATCH")
    balanced = parse_int(row.get("n_wgs_balanced_gene_snps"))
    if balanced is None or balanced <= 0:
        reasons.append("NO_WGS_BALANCED_SNP")
    conflict = parse_float(row.get("fragment_hap_conflict_fraction"))
    if conflict is None:
        reasons.append("FRAGMENT_CONFLICT_NOT_ASSESSED")
    elif conflict > thresholds["max_fragment_hap_conflict_fraction"]:
        reasons.append("FRAGMENT_CONFLICT")
    if not row.get("gene_id") or not row.get("haplotype_block_id"):
        reasons.append("MISSING_BLOCK_OR_GENE_ID")
    return ("FAIL", reasons) if reasons else ("PASS", [])


def asm_readiness(row):
    status = row.get("wgs_phase_status", "")
    phase_set = row.get("wgs_phase_set", "").strip()
    phase_set_present = bool(phase_set) and phase_set.upper() not in {
        "NA", "NAN", "NONE", "."
    }
    if status == "CONCORDANT" and phase_set_present:
        return "READY_SINGLE_WGS_PS"
    if status == "MULTI_PHASE_SET":
        return "MULTI_PS_SPLIT_REQUIRED"
    if status == "DISCORDANT":
        return "ORIENTATION_DISCORDANT"
    return "PHASE_INSUFFICIENT"


def ase_direction(row):
    a_count = parse_int(row.get("ase_a_count"))
    b_count = parse_int(row.get("ase_b_count"))
    if a_count is None or b_count is None:
        return "NOT_ASSESSED"
    if a_count > b_count:
        return "HAP_A"
    if b_count > a_count:
        return "HAP_B"
    return "TIE"


def derive_evidence_class(technical_status, signal, support, phase):
    if technical_status == "EXCLUDED":
        return "EXCLUDED_TISSUE"
    if technical_status != "PASS":
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


def reclassification_reasons(row, signal, shape, support, phase, technical_reasons):
    reasons = list(technical_reasons)
    if signal != "NO_SELECTED_SIGNAL":
        reasons.append(signal)
    if shape in {"MONOALLELIC_LIKE", "LOW_MINOR"}:
        reasons.append(shape)
    if support:
        reasons.append(support)
    if phase and phase != "PASS":
        reasons.append("LOCAL_PHASE_{}".format(phase))
    old_reasons = set(
        item for item in row.get("exclusion_reasons", "").split(",") if item
    )
    if old_reasons.intersection({"EXTREME_MINOR", "LOW_MINOR"}):
        reasons.append("V025_MINOR_HARD_FAILURE_REMOVED")
    return ",".join(dict.fromkeys(reasons))


def validate_numeric_invariants(row, row_number):
    a_count = parse_int(row.get("ase_a_count"))
    b_count = parse_int(row.get("ase_b_count"))
    total = parse_int(row.get("ase_total_count"))
    if a_count is None or b_count is None or total is None:
        raise ValueError("row {} has invalid ASE counts".format(row_number))
    if a_count < 0 or b_count < 0 or total < 0 or a_count + b_count != total:
        raise ValueError("row {} violates ASE count conservation".format(row_number))
    observed = parse_float(row.get("major_haplotype_fraction"))
    expected = max(a_count, b_count) / float(total) if total else None
    if observed is not None and expected is not None and abs(observed - expected) > 1e-9:
        raise ValueError(
            "row {} has inconsistent major_haplotype_fraction".format(row_number)
        )


def classify_row(row, thresholds, assignment_index):
    technical_status, technical_reasons = technical_qc(row, thresholds)
    signal = statistical_signal(row, thresholds)
    shape = imbalance_shape(row)
    support = row.get("snp_support_status", "") or "SNP_SUPPORT_UNRESOLVED"
    phase = row.get("local_phase_status", "") or "NOT_ASSESSED"
    evidence = derive_evidence_class(technical_status, signal, support, phase)
    classified = dict(row)
    classified.update(
        {
            "candidate_class_v025": row.get("candidate_class", ""),
            "is_primary_candidate_v025": row.get("is_primary_candidate", ""),
            "v026_technical_qc_status": technical_status,
            "v026_technical_qc_reasons": ",".join(technical_reasons),
            "v026_statistical_signal": signal,
            "v026_imbalance_shape": shape,
            "v026_snp_support_class": support,
            "v026_phase_qc_status": phase,
            "v026_assignment_dependency": assignment_dependency(
                row.get("measurement_id", ""), assignment_index
            ),
            "v026_asm_integration_readiness": asm_readiness(row),
            "v026_ase_direction": ase_direction(row),
            "v026_evidence_class": evidence,
            "v026_is_fdr_strong_signal": bool_text(
                signal == "FDR_STRONG_EFFECT"
            ),
            "v026_is_core_ase": bool_text(evidence in CORE_CLASSES),
            "v026_reclassification_reasons": reclassification_reasons(
                row,
                signal,
                shape,
                support,
                phase,
                technical_reasons,
            ),
        }
    )
    return classified


def create_gene_db(path):
    connection = sqlite3.connect(str(path))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute(
        """
        CREATE TABLE block_gene (
            sample TEXT NOT NULL,
            tissue TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            gene_name TEXT,
            gene_biotype TEXT,
            block_id TEXT NOT NULL,
            measurement_id TEXT,
            total_count INTEGER,
            qvalue REAL,
            major_fraction REAL,
            test_status TEXT,
            evidence_class TEXT,
            statistical_signal TEXT,
            imbalance_shape TEXT,
            snp_support TEXT,
            phase_status TEXT,
            assignment_dependency TEXT,
            asm_readiness TEXT,
            ase_direction TEXT,
            technical_status TEXT
        )
        """
    )
    return connection


def insert_gene_row(cursor, row):
    cursor.execute(
        """
        INSERT INTO block_gene VALUES (
            ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
        )
        """,
        (
            row["sample"],
            row["tissue"],
            row["gene_id"],
            row.get("gene_name", ""),
            row.get("gene_biotype", ""),
            row["haplotype_block_id"],
            row.get("measurement_id", ""),
            parse_int(row.get("ase_total_count")),
            parse_float(row.get("ase_q_bh")),
            parse_float(row.get("major_haplotype_fraction")),
            row.get("ase_test_status", ""),
            row["v026_evidence_class"],
            row["v026_statistical_signal"],
            row["v026_imbalance_shape"],
            row["v026_snp_support_class"],
            row["v026_phase_qc_status"],
            row["v026_assignment_dependency"],
            row["v026_asm_integration_readiness"],
            row["v026_ase_direction"],
            row["v026_technical_qc_status"],
        ),
    )


def top_sort_key(row):
    qvalue = row[8] if row[8] is not None else float("inf")
    total = row[7] if row[7] is not None else -1
    fraction = row[9] if row[9] is not None else -1.0
    evidence = row[11]
    return (
        CLASS_PRIORITY.get(evidence, 999),
        qvalue,
        -total,
        -fraction,
        row[5],
    )


def gene_status(rows):
    classes = {row[11] for row in rows}
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
    if not any(row[10] == "TESTED" and row[19] == "PASS" for row in rows):
        return "UNTESTABLE"
    return "NOT_ASE_GENE"


def summarize_gene(rows):
    first = rows[0]
    counts = Counter(row[11] for row in rows)
    shape_counts = Counter(row[13] for row in rows)
    blocks = sorted({row[5] for row in rows if (row[7] or 0) > 0})
    core_blocks = sorted({row[5] for row in rows if row[11] in CORE_CLASSES})
    supported_blocks = sorted(
        {row[5] for row in rows if row[11] in SUPPORTED_CLASSES}
    )
    sensitivity_blocks = sorted(
        {row[5] for row in rows if row[11] in SENSITIVITY_CLASSES}
    )
    top = min(rows, key=top_sort_key)
    phase_scope = "SINGLE_BLOCK" if len(blocks) <= 1 else "MULTI_BLOCK_UNRESOLVED"
    return {
        "sample": first[0],
        "tissue": first[1],
        "gene_id": first[2],
        "gene_name": first[3],
        "gene_biotype": first[4],
        "n_informative_haplotype_blocks": len(blocks),
        "n_testable_haplotype_blocks": sum(row[10] == "TESTED" for row in rows),
        "n_fdr_strong_blocks": sum(
            row[12] == "FDR_STRONG_EFFECT" for row in rows
        ),
        "n_core_robust_blocks": counts["ASE_CORE_ROBUST"],
        "n_core_supported_blocks": counts["ASE_CORE_SUPPORTED"],
        "n_supported_single_snp_blocks": counts["ASE_SUPPORTED_SINGLE_SNP"],
        "n_phase_unresolved_blocks": counts[
            "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED"
        ],
        "n_support_unresolved_blocks": counts["ASE_SUPPORT_UNRESOLVED"],
        "n_fdr_moderate_effect_blocks": counts["ASE_FDR_MODERATE_EFFECT"],
        "n_near_fdr_blocks": counts["ASE_NEAR_FDR"],
        "n_nominal_directional_blocks": counts["ASE_NOMINAL_DIRECTIONAL"],
        "n_statistical_small_effect_blocks": counts[
            "ASE_STATISTICAL_SMALL_EFFECT"
        ],
        "n_monoallelic_like_blocks": shape_counts["MONOALLELIC_LIKE"],
        "n_low_minor_blocks": shape_counts["LOW_MINOR"],
        "n_two_haplotype_supported_blocks": shape_counts[
            "TWO_HAPLOTYPE_SUPPORTED"
        ],
        "gene_phase_scope": phase_scope,
        "gene_ase_status_v026": gene_status(rows),
        "top_haplotype_block_id": top[5],
        "top_measurement_id": top[6],
        "top_evidence_class_v026": top[11],
        "top_statistical_signal_v026": top[12],
        "top_imbalance_shape_v026": top[13],
        "top_snp_support_class_v026": top[14],
        "top_phase_qc_status_v026": top[15],
        "top_assignment_dependency_v026": top[16],
        "top_asm_integration_readiness_v026": top[17],
        "top_ase_direction_v026": top[18],
        "top_ase_q_bh": "" if top[8] is None else top[8],
        "top_major_haplotype_fraction": "" if top[9] is None else top[9],
        "top_total_count": "" if top[7] is None else top[7],
        "core_haplotype_block_ids": ",".join(core_blocks),
        "supported_haplotype_block_ids": ",".join(supported_blocks),
        "sensitivity_haplotype_block_ids": ",".join(sensitivity_blocks),
    }


def write_gene_outputs(connection, output_dir):
    query = """
        SELECT sample,tissue,gene_id,gene_name,gene_biotype,block_id,
               measurement_id,total_count,qvalue,major_fraction,test_status,
               evidence_class,statistical_signal,imbalance_shape,snp_support,
               phase_status,assignment_dependency,asm_readiness,ase_direction,
               technical_status
        FROM block_gene
        ORDER BY sample,tissue,gene_id,gene_name,gene_biotype,block_id
    """
    summary_path = output_dir / "ase_gene_summary.tsv"
    core_path = output_dir / "ase_core_genes.tsv"
    supported_path = output_dir / "ase_supported_genes.tsv"
    sensitivity_path = output_dir / "ase_sensitivity_genes.tsv"
    status_counts = Counter()
    n_genes = 0
    with open(summary_path, "w", encoding="utf-8", newline="") as summary_handle, \
            open(core_path, "w", encoding="utf-8", newline="") as core_handle, \
            open(supported_path, "w", encoding="utf-8", newline="") as supported_handle, \
            open(sensitivity_path, "w", encoding="utf-8", newline="") as sensitivity_handle:
        summary_writer = csv.DictWriter(summary_handle, fieldnames=GENE_FIELDS, delimiter="\t")
        core_writer = csv.DictWriter(core_handle, fieldnames=GENE_FIELDS, delimiter="\t")
        supported_writer = csv.DictWriter(
            supported_handle, fieldnames=GENE_FIELDS, delimiter="\t"
        )
        sensitivity_writer = csv.DictWriter(
            sensitivity_handle, fieldnames=GENE_FIELDS, delimiter="\t"
        )
        for writer in (summary_writer, core_writer, supported_writer, sensitivity_writer):
            writer.writeheader()

        current_key = None
        group = []
        for row in connection.execute(query):
            key = row[:5]
            if current_key is not None and key != current_key:
                summary = summarize_gene(group)
                summary_writer.writerow(summary)
                status = summary["gene_ase_status_v026"]
                status_counts[status] += 1
                n_genes += 1
                if status == "ASE_GENE_CORE":
                    core_writer.writerow(summary)
                elif status in {
                    "ASE_GENE_SUPPORTED_SINGLE_SNP",
                    "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED",
                }:
                    supported_writer.writerow(summary)
                elif status not in {"NOT_ASE_GENE", "UNTESTABLE"}:
                    sensitivity_writer.writerow(summary)
                group = []
            current_key = key
            group.append(row)
        if group:
            summary = summarize_gene(group)
            summary_writer.writerow(summary)
            status = summary["gene_ase_status_v026"]
            status_counts[status] += 1
            n_genes += 1
            if status == "ASE_GENE_CORE":
                core_writer.writerow(summary)
            elif status in {
                "ASE_GENE_SUPPORTED_SINGLE_SNP",
                "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED",
            }:
                supported_writer.writerow(summary)
            elif status not in {"NOT_ASE_GENE", "UNTESTABLE"}:
                sensitivity_writer.writerow(summary)
    return n_genes, status_counts


def run_reclassification(args):
    started_at = utc_now()
    thresholds = dict(DEFAULT_THRESHOLDS)
    thresholds.update(
        {
            "fdr": args.fdr,
            "near_fdr": args.near_fdr,
            "nominal_p": args.nominal_p,
            "strong_major_fraction": args.strong_major_fraction,
            "moderate_major_fraction": args.moderate_major_fraction,
            "max_fragment_hap_conflict_fraction": args.max_fragment_hap_conflict_fraction,
        }
    )
    if thresholds["moderate_major_fraction"] >= thresholds["strong_major_fraction"]:
        raise ValueError("moderate major fraction must be below strong fraction")
    if thresholds["near_fdr"] < thresholds["fdr"]:
        raise ValueError("near-FDR cutoff must be greater than or equal to FDR")

    block_path = resolve_input_file(
        args.input_dir,
        args.block_results,
        ["ase_haplotype_block_results.tsv", "ase_haplotype_block_results.tsv.gz"],
    )
    validation_path = resolve_input_file(
        args.input_dir,
        args.validation_json,
        ["ase_v025_validation.json"],
        required=not args.allow_unvalidated_input,
    )
    qc_path = resolve_input_file(
        args.input_dir,
        args.qc_json,
        ["ase_v025_qc.json"],
        required=not args.allow_unvalidated_input,
    )
    assignment_path = resolve_input_file(
        args.input_dir,
        args.assignment_audit,
        [os.path.join("audit", "ase_measurement_unique.tsv"), os.path.join("audit", "ase_measurement_unique.tsv.gz")],
        required=False,
    )
    source_validation, source_qc, warnings = validate_source(
        validation_path, qc_path, args.allow_unvalidated_input
    )
    assignment_index, assignment_warnings = load_assignment_audit(assignment_path)
    warnings.extend(assignment_warnings[:100])
    if assignment_path is None:
        warnings.append("ASSIGNMENT_DEPENDENCY_NOT_ASSESSED")

    output_dir = prepare_output_dir(args.out_dir, args.force)
    db_fd, db_name = tempfile.mkstemp(
        prefix="ase_v026_gene_", suffix=".sqlite", dir=str(output_dir)
    )
    os.close(db_fd)
    connection = create_gene_db(db_name)
    cursor = connection.cursor()

    block_output = output_dir / "ase_haplotype_block_reclassified.tsv"
    core_output = output_dir / "ase_core_haplotype_blocks.tsv"
    supported_output = output_dir / "ase_supported_haplotype_blocks.tsv"
    sensitivity_output = output_dir / "ase_sensitivity_haplotype_blocks.tsv"
    transition_output = output_dir / "audit" / "ase_v025_to_v026_transition.tsv"

    class_counts = Counter()
    signal_counts = Counter()
    shape_counts = Counter()
    transition_counts = Counter()
    sample_counts = Counter()
    tissue_counts = Counter()
    n_rows = 0

    try:
        with open_text(block_path, "rt") as input_handle, \
                open(block_output, "w", encoding="utf-8", newline="") as block_handle, \
                open(core_output, "w", encoding="utf-8", newline="") as core_handle, \
                open(supported_output, "w", encoding="utf-8", newline="") as supported_handle, \
                open(sensitivity_output, "w", encoding="utf-8", newline="") as sensitivity_handle, \
                open(transition_output, "w", encoding="utf-8", newline="") as transition_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            source_fields = reader.fieldnames or []
            missing = REQUIRED_BLOCK_FIELDS.difference(source_fields)
            if missing:
                raise ValueError(
                    "block input missing fields: {}".format(
                        ",".join(sorted(missing))
                    )
                )
            output_fields = list(source_fields)
            for field in V026_FIELDS:
                if field not in output_fields:
                    output_fields.append(field)
            writers = {
                "all": csv.DictWriter(block_handle, fieldnames=output_fields, delimiter="\t"),
                "core": csv.DictWriter(core_handle, fieldnames=output_fields, delimiter="\t"),
                "supported": csv.DictWriter(supported_handle, fieldnames=output_fields, delimiter="\t"),
                "sensitivity": csv.DictWriter(sensitivity_handle, fieldnames=output_fields, delimiter="\t"),
            }
            transition_writer = csv.DictWriter(
                transition_handle, fieldnames=TRANSITION_FIELDS, delimiter="\t"
            )
            for writer in writers.values():
                writer.writeheader()
            transition_writer.writeheader()

            for row_number, row in enumerate(reader, start=2):
                validate_numeric_invariants(row, row_number)
                if args.sample and row.get("sample") != args.sample:
                    raise ValueError(
                        "row {} sample {} does not match --sample {}".format(
                            row_number, row.get("sample"), args.sample
                        )
                    )
                classified = classify_row(row, thresholds, assignment_index)
                writers["all"].writerow(classified)
                evidence = classified["v026_evidence_class"]
                if evidence in CORE_CLASSES:
                    writers["core"].writerow(classified)
                elif evidence in SUPPORTED_CLASSES:
                    writers["supported"].writerow(classified)
                elif evidence in SENSITIVITY_CLASSES:
                    writers["sensitivity"].writerow(classified)

                transition_writer.writerow(
                    {
                        "sample": row["sample"],
                        "tissue": row["tissue"],
                        "gene_id": row["gene_id"],
                        "gene_name": row.get("gene_name", ""),
                        "haplotype_block_id": row["haplotype_block_id"],
                        "measurement_id": row["measurement_id"],
                        "candidate_class_v025": row.get("candidate_class", ""),
                        "exclusion_reasons_v025": row.get("exclusion_reasons", ""),
                        "v026_evidence_class": evidence,
                        "v026_statistical_signal": classified[
                            "v026_statistical_signal"
                        ],
                        "v026_imbalance_shape": classified[
                            "v026_imbalance_shape"
                        ],
                        "v026_snp_support_class": classified[
                            "v026_snp_support_class"
                        ],
                        "v026_phase_qc_status": classified[
                            "v026_phase_qc_status"
                        ],
                        "v026_technical_qc_status": classified[
                            "v026_technical_qc_status"
                        ],
                        "v026_reclassification_reasons": classified[
                            "v026_reclassification_reasons"
                        ],
                    }
                )
                if classified.get("analysis_set") == "MAIN":
                    insert_gene_row(cursor, classified)
                n_rows += 1
                class_counts[evidence] += 1
                signal_counts[classified["v026_statistical_signal"]] += 1
                shape_counts[classified["v026_imbalance_shape"]] += 1
                transition_counts[(row.get("candidate_class", ""), evidence)] += 1
                sample_counts[row["sample"]] += 1
                tissue_counts[(row["sample"], row["tissue"])] += 1
                if n_rows % 10000 == 0:
                    connection.commit()
            connection.commit()

        expected_source_rows = None
        if source_validation:
            expected_source_rows = source_validation.get("row_counts", {}).get("blocks")
        if expected_source_rows is not None and n_rows != expected_source_rows:
            raise ValueError(
                "source row count mismatch: observed={} validated={}".format(
                    n_rows, expected_source_rows
                )
            )

        n_genes, gene_status_counts = write_gene_outputs(connection, output_dir)
    finally:
        connection.close()
        try:
            os.remove(db_name)
        except FileNotFoundError:
            pass

    transition_count_path = output_dir / "audit" / "ase_reclassification_reason_counts.tsv"
    with open(transition_count_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["candidate_class_v025", "v026_evidence_class", "n_rows"])
        for (old_class, new_class), count in sorted(transition_counts.items()):
            writer.writerow([old_class, new_class, count])

    metadata = {
        "version": VERSION,
        "started_at": started_at,
        "completed_at": utc_now(),
        "source_block_results": str(block_path.resolve()),
        "source_block_results_sha256": sha256_file(block_path),
        "source_validation_json": str(validation_path.resolve()) if validation_path else None,
        "source_qc_json": str(qc_path.resolve()) if qc_path else None,
        "assignment_audit": str(assignment_path.resolve()) if assignment_path else None,
        "source_validation_status": (
            source_validation.get("status") if source_validation else "NOT_PROVIDED"
        ),
        "source_qc_loaded": source_qc is not None,
        "thresholds": thresholds,
        "n_block_gene_rows": n_rows,
        "n_gene_summary_rows": n_genes,
        "n_assignment_measurements": len(assignment_index),
        "block_evidence_class_counts": dict(sorted(class_counts.items())),
        "statistical_signal_counts": dict(sorted(signal_counts.items())),
        "imbalance_shape_counts": dict(sorted(shape_counts.items())),
        "gene_status_counts": dict(sorted(gene_status_counts.items())),
        "sample_row_counts": dict(sorted(sample_counts.items())),
        "sample_tissue_row_counts": {
            "{}|{}".format(sample, tissue): count
            for (sample, tissue), count in sorted(tissue_counts.items())
        },
        "warnings": warnings[:100],
        "command": sys.argv,
    }
    with open(output_dir / "ase_v026_qc.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def build_parser():
    parser = argparse.ArgumentParser(
        description="Reclassify ASE v0.2.5 block results as ASE v0.2.6"
    )
    parser.add_argument("--input-dir")
    parser.add_argument("--block-results")
    parser.add_argument("--validation-json")
    parser.add_argument("--qc-json")
    parser.add_argument("--assignment-audit")
    parser.add_argument("--sample")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--allow-unvalidated-input", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--near-fdr", type=float, default=0.10)
    parser.add_argument("--nominal-p", type=float, default=0.05)
    parser.add_argument("--strong-major-fraction", type=float, default=0.65)
    parser.add_argument("--moderate-major-fraction", type=float, default=0.60)
    parser.add_argument(
        "--max-fragment-hap-conflict-fraction", type=float, default=0.10
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        metadata = run_reclassification(args)
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        return 2
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
