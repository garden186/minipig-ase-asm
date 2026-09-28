#!/usr/bin/env python3
"""Build cohort recurrence tables from validated ASE v0.2.6 gene statuses."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


VERSION = "ase_candidates_v0.2.6_cohort_v0.3.0"

DEFAULT_SAMPLES = [
    "0235", "0242", "0276", "0309", "0326",
    "0349", "0355", "0377", "0384", "0385",
]
DEFAULT_TISSUE_ALIASES = {"L-blood": "Blood", "blood": "Blood"}

CORE_STATUS = "ASE_GENE_CORE"
SUPPORTED_STATUSES = {
    "ASE_GENE_SUPPORTED_SINGLE_SNP",
    "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED",
}
SENSITIVITY_STATUSES = {
    "ASE_GENE_FDR_MODERATE",
    "ASE_GENE_NEAR_FDR",
    "ASE_GENE_NOMINAL",
    "ASE_GENE_STATISTICAL_SMALL_EFFECT",
}
NEGATIVE_STATUSES = {"NOT_ASE_GENE", "UNTESTABLE"}
ALLOWED_STATUSES = (
    {CORE_STATUS} | SUPPORTED_STATUSES | SENSITIVITY_STATUSES | NEGATIVE_STATUSES
)
SELECTED_STATUSES = {CORE_STATUS} | SUPPORTED_STATUSES | SENSITIVITY_STATUSES

GENE_HEADER = [
    "sample", "tissue", "gene_id", "gene_name", "gene_biotype",
    "n_informative_haplotype_blocks", "n_testable_haplotype_blocks",
    "n_fdr_strong_blocks", "n_core_robust_blocks",
    "n_core_supported_blocks", "n_supported_single_snp_blocks",
    "n_phase_unresolved_blocks", "n_support_unresolved_blocks",
    "n_fdr_moderate_effect_blocks", "n_near_fdr_blocks",
    "n_nominal_directional_blocks", "n_statistical_small_effect_blocks",
    "n_monoallelic_like_blocks", "n_low_minor_blocks",
    "n_two_haplotype_supported_blocks", "gene_phase_scope",
    "gene_ase_status_v026", "top_haplotype_block_id", "top_measurement_id",
    "top_evidence_class_v026", "top_statistical_signal_v026",
    "top_imbalance_shape_v026", "top_snp_support_class_v026",
    "top_phase_qc_status_v026", "top_assignment_dependency_v026",
    "top_asm_integration_readiness_v026", "top_ase_direction_v026",
    "top_ase_q_bh", "top_major_haplotype_fraction", "top_total_count",
    "core_haplotype_block_ids", "supported_haplotype_block_ids",
    "sensitivity_haplotype_block_ids",
]

RECURRENCE_HEADER = [
    "gene_id", "gene_name", "gene_biotype", "tissue",
    "n_samples_testable", "n_samples_core", "samples_core",
    "expected_core_under_null", "recurrence_p", "recurrence_q_bh",
    "is_statistically_recurrent_core",
]

STATUS_COUNT_HEADER = [
    "gene_id", "gene_name", "gene_biotype", "tissue",
    "n_samples_testable",
    "n_samples_core", "samples_core",
    "n_samples_supported_single_snp", "samples_supported_single_snp",
    "n_samples_phase_or_support_unresolved",
    "samples_phase_or_support_unresolved",
    "n_samples_fdr_moderate", "samples_fdr_moderate",
    "n_samples_near_fdr", "samples_near_fdr",
    "n_samples_nominal", "samples_nominal",
    "n_samples_statistical_small_effect",
    "samples_statistical_small_effect",
]

GENE_OVERVIEW_HEADER = [
    "gene_id", "gene_name", "gene_biotype",
    "n_selected_measurements", "n_core_measurements",
    "n_supported_single_snp_measurements",
    "n_phase_or_support_unresolved_measurements",
    "n_fdr_moderate_measurements", "n_near_fdr_measurements",
    "n_nominal_measurements", "n_statistical_small_effect_measurements",
    "n_samples_any_selected_tissue", "samples_any_selected_tissue",
    "n_tissues_with_selected_ase", "tissues_with_selected_ase",
    "n_samples_any_core_tissue", "samples_any_core_tissue",
    "n_tissues_with_core_ase", "tissues_with_core_ase",
    "n_tested_gene_tissue_pairs",
    "n_statistically_recurrent_core_gene_tissue_pairs",
    "min_core_recurrence_q_bh", "statistically_recurrent_core_tissues",
]

BURDEN_HEADER = [
    "sample", "tissue", "n_testable_genes", "n_core_genes",
    "core_null_probability", "n_supported_single_snp_genes",
    "n_phase_or_support_unresolved_genes", "n_fdr_moderate_genes",
    "n_near_fdr_genes", "n_nominal_genes",
    "n_statistical_small_effect_genes",
]

TISSUE_SUMMARY_HEADER = [
    "tissue", "n_samples_available", "sample_ids_available",
    "n_testable_gene_measurements", "n_core_gene_measurements",
    "n_supported_single_snp_gene_measurements",
    "n_phase_or_support_unresolved_gene_measurements",
    "n_fdr_moderate_gene_measurements", "n_near_fdr_gene_measurements",
    "n_nominal_gene_measurements",
    "n_statistical_small_effect_gene_measurements",
    "n_unique_core_genes", "n_recurrence_hypotheses",
    "n_statistically_recurrent_core_genes",
]

OUTPUT_SPECS = [
    ("core_gene_measurements.tsv", GENE_HEADER),
    ("supported_gene_measurements.tsv", GENE_HEADER),
    ("sensitivity_gene_measurements.tsv", GENE_HEADER),
    ("gene_tissue_core_recurrence_test_results.tsv", RECURRENCE_HEADER),
    ("statistically_recurrent_core_gene_tissue.tsv", RECURRENCE_HEADER),
    ("gene_tissue_status_counts.tsv", STATUS_COUNT_HEADER),
    ("gene_overview.tsv", GENE_OVERVIEW_HEADER),
    ("sample_tissue_burden.tsv", BURDEN_HEADER),
    ("tissue_summary.tsv", TISSUE_SUMMARY_HEADER),
]


class CohortError(RuntimeError):
    pass


def parse_csv_list(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def parse_tissue_aliases(text: str) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for item in parse_csv_list(text):
        if "=" not in item:
            raise CohortError(
                f"invalid tissue alias {item!r}; expected source=canonical"
            )
        source, target = (part.strip() for part in item.split("=", 1))
        if not source or not target or source == target:
            raise CohortError(f"invalid tissue alias: {item!r}")
        if source in aliases and aliases[source] != target:
            raise CohortError(f"conflicting tissue aliases for {source!r}")
        aliases[source] = target
    chained = sorted(set(aliases.values()) & set(aliases))
    if chained:
        raise CohortError(f"chained tissue aliases are not allowed: {chained}")
    return aliases


def normalize_tissue(tissue: str, aliases: Mapping[str, str]) -> str:
    return aliases.get(tissue, tissue)


def normalize_row(row: Mapping[str, str], aliases: Mapping[str, str]) -> Dict[str, str]:
    result = dict(row)
    result["tissue"] = normalize_tissue(result["tissue"], aliases)
    return result


def row_key(row: Mapping[str, str]) -> Tuple[str, str, str]:
    return row["sample"], row["tissue"], row["gene_id"]


def read_tsv(path: Path, expected_header: Sequence[str]) -> List[Dict[str, str]]:
    if not path.is_file():
        raise CohortError(f"missing input file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(expected_header):
            raise CohortError(
                f"unexpected header in {path}; "
                f"expected={list(expected_header)}, observed={reader.fieldnames}"
            )
        return [dict(row) for row in reader]


def write_tsv(
    path: Path,
    header: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(header),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in header})


def to_int(value: object, field: str, context: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise CohortError(f"invalid integer {field}={value!r} at {context}") from exc


def fmt_fraction(value: float) -> str:
    return f"{value:.6f}"


def fmt_probability(value: float) -> str:
    return f"{value:.12g}"


def poisson_binomial_tail(probabilities: Sequence[float], observed: int) -> float:
    if observed <= 0:
        return 1.0
    if observed > len(probabilities):
        return 0.0
    pmf = [1.0]
    for probability in probabilities:
        if not 0.0 <= probability <= 1.0:
            raise CohortError(f"invalid null probability: {probability}")
        updated = [0.0] * (len(pmf) + 1)
        for count, mass in enumerate(pmf):
            updated[count] += mass * (1.0 - probability)
            updated[count + 1] += mass * probability
        pmf = updated
    return min(1.0, max(0.0, math.fsum(pmf[observed:])))


def benjamini_hochberg(
    p_values: Mapping[Tuple[str, str], float],
) -> Dict[Tuple[str, str], float]:
    ordered = sorted(p_values.items(), key=lambda item: (item[1], item[0]))
    total = len(ordered)
    adjusted: Dict[Tuple[str, str], float] = {}
    running = 1.0
    for rank in range(total, 0, -1):
        key, p_value = ordered[rank - 1]
        running = min(running, min(1.0, p_value * total / rank))
        adjusted[key] = running
    return adjusted


def canonical_annotation(
    gene_id: str,
    annotations: Mapping[str, Mapping[str, set]],
    field: str,
) -> str:
    values = {
        value
        for value in annotations.get(gene_id, {}).get(field, set())
        if value not in {"", "."}
    }
    if len(values) > 1:
        raise CohortError(
            f"inconsistent {field} for gene_id={gene_id}: {sorted(values)}"
        )
    return next(iter(values)) if values else ""


def read_json(path: Path) -> Dict[str, object]:
    if not path.is_file():
        raise CohortError(f"missing JSON input: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_source_jsons(
    sample_dir: Path,
    sample: str,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    validation = read_json(sample_dir / "ase_v026_validation.json")
    qc = read_json(sample_dir / "ase_v026_qc.json")
    if validation.get("status") != "PASS" or int(validation.get("n_errors", -1)) != 0:
        raise CohortError(f"{sample}: ASE v0.2.6 validation did not PASS")
    if str(validation.get("version")) != "0.2.6":
        raise CohortError(f"{sample}: unexpected validation version")
    if str(qc.get("version")) != "0.2.6":
        raise CohortError(f"{sample}: unexpected QC version")
    if qc.get("source_validation_status") != "PASS":
        raise CohortError(f"{sample}: source v0.2.5 validation was not PASS")
    if qc.get("source_qc_loaded") is not True:
        raise CohortError(f"{sample}: source v0.2.5 QC was not loaded")
    return validation, qc


def validate_gene_row(row: Mapping[str, str], sample: str) -> None:
    context = "|".join(row_key(row))
    if row["sample"] != sample:
        raise CohortError(
            f"sample directory={sample}, gene row sample={row['sample']} at {context}"
        )
    status = row["gene_ase_status_v026"]
    if status not in ALLOWED_STATUSES:
        raise CohortError(f"invalid v0.2.6 gene status {status!r} at {context}")
    numeric_fields = [
        "n_informative_haplotype_blocks", "n_testable_haplotype_blocks",
        "n_fdr_strong_blocks", "n_core_robust_blocks",
        "n_core_supported_blocks", "n_supported_single_snp_blocks",
        "n_phase_unresolved_blocks", "n_support_unresolved_blocks",
        "n_fdr_moderate_effect_blocks", "n_near_fdr_blocks",
        "n_nominal_directional_blocks", "n_statistical_small_effect_blocks",
        "n_monoallelic_like_blocks", "n_low_minor_blocks",
        "n_two_haplotype_supported_blocks",
    ]
    parsed = {field: to_int(row[field], field, context) for field in numeric_fields}
    if any(value < 0 for value in parsed.values()):
        raise CohortError(f"negative gene block count at {context}")
    if status in SELECTED_STATUSES and parsed["n_testable_haplotype_blocks"] < 1:
        raise CohortError(f"selected ASE gene has no testable block at {context}")
    if status == CORE_STATUS and (
        parsed["n_core_robust_blocks"] + parsed["n_core_supported_blocks"] < 1
    ):
        raise CohortError(f"core gene has no core block at {context}")
    if status == "ASE_GENE_SUPPORTED_SINGLE_SNP" and parsed[
        "n_supported_single_snp_blocks"
    ] < 1:
        raise CohortError(f"single-SNP supported gene lacks supporting block at {context}")
    if status == "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED" and (
        parsed["n_phase_unresolved_blocks"] + parsed["n_support_unresolved_blocks"] < 1
    ):
        raise CohortError(f"phase/support unresolved gene lacks matching block at {context}")


def exact_subset(
    summary: Sequence[Mapping[str, str]],
    subset: Sequence[Mapping[str, str]],
    statuses: set,
    label: str,
    sample: str,
) -> None:
    expected = {row_key(row): dict(row) for row in summary if row["gene_ase_status_v026"] in statuses}
    observed = {row_key(row): dict(row) for row in subset}
    if len(observed) != len(subset):
        raise CohortError(f"{sample}: duplicate rows in {label}")
    if expected != observed:
        missing = sorted(set(expected) - set(observed))[:10]
        extra = sorted(set(observed) - set(expected))[:10]
        changed = sorted(
            key
            for key in set(expected) & set(observed)
            if expected[key] != observed[key]
        )[:10]
        raise CohortError(
            f"{sample}: {label} is not exact status subset; "
            f"missing={missing}, extra={extra}, changed={changed}"
        )


def load_inputs(
    input_root: Path,
    samples: Sequence[str],
    excluded_tissues: set,
    tissue_aliases: Mapping[str, str],
) -> Dict[str, object]:
    core_rows: List[Dict[str, str]] = []
    supported_rows: List[Dict[str, str]] = []
    sensitivity_rows: List[Dict[str, str]] = []
    selected_rows: List[Dict[str, str]] = []
    sample_tissue_combinations = set()
    normalized_tissue_samples: MutableMapping[Tuple[str, str], set] = defaultdict(set)
    annotations: MutableMapping[str, Dict[str, set]] = defaultdict(
        lambda: {"gene_name": set(), "gene_biotype": set()}
    )
    testable_samples: MutableMapping[Tuple[str, str], set] = defaultdict(set)
    status_samples: MutableMapping[Tuple[str, str, str], set] = defaultdict(set)
    sample_tissue_counts: MutableMapping[Tuple[str, str], Counter] = defaultdict(Counter)
    available_samples_by_tissue: MutableMapping[str, set] = defaultdict(set)
    validations: Dict[str, Dict[str, object]] = {}
    qcs: Dict[str, Dict[str, object]] = {}
    threshold_signatures = set()
    seen_normalized_keys = set()

    for sample in samples:
        sample_dir = input_root / sample
        validation, qc = validate_source_jsons(sample_dir, sample)
        validations[sample] = validation
        qcs[sample] = qc
        threshold_signatures.add(json.dumps(qc.get("thresholds", {}), sort_keys=True))

        summary_source = read_tsv(sample_dir / "ase_gene_summary.tsv", GENE_HEADER)
        core_source = read_tsv(sample_dir / "ase_core_genes.tsv", GENE_HEADER)
        supported_source = read_tsv(sample_dir / "ase_supported_genes.tsv", GENE_HEADER)
        sensitivity_source = read_tsv(sample_dir / "ase_sensitivity_genes.tsv", GENE_HEADER)

        observed_status = Counter()
        source_keys = set()
        for row in summary_source:
            validate_gene_row(row, sample)
            key = row_key(row)
            if key in source_keys:
                raise CohortError(f"{sample}: duplicate gene summary row {key}")
            source_keys.add(key)
            observed_status[row["gene_ase_status_v026"]] += 1

        if int(validation.get("n_gene_summary_rows", -1)) != len(summary_source):
            raise CohortError(f"{sample}: gene summary row count differs from validation")
        expected_status = {
            str(key): int(value)
            for key, value in validation.get("gene_status_counts", {}).items()
        }
        if dict(sorted(observed_status.items())) != dict(sorted(expected_status.items())):
            raise CohortError(f"{sample}: gene status counts differ from validation")
        qc_status = {
            str(key): int(value)
            for key, value in qc.get("gene_status_counts", {}).items()
        }
        if dict(sorted(observed_status.items())) != dict(sorted(qc_status.items())):
            raise CohortError(f"{sample}: gene status counts differ from QC")

        exact_subset(summary_source, core_source, {CORE_STATUS}, "ase_core_genes.tsv", sample)
        exact_subset(
            summary_source,
            supported_source,
            SUPPORTED_STATUSES,
            "ase_supported_genes.tsv",
            sample,
        )
        exact_subset(
            summary_source,
            sensitivity_source,
            SENSITIVITY_STATUSES,
            "ase_sensitivity_genes.tsv",
            sample,
        )

        for source_row in summary_source:
            source_tissue = source_row["tissue"]
            canonical_tissue = normalize_tissue(source_tissue, tissue_aliases)
            if source_tissue != canonical_tissue:
                normalized_tissue_samples[(source_tissue, canonical_tissue)].add(sample)
            row = normalize_row(source_row, tissue_aliases)
            if row["tissue"] in excluded_tissues:
                continue
            key = row_key(row)
            if key in seen_normalized_keys:
                raise CohortError(f"duplicate gene after tissue normalization: {key}")
            seen_normalized_keys.add(key)
            sample_tissue = (sample, row["tissue"])
            sample_tissue_combinations.add(sample_tissue)
            available_samples_by_tissue[row["tissue"]].add(sample)
            counts = sample_tissue_counts[sample_tissue]
            status = row["gene_ase_status_v026"]
            counts[status] += 1
            n_testable = to_int(
                row["n_testable_haplotype_blocks"],
                "n_testable_haplotype_blocks",
                "|".join(key),
            )
            if n_testable >= 1:
                counts["TESTABLE"] += 1
                gene_tissue = (row["gene_id"], row["tissue"])
                testable_samples[gene_tissue].add(sample)
                annotations[row["gene_id"]]["gene_name"].add(row["gene_name"])
                annotations[row["gene_id"]]["gene_biotype"].add(row["gene_biotype"])
            if status in SELECTED_STATUSES:
                status_samples[(row["gene_id"], row["tissue"], status)].add(sample)
                selected_rows.append(row)
                if status == CORE_STATUS:
                    core_rows.append(row)
                elif status in SUPPORTED_STATUSES:
                    supported_rows.append(row)
                else:
                    sensitivity_rows.append(row)

    if len(threshold_signatures) != 1:
        raise CohortError("v0.2.6 threshold sets differ among samples")

    for row in core_rows:
        if to_int(
            row["n_testable_haplotype_blocks"],
            "n_testable_haplotype_blocks",
            "|".join(row_key(row)),
        ) < 1:
            raise CohortError(f"core gene is not testable: {row_key(row)}")

    return {
        "core_rows": core_rows,
        "supported_rows": supported_rows,
        "sensitivity_rows": sensitivity_rows,
        "selected_rows": selected_rows,
        "sample_tissue_combinations": sample_tissue_combinations,
        "normalized_tissue_samples": normalized_tissue_samples,
        "annotations": annotations,
        "testable_samples": testable_samples,
        "status_samples": status_samples,
        "sample_tissue_counts": sample_tissue_counts,
        "available_samples_by_tissue": available_samples_by_tissue,
        "validations": validations,
        "qcs": qcs,
        "thresholds": json.loads(next(iter(threshold_signatures))),
    }


def status_sample_set(
    status_samples: Mapping[Tuple[str, str, str], set],
    gene_id: str,
    tissue: str,
    status: str,
) -> set:
    return status_samples.get((gene_id, tissue, status), set())


def build_recurrence(
    loaded: Mapping[str, object],
    q_threshold: float,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Dict[Tuple[str, str], Dict[str, object]]]:
    if not 0.0 < q_threshold <= 1.0:
        raise CohortError(f"recurrence q threshold must be in (0,1]: {q_threshold}")
    testable_samples = loaded["testable_samples"]
    status_samples = loaded["status_samples"]
    sample_tissue_counts = loaded["sample_tissue_counts"]
    annotations = loaded["annotations"]

    null_by_sample_tissue: Dict[Tuple[str, str], Dict[str, object]] = {}
    for sample_tissue, counts in sorted(sample_tissue_counts.items()):
        n_testable = counts["TESTABLE"]
        n_core = counts[CORE_STATUS]
        if n_testable < 1:
            continue
        if n_core > n_testable:
            raise CohortError(f"core burden exceeds testable burden: {sample_tissue}")
        null_by_sample_tissue[sample_tissue] = {
            "n_testable_genes": n_testable,
            "n_core_genes": n_core,
            "core_probability": n_core / n_testable,
        }

    raw: Dict[Tuple[str, str], Dict[str, object]] = {}
    p_values: Dict[Tuple[str, str], float] = {}
    for (gene_id, tissue), testable_set in sorted(testable_samples.items()):
        samples_testable = sorted(testable_set)
        if len(samples_testable) < 2:
            continue
        core_set = status_sample_set(status_samples, gene_id, tissue, CORE_STATUS)
        if core_set - testable_set:
            raise CohortError(f"core sample is not testable for {gene_id}|{tissue}")
        probabilities = [
            float(null_by_sample_tissue[(sample, tissue)]["core_probability"])
            for sample in samples_testable
        ]
        observed = len(core_set)
        p_value = poisson_binomial_tail(probabilities, observed)
        key = (gene_id, tissue)
        p_values[key] = p_value
        raw[key] = {
            "n_samples_testable": len(samples_testable),
            "n_samples_core": observed,
            "samples_core": ",".join(sorted(core_set)),
            "expected_core_under_null": math.fsum(probabilities),
            "recurrence_p": p_value,
        }

    q_values = benjamini_hochberg(p_values)
    rows: List[Dict[str, object]] = []
    recurrent: List[Dict[str, object]] = []
    for key in sorted(raw):
        gene_id, tissue = key
        entry = raw[key]
        q_value = q_values[key]
        flag = int(entry["n_samples_core"] >= 2 and q_value <= q_threshold)
        row = {
            "gene_id": gene_id,
            "gene_name": canonical_annotation(gene_id, annotations, "gene_name"),
            "gene_biotype": canonical_annotation(gene_id, annotations, "gene_biotype"),
            "tissue": tissue,
            "n_samples_testable": entry["n_samples_testable"],
            "n_samples_core": entry["n_samples_core"],
            "samples_core": entry["samples_core"],
            "expected_core_under_null": fmt_fraction(
                float(entry["expected_core_under_null"])
            ),
            "recurrence_p": fmt_probability(float(entry["recurrence_p"])),
            "recurrence_q_bh": fmt_probability(q_value),
            "is_statistically_recurrent_core": flag,
        }
        rows.append(row)
        if flag:
            recurrent.append(dict(row))
    return rows, recurrent, null_by_sample_tissue


def build_status_count_rows(loaded: Mapping[str, object]) -> List[Dict[str, object]]:
    testable_samples = loaded["testable_samples"]
    status_samples = loaded["status_samples"]
    annotations = loaded["annotations"]
    output = []
    for (gene_id, tissue), testable_set in sorted(testable_samples.items()):
        sets = {
            status: status_sample_set(status_samples, gene_id, tissue, status)
            for status in SELECTED_STATUSES
        }
        if not any(sets.values()):
            continue
        def values(status: str) -> Tuple[int, str]:
            sample_set = sets[status]
            return len(sample_set), ",".join(sorted(sample_set))
        core_n, core_s = values(CORE_STATUS)
        single_n, single_s = values("ASE_GENE_SUPPORTED_SINGLE_SNP")
        unresolved_n, unresolved_s = values("ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED")
        moderate_n, moderate_s = values("ASE_GENE_FDR_MODERATE")
        near_n, near_s = values("ASE_GENE_NEAR_FDR")
        nominal_n, nominal_s = values("ASE_GENE_NOMINAL")
        small_n, small_s = values("ASE_GENE_STATISTICAL_SMALL_EFFECT")
        output.append({
            "gene_id": gene_id,
            "gene_name": canonical_annotation(gene_id, annotations, "gene_name"),
            "gene_biotype": canonical_annotation(gene_id, annotations, "gene_biotype"),
            "tissue": tissue,
            "n_samples_testable": len(testable_set),
            "n_samples_core": core_n,
            "samples_core": core_s,
            "n_samples_supported_single_snp": single_n,
            "samples_supported_single_snp": single_s,
            "n_samples_phase_or_support_unresolved": unresolved_n,
            "samples_phase_or_support_unresolved": unresolved_s,
            "n_samples_fdr_moderate": moderate_n,
            "samples_fdr_moderate": moderate_s,
            "n_samples_near_fdr": near_n,
            "samples_near_fdr": near_s,
            "n_samples_nominal": nominal_n,
            "samples_nominal": nominal_s,
            "n_samples_statistical_small_effect": small_n,
            "samples_statistical_small_effect": small_s,
        })
    return output


def build_gene_overview(
    loaded: Mapping[str, object],
    recurrence_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    annotations = loaded["annotations"]
    by_gene: MutableMapping[str, List[Mapping[str, str]]] = defaultdict(list)
    recurrence_by_gene: MutableMapping[str, List[Mapping[str, object]]] = defaultdict(list)
    for row in loaded["selected_rows"]:
        by_gene[row["gene_id"]].append(row)
    for row in recurrence_rows:
        recurrence_by_gene[str(row["gene_id"])].append(row)
    output = []
    for gene_id, rows in sorted(by_gene.items()):
        statuses = Counter(row["gene_ase_status_v026"] for row in rows)
        core_rows = [row for row in rows if row["gene_ase_status_v026"] == CORE_STATUS]
        recurrence = recurrence_by_gene.get(gene_id, [])
        recurrent = [
            row
            for row in recurrence
            if int(row["is_statistically_recurrent_core"]) == 1
        ]
        q_values = [float(row["recurrence_q_bh"]) for row in recurrence]
        samples_selected = sorted({row["sample"] for row in rows})
        tissues_selected = sorted({row["tissue"] for row in rows})
        samples_core = sorted({row["sample"] for row in core_rows})
        tissues_core = sorted({row["tissue"] for row in core_rows})
        output.append({
            "gene_id": gene_id,
            "gene_name": canonical_annotation(gene_id, annotations, "gene_name"),
            "gene_biotype": canonical_annotation(gene_id, annotations, "gene_biotype"),
            "n_selected_measurements": len(rows),
            "n_core_measurements": statuses[CORE_STATUS],
            "n_supported_single_snp_measurements": statuses[
                "ASE_GENE_SUPPORTED_SINGLE_SNP"
            ],
            "n_phase_or_support_unresolved_measurements": statuses[
                "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED"
            ],
            "n_fdr_moderate_measurements": statuses["ASE_GENE_FDR_MODERATE"],
            "n_near_fdr_measurements": statuses["ASE_GENE_NEAR_FDR"],
            "n_nominal_measurements": statuses["ASE_GENE_NOMINAL"],
            "n_statistical_small_effect_measurements": statuses[
                "ASE_GENE_STATISTICAL_SMALL_EFFECT"
            ],
            "n_samples_any_selected_tissue": len(samples_selected),
            "samples_any_selected_tissue": ",".join(samples_selected),
            "n_tissues_with_selected_ase": len(tissues_selected),
            "tissues_with_selected_ase": ",".join(tissues_selected),
            "n_samples_any_core_tissue": len(samples_core),
            "samples_any_core_tissue": ",".join(samples_core),
            "n_tissues_with_core_ase": len(tissues_core),
            "tissues_with_core_ase": ",".join(tissues_core),
            "n_tested_gene_tissue_pairs": len(recurrence),
            "n_statistically_recurrent_core_gene_tissue_pairs": len(recurrent),
            "min_core_recurrence_q_bh": (
                fmt_probability(min(q_values)) if q_values else ""
            ),
            "statistically_recurrent_core_tissues": ",".join(
                sorted(str(row["tissue"]) for row in recurrent)
            ),
        })
    return output


def build_burden_rows(
    loaded: Mapping[str, object],
    null_by_sample_tissue: Mapping[Tuple[str, str], Mapping[str, object]],
) -> List[Dict[str, object]]:
    output = []
    for sample_tissue in sorted(loaded["sample_tissue_combinations"]):
        sample, tissue = sample_tissue
        counts = loaded["sample_tissue_counts"][sample_tissue]
        null = null_by_sample_tissue.get(sample_tissue)
        output.append({
            "sample": sample,
            "tissue": tissue,
            "n_testable_genes": counts["TESTABLE"],
            "n_core_genes": counts[CORE_STATUS],
            "core_null_probability": (
                fmt_fraction(float(null["core_probability"])) if null else ""
            ),
            "n_supported_single_snp_genes": counts[
                "ASE_GENE_SUPPORTED_SINGLE_SNP"
            ],
            "n_phase_or_support_unresolved_genes": counts[
                "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED"
            ],
            "n_fdr_moderate_genes": counts["ASE_GENE_FDR_MODERATE"],
            "n_near_fdr_genes": counts["ASE_GENE_NEAR_FDR"],
            "n_nominal_genes": counts["ASE_GENE_NOMINAL"],
            "n_statistical_small_effect_genes": counts[
                "ASE_GENE_STATISTICAL_SMALL_EFFECT"
            ],
        })
    return output


def build_tissue_summary(
    loaded: Mapping[str, object],
    recurrence_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    recurrence_by_tissue: MutableMapping[str, List[Mapping[str, object]]] = defaultdict(list)
    for row in recurrence_rows:
        recurrence_by_tissue[str(row["tissue"])].append(row)
    output = []
    for tissue, samples in sorted(loaded["available_samples_by_tissue"].items()):
        counts = Counter()
        core_genes = set()
        for sample in samples:
            current = loaded["sample_tissue_counts"].get((sample, tissue), Counter())
            counts.update(current)
        for row in loaded["core_rows"]:
            if row["tissue"] == tissue:
                core_genes.add(row["gene_id"])
        recurrence = recurrence_by_tissue.get(tissue, [])
        output.append({
            "tissue": tissue,
            "n_samples_available": len(samples),
            "sample_ids_available": ",".join(sorted(samples)),
            "n_testable_gene_measurements": counts["TESTABLE"],
            "n_core_gene_measurements": counts[CORE_STATUS],
            "n_supported_single_snp_gene_measurements": counts[
                "ASE_GENE_SUPPORTED_SINGLE_SNP"
            ],
            "n_phase_or_support_unresolved_gene_measurements": counts[
                "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED"
            ],
            "n_fdr_moderate_gene_measurements": counts["ASE_GENE_FDR_MODERATE"],
            "n_near_fdr_gene_measurements": counts["ASE_GENE_NEAR_FDR"],
            "n_nominal_gene_measurements": counts["ASE_GENE_NOMINAL"],
            "n_statistical_small_effect_gene_measurements": counts[
                "ASE_GENE_STATISTICAL_SMALL_EFFECT"
            ],
            "n_unique_core_genes": len(core_genes),
            "n_recurrence_hypotheses": len(recurrence),
            "n_statistically_recurrent_core_genes": sum(
                int(row["is_statistically_recurrent_core"]) == 1
                for row in recurrence
            ),
        })
    return output


def build_expected(
    input_root: Path,
    samples: Sequence[str],
    excluded_tissues: set,
    tissue_aliases: Mapping[str, str] | None = None,
    expect_sample_tissue_combinations: int = 0,
    expect_tissues: int = 0,
    recurrence_q_threshold: float = 0.05,
) -> Tuple[Dict[str, List[Dict[str, object]]], Dict[str, object]]:
    tissue_aliases = dict(tissue_aliases or {})
    excluded_tissues = {
        normalize_tissue(tissue, tissue_aliases) for tissue in excluded_tissues
    }
    loaded = load_inputs(input_root, samples, excluded_tissues, tissue_aliases)
    observed_combinations = len(loaded["sample_tissue_combinations"])
    if expect_sample_tissue_combinations and observed_combinations != expect_sample_tissue_combinations:
        raise CohortError(
            "nonexcluded sample-tissue denominator mismatch: "
            f"observed={observed_combinations}, expected={expect_sample_tissue_combinations}"
        )
    if expect_tissues and len(loaded["available_samples_by_tissue"]) != expect_tissues:
        raise CohortError(
            "nonexcluded tissue count mismatch: "
            f"observed={len(loaded['available_samples_by_tissue'])}, expected={expect_tissues}"
        )

    recurrence_rows, recurrent_rows, null_by_sample_tissue = build_recurrence(
        loaded, recurrence_q_threshold
    )
    status_count_rows = build_status_count_rows(loaded)
    gene_overview_rows = build_gene_overview(loaded, recurrence_rows)
    burden_rows = build_burden_rows(loaded, null_by_sample_tissue)
    tissue_summary_rows = build_tissue_summary(loaded, recurrence_rows)
    for key in ("core_rows", "supported_rows", "sensitivity_rows"):
        loaded[key].sort(key=lambda row: (row["sample"], row["tissue"], row["gene_id"]))

    outputs = {
        "core_gene_measurements.tsv": loaded["core_rows"],
        "supported_gene_measurements.tsv": loaded["supported_rows"],
        "sensitivity_gene_measurements.tsv": loaded["sensitivity_rows"],
        "gene_tissue_core_recurrence_test_results.tsv": recurrence_rows,
        "statistically_recurrent_core_gene_tissue.tsv": recurrent_rows,
        "gene_tissue_status_counts.tsv": status_count_rows,
        "gene_overview.tsv": gene_overview_rows,
        "sample_tissue_burden.tsv": burden_rows,
        "tissue_summary.tsv": tissue_summary_rows,
    }
    unique_core_genes = {row["gene_id"] for row in loaded["core_rows"]}
    recurrent_genes = {row["gene_id"] for row in recurrent_rows}
    metadata = {
        "version": VERSION,
        "input_root": str(input_root.resolve()),
        "samples": list(samples),
        "excluded_tissues": sorted(excluded_tissues),
        "tissue_aliases": dict(sorted(tissue_aliases.items())),
        "tissue_normalization": {
            f"{source}->{target}": {
                "n_samples": len(sample_set),
                "samples": sorted(sample_set),
            }
            for (source, target), sample_set in sorted(
                loaded["normalized_tissue_samples"].items()
            )
        },
        "n_sample_tissue_combinations": observed_combinations,
        "n_tissues": len(loaded["available_samples_by_tissue"]),
        "n_unique_core_gene_ids": len(unique_core_genes),
        "n_core_gene_measurements": len(loaded["core_rows"]),
        "n_supported_gene_measurements": len(loaded["supported_rows"]),
        "n_sensitivity_gene_measurements": len(loaded["sensitivity_rows"]),
        "n_core_gene_tissue_pairs": len({
            (row["gene_id"], row["tissue"]) for row in loaded["core_rows"]
        }),
        "n_statistical_recurrence_hypotheses": len(recurrence_rows),
        "n_statistically_recurrent_core_gene_tissue_pairs": len(recurrent_rows),
        "n_genes_with_statistically_recurrent_core": len(recurrent_genes),
        "row_counts": {name: len(rows) for name, rows in outputs.items()},
        "gene_status_measurement_counts": dict(sorted(Counter(
            row["gene_ase_status_v026"] for row in loaded["selected_rows"]
        ).items())),
        "v026_thresholds": loaded["thresholds"],
        "statistical_recurrence": {
            "main_label": CORE_STATUS,
            "null_model": (
                "sample-wise core-label permutation among testable genes; "
                "exact marginal Poisson-binomial tail"
            ),
            "testable_definition": "n_testable_haplotype_blocks >= 1",
            "minimum_testable_samples": 2,
            "minimum_core_samples_for_recurrence_call": 2,
            "bh_scope": (
                "global across all gene-tissue hypotheses with at least "
                "2 testable samples"
            ),
            "q_threshold": recurrence_q_threshold,
            "supported_and_sensitivity_policy": (
                "descriptive only; excluded from the core recurrence null burden"
            ),
            "sample_tissue_null_rates": {
                f"{sample}|{tissue}": {
                    "n_testable_genes": int(entry["n_testable_genes"]),
                    "n_core_genes": int(entry["n_core_genes"]),
                    "core_probability": float(entry["core_probability"]),
                }
                for (sample, tissue), entry in sorted(null_by_sample_tissue.items())
            },
        },
        "source_validation": {
            sample: {
                "status": validation.get("status"),
                "n_errors": validation.get("n_errors"),
                "n_warnings": validation.get("n_warnings"),
                "n_gene_summary_rows": validation.get("n_gene_summary_rows"),
            }
            for sample, validation in loaded["validations"].items()
        },
        "source_qc": {
            sample: {
                "source_validation_status": qc.get("source_validation_status"),
                "source_qc_loaded": qc.get("source_qc_loaded"),
                "source_block_results_sha256": qc.get("source_block_results_sha256"),
            }
            for sample, qc in loaded["qcs"].items()
        },
    }
    return outputs, metadata


def make_backup_path(out_dir: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = out_dir.with_name(f"{out_dir.name}.backup.{stamp}")
    serial = 1
    while candidate.exists():
        candidate = out_dir.with_name(f"{out_dir.name}.backup.{stamp}.{serial}")
        serial += 1
    return candidate


def materialize_outputs(
    outputs: Mapping[str, Sequence[Mapping[str, object]]],
    metadata: Mapping[str, object],
    out_dir: Path,
    project_tmp: Path,
    force: bool,
) -> Path | None:
    project_tmp.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="ase_v026_cohort.", dir=project_tmp))
    try:
        for name, header in OUTPUT_SPECS:
            write_tsv(staging / name, header, outputs[name])
        with (staging / "cohort_qc.json").open("w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")
        backup = None
        if out_dir.exists():
            if not force:
                raise CohortError(
                    f"output directory already exists: {out_dir}; use --force 1"
                )
            backup = make_backup_path(out_dir)
            out_dir.rename(backup)
        out_dir.parent.mkdir(parents=True, exist_ok=True)
        staging.rename(out_dir)
        return backup
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build ASE v0.2.6 core cohort recurrence tables"
    )
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--project-tmp", required=True, type=Path)
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
    parser.add_argument("--force", type=int, choices=(0, 1), default=0)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    samples = parse_csv_list(args.samples)
    if len(samples) != len(set(samples)):
        raise CohortError(f"duplicate samples: {samples}")
    tissue_aliases = parse_tissue_aliases(args.tissue_aliases)
    excluded = {
        normalize_tissue(tissue, tissue_aliases)
        for tissue in parse_csv_list(args.exclude_tissues)
    }
    outputs, metadata = build_expected(
        args.input_root,
        samples,
        excluded,
        tissue_aliases=tissue_aliases,
        expect_sample_tissue_combinations=args.expect_sample_tissue_combinations,
        expect_tissues=args.expect_tissues,
        recurrence_q_threshold=args.recurrence_q_threshold,
    )
    backup = materialize_outputs(
        outputs,
        metadata,
        args.out_dir,
        args.project_tmp,
        force=bool(args.force),
    )
    print(f"[PASS] built ASE v0.2.6 cohort: {args.out_dir}")
    if backup:
        print(f"[BACKUP] {backup}")
    print(
        "[COUNTS] "
        f"core_measurements={metadata['n_core_gene_measurements']} "
        f"unique_core_genes={metadata['n_unique_core_gene_ids']} "
        f"hypotheses={metadata['n_statistical_recurrence_hypotheses']} "
        f"recurrent_pairs={metadata['n_statistically_recurrent_core_gene_tissue_pairs']}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CohortError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
