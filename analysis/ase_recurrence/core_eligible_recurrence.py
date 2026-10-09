#!/usr/bin/env python3
"""Outcome-independent high-confidence-eligible ASE recurrence analysis.

The analysis separates three nested concepts:

1. binomial-testable block: total haplotype-assigned fragments >= threshold;
2. high-confidence-eligible block: testable plus technical PASS, at least two
   contributing WGS-balanced SNPs, and local-phase PASS;
3. high-confidence/core ASE block: eligible plus the v0.2.6 statistical,
   effect-size, and multi-SNP support criteria.

Recurrence is evaluated at gene x tissue level among eligible animal
opportunities.  The primary null uses leave-one-gene-out sample x tissue core
burdens.  A separate opportunity-stratified null and matched-label permutation
diagnostic account for block number, fragment support, and SNP support.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import random
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "0.1.0"
CORE_SUPPORT_CLASSES = {"MULTI_SNP_LOO_ROBUST", "MULTI_SNP_SUPPORTED"}
DEFAULT_THRESHOLDS = (15, 20, 30)
REQUIRED_COLUMNS = {
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
    "ase_p_exact",
    "ase_q_bh",
    "n_gene_variants_with_reads",
    "snp_support_status",
    "local_phase_status",
}


class AnalysisError(RuntimeError):
    pass


@dataclass
class GeneMeasurement:
    sample: str
    tissue: str
    gene_id: str
    gene_name: str
    gene_biotype: str
    n_tested_blocks: int = 0
    n_eligible_blocks: int = 0
    n_core_blocks: int = 0
    max_eligible_total_fragments: int = 0
    max_eligible_contributing_snps: int = 0

    @property
    def eligible(self) -> bool:
        return self.n_eligible_blocks > 0

    @property
    def core(self) -> bool:
        return self.n_core_blocks > 0


def normalize_tissue(tissue: str) -> str:
    aliases = {
        "blood": "Blood",
        "L-blood": "Blood",
        "Blood": "Blood",
        "Tenderlo": "Tenderloin",
        "tenderloin": "Tenderloin",
        "Tenderloin": "Tenderloin",
    }
    return aliases.get(tissue, tissue)


def to_int(value: str, field: str, context: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        raise AnalysisError(f"invalid integer {field}={value!r} at {context}")


def to_float(value: str, field: str, context: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        raise AnalysisError(f"invalid float {field}={value!r} at {context}")
    if not math.isfinite(parsed):
        raise AnalysisError(f"non-finite {field}={value!r} at {context}")
    return parsed


def fmt(value: float) -> str:
    return format(value, ".12g")


def open_source_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def infer_sample(path: Path) -> str:
    if path.name.startswith("ase_haplotype_block_results.tsv"):
        return path.parent.name
    return path.name.split(".", 1)[0]


def discover_source_files(input_dir: Path) -> List[Path]:
    patterns = (
        "*.ase_haplotype_block_results.tsv.gz",
        "*.ase_haplotype_block_results.tsv",
        "*/ase_haplotype_block_results.tsv.gz",
        "*/ase_haplotype_block_results.tsv",
    )
    by_sample: Dict[str, List[Path]] = defaultdict(list)
    for pattern in patterns:
        for path in input_dir.glob(pattern):
            by_sample[infer_sample(path)].append(path)
    duplicates = {sample: paths for sample, paths in by_sample.items() if len(paths) != 1}
    if duplicates:
        detail = "; ".join(
            f"{sample}={','.join(str(path) for path in paths)}"
            for sample, paths in sorted(duplicates.items())
        )
        raise AnalysisError(f"ambiguous source files: {detail}")
    return [by_sample[sample][0] for sample in sorted(by_sample)]


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        raw = path.open("wb")
        gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        handle = io.TextIOWrapper(gz, encoding="utf-8", newline="")
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    finally:
        handle.close()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_values(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    adjusted = [1.0] * len(values)
    running = 1.0
    for reverse_index in range(len(order) - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, min(1.0, values[index] * len(values) / rank))
        adjusted[index] = running
    return adjusted


def poisson_binomial_tail(probabilities: Sequence[float], observed: int) -> float:
    if observed <= 0:
        return 1.0
    if observed > len(probabilities):
        return 0.0
    pmf = [1.0]
    for probability in probabilities:
        if not 0.0 <= probability <= 1.0:
            raise AnalysisError(f"invalid null probability: {probability}")
        updated = [0.0] * (len(pmf) + 1)
        for count, mass in enumerate(pmf):
            updated[count] += mass * (1.0 - probability)
            updated[count + 1] += mass * probability
        pmf = updated
    return min(1.0, max(0.0, math.fsum(pmf[observed:])))


def technical_pass(row: Mapping[str, str], max_conflict: float) -> bool:
    context = row.get("measurement_id", "")
    if row.get("analysis_set") != "MAIN" or row.get("count_validation") != "MATCH":
        return False
    if not row.get("gene_id") or not row.get("haplotype_block_id"):
        return False
    if to_int(row.get("n_wgs_balanced_gene_snps", ""), "n_wgs_balanced_gene_snps", context) <= 0:
        return False
    conflict = row.get("fragment_hap_conflict_fraction", "")
    if conflict == "":
        return False
    return to_float(conflict, "fragment_hap_conflict_fraction", context) <= max_conflict


def eligibility_pass(row: Mapping[str, str], threshold: int, max_conflict: float) -> bool:
    context = row.get("measurement_id", "")
    return (
        to_int(row.get("ase_total_count", ""), "ase_total_count", context) >= threshold
        and technical_pass(row, max_conflict)
        and to_int(
            row.get("n_gene_variants_with_reads", ""),
            "n_gene_variants_with_reads",
            context,
        )
        >= 2
        and row.get("local_phase_status") == "PASS"
    )


def core_pass(
    row: Mapping[str, str],
    q_value: float,
    threshold: int,
    max_conflict: float,
    fdr: float,
    min_major_fraction: float,
) -> bool:
    if not eligibility_pass(row, threshold, max_conflict):
        return False
    context = row.get("measurement_id", "")
    return (
        q_value <= fdr
        and to_float(
            row.get("major_haplotype_fraction", ""),
            "major_haplotype_fraction",
            context,
        )
        >= min_major_fraction
        and row.get("snp_support_status") in CORE_SUPPORT_CLASSES
    )


def validate_header(fieldnames: Sequence[str] | None, path: Path) -> None:
    if fieldnames is None:
        raise AnalysisError(f"empty input: {path}")
    missing = sorted(REQUIRED_COLUMNS.difference(fieldnames))
    if missing:
        raise AnalysisError(f"{path}: missing columns: {','.join(missing)}")


def build_q_maps(
    path: Path,
    sample: str,
    thresholds: Sequence[int],
) -> Tuple[Dict[int, Dict[Tuple[str, str], float]], Dict[str, object]]:
    by_threshold: Dict[int, MutableMapping[Tuple[str, str], Dict[str, Tuple[int, int, float]]]] = {
        threshold: defaultdict(dict) for threshold in thresholds
    }
    n_rows = 0
    tissue_units: Set[Tuple[str, str]] = set()
    with open_source_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(reader.fieldnames, path)
        for row in reader:
            n_rows += 1
            if row["sample"] != sample:
                raise AnalysisError(f"{path}: row sample={row['sample']} expected={sample}")
            if row["analysis_set"] != "MAIN":
                continue
            tissue_units.add((sample, normalize_tissue(row["tissue"])))
            context = row["measurement_id"]
            total = to_int(row["ase_total_count"], "ase_total_count", context)
            a_count = to_int(row["ase_a_count"], "ase_a_count", context)
            b_count = to_int(row["ase_b_count"], "ase_b_count", context)
            if a_count + b_count != total:
                raise AnalysisError(f"count conservation failure at {context}")
            for threshold in thresholds:
                if total < threshold:
                    continue
                p_value = to_float(row["ase_p_exact"], "ase_p_exact", context)
                family = (sample, row["tissue"])
                previous = by_threshold[threshold][family].get(context)
                signature = (a_count, b_count, p_value)
                if previous is not None and previous != signature:
                    raise AnalysisError(f"inconsistent duplicate measurement {context}")
                by_threshold[threshold][family][context] = signature

    q_maps: Dict[int, Dict[Tuple[str, str], float]] = {}
    family_counts: Dict[str, Dict[str, int]] = {}
    for threshold in thresholds:
        q_map: Dict[Tuple[str, str], float] = {}
        threshold_counts: Dict[str, int] = {}
        for family, measurements in sorted(by_threshold[threshold].items()):
            ids = sorted(measurements)
            q_values = bh_values([measurements[mid][2] for mid in ids])
            for mid, q_value in zip(ids, q_values):
                q_map[(family[1], mid)] = q_value
            threshold_counts[family[1]] = len(ids)
        q_maps[threshold] = q_map
        family_counts[str(threshold)] = threshold_counts
    return q_maps, {
        "n_block_gene_rows": n_rows,
        "tissue_units": sorted("|".join(item) for item in tissue_units),
        "bh_family_measurement_counts": family_counts,
    }


def add_block_to_measurement(
    measurement: GeneMeasurement,
    row: Mapping[str, str],
    eligible: bool,
    core: bool,
) -> None:
    measurement.n_tested_blocks += 1
    if not eligible:
        return
    measurement.n_eligible_blocks += 1
    context = row["measurement_id"]
    total = to_int(row["ase_total_count"], "ase_total_count", context)
    n_snps = to_int(
        row["n_gene_variants_with_reads"], "n_gene_variants_with_reads", context
    )
    measurement.max_eligible_total_fragments = max(
        measurement.max_eligible_total_fragments, total
    )
    measurement.max_eligible_contributing_snps = max(
        measurement.max_eligible_contributing_snps, n_snps
    )
    if core:
        measurement.n_core_blocks += 1


def load_sample_measurements(
    path: Path,
    sample: str,
    thresholds: Sequence[int],
    q_maps: Mapping[int, Mapping[Tuple[str, str], float]],
    max_conflict: float,
    fdr: float,
    min_major_fraction: float,
) -> Tuple[Dict[int, Dict[Tuple[str, str, str], GeneMeasurement]], Dict[str, object]]:
    sample_measurements: Dict[int, Dict[Tuple[str, str, str], GeneMeasurement]] = {
        threshold: {} for threshold in thresholds
    }
    q_comparisons = 0
    q_mismatches = 0
    q_max_abs_difference = 0.0
    with open_source_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(reader.fieldnames, path)
        for row in reader:
            if row["analysis_set"] != "MAIN":
                continue
            context = row["measurement_id"]
            total = to_int(row["ase_total_count"], "ase_total_count", context)
            tissue = normalize_tissue(row["tissue"])
            key = (sample, tissue, row["gene_id"])
            for threshold in thresholds:
                if total < threshold:
                    continue
                measurement = sample_measurements[threshold].get(key)
                if measurement is None:
                    measurement = GeneMeasurement(
                        sample=sample,
                        tissue=tissue,
                        gene_id=row["gene_id"],
                        gene_name=row["gene_name"],
                        gene_biotype=row["gene_biotype"],
                    )
                    sample_measurements[threshold][key] = measurement
                elif (
                    measurement.gene_name != row["gene_name"]
                    or measurement.gene_biotype != row["gene_biotype"]
                ):
                    raise AnalysisError(f"inconsistent gene annotation at {context}")
                q_key = (row["tissue"], context)
                if q_key not in q_maps[threshold]:
                    raise AnalysisError(f"missing recomputed q-value at {context}, threshold={threshold}")
                q_value = q_maps[threshold][q_key]
                is_eligible = eligibility_pass(row, threshold, max_conflict)
                is_core = core_pass(
                    row,
                    q_value,
                    threshold,
                    max_conflict,
                    fdr,
                    min_major_fraction,
                )
                add_block_to_measurement(measurement, row, is_eligible, is_core)
                if threshold == 15:
                    source_q = to_float(row["ase_q_bh"], "ase_q_bh", context)
                    difference = abs(source_q - q_value)
                    q_comparisons += 1
                    q_max_abs_difference = max(q_max_abs_difference, difference)
                    if difference > 1e-9:
                        q_mismatches += 1
    return sample_measurements, {
        "q_comparisons_threshold15": q_comparisons,
        "q_mismatches_threshold15": q_mismatches,
        "q_max_abs_difference_threshold15": q_max_abs_difference,
    }


def read_canonical_core(path: Path) -> Set[Tuple[str, str, str]]:
    result: Set[Tuple[str, str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["sample"], normalize_tissue(row["tissue"]), row["gene_id"])
            if key in result:
                raise AnalysisError(f"duplicate canonical core key: {key}")
            result.add(key)
    return result


def build_recurrence(
    eligible_measurements: Mapping[Tuple[str, str, str], GeneMeasurement],
    tissue_animals: Mapping[str, Set[str]],
    fdr: float,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    eligible_by_gt: MutableMapping[Tuple[str, str], Set[str]] = defaultdict(set)
    core_by_gt: MutableMapping[Tuple[str, str], Set[str]] = defaultdict(set)
    st_eligible: Counter = Counter()
    st_core: Counter = Counter()
    annotations: Dict[str, Tuple[str, str]] = {}
    for (sample, tissue, gene_id), measurement in eligible_measurements.items():
        if not measurement.eligible:
            continue
        eligible_by_gt[(gene_id, tissue)].add(sample)
        st_eligible[(sample, tissue)] += 1
        annotations[gene_id] = (measurement.gene_name, measurement.gene_biotype)
        if measurement.core:
            core_by_gt[(gene_id, tissue)].add(sample)
            st_core[(sample, tissue)] += 1

    raw: List[Dict[str, object]] = []
    p_inclusive: List[float] = []
    p_loo: List[float] = []
    for (gene_id, tissue), eligible_samples in sorted(eligible_by_gt.items()):
        if len(eligible_samples) < 2:
            continue
        core_samples = core_by_gt.get((gene_id, tissue), set()) & eligible_samples
        inclusive_probabilities: List[float] = []
        loo_probabilities: List[float] = []
        for sample in sorted(eligible_samples):
            n_eligible = st_eligible[(sample, tissue)]
            n_core = st_core[(sample, tissue)]
            y_value = int(sample in core_samples)
            inclusive_probabilities.append(n_core / n_eligible)
            if n_eligible <= 1:
                raise AnalysisError(f"cannot compute leave-one-gene-out null: {sample}|{tissue}")
            loo_probabilities.append((n_core - y_value) / (n_eligible - 1))
        observed = len(core_samples)
        inclusive_p = poisson_binomial_tail(inclusive_probabilities, observed)
        loo_p = poisson_binomial_tail(loo_probabilities, observed)
        name, biotype = annotations[gene_id]
        raw.append(
            {
                "gene_id": gene_id,
                "gene_name": name,
                "gene_biotype": biotype,
                "tissue": tissue,
                "n_tissue_animals_total": len(tissue_animals[tissue]),
                "n_samples_eligible": len(eligible_samples),
                "samples_eligible": ",".join(sorted(eligible_samples)),
                "n_samples_core": observed,
                "samples_core": ",".join(sorted(core_samples)),
                "eligible_fraction_of_tissue_animals": len(eligible_samples)
                / len(tissue_animals[tissue]),
                "expected_core_inclusive": math.fsum(inclusive_probabilities),
                "expected_core_leave_one_gene_out": math.fsum(loo_probabilities),
                "recurrence_p_inclusive": inclusive_p,
                "recurrence_p_leave_one_gene_out": loo_p,
            }
        )
        p_inclusive.append(inclusive_p)
        p_loo.append(loo_p)

    q_inclusive = bh_values(p_inclusive)
    q_loo = bh_values(p_loo)
    recurrent: List[Dict[str, object]] = []
    for row, qi, ql in zip(raw, q_inclusive, q_loo):
        row["recurrence_q_bh_inclusive"] = qi
        row["recurrence_q_bh_leave_one_gene_out"] = ql
        flag = int(int(row["n_samples_core"]) >= 2 and ql <= fdr)
        row["is_statistically_recurrent_core"] = flag
        row["replication_support"] = (
            "AT_LEAST_THREE_CORE_ANIMALS"
            if int(row["n_samples_core"]) >= 3
            else "TWO_CORE_ANIMALS"
            if int(row["n_samples_core"]) == 2
            else "FEWER_THAN_TWO_CORE_ANIMALS"
        )
        if flag:
            recurrent.append(dict(row))

    burden_rows: List[Dict[str, object]] = []
    for sample, tissue in sorted(st_eligible):
        burden_rows.append(
            {
                "sample": sample,
                "tissue": tissue,
                "n_tissue_animals_total": len(tissue_animals[tissue]),
                "n_eligible_genes": st_eligible[(sample, tissue)],
                "n_core_genes": st_core[(sample, tissue)],
                "core_probability_inclusive": st_core[(sample, tissue)]
                / st_eligible[(sample, tissue)],
            }
        )

    gene_tissues: MutableMapping[str, Set[str]] = defaultdict(set)
    for row in recurrent:
        gene_tissues[str(row["gene_id"])].add(str(row["tissue"]))
    summary = {
        "n_eligible_gene_measurements": len(eligible_measurements),
        "n_core_gene_measurements": sum(m.core for m in eligible_measurements.values()),
        "n_recurrence_hypotheses": len(raw),
        "n_recurrent_gene_tissue_pairs": len(recurrent),
        "n_recurrent_genes": len(gene_tissues),
        "n_recurrent_in_one_tissue": sum(len(value) == 1 for value in gene_tissues.values()),
        "n_recurrent_in_multiple_tissues": sum(len(value) > 1 for value in gene_tissues.values()),
        "n_recurrent_pairs_with_two_core_animals": sum(
            int(row["n_samples_core"]) == 2 for row in recurrent
        ),
        "n_recurrent_pairs_with_at_least_three_core_animals": sum(
            int(row["n_samples_core"]) >= 3 for row in recurrent
        ),
    }
    return raw, recurrent, burden_rows, summary


def block_bin(value: int) -> str:
    return "1" if value == 1 else "2" if value == 2 else "3_PLUS"


def coverage_bin(value: int) -> str:
    if value < 30:
        return "15_29"
    if value < 60:
        return "30_59"
    if value < 120:
        return "60_119"
    if value < 240:
        return "120_239"
    return "240_PLUS"


def snp_bin(value: int) -> str:
    if value == 2:
        return "2"
    if value == 3:
        return "3"
    if value <= 5:
        return "4_5"
    return "6_PLUS"


def opportunity_stratum(measurement: GeneMeasurement) -> Tuple[str, str, str]:
    return (
        block_bin(measurement.n_eligible_blocks),
        coverage_bin(measurement.max_eligible_total_fragments),
        snp_bin(measurement.max_eligible_contributing_snps),
    )


def build_matched_analysis(
    measurements: Mapping[Tuple[str, str, str], GeneMeasurement],
    tissue_animals: Mapping[str, Set[str]],
    fdr: float,
    n_permutations: int,
    seed: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    groups: MutableMapping[
        Tuple[str, str, str, str, str], List[Tuple[str, str, str]]
    ] = defaultdict(list)
    for key, measurement in measurements.items():
        bbin, cbin, sbin = opportunity_stratum(measurement)
        groups[(measurement.sample, measurement.tissue, bbin, cbin, sbin)].append(key)

    group_probability: Dict[Tuple[str, str, str], float] = {}
    strata_rows: List[Dict[str, object]] = []
    for group_key, members in sorted(groups.items()):
        sample, tissue, bbin, cbin, sbin = group_key
        n_core = sum(measurements[key].core for key in members)
        probability = n_core / len(members)
        for key in members:
            group_probability[key] = probability
        strata_rows.append(
            {
                "sample": sample,
                "tissue": tissue,
                "eligible_block_bin": bbin,
                "max_fragment_bin": cbin,
                "max_contributing_snp_bin": sbin,
                "n_eligible_gene_measurements": len(members),
                "n_core_gene_measurements": n_core,
                "core_probability": probability,
            }
        )

    eligible_by_gt: MutableMapping[Tuple[str, str], Set[str]] = defaultdict(set)
    core_by_gt: MutableMapping[Tuple[str, str], Set[str]] = defaultdict(set)
    annotations: Dict[str, Tuple[str, str]] = {}
    for key, measurement in measurements.items():
        sample, tissue, gene_id = key
        eligible_by_gt[(gene_id, tissue)].add(sample)
        annotations[gene_id] = (measurement.gene_name, measurement.gene_biotype)
        if measurement.core:
            core_by_gt[(gene_id, tissue)].add(sample)

    raw: List[Dict[str, object]] = []
    p_values: List[float] = []
    hypothesis_keys: List[Tuple[str, str]] = []
    for gene_tissue, eligible_samples in sorted(eligible_by_gt.items()):
        if len(eligible_samples) < 2:
            continue
        gene_id, tissue = gene_tissue
        core_samples = core_by_gt.get(gene_tissue, set()) & eligible_samples
        probabilities = [
            group_probability[(sample, tissue, gene_id)]
            for sample in sorted(eligible_samples)
        ]
        p_value = poisson_binomial_tail(probabilities, len(core_samples))
        name, biotype = annotations[gene_id]
        raw.append(
            {
                "gene_id": gene_id,
                "gene_name": name,
                "gene_biotype": biotype,
                "tissue": tissue,
                "n_tissue_animals_total": len(tissue_animals[tissue]),
                "n_samples_eligible": len(eligible_samples),
                "n_samples_core": len(core_samples),
                "samples_core": ",".join(sorted(core_samples)),
                "expected_core_opportunity_matched": math.fsum(probabilities),
                "recurrence_p_opportunity_matched": p_value,
            }
        )
        p_values.append(p_value)
        hypothesis_keys.append(gene_tissue)

    q_values = bh_values(p_values)
    recurrent: List[Dict[str, object]] = []
    for row, q_value in zip(raw, q_values):
        row["recurrence_q_bh_opportunity_matched"] = q_value
        flag = int(int(row["n_samples_core"]) >= 2 and q_value <= fdr)
        row["is_statistically_recurrent_core"] = flag
        row["matched_permutation_p_count_tail"] = ""
        if flag:
            recurrent.append(dict(row))

    permutation_summary: Dict[str, object] = {
        "n_permutations": n_permutations,
        "seed": seed,
        "role": "diagnostic_count-tail_calibration_not_primary_fdr",
    }
    if n_permutations > 0:
        try:
            import numpy as np
        except ImportError as exc:
            raise AnalysisError("numpy is required when --n-permutations > 0") from exc

        hypothesis_index = {key: index for index, key in enumerate(hypothesis_keys)}
        observed = np.asarray([int(row["n_samples_core"]) for row in raw], dtype=np.uint8)
        exceed = np.zeros(len(raw), dtype=np.int64)
        ge2_distribution: List[int] = []
        ge3_distribution: List[int] = []
        rng = np.random.default_rng(seed)
        permutation_groups = []
        for members in groups.values():
            n_core = sum(measurements[key].core for key in members)
            member_indices = np.asarray(
                [hypothesis_index.get((key[2], key[1]), -1) for key in members],
                dtype=np.int64,
            )
            permutation_groups.append((member_indices, n_core))

        completed = 0
        batch_size = 50
        while completed < n_permutations:
            batch = min(batch_size, n_permutations - completed)
            perm_counts = np.zeros((batch, len(raw)), dtype=np.uint8)
            batch_rows = np.arange(batch)[:, None]
            for member_indices, n_core in permutation_groups:
                n_members = len(member_indices)
                if n_core == 0:
                    continue
                if n_core == n_members:
                    selected = np.broadcast_to(member_indices, (batch, n_members))
                else:
                    scores = rng.random((batch, n_members))
                    selected_positions = np.argpartition(
                        scores, n_core - 1, axis=1
                    )[:, :n_core]
                    selected = member_indices[selected_positions]
                row_indices = np.broadcast_to(batch_rows, selected.shape)
                mask = selected >= 0
                np.add.at(perm_counts, (row_indices[mask], selected[mask]), 1)
            exceed += np.sum(perm_counts >= observed[None, :], axis=0)
            ge2_distribution.extend(np.sum(perm_counts >= 2, axis=1).astype(int).tolist())
            ge3_distribution.extend(np.sum(perm_counts >= 3, axis=1).astype(int).tolist())
            completed += batch

        empirical = (exceed + 1) / (n_permutations + 1)
        for row, value in zip(raw, empirical.tolist()):
            row["matched_permutation_p_count_tail"] = value
        recurrent_keys = {
            (str(row["gene_id"]), str(row["tissue"])) for row in recurrent
        }
        recurrent = [
            dict(row)
            for row in raw
            if (str(row["gene_id"]), str(row["tissue"])) in recurrent_keys
        ]
        permutation_summary.update(
            {
                "mean_null_pairs_with_at_least_two_core_animals": float(
                    sum(ge2_distribution) / len(ge2_distribution)
                ),
                "min_null_pairs_with_at_least_two_core_animals": min(ge2_distribution),
                "max_null_pairs_with_at_least_two_core_animals": max(ge2_distribution),
                "mean_null_pairs_with_at_least_three_core_animals": float(
                    sum(ge3_distribution) / len(ge3_distribution)
                ),
                "min_null_pairs_with_at_least_three_core_animals": min(ge3_distribution),
                "max_null_pairs_with_at_least_three_core_animals": max(ge3_distribution),
            }
        )

    gene_tissues: MutableMapping[str, Set[str]] = defaultdict(set)
    for row in recurrent:
        gene_tissues[str(row["gene_id"])].add(str(row["tissue"]))
    summary = {
        "n_opportunity_strata": len(strata_rows),
        "n_recurrence_hypotheses": len(raw),
        "n_recurrent_gene_tissue_pairs": len(recurrent),
        "n_recurrent_genes": len(gene_tissues),
        "n_recurrent_in_one_tissue": sum(len(value) == 1 for value in gene_tissues.values()),
        "n_recurrent_in_multiple_tissues": sum(len(value) > 1 for value in gene_tissues.values()),
        "stratum_size_min": min(int(row["n_eligible_gene_measurements"]) for row in strata_rows),
        "stratum_size_median": sorted(
            int(row["n_eligible_gene_measurements"]) for row in strata_rows
        )[len(strata_rows) // 2],
        "stratum_size_max": max(int(row["n_eligible_gene_measurements"]) for row in strata_rows),
        "permutation": permutation_summary,
    }
    return raw, recurrent, strata_rows, summary


def measurement_rows(
    measurements: Mapping[Tuple[str, str, str], GeneMeasurement]
) -> Iterable[Dict[str, object]]:
    for key in sorted(measurements):
        measurement = measurements[key]
        row = asdict(measurement)
        row["is_core_ase"] = int(measurement.core)
        yield row


def serialize_rows(rows: Iterable[Mapping[str, object]]) -> Iterable[Dict[str, object]]:
    for source in rows:
        row = dict(source)
        for key, value in list(row.items()):
            if isinstance(value, float):
                row[key] = fmt(value)
        yield row


def read_recurrent_keys(path: Path) -> Set[Tuple[str, str]]:
    result: Set[Tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result.add((row["gene_id"], normalize_tissue(row["tissue"])))
    return result


def output_comparison(
    out_path: Path,
    canonical_keys: Set[Tuple[str, str]],
    refined_rows: Sequence[Mapping[str, object]],
) -> Dict[str, int]:
    refined_keys = {(str(row["gene_id"]), str(row["tissue"])) for row in refined_rows}
    rows = []
    for gene_id, tissue in sorted(canonical_keys | refined_keys):
        in_canonical = (gene_id, tissue) in canonical_keys
        in_refined = (gene_id, tissue) in refined_keys
        status = (
            "RETAINED" if in_canonical and in_refined else
            "LOST_AFTER_ELIGIBILITY_REFINEMENT" if in_canonical else
            "NEW_AFTER_ELIGIBILITY_REFINEMENT"
        )
        rows.append(
            {
                "gene_id": gene_id,
                "tissue": tissue,
                "canonical_recurrent": int(in_canonical),
                "refined_recurrent": int(in_refined),
                "comparison_status": status,
            }
        )
    write_tsv(
        out_path,
        ["gene_id", "tissue", "canonical_recurrent", "refined_recurrent", "comparison_status"],
        rows,
    )
    counts = Counter(row["comparison_status"] for row in rows)
    return dict(sorted(counts.items()))


MEASUREMENT_FIELDS = [
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "n_tested_blocks",
    "n_eligible_blocks",
    "n_core_blocks",
    "max_eligible_total_fragments",
    "max_eligible_contributing_snps",
    "is_core_ase",
]

RECURRENCE_FIELDS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "tissue",
    "n_tissue_animals_total",
    "n_samples_eligible",
    "samples_eligible",
    "n_samples_core",
    "samples_core",
    "eligible_fraction_of_tissue_animals",
    "expected_core_inclusive",
    "expected_core_leave_one_gene_out",
    "recurrence_p_inclusive",
    "recurrence_p_leave_one_gene_out",
    "recurrence_q_bh_inclusive",
    "recurrence_q_bh_leave_one_gene_out",
    "is_statistically_recurrent_core",
    "replication_support",
]

BURDEN_FIELDS = [
    "sample",
    "tissue",
    "n_tissue_animals_total",
    "n_eligible_genes",
    "n_core_genes",
    "core_probability_inclusive",
]

MATCHED_FIELDS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "tissue",
    "n_tissue_animals_total",
    "n_samples_eligible",
    "n_samples_core",
    "samples_core",
    "expected_core_opportunity_matched",
    "recurrence_p_opportunity_matched",
    "recurrence_q_bh_opportunity_matched",
    "matched_permutation_p_count_tail",
    "is_statistically_recurrent_core",
]

STRATA_FIELDS = [
    "sample",
    "tissue",
    "eligible_block_bin",
    "max_fragment_bin",
    "max_contributing_snp_bin",
    "n_eligible_gene_measurements",
    "n_core_gene_measurements",
    "core_probability",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--canonical-core", type=Path, required=True)
    parser.add_argument("--canonical-recurrent", type=Path, required=True)
    parser.add_argument("--thresholds", default="15,20,30")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-major-fraction", type=float, default=0.65)
    parser.add_argument("--max-fragment-conflict", type=float, default=0.10)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=260816)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    thresholds = tuple(sorted({int(value) for value in args.thresholds.split(",")}))
    if not thresholds or min(thresholds) < 15:
        raise AnalysisError("this implementation requires fragment thresholds >=15")
    if 15 not in thresholds:
        raise AnalysisError("threshold 15 is required for canonical validation")
    files = discover_source_files(args.input_dir)
    if len(files) != 10:
        raise AnalysisError(f"expected 10 source files, found {len(files)}")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.force:
        raise AnalysisError(f"output directory is not empty: {args.out_dir}; use --force")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    measurements_by_threshold: Dict[
        int, Dict[Tuple[str, str, str], GeneMeasurement]
    ] = {threshold: {} for threshold in thresholds}
    tissue_animals: MutableMapping[str, Set[str]] = defaultdict(set)
    input_qc: Dict[str, object] = {}
    q_qc = Counter()
    q_max_difference = 0.0
    source_files = []

    for path in files:
        sample = infer_sample(path)
        print(f"[{sample}] first pass: threshold-specific BH", flush=True)
        q_maps, first_qc = build_q_maps(path, sample, thresholds)
        for unit in first_qc["tissue_units"]:
            unit_sample, tissue = str(unit).split("|", 1)
            tissue_animals[tissue].add(unit_sample)
        print(f"[{sample}] second pass: eligibility/core roll-up", flush=True)
        sample_measurements, second_qc = load_sample_measurements(
            path,
            sample,
            thresholds,
            q_maps,
            args.max_fragment_conflict,
            args.fdr,
            args.min_major_fraction,
        )
        for threshold in thresholds:
            eligible_subset = {
                key: value
                for key, value in sample_measurements[threshold].items()
                if value.eligible
            }
            overlap = set(measurements_by_threshold[threshold]).intersection(eligible_subset)
            if overlap:
                raise AnalysisError(f"duplicate cohort measurement keys: {sorted(overlap)[:3]}")
            measurements_by_threshold[threshold].update(eligible_subset)
            print(
                f"[{sample}] threshold={threshold} eligible={len(eligible_subset)} "
                f"core={sum(value.core for value in eligible_subset.values())}",
                flush=True,
            )
        q_qc["comparisons"] += int(second_qc["q_comparisons_threshold15"])
        q_qc["mismatches"] += int(second_qc["q_mismatches_threshold15"])
        q_max_difference = max(
            q_max_difference, float(second_qc["q_max_abs_difference_threshold15"])
        )
        input_qc[sample] = {**first_qc, **second_qc}
        source_files.append(
            {
                "sample": sample,
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

    canonical_core = read_canonical_core(args.canonical_core)
    reconstructed_core = {
        key for key, measurement in measurements_by_threshold[15].items() if measurement.core
    }
    missing_core = canonical_core - reconstructed_core
    extra_core = reconstructed_core - canonical_core
    if missing_core or extra_core:
        raise AnalysisError(
            f"threshold15 core mismatch: missing={len(missing_core)} extra={len(extra_core)}"
        )
    if q_qc["mismatches"]:
        raise AnalysisError(f"threshold15 BH mismatches: {q_qc['mismatches']}")

    threshold_summaries: Dict[str, object] = {}
    threshold_recurrent: Dict[int, List[Dict[str, object]]] = {}
    for threshold in thresholds:
        print(f"[threshold={threshold}] recurrence", flush=True)
        measurements = measurements_by_threshold[threshold]
        recurrence, recurrent, burden, summary = build_recurrence(
            measurements, tissue_animals, args.fdr
        )
        threshold_recurrent[threshold] = recurrent
        threshold_summaries[str(threshold)] = summary
        prefix = f"threshold_{threshold}"
        write_tsv(
            args.out_dir / f"{prefix}_eligible_gene_measurements.tsv.gz",
            MEASUREMENT_FIELDS,
            measurement_rows(measurements),
        )
        write_tsv(
            args.out_dir / f"{prefix}_sample_tissue_burden.tsv",
            BURDEN_FIELDS,
            serialize_rows(burden),
        )
        write_tsv(
            args.out_dir / f"{prefix}_recurrence_all.tsv.gz",
            RECURRENCE_FIELDS,
            serialize_rows(recurrence),
        )
        write_tsv(
            args.out_dir / f"{prefix}_statistically_recurrent.tsv",
            RECURRENCE_FIELDS,
            serialize_rows(recurrent),
        )

    print("[threshold=15] opportunity-matched analysis", flush=True)
    matched, matched_recurrent, strata, matched_summary = build_matched_analysis(
        measurements_by_threshold[15],
        tissue_animals,
        args.fdr,
        args.n_permutations,
        args.seed,
    )
    write_tsv(
        args.out_dir / "opportunity_strata.tsv.gz",
        STRATA_FIELDS,
        serialize_rows(strata),
    )
    write_tsv(
        args.out_dir / "opportunity_matched_recurrence_all.tsv.gz",
        MATCHED_FIELDS,
        serialize_rows(matched),
    )
    write_tsv(
        args.out_dir / "opportunity_matched_statistically_recurrent.tsv",
        MATCHED_FIELDS,
        serialize_rows(matched_recurrent),
    )

    canonical_recurrent = read_recurrent_keys(args.canonical_recurrent)
    comparison_counts = output_comparison(
        args.out_dir / "comparison_with_canonical_recurrence.tsv",
        canonical_recurrent,
        threshold_recurrent[15],
    )
    threshold_summary_rows = []
    for threshold in thresholds:
        summary = threshold_summaries[str(threshold)]
        threshold_summary_rows.append({"fragment_threshold": threshold, **summary})
    threshold_fields = ["fragment_threshold"] + list(threshold_summary_rows[0].keys())[1:]
    write_tsv(
        args.out_dir / "threshold_summary.tsv",
        threshold_fields,
        threshold_summary_rows,
    )

    qc = {
        "status": "PASS",
        "version": VERSION,
        "definition": {
            "binomial_testable": "analysis_set MAIN and total fragments >= threshold",
            "high_confidence_eligible": (
                "binomial-testable, technical PASS, >=2 contributing SNPs, local phase PASS"
            ),
            "core": (
                "eligible, threshold-specific within-sample-tissue BH q<=0.05, "
                "major fraction>=0.65, multi-SNP supported/LOO-robust"
            ),
            "primary_recurrence_null": "leave-one-gene-out sample-tissue core/eligible burden",
            "opportunity_strata": (
                "eligible-block count (1,2,3+), max fragments "
                "(15-29,30-59,60-119,120-239,240+), max contributing SNPs "
                "(2,3,4-5,6+)"
            ),
        },
        "parameters": {
            "thresholds": thresholds,
            "fdr": args.fdr,
            "min_major_fraction": args.min_major_fraction,
            "max_fragment_conflict": args.max_fragment_conflict,
            "n_permutations": args.n_permutations,
            "seed": args.seed,
        },
        "source_files": source_files,
        "input_qc": input_qc,
        "threshold15_q_recalculation": {
            "n_comparisons": q_qc["comparisons"],
            "n_mismatches_tolerance_1e-9": q_qc["mismatches"],
            "max_abs_difference": q_max_difference,
        },
        "canonical_core_crosscheck": {
            "canonical": len(canonical_core),
            "reconstructed": len(reconstructed_core),
            "missing": len(missing_core),
            "extra": len(extra_core),
        },
        "tissue_animal_counts": {
            tissue: len(samples) for tissue, samples in sorted(tissue_animals.items())
        },
        "threshold_summaries": threshold_summaries,
        "opportunity_matched_summary": matched_summary,
        "canonical_recurrence_comparison": comparison_counts,
    }
    with (args.out_dir / "analysis_qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": "PASS", "thresholds": threshold_summaries, "matched": matched_summary}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
