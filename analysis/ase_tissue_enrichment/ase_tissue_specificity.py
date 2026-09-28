#!/usr/bin/env python3
"""Formal tissue-enrichment tests for high-confidence ASE recurrence.

For each gene, the binary core-ASE labels are compared across tissues while:

* using only outcome-independent high-confidence-eligible measurements;
* conditioning on each animal's number of core-positive eligible tissues;
* weighting tissue assignment by leave-one-gene-out sample-tissue core burden;
* testing every evaluable gene-tissue contrast before looking at recurrence;
* applying global BH correction to all evaluable gene-tissue contrasts.

The formal result is called tissue-enriched ASE, rather than absolute tissue
specificity, because lack of ASE in other tissues may still reflect incomplete
measurement opportunity.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "0.1.0"
DEFAULT_THRESHOLDS = (15, 20, 30)
PROFILE_LR_95 = 3.841458820694124


class SpecificityError(RuntimeError):
    pass


@dataclass
class GeneData:
    gene_id: str
    gene_name: str
    gene_biotype: str
    by_sample: MutableMapping[str, Dict[str, int]]


def normalize_tissue(tissue: str) -> str:
    return {
        "blood": "Blood",
        "L-blood": "Blood",
        "Blood": "Blood",
        "Tenderlo": "Tenderloin",
        "tenderloin": "Tenderloin",
        "Tenderloin": "Tenderloin",
    }.get(tissue, tissue)


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return format(value, ".12g")
    return str(value)


def number(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    if value == "Inf":
        return math.inf
    if value == "-Inf":
        return -math.inf
    return float(value)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        raw = path.open("wb")
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in rows:
            writer.writerow({field: fmt(source.get(field)) for field in fields})
    finally:
        handle.close()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_values(values: Sequence[float]) -> List[float]:
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
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
            raise SpecificityError(f"invalid Bernoulli probability: {probability}")
        updated = [0.0] * (len(pmf) + 1)
        for count, mass in enumerate(pmf):
            updated[count] += mass * (1.0 - probability)
            updated[count + 1] += mass * probability
        pmf = updated
    return min(1.0, max(0.0, math.fsum(pmf[observed:])))


def elementary_symmetric(weights: Sequence[float], order: int) -> float:
    if order < 0 or order > len(weights):
        return 0.0
    coefficients = [0.0] * (order + 1)
    coefficients[0] = 1.0
    for weight in weights:
        for index in range(order, 0, -1):
            coefficients[index] += weight * coefficients[index - 1]
    return coefficients[order]


def conditional_marginals(
    tissues: Sequence[str],
    odds: Sequence[float],
    n_positive: int,
) -> Dict[str, float]:
    """P(each tissue is positive | total positives) under weighted Bernoulli null."""
    if len(tissues) != len(odds) or len(set(tissues)) != len(tissues):
        raise SpecificityError("invalid tissue/odds vectors")
    if not 0 < n_positive < len(tissues):
        raise SpecificityError("conditional marginals require 0 < positives < tissues")
    if any(not math.isfinite(weight) or weight <= 0.0 for weight in odds):
        raise SpecificityError(f"conditional odds must be positive and finite: {odds}")
    denominator = elementary_symmetric(odds, n_positive)
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise SpecificityError("invalid conditional denominator")
    result = {}
    for index, tissue in enumerate(tissues):
        others = list(odds[:index]) + list(odds[index + 1 :])
        numerator = odds[index] * elementary_symmetric(others, n_positive - 1)
        result[tissue] = min(1.0, max(0.0, numerator / denominator))
    if not math.isclose(math.fsum(result.values()), n_positive, rel_tol=1e-9, abs_tol=1e-9):
        raise SpecificityError("conditional marginal probabilities do not sum to m")
    return result


def logistic_from_offset(logit0: float, beta: float) -> float:
    eta = max(-745.0, min(709.0, logit0 + beta))
    if eta >= 0.0:
        return 1.0 / (1.0 + math.exp(-eta))
    exp_eta = math.exp(eta)
    return exp_eta / (1.0 + exp_eta)


def conditional_log_likelihood(
    outcomes: Sequence[int],
    baseline_probabilities: Sequence[float],
    beta: float,
) -> float:
    likelihood = 0.0
    for outcome, probability0 in zip(outcomes, baseline_probabilities):
        logit0 = math.log(probability0 / (1.0 - probability0))
        probability = logistic_from_offset(logit0, beta)
        if outcome:
            likelihood += math.log(max(probability, 1e-323))
        else:
            likelihood += math.log(max(1.0 - probability, 1e-323))
    return likelihood


def solve_beta(
    baseline_probabilities: Sequence[float],
    target_sum: int,
    lower: float = -60.0,
    upper: float = 60.0,
) -> float:
    for _ in range(200):
        midpoint = (lower + upper) / 2.0
        fitted = math.fsum(
            logistic_from_offset(math.log(p / (1.0 - p)), midpoint)
            for p in baseline_probabilities
        )
        if fitted < target_sum:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def profile_log_odds_ratio(
    outcomes: Sequence[int],
    baseline_probabilities: Sequence[float],
) -> Tuple[float, float, float]:
    """Conditional common odds ratio and two-sided 95% profile-likelihood CI."""
    if not outcomes or len(outcomes) != len(baseline_probabilities):
        return math.nan, math.nan, math.nan
    if any(outcome not in (0, 1) for outcome in outcomes):
        raise SpecificityError("conditional outcomes must be binary")
    if any(not 0.0 < probability < 1.0 for probability in baseline_probabilities):
        raise SpecificityError("profile likelihood requires probabilities in (0,1)")
    observed = sum(outcomes)
    if observed == 0:
        beta_hat = -math.inf
        max_likelihood = 0.0
    elif observed == len(outcomes):
        beta_hat = math.inf
        max_likelihood = 0.0
    else:
        beta_hat = solve_beta(baseline_probabilities, observed)
        max_likelihood = conditional_log_likelihood(outcomes, baseline_probabilities, beta_hat)

    def objective(beta: float) -> float:
        return 2.0 * (
            max_likelihood
            - conditional_log_likelihood(outcomes, baseline_probabilities, beta)
        ) - PROFILE_LR_95

    def root(left: float, right: float) -> float:
        left_value = objective(left)
        right_value = objective(right)
        if left_value * right_value > 0.0:
            raise SpecificityError("profile-likelihood CI root is not bracketed")
        for _ in range(200):
            midpoint = (left + right) / 2.0
            mid_value = objective(midpoint)
            if left_value * mid_value <= 0.0:
                right = midpoint
                right_value = mid_value
            else:
                left = midpoint
                left_value = mid_value
        return (left + right) / 2.0

    if math.isinf(beta_hat) and beta_hat < 0:
        lower_beta = -math.inf
        upper_beta = root(-60.0, 60.0)
    elif math.isinf(beta_hat):
        lower_beta = root(-60.0, 60.0)
        upper_beta = math.inf
    else:
        lower_beta = root(-60.0, beta_hat)
        upper_beta = root(beta_hat, 60.0)

    def exponentiate(beta: float) -> float:
        if beta == math.inf:
            return math.inf
        if beta == -math.inf:
            return 0.0
        return math.exp(beta)

    return exponentiate(beta_hat), exponentiate(lower_beta), exponentiate(upper_beta)


def load_burden(path: Path) -> Dict[Tuple[str, str], Tuple[int, int]]:
    result = {}
    for row in read_tsv(path):
        key = (row["sample"], normalize_tissue(row["tissue"]))
        if key in result:
            raise SpecificityError(f"duplicate burden key: {key}")
        n_eligible = int(row["n_eligible_genes"])
        n_core = int(row["n_core_genes"])
        if not 0 <= n_core <= n_eligible or n_eligible <= 1:
            raise SpecificityError(f"invalid burden at {key}: {n_core}/{n_eligible}")
        result[key] = (n_eligible, n_core)
    return result


def load_measurements(path: Path) -> Tuple[Dict[str, GeneData], Dict[str, int]]:
    genes: Dict[str, GeneData] = {}
    n_rows = 0
    n_core = 0
    for row in read_tsv(path):
        n_rows += 1
        if int(row["n_eligible_blocks"]) < 1:
            raise SpecificityError(f"non-eligible row in eligible input: {row}")
        gene_id = row["gene_id"]
        gene = genes.get(gene_id)
        if gene is None:
            gene = GeneData(
                gene_id=gene_id,
                gene_name=row["gene_name"],
                gene_biotype=row["gene_biotype"],
                by_sample=defaultdict(dict),
            )
            genes[gene_id] = gene
        elif gene.gene_name != row["gene_name"] or gene.gene_biotype != row["gene_biotype"]:
            raise SpecificityError(f"inconsistent annotation for {gene_id}")
        sample = row["sample"]
        tissue = normalize_tissue(row["tissue"])
        if tissue in gene.by_sample[sample]:
            raise SpecificityError(f"duplicate measurement: {sample}|{tissue}|{gene_id}")
        core = int(row["is_core_ase"])
        if core not in (0, 1):
            raise SpecificityError(f"invalid core flag: {row['is_core_ase']}")
        gene.by_sample[sample][tissue] = core
        n_core += core
    return genes, {"n_eligible_gene_measurements": n_rows, "n_core_gene_measurements": n_core}


def load_recurrence(
    path: Path,
) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], Dict[str, Set[str]]]:
    result = {}
    tissues_by_gene: MutableMapping[str, Set[str]] = defaultdict(set)
    for row in read_tsv(path):
        tissue = normalize_tissue(row["tissue"])
        key = (row["gene_id"], tissue)
        if key in result:
            raise SpecificityError(f"duplicate recurrent pair: {key}")
        copy = dict(row)
        copy["tissue"] = tissue
        result[key] = copy
        tissues_by_gene[row["gene_id"]].add(tissue)
    return result, dict(tissues_by_gene)


def leave_one_gene_out_probability(
    burden: Mapping[Tuple[str, str], Tuple[int, int]],
    sample: str,
    tissue: str,
    outcome: int,
) -> float:
    if (sample, tissue) not in burden:
        raise SpecificityError(f"missing burden: {sample}|{tissue}")
    n_eligible, n_core = burden[(sample, tissue)]
    probability = (n_core - outcome) / (n_eligible - 1)
    if not 0.0 < probability < 1.0:
        raise SpecificityError(
            f"boundary leave-one-gene-out probability at {sample}|{tissue}: {probability}"
        )
    return probability


ALL_CONTRAST_FIELDS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "target_tissue",
    "n_gene_eligible_tissues",
    "gene_eligible_tissues",
    "n_target_eligible_animals_total",
    "n_target_core_animals_total",
    "target_core_fraction_total",
    "n_other_eligible_measurements_total",
    "n_other_core_measurements_total",
    "other_core_fraction_total",
    "n_informative_animals",
    "informative_animals",
    "n_target_core_informative",
    "expected_target_core_conditional_null",
    "conditional_core_excess",
    "conditional_observed_to_expected_ratio",
    "conditional_common_odds_ratio",
    "conditional_common_odds_ratio_CI95_lower",
    "conditional_common_odds_ratio_CI95_upper",
    "tissue_enrichment_exact_p_one_sided",
    "tissue_enrichment_q_global_BH",
    "is_contrast_evaluable",
    "gene_omnibus_p_bonferroni",
    "gene_omnibus_q_global_BH",
    "gene_omnibus_fdr_pass",
    "is_primary_recurrent",
    "n_primary_recurrent_tissues_for_gene",
    "primary_recurrent_tissues_for_gene",
    "primary_recurrence_pattern",
    "is_formally_tissue_enriched_ase",
    "is_formally_tissue_restricted_recurrence",
    "formal_status",
]

OMNIBUS_FIELDS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "n_gene_eligible_tissues",
    "gene_eligible_tissues",
    "n_evaluable_tissue_contrasts",
    "minimum_tissue_contrast_p",
    "minimum_p_target_tissue",
    "gene_omnibus_p_bonferroni",
    "gene_omnibus_q_global_BH",
    "gene_omnibus_fdr_pass",
    "n_primary_recurrent_tissues_for_gene",
    "primary_recurrent_tissues_for_gene",
    "n_formally_tissue_enriched_tissues",
    "formally_tissue_enriched_tissues",
]

RECURRENT_STATUS_FIELDS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "tissue",
    "n_samples_eligible_recurrence",
    "n_samples_core_recurrence",
    "recurrence_q_bh_leave_one_gene_out",
    "n_primary_recurrent_tissues_for_gene",
    "primary_recurrence_pattern",
    "specificity_test_status",
    "n_informative_animals",
    "n_target_core_informative",
    "expected_target_core_conditional_null",
    "conditional_common_odds_ratio",
    "conditional_common_odds_ratio_CI95_lower",
    "conditional_common_odds_ratio_CI95_upper",
    "tissue_enrichment_exact_p_one_sided",
    "tissue_enrichment_q_global_BH",
    "gene_omnibus_q_global_BH",
    "is_formally_tissue_enriched_ase",
    "is_formally_tissue_restricted_recurrence",
    "formal_status",
]


def analyze_threshold(
    threshold: int,
    measurement_path: Path,
    burden_path: Path,
    recurrent_path: Path,
    out_dir: Path,
    min_informative_animals: int,
    min_informative_target_core: int,
    fdr: float,
) -> Tuple[Dict[str, object], Set[Tuple[str, str]], Dict[str, object]]:
    burden = load_burden(burden_path)
    genes, input_counts = load_measurements(measurement_path)
    recurrent, recurrent_tissues = load_recurrence(recurrent_path)
    contrast_rows: List[Dict[str, object]] = []
    evaluable_indices: List[int] = []
    gene_contrast_indices: MutableMapping[str, List[int]] = defaultdict(list)

    for gene_id in sorted(genes):
        gene = genes[gene_id]
        tissue_eligible = Counter()
        tissue_core = Counter()
        total_eligible = 0
        total_core = 0
        all_tissues: Set[str] = set()
        informative: MutableMapping[str, List[Tuple[str, int, float]]] = defaultdict(list)
        for sample, cells in sorted(gene.by_sample.items()):
            for tissue, outcome in cells.items():
                tissue_eligible[tissue] += 1
                tissue_core[tissue] += outcome
                total_eligible += 1
                total_core += outcome
                all_tissues.add(tissue)
            n_tissues = len(cells)
            n_positive = sum(cells.values())
            if n_tissues < 2 or n_positive == 0 or n_positive == n_tissues:
                continue
            tissues = sorted(cells)
            probabilities = [
                leave_one_gene_out_probability(burden, sample, tissue, cells[tissue])
                for tissue in tissues
            ]
            odds = [probability / (1.0 - probability) for probability in probabilities]
            marginals = conditional_marginals(tissues, odds, n_positive)
            for tissue in tissues:
                informative[tissue].append((sample, cells[tissue], marginals[tissue]))

        if len(all_tissues) < 2:
            continue
        eligible_tissue_text = ",".join(sorted(all_tissues))
        gene_recurrent_tissues = sorted(recurrent_tissues.get(gene_id, set()))
        for target_tissue in sorted(all_tissues):
            entries = informative.get(target_tissue, [])
            outcomes = [entry[1] for entry in entries]
            probabilities = [entry[2] for entry in entries]
            observed = sum(outcomes)
            expected = math.fsum(probabilities)
            is_evaluable = len(entries) >= min_informative_animals
            exact_p = poisson_binomial_tail(probabilities, observed) if is_evaluable else math.nan
            if entries:
                odds_ratio, ci_lower, ci_upper = profile_log_odds_ratio(outcomes, probabilities)
            else:
                odds_ratio, ci_lower, ci_upper = math.nan, math.nan, math.nan
            other_eligible = total_eligible - tissue_eligible[target_tissue]
            other_core = total_core - tissue_core[target_tissue]
            row: Dict[str, object] = {
                "gene_id": gene_id,
                "gene_name": gene.gene_name,
                "gene_biotype": gene.gene_biotype,
                "target_tissue": target_tissue,
                "n_gene_eligible_tissues": len(all_tissues),
                "gene_eligible_tissues": eligible_tissue_text,
                "n_target_eligible_animals_total": tissue_eligible[target_tissue],
                "n_target_core_animals_total": tissue_core[target_tissue],
                "target_core_fraction_total": tissue_core[target_tissue]
                / tissue_eligible[target_tissue],
                "n_other_eligible_measurements_total": other_eligible,
                "n_other_core_measurements_total": other_core,
                "other_core_fraction_total": other_core / other_eligible if other_eligible else math.nan,
                "n_informative_animals": len(entries),
                "informative_animals": ",".join(entry[0] for entry in entries),
                "n_target_core_informative": observed,
                "expected_target_core_conditional_null": expected,
                "conditional_core_excess": observed - expected,
                "conditional_observed_to_expected_ratio": observed / expected
                if expected > 0.0
                else math.inf,
                "conditional_common_odds_ratio": odds_ratio,
                "conditional_common_odds_ratio_CI95_lower": ci_lower,
                "conditional_common_odds_ratio_CI95_upper": ci_upper,
                "tissue_enrichment_exact_p_one_sided": exact_p,
                "tissue_enrichment_q_global_BH": math.nan,
                "is_contrast_evaluable": int(is_evaluable),
                "gene_omnibus_p_bonferroni": math.nan,
                "gene_omnibus_q_global_BH": math.nan,
                "gene_omnibus_fdr_pass": 0,
                "is_primary_recurrent": int((gene_id, target_tissue) in recurrent),
                "n_primary_recurrent_tissues_for_gene": len(gene_recurrent_tissues),
                "primary_recurrent_tissues_for_gene": ",".join(gene_recurrent_tissues),
                "primary_recurrence_pattern": (
                    "ONE_RECURRENT_TISSUE"
                    if len(gene_recurrent_tissues) == 1
                    else "MULTIPLE_RECURRENT_TISSUES"
                    if len(gene_recurrent_tissues) > 1
                    else "NOT_PRIMARY_RECURRENT"
                ),
                "is_formally_tissue_enriched_ase": 0,
                "is_formally_tissue_restricted_recurrence": 0,
                "formal_status": "PENDING",
            }
            contrast_rows.append(row)
            index = len(contrast_rows) - 1
            if is_evaluable:
                evaluable_indices.append(index)
                gene_contrast_indices[gene_id].append(index)

    contrast_q = bh_values(
        [float(contrast_rows[index]["tissue_enrichment_exact_p_one_sided"]) for index in evaluable_indices]
    )
    for index, q_value in zip(evaluable_indices, contrast_q):
        contrast_rows[index]["tissue_enrichment_q_global_BH"] = q_value

    omnibus_rows: List[Dict[str, object]] = []
    for gene_id in sorted(gene_contrast_indices):
        indices = gene_contrast_indices[gene_id]
        p_values = [float(contrast_rows[index]["tissue_enrichment_exact_p_one_sided"]) for index in indices]
        minimum_index = min(indices, key=lambda index: (
            float(contrast_rows[index]["tissue_enrichment_exact_p_one_sided"]),
            str(contrast_rows[index]["target_tissue"]),
        ))
        omnibus_p = min(1.0, min(p_values) * len(p_values))
        source = contrast_rows[indices[0]]
        omnibus_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": source["gene_name"],
                "gene_biotype": source["gene_biotype"],
                "n_gene_eligible_tissues": source["n_gene_eligible_tissues"],
                "gene_eligible_tissues": source["gene_eligible_tissues"],
                "n_evaluable_tissue_contrasts": len(indices),
                "minimum_tissue_contrast_p": min(p_values),
                "minimum_p_target_tissue": contrast_rows[minimum_index]["target_tissue"],
                "gene_omnibus_p_bonferroni": omnibus_p,
                "gene_omnibus_q_global_BH": math.nan,
                "gene_omnibus_fdr_pass": 0,
                "n_primary_recurrent_tissues_for_gene": source[
                    "n_primary_recurrent_tissues_for_gene"
                ],
                "primary_recurrent_tissues_for_gene": source[
                    "primary_recurrent_tissues_for_gene"
                ],
                "n_formally_tissue_enriched_tissues": 0,
                "formally_tissue_enriched_tissues": "",
            }
        )

    omnibus_q = bh_values([float(row["gene_omnibus_p_bonferroni"]) for row in omnibus_rows])
    omnibus_by_gene = {}
    for row, q_value in zip(omnibus_rows, omnibus_q):
        row["gene_omnibus_q_global_BH"] = q_value
        row["gene_omnibus_fdr_pass"] = int(q_value <= fdr)
        omnibus_by_gene[str(row["gene_id"])] = row

    formal_rows: List[Dict[str, object]] = []
    formal_by_gene: MutableMapping[str, List[str]] = defaultdict(list)
    contrast_lookup: Dict[Tuple[str, str], Dict[str, object]] = {}
    for row in contrast_rows:
        gene_id = str(row["gene_id"])
        omnibus = omnibus_by_gene.get(gene_id)
        if omnibus is not None:
            row["gene_omnibus_p_bonferroni"] = omnibus["gene_omnibus_p_bonferroni"]
            row["gene_omnibus_q_global_BH"] = omnibus["gene_omnibus_q_global_BH"]
            row["gene_omnibus_fdr_pass"] = omnibus["gene_omnibus_fdr_pass"]
        failures = []
        if not int(row["is_contrast_evaluable"]):
            failures.append("MIN_INFORMATIVE_ANIMALS")
        else:
            if float(row["tissue_enrichment_q_global_BH"]) > fdr:
                failures.append("TARGET_CONTRAST_FDR")
            if omnibus is None or int(omnibus["gene_omnibus_fdr_pass"]) != 1:
                failures.append("GENE_OMNIBUS_FDR")
            if int(row["n_target_core_informative"]) < min_informative_target_core:
                failures.append("MIN_INFORMATIVE_TARGET_CORE")
            if float(row["conditional_core_excess"]) <= 0.0:
                failures.append("TARGET_NOT_ENRICHED")
        if not int(row["is_primary_recurrent"]):
            failures.append("NOT_PRIMARY_RECURRENT")
        formal = int(not failures)
        row["is_formally_tissue_enriched_ase"] = formal
        row["is_formally_tissue_restricted_recurrence"] = int(
            formal and int(row["n_primary_recurrent_tissues_for_gene"]) == 1
        )
        row["formal_status"] = "PASS" if formal else "FAIL:" + ",".join(failures)
        contrast_lookup[(gene_id, str(row["target_tissue"]))] = row
        if formal:
            formal_rows.append(dict(row))
            formal_by_gene[gene_id].append(str(row["target_tissue"]))

    for row in omnibus_rows:
        tissues = sorted(formal_by_gene.get(str(row["gene_id"]), []))
        row["n_formally_tissue_enriched_tissues"] = len(tissues)
        row["formally_tissue_enriched_tissues"] = ",".join(tissues)

    recurrent_status = []
    for key, recurrence_row in sorted(recurrent.items()):
        gene_id, tissue = key
        contrast = contrast_lookup.get(key)
        n_recurrent_tissues = len(recurrent_tissues[gene_id])
        status: Dict[str, object] = {
            "gene_id": gene_id,
            "gene_name": recurrence_row["gene_name"],
            "gene_biotype": recurrence_row["gene_biotype"],
            "tissue": tissue,
            "n_samples_eligible_recurrence": recurrence_row["n_samples_eligible"],
            "n_samples_core_recurrence": recurrence_row["n_samples_core"],
            "recurrence_q_bh_leave_one_gene_out": recurrence_row[
                "recurrence_q_bh_leave_one_gene_out"
            ],
            "n_primary_recurrent_tissues_for_gene": n_recurrent_tissues,
            "primary_recurrence_pattern": (
                "ONE_RECURRENT_TISSUE" if n_recurrent_tissues == 1 else "MULTIPLE_RECURRENT_TISSUES"
            ),
        }
        if contrast is None:
            status.update(
                {
                    "specificity_test_status": "NOT_EVALUABLE_NO_CROSS_TISSUE_OPPORTUNITY",
                    "n_informative_animals": 0,
                    "n_target_core_informative": 0,
                    "expected_target_core_conditional_null": math.nan,
                    "conditional_common_odds_ratio": math.nan,
                    "conditional_common_odds_ratio_CI95_lower": math.nan,
                    "conditional_common_odds_ratio_CI95_upper": math.nan,
                    "tissue_enrichment_exact_p_one_sided": math.nan,
                    "tissue_enrichment_q_global_BH": math.nan,
                    "gene_omnibus_q_global_BH": math.nan,
                    "is_formally_tissue_enriched_ase": 0,
                    "is_formally_tissue_restricted_recurrence": 0,
                    "formal_status": "FAIL:NO_CROSS_TISSUE_OPPORTUNITY",
                }
            )
        else:
            status.update(
                {
                    "specificity_test_status": (
                        "TESTED" if int(contrast["is_contrast_evaluable"]) else "NOT_EVALUABLE_MIN_INFORMATIVE_ANIMALS"
                    ),
                    **{
                        field: contrast[field]
                        for field in RECURRENT_STATUS_FIELDS
                        if field in contrast
                        and field
                        not in {
                            "gene_id",
                            "gene_name",
                            "gene_biotype",
                            "tissue",
                            "n_primary_recurrent_tissues_for_gene",
                            "primary_recurrence_pattern",
                        }
                    },
                }
            )
        recurrent_status.append(status)

    prefix = f"threshold_{threshold}"
    write_tsv(out_dir / f"{prefix}_tissue_contrasts_all.tsv.gz", ALL_CONTRAST_FIELDS, contrast_rows)
    write_tsv(out_dir / f"{prefix}_gene_omnibus.tsv.gz", OMNIBUS_FIELDS, omnibus_rows)
    write_tsv(out_dir / f"{prefix}_formal_tissue_enriched_ase.tsv", ALL_CONTRAST_FIELDS, formal_rows)
    write_tsv(
        out_dir / f"{prefix}_recurrent_specificity_status.tsv",
        RECURRENT_STATUS_FIELDS,
        recurrent_status,
    )

    tested_contrasts = [row for row in contrast_rows if int(row["is_contrast_evaluable"])]
    tested_recurrent = [
        row
        for row in recurrent_status
        if row["specificity_test_status"] == "TESTED"
    ]
    formal_genes = {str(row["gene_id"]) for row in formal_rows}
    summary: Dict[str, object] = {
        "fragment_threshold": threshold,
        **input_counts,
        "n_unique_eligible_genes": len(genes),
        "n_cross_tissue_candidate_contrasts": len(contrast_rows),
        "n_tissue_contrasts_tested": len(tested_contrasts),
        "n_genes_omnibus_tested": len(omnibus_rows),
        "n_genes_omnibus_fdr_pass": sum(int(row["gene_omnibus_fdr_pass"]) for row in omnibus_rows),
        "n_tissue_contrasts_global_fdr_pass": sum(
            float(row["tissue_enrichment_q_global_BH"]) <= fdr for row in tested_contrasts
        ),
        "n_primary_recurrent_pairs": len(recurrent),
        "n_primary_recurrent_pairs_specificity_tested": len(tested_recurrent),
        "n_primary_recurrent_pairs_not_evaluable": len(recurrent) - len(tested_recurrent),
        "n_formal_tissue_enriched_recurrent_pairs": len(formal_rows),
        "n_formal_tissue_enriched_genes": len(formal_genes),
        "n_formal_tissue_restricted_recurrence_pairs": sum(
            int(row["is_formally_tissue_restricted_recurrence"]) for row in formal_rows
        ),
    }
    extra = {
        "formal_keys": {(str(row["gene_id"]), str(row["target_tissue"])) for row in formal_rows},
        "formal_rows": formal_rows,
        "recurrent_status": recurrent_status,
    }
    return summary, extra["formal_keys"], extra


SUMMARY_FIELDS = [
    "fragment_threshold",
    "n_eligible_gene_measurements",
    "n_core_gene_measurements",
    "n_unique_eligible_genes",
    "n_cross_tissue_candidate_contrasts",
    "n_tissue_contrasts_tested",
    "n_genes_omnibus_tested",
    "n_genes_omnibus_fdr_pass",
    "n_tissue_contrasts_global_fdr_pass",
    "n_primary_recurrent_pairs",
    "n_primary_recurrent_pairs_specificity_tested",
    "n_primary_recurrent_pairs_not_evaluable",
    "n_formal_tissue_enriched_recurrent_pairs",
    "n_formal_tissue_enriched_genes",
    "n_formal_tissue_restricted_recurrence_pairs",
    "n_threshold15_formal_pairs_retained",
    "n_formal_pairs_new_vs_threshold15",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recurrence-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--thresholds", default="15,20,30")
    parser.add_argument("--min-informative-animals", type=int, default=3)
    parser.add_argument("--min-informative-target-core", type=int, default=2)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    thresholds = tuple(sorted({int(value) for value in args.thresholds.split(",")}))
    if 15 not in thresholds or min(thresholds) < 15:
        raise SpecificityError("threshold 15 is required and thresholds must be >=15")
    if args.min_informative_animals < 2:
        raise SpecificityError("min-informative-animals must be >=2")
    if args.min_informative_target_core < 1:
        raise SpecificityError("min-informative-target-core must be positive")
    if not 0.0 < args.fdr <= 1.0:
        raise SpecificityError("fdr must be in (0,1]")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.force:
        raise SpecificityError(f"output directory is not empty: {args.out_dir}; use --force")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    formal_keys_by_threshold: Dict[int, Set[Tuple[str, str]]] = {}
    extras_by_threshold: Dict[int, Dict[str, object]] = {}
    source_files = []
    for threshold in thresholds:
        measurement_path = args.recurrence_dir / f"threshold_{threshold}_eligible_gene_measurements.tsv.gz"
        burden_path = args.recurrence_dir / f"threshold_{threshold}_sample_tissue_burden.tsv"
        recurrent_path = args.recurrence_dir / f"threshold_{threshold}_statistically_recurrent.tsv"
        for path in (measurement_path, burden_path, recurrent_path):
            if not path.is_file():
                raise SpecificityError(f"missing input: {path}")
            source_files.append(
                {
                    "fragment_threshold": threshold,
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
        print(f"[threshold={threshold}] loading and testing tissue contrasts", flush=True)
        summary, formal_keys, extra = analyze_threshold(
            threshold,
            measurement_path,
            burden_path,
            recurrent_path,
            args.out_dir,
            args.min_informative_animals,
            args.min_informative_target_core,
            args.fdr,
        )
        summaries.append(summary)
        formal_keys_by_threshold[threshold] = formal_keys
        extras_by_threshold[threshold] = extra
        print(json.dumps(summary, sort_keys=True), flush=True)

    primary_keys = formal_keys_by_threshold[15]
    for summary in summaries:
        threshold = int(summary["fragment_threshold"])
        keys = formal_keys_by_threshold[threshold]
        summary["n_threshold15_formal_pairs_retained"] = len(primary_keys & keys)
        summary["n_formal_pairs_new_vs_threshold15"] = len(keys - primary_keys)
    write_tsv(args.out_dir / "threshold_summary.tsv", SUMMARY_FIELDS, summaries)

    sensitivity_rows = []
    for row in extras_by_threshold[15]["formal_rows"]:
        key = (str(row["gene_id"]), str(row["target_tissue"]))
        sensitivity_rows.append(
            {
                **row,
                "formal_at_threshold_20": int(key in formal_keys_by_threshold.get(20, set())),
                "formal_at_threshold_30": int(key in formal_keys_by_threshold.get(30, set())),
                "formal_sensitivity_support": (
                    "THRESHOLD_15_20_30"
                    if key in formal_keys_by_threshold.get(20, set())
                    and key in formal_keys_by_threshold.get(30, set())
                    else "THRESHOLD_15_20"
                    if key in formal_keys_by_threshold.get(20, set())
                    else "THRESHOLD_15_30"
                    if key in formal_keys_by_threshold.get(30, set())
                    else "THRESHOLD_15_ONLY"
                ),
            }
        )
    write_tsv(
        args.out_dir / "primary_formal_tissue_enriched_with_sensitivity.tsv",
        ALL_CONTRAST_FIELDS
        + ["formal_at_threshold_20", "formal_at_threshold_30", "formal_sensitivity_support"],
        sensitivity_rows,
    )

    metadata = {
        "status": "PASS",
        "version": VERSION,
        "method": {
            "outcome": "binary core ASE among high-confidence-eligible animal-tissue-gene measurements",
            "repeated_measure_control": (
                "condition on each animal's total number of core-positive eligible tissues for the gene"
            ),
            "background": (
                "leave-one-gene-out sample-tissue core/eligible probability converted to conditional odds"
            ),
            "target_test": (
                "one-sided exact Poisson-binomial tail for target-tissue positives across informative animals"
            ),
            "target_fdr_scope": "global BH across all evaluable gene-tissue contrasts within threshold",
            "gene_omnibus": (
                "Bonferroni union test from the minimum pre-specified tissue-contrast p-value, "
                "then global BH across evaluable genes"
            ),
            "effect_size": "conditional common odds ratio with two-sided 95% profile-likelihood interval",
            "formal_call": (
                "primary recurrent pair; target global q<=FDR; gene omnibus global q<=FDR; "
                "minimum informative animals and target-positive informative animals; positive conditional excess"
            ),
        },
        "parameters": {
            "thresholds": thresholds,
            "min_informative_animals": args.min_informative_animals,
            "min_informative_target_core": args.min_informative_target_core,
            "fdr": args.fdr,
            "profile_likelihood_CI": 0.95,
        },
        "source_files": source_files,
        "threshold_summaries": {str(row["fragment_threshold"]): row for row in summaries},
        "interpretation": {
            "preferred_term": "formally tissue-enriched ASE",
            "not_implied": (
                "absolute tissue exclusivity, parental-origin imprinting direction, or independence from kinship"
            ),
        },
    }
    (args.out_dir / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "summaries": summaries}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SpecificityError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
