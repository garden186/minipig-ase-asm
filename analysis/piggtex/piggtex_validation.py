#!/usr/bin/env python3
"""Matched-tissue PigGTEx validation of miniature-pig ASE results.

The external validation has two deliberately separate components.

1. PigGTEx ASE is a significant-only table. Its presence is reported as
   descriptive matched-tissue support; absence is never treated as non-ASE.
2. PigGTEx permutation cis-eQTL tables contain both tested eGenes and tested
   non-eGenes. They are used for formal tissue-stratified enrichment tests,
   with the internal recurrence hypothesis universe as the background.

Primary enrichment uses a one-sided Fisher exact test. A deterministic exact
stratified sensitivity conditions on internal ASE measurement opportunity:
eligible-animal count, median fragment support, median contributing-SNP count,
and median eligible-block count. Thresholds 20 and 30 and the existing
opportunity-matched recurrence labels are additional sensitivity analyses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


VERSION = "0.1.0"
THRESHOLDS = (15, 20, 30)
TISSUE_MAP = {
    "Backfat": "Adipose",
    "Blood": "Blood",
    "Brain": "Brain",
    "Heart": "Heart",
    "Kidney": "Kidney",
    "Liver": "Liver",
    "Loin": "Muscle",
    "Lung": "Lung",
    "Lymph": "Lymph_node",
    "Spleen": "Spleen",
    "Tenderloin": "Muscle",
}
TISSUE_ORDER = tuple(TISSUE_MAP)
FORMAL_EXPECTED = 8


class ValidationError(RuntimeError):
    pass


@dataclass
class Opportunity:
    fragments: List[int]
    snps: List[int]
    blocks: List[int]


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
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for source in rows:
            row = {}
            for field in fields:
                value = source.get(field, "")
                if isinstance(value, float):
                    if math.isnan(value):
                        value = "NA"
                    elif math.isinf(value):
                        value = "Inf" if value > 0 else "-Inf"
                    else:
                        value = format(value, ".12g")
                row[field] = value
            writer.writerow(row)
    finally:
        handle.close()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def normalize_tissue(tissue: str) -> str:
    return {
        "blood": "Blood",
        "L-blood": "Blood",
        "Tenderlo": "Tenderloin",
    }.get(tissue, tissue)


def median(values: Sequence[int]) -> float:
    if not values:
        raise ValidationError("median requested for an empty opportunity vector")
    return float(statistics.median(values))


def bh_values(values: Sequence[float]) -> List[float]:
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for reverse_index in range(total - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, values[index] * total / rank, 1.0)
        adjusted[index] = running
    return adjusted


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_pmf(population: int, successes: int, draws: int) -> Dict[int, float]:
    lower = max(0, draws - (population - successes))
    upper = min(draws, successes)
    values = {}
    logs = []
    for observed in range(lower, upper + 1):
        value = (
            log_choose(successes, observed)
            + log_choose(population - successes, draws - observed)
            - log_choose(population, draws)
        )
        logs.append((observed, value))
    maximum = max(value for _, value in logs)
    scale = math.fsum(math.exp(value - maximum) for _, value in logs)
    for observed, value in logs:
        values[observed] = math.exp(value - maximum) / scale
    return values


def fisher_upper_tail(a: int, b: int, c: int, d: int) -> float:
    population = a + b + c + d
    successes = a + c
    draws = a + b
    distribution = hypergeom_pmf(population, successes, draws)
    return min(1.0, math.fsum(value for key, value in distribution.items() if key >= a))


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> Tuple[float, float, float]:
    raw = math.inf if b * c == 0 and a * d > 0 else 0.0 if a * d == 0 else a * d / (b * c)
    corrected = [float(a), float(b), float(c), float(d)]
    if min(corrected) == 0.0:
        corrected = [value + 0.5 for value in corrected]
    ca, cb, cc, cd = corrected
    estimate = ca * cd / (cb * cc)
    standard_error = math.sqrt(1.0 / ca + 1.0 / cb + 1.0 / cc + 1.0 / cd)
    lower = math.exp(math.log(estimate) - 1.959963984540054 * standard_error)
    upper = math.exp(math.log(estimate) + 1.959963984540054 * standard_error)
    return raw, lower, upper


def convolve(left: Mapping[int, float], right: Mapping[int, float]) -> Dict[int, float]:
    result: MutableMapping[int, float] = defaultdict(float)
    for left_count, left_mass in left.items():
        for right_count, right_mass in right.items():
            result[left_count + right_count] += left_mass * right_mass
    total = math.fsum(result.values())
    return {key: value / total for key, value in result.items()}


def stratified_exact_tail(strata: Sequence[Tuple[int, int, int, int]]) -> float:
    distribution: Dict[int, float] = {0: 1.0}
    observed_total = 0
    for a, b, c, d in strata:
        population = a + b + c + d
        successes = a + c
        draws = a + b
        if draws == 0:
            continue
        observed_total += a
        distribution = convolve(
            distribution, hypergeom_pmf(population, successes, draws)
        )
    return min(
        1.0,
        math.fsum(
            probability
            for count, probability in distribution.items()
            if count >= observed_total
        ),
    )


def mantel_haenszel_or(strata: Sequence[Tuple[int, int, int, int]]) -> float:
    numerator = 0.0
    denominator = 0.0
    for a, b, c, d in strata:
        total = a + b + c + d
        if total == 0:
            continue
        numerator += a * d / total
        denominator += b * c / total
    if denominator == 0.0:
        return math.inf if numerator > 0.0 else math.nan
    return numerator / denominator


def eligible_animal_bin(value: int) -> str:
    if value == 2:
        return "2"
    if value <= 4:
        return "3_4"
    if value <= 6:
        return "5_6"
    return "7_PLUS"


def fragment_bin(value: float) -> str:
    if value < 30:
        return "15_29"
    if value < 60:
        return "30_59"
    if value < 120:
        return "60_119"
    return "120_PLUS"


def snp_bin(value: float) -> str:
    if value <= 2:
        return "2"
    if value <= 3:
        return "3"
    if value <= 5:
        return "4_5"
    return "6_PLUS"


def block_bin(value: float) -> str:
    if value <= 1:
        return "1"
    if value <= 2:
        return "2"
    return "3_PLUS"


def opportunity_stratum(n_eligible: int, record: Opportunity) -> str:
    return "|".join(
        (
            eligible_animal_bin(n_eligible),
            fragment_bin(median(record.fragments)),
            snp_bin(median(record.snps)),
            block_bin(median(record.blocks)),
        )
    )


def parse_gtf_gene_names(path: Path) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    gene_to_name: Dict[str, str] = {}
    name_to_genes: MutableMapping[str, List[str]] = defaultdict(list)
    attribute_pattern = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(attribute_pattern.findall(fields[8]))
            gene_id = strip_version(attributes.get("gene_id", ""))
            gene_name = attributes.get("gene_name", "")
            if not gene_id:
                continue
            gene_to_name[gene_id] = gene_name
            if gene_name:
                name_to_genes[gene_name].append(gene_id)
    return gene_to_name, dict(name_to_genes)


def load_internal_annotations(recurrence_dir: Path) -> Dict[str, Tuple[str, str]]:
    result: Dict[str, Tuple[str, str]] = {}
    path = recurrence_dir / "threshold_15_recurrence_all.tsv.gz"
    for row in read_tsv(path):
        gene_id = strip_version(row["gene_id"])
        annotation = (row["gene_name"], row["gene_biotype"])
        if gene_id in result and result[gene_id] != annotation:
            raise ValidationError(f"inconsistent internal annotation for {gene_id}")
        result[gene_id] = annotation
    return result


def load_eqtl_tables(
    root: Path,
) -> Tuple[Dict[str, Dict[str, Dict[str, str]]], Dict[str, object]]:
    tables: Dict[str, Dict[str, Dict[str, str]]] = {}
    summary = {}
    for tissue in sorted(set(TISSUE_MAP.values())):
        path = root / f"{tissue}.cis_qtl_fdr0.05.txt.gz"
        if not path.is_file():
            raise ValidationError(f"missing PigGTEx eQTL file: {path}")
        records = {}
        n_egenes = 0
        for row in read_tsv(path):
            gene_id = strip_version(row["phenotype_id"])
            if gene_id in records:
                raise ValidationError(f"duplicate phenotype in {path}: {gene_id}")
            if row["is_eGene"] not in {"TRUE", "FALSE"}:
                raise ValidationError(f"invalid is_eGene in {path}: {row['is_eGene']}")
            records[gene_id] = row
            n_egenes += row["is_eGene"] == "TRUE"
        tables[tissue] = records
        summary[tissue] = {
            "path": str(path.resolve()),
            "sha256": sha256(path),
            "n_tested_genes": len(records),
            "n_egenes": n_egenes,
        }
    return tables, summary


def build_crosswalk(
    annotations: Mapping[str, Tuple[str, str]],
    gtf_gene_to_name: Mapping[str, str],
    gtf_name_to_genes: Mapping[str, Sequence[str]],
    external_gene_ids: set[str],
) -> Tuple[Dict[str, str], List[Dict[str, object]]]:
    mapping = {}
    rows = []
    for gene_id in sorted(annotations):
        gene_name, biotype = annotations[gene_id]
        external_id = ""
        status = "UNMAPPED"
        if gene_id in external_gene_ids or gene_id in gtf_gene_to_name:
            external_id = gene_id
            status = "DIRECT_ID"
        elif gene_name:
            candidates = sorted(
                gene
                for gene in set(gtf_name_to_genes.get(gene_name, []))
                if gene in external_gene_ids
            )
            if len(candidates) == 1:
                external_id = candidates[0]
                status = "UNIQUE_GENE_NAME"
            elif len(candidates) > 1:
                status = "AMBIGUOUS_GENE_NAME"
        if external_id:
            mapping[gene_id] = external_id
        rows.append(
            {
                "internal_gene_id": gene_id,
                "internal_gene_name": gene_name,
                "internal_gene_biotype": biotype,
                "piggtex_gene_id": external_id,
                "piggtex_gene_name_v100": gtf_gene_to_name.get(external_id, ""),
                "mapping_status": status,
            }
        )
    return mapping, rows


def load_formal(path: Path) -> List[Dict[str, str]]:
    rows = list(read_tsv(path))
    if len(rows) != FORMAL_EXPECTED:
        raise ValidationError(f"expected {FORMAL_EXPECTED} formal pairs, found {len(rows)}")
    return rows


def load_ase_significant_summary(
    path: Path,
    wanted_gene_ids: set[str],
    formal_external_keys: set[Tuple[str, str]],
) -> Tuple[Dict[Tuple[str, str], Dict[str, object]], Dict[str, object]]:
    summaries: MutableMapping[Tuple[str, str], Dict[str, object]] = defaultdict(
        lambda: {
            "n_rows": 0,
            "min_fdr": 1.0,
            "formal_samples": set(),
            "formal_breeds": set(),
            "formal_effects": [],
            "formal_totals": [],
        }
    )
    total_rows = 0
    samples = set()
    tissues = set()
    genes = set()
    minimum_fdr = 1.0
    maximum_fdr = 0.0
    malformed = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = next(handle).rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 13:
                malformed += 1
                continue
            total_rows += 1
            gene_id = strip_version(fields[3])
            fdr = float(fields[9])
            samples.add(fields[10])
            tissues.add(fields[12])
            genes.add(gene_id)
            minimum_fdr = min(minimum_fdr, fdr)
            maximum_fdr = max(maximum_fdr, fdr)
            if gene_id not in wanted_gene_ids or fields[12] not in set(TISSUE_MAP.values()):
                continue
            key = (gene_id, fields[12])
            summary = summaries[key]
            summary["n_rows"] += 1
            summary["min_fdr"] = min(float(summary["min_fdr"]), fdr)
            if key in formal_external_keys:
                summary["formal_samples"].add(fields[10])
                summary["formal_breeds"].add(fields[11])
                summary["formal_effects"].append(abs(float(fields[7])))
                summary["formal_totals"].append(int(fields[6]))
    if malformed:
        raise ValidationError(f"malformed PigGTEx ASE rows: {malformed}")
    if maximum_fdr > 0.05 + 1e-12:
        raise ValidationError(f"PigGTEx ASE table is not significant-only: max FDR={maximum_fdr}")
    qc = {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "header_field_count": len(header),
        "data_field_count": 13,
        "header_known_issue": "p-value and FDR labels are concatenated in the source header",
        "n_rows": total_rows,
        "n_unique_genes": len(genes),
        "n_unique_samples": len(samples),
        "n_tissues": len(tissues),
        "minimum_fdr": minimum_fdr,
        "maximum_fdr": maximum_fdr,
        "significant_only": True,
    }
    return dict(summaries), qc


def load_opportunity(path: Path) -> Dict[Tuple[str, str], Opportunity]:
    result: MutableMapping[Tuple[str, str], Opportunity] = {}
    for row in read_tsv(path):
        key = (strip_version(row["gene_id"]), normalize_tissue(row["tissue"]))
        record = result.get(key)
        if record is None:
            record = Opportunity([], [], [])
            result[key] = record
        record.fragments.append(int(row["max_eligible_total_fragments"]))
        record.snps.append(int(row["max_eligible_contributing_snps"]))
        record.blocks.append(int(row["n_eligible_blocks"]))
    return dict(result)


def load_recurrence_rows(path: Path) -> List[Dict[str, str]]:
    rows = []
    for row in read_tsv(path):
        copy = dict(row)
        copy["gene_id"] = strip_version(copy["gene_id"])
        copy["tissue"] = normalize_tissue(copy["tissue"])
        rows.append(copy)
    return rows


def load_matched_labels(path: Path) -> Dict[Tuple[str, str], int]:
    result = {}
    for row in read_tsv(path):
        key = (strip_version(row["gene_id"]), normalize_tissue(row["tissue"]))
        result[key] = int(row["is_statistically_recurrent_core"])
    return result


UNIVERSE_FIELDS = [
    "analysis_mode",
    "fragment_threshold",
    "internal_gene_id",
    "internal_gene_name",
    "internal_tissue",
    "piggtex_gene_id",
    "piggtex_tissue",
    "mapping_status",
    "n_samples_eligible",
    "n_samples_core",
    "is_internal_recurrent",
    "median_max_fragments",
    "median_max_contributing_snps",
    "median_eligible_blocks",
    "opportunity_stratum",
    "piggtex_eqtl_tested",
    "piggtex_is_eGene",
    "piggtex_lead_variant",
    "piggtex_eqtl_slope",
    "piggtex_eqtl_pval_adj_BH",
    "piggtex_ase_observed_significant_only",
    "piggtex_ase_significant_rows",
]


STRATA_FIELDS = [
    "analysis_mode",
    "fragment_threshold",
    "internal_tissue",
    "piggtex_tissue",
    "opportunity_stratum",
    "recurrent_eGene",
    "recurrent_non_eGene",
    "nonrecurrent_eGene",
    "nonrecurrent_non_eGene",
]


ENRICHMENT_FIELDS = [
    "analysis_mode",
    "fragment_threshold",
    "internal_tissue",
    "piggtex_tissue",
    "n_internal_hypotheses",
    "n_gene_mapped",
    "n_double_testable",
    "n_recurrent_double_testable",
    "n_nonrecurrent_double_testable",
    "recurrent_eGene",
    "recurrent_non_eGene",
    "nonrecurrent_eGene",
    "nonrecurrent_non_eGene",
    "recurrent_eGene_fraction",
    "nonrecurrent_eGene_fraction",
    "fisher_odds_ratio",
    "fisher_odds_ratio_CI95_lower",
    "fisher_odds_ratio_CI95_upper",
    "fisher_p_one_sided",
    "fisher_q_BH_across_tissues",
    "n_opportunity_strata",
    "n_informative_opportunity_strata",
    "mantel_haenszel_odds_ratio",
    "opportunity_stratified_exact_p_one_sided",
    "opportunity_stratified_q_BH_across_tissues",
    "positive_direction",
    "primary_fdr_pass",
    "opportunity_stratified_fdr_pass",
]


def analyze_mode(
    mode: str,
    threshold: int,
    recurrence_rows: Sequence[Mapping[str, str]],
    opportunities: Mapping[Tuple[str, str], Opportunity],
    mapping: Mapping[str, str],
    mapping_status: Mapping[str, str],
    eqtl_tables: Mapping[str, Mapping[str, Mapping[str, str]]],
    ase_summary: Mapping[Tuple[str, str], Mapping[str, object]],
    matched_labels: Mapping[Tuple[str, str], int] | None = None,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    universe = []
    strata_output = []
    enrichments = []
    rows_by_tissue: MutableMapping[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in recurrence_rows:
        rows_by_tissue[str(row["tissue"])].append(row)

    for tissue in TISSUE_ORDER:
        external_tissue = TISSUE_MAP[tissue]
        eqtl = eqtl_tables[external_tissue]
        internal_rows = rows_by_tissue.get(tissue, [])
        table = [0, 0, 0, 0]
        mapped_count = 0
        strata_counts: MutableMapping[str, List[int]] = defaultdict(lambda: [0, 0, 0, 0])
        for row in internal_rows:
            gene_id = str(row["gene_id"])
            external_gene = mapping.get(gene_id, "")
            if external_gene:
                mapped_count += 1
            if not external_gene or external_gene not in eqtl:
                continue
            opportunity = opportunities.get((gene_id, tissue))
            if opportunity is None:
                raise ValidationError(f"missing opportunity for {threshold}|{gene_id}|{tissue}")
            n_eligible = int(row["n_samples_eligible"])
            if len(opportunity.fragments) != n_eligible:
                raise ValidationError(
                    f"eligible-animal mismatch for {threshold}|{gene_id}|{tissue}: "
                    f"{len(opportunity.fragments)} != {n_eligible}"
                )
            recurrent = int(row["is_statistically_recurrent_core"])
            if matched_labels is not None:
                recurrent = int(matched_labels.get((gene_id, tissue), 0))
            is_egene = int(eqtl[external_gene]["is_eGene"] == "TRUE")
            index = 0 if recurrent and is_egene else 1 if recurrent else 2 if is_egene else 3
            table[index] += 1
            stratum = opportunity_stratum(n_eligible, opportunity)
            strata_counts[stratum][index] += 1
            external_ase = ase_summary.get((external_gene, external_tissue), {})
            universe.append(
                {
                    "analysis_mode": mode,
                    "fragment_threshold": threshold,
                    "internal_gene_id": gene_id,
                    "internal_gene_name": row["gene_name"],
                    "internal_tissue": tissue,
                    "piggtex_gene_id": external_gene,
                    "piggtex_tissue": external_tissue,
                    "mapping_status": mapping_status.get(gene_id, "UNMAPPED"),
                    "n_samples_eligible": n_eligible,
                    "n_samples_core": int(row["n_samples_core"]),
                    "is_internal_recurrent": recurrent,
                    "median_max_fragments": median(opportunity.fragments),
                    "median_max_contributing_snps": median(opportunity.snps),
                    "median_eligible_blocks": median(opportunity.blocks),
                    "opportunity_stratum": stratum,
                    "piggtex_eqtl_tested": 1,
                    "piggtex_is_eGene": is_egene,
                    "piggtex_lead_variant": eqtl[external_gene]["variant_id"],
                    "piggtex_eqtl_slope": float(eqtl[external_gene]["slope"]),
                    "piggtex_eqtl_pval_adj_BH": float(eqtl[external_gene]["pval_adj_BH"]),
                    "piggtex_ase_observed_significant_only": int(bool(external_ase)),
                    "piggtex_ase_significant_rows": int(external_ase.get("n_rows", 0)),
                }
            )

        a, b, c, d = table
        if a + b == 0 or c + d == 0:
            raise ValidationError(f"non-informative recurrent table for {mode}|{tissue}: {table}")
        fisher_p = fisher_upper_tail(a, b, c, d)
        odds_ratio, ci_lower, ci_upper = odds_ratio_ci(a, b, c, d)
        strata = [tuple(values) for values in strata_counts.values()]
        stratified_p = stratified_exact_tail(strata)
        mh_or = mantel_haenszel_or(strata)
        for stratum, values in sorted(strata_counts.items()):
            strata_output.append(
                {
                    "analysis_mode": mode,
                    "fragment_threshold": threshold,
                    "internal_tissue": tissue,
                    "piggtex_tissue": external_tissue,
                    "opportunity_stratum": stratum,
                    "recurrent_eGene": values[0],
                    "recurrent_non_eGene": values[1],
                    "nonrecurrent_eGene": values[2],
                    "nonrecurrent_non_eGene": values[3],
                }
            )
        enrichments.append(
            {
                "analysis_mode": mode,
                "fragment_threshold": threshold,
                "internal_tissue": tissue,
                "piggtex_tissue": external_tissue,
                "n_internal_hypotheses": len(internal_rows),
                "n_gene_mapped": mapped_count,
                "n_double_testable": sum(table),
                "n_recurrent_double_testable": a + b,
                "n_nonrecurrent_double_testable": c + d,
                "recurrent_eGene": a,
                "recurrent_non_eGene": b,
                "nonrecurrent_eGene": c,
                "nonrecurrent_non_eGene": d,
                "recurrent_eGene_fraction": a / (a + b),
                "nonrecurrent_eGene_fraction": c / (c + d),
                "fisher_odds_ratio": odds_ratio,
                "fisher_odds_ratio_CI95_lower": ci_lower,
                "fisher_odds_ratio_CI95_upper": ci_upper,
                "fisher_p_one_sided": fisher_p,
                "fisher_q_BH_across_tissues": math.nan,
                "n_opportunity_strata": len(strata),
                "n_informative_opportunity_strata": sum(
                    int((x[0] + x[1]) > 0 and (x[2] + x[3]) > 0) for x in strata
                ),
                "mantel_haenszel_odds_ratio": mh_or,
                "opportunity_stratified_exact_p_one_sided": stratified_p,
                "opportunity_stratified_q_BH_across_tissues": math.nan,
                "positive_direction": int(a / (a + b) > c / (c + d)),
                "primary_fdr_pass": 0,
                "opportunity_stratified_fdr_pass": 0,
            }
        )
    fisher_q = bh_values([float(row["fisher_p_one_sided"]) for row in enrichments])
    stratified_q = bh_values(
        [float(row["opportunity_stratified_exact_p_one_sided"]) for row in enrichments]
    )
    for row, q_fisher, q_stratified in zip(enrichments, fisher_q, stratified_q):
        row["fisher_q_BH_across_tissues"] = q_fisher
        row["opportunity_stratified_q_BH_across_tissues"] = q_stratified
        row["primary_fdr_pass"] = int(q_fisher <= 0.05 and row["positive_direction"])
        row["opportunity_stratified_fdr_pass"] = int(
            q_stratified <= 0.05 and row["positive_direction"]
        )
    return universe, strata_output, enrichments


FORMAL_FIELDS = [
    "internal_gene_id",
    "gene_name",
    "internal_target_tissue",
    "piggtex_gene_id",
    "piggtex_target_tissue",
    "mapping_status",
    "formal_sensitivity_support",
    "target_core_animals",
    "target_eligible_animals",
    "other_core_measurements",
    "other_eligible_measurements",
    "internal_tissue_enrichment_q",
    "internal_gene_omnibus_q",
    "piggtex_ase_observed_significant_only",
    "piggtex_ase_significant_samples",
    "piggtex_ase_breeds",
    "piggtex_ase_median_abs_log2_aFC",
    "piggtex_ase_median_total_count",
    "piggtex_eqtl_tested",
    "piggtex_is_eGene",
    "piggtex_lead_variant",
    "piggtex_tss_distance",
    "piggtex_eqtl_slope",
    "piggtex_eqtl_pval_adj_BH",
    "external_support_class",
]


def formal_external_rows(
    formal: Sequence[Mapping[str, str]],
    mapping: Mapping[str, str],
    mapping_status: Mapping[str, str],
    eqtl_tables: Mapping[str, Mapping[str, Mapping[str, str]]],
    ase_summary: Mapping[Tuple[str, str], Mapping[str, object]],
) -> List[Dict[str, object]]:
    result = []
    for row in formal:
        internal_gene = strip_version(row["gene_id"])
        external_gene = mapping.get(internal_gene, "")
        internal_tissue = normalize_tissue(row["target_tissue"])
        external_tissue = TISSUE_MAP[internal_tissue]
        ase = ase_summary.get((external_gene, external_tissue), {})
        eqtl = eqtl_tables[external_tissue].get(external_gene)
        ase_observed = bool(ase)
        is_egene = bool(eqtl and eqtl["is_eGene"] == "TRUE")
        support = (
            "MATCHED_ASE_AND_CIS_EGENE"
            if ase_observed and is_egene
            else "MATCHED_ASE_ONLY"
            if ase_observed
            else "CIS_EGENE_ONLY"
            if is_egene
            else "NO_EXTERNAL_SUPPORT_OBSERVED"
        )
        effects = list(ase.get("formal_effects", []))
        totals = list(ase.get("formal_totals", []))
        result.append(
            {
                "internal_gene_id": internal_gene,
                "gene_name": row["gene_name"],
                "internal_target_tissue": internal_tissue,
                "piggtex_gene_id": external_gene,
                "piggtex_target_tissue": external_tissue,
                "mapping_status": mapping_status.get(internal_gene, "UNMAPPED"),
                "formal_sensitivity_support": row["formal_sensitivity_support"],
                "target_core_animals": int(row["n_target_core_animals_total"]),
                "target_eligible_animals": int(row["n_target_eligible_animals_total"]),
                "other_core_measurements": int(row["n_other_core_measurements_total"]),
                "other_eligible_measurements": int(row["n_other_eligible_measurements_total"]),
                "internal_tissue_enrichment_q": float(row["tissue_enrichment_q_global_BH"]),
                "internal_gene_omnibus_q": float(row["gene_omnibus_q_global_BH"]),
                "piggtex_ase_observed_significant_only": int(ase_observed),
                "piggtex_ase_significant_samples": len(ase.get("formal_samples", set())),
                "piggtex_ase_breeds": len(ase.get("formal_breeds", set())),
                "piggtex_ase_median_abs_log2_aFC": statistics.median(effects) if effects else math.nan,
                "piggtex_ase_median_total_count": statistics.median(totals) if totals else math.nan,
                "piggtex_eqtl_tested": int(eqtl is not None),
                "piggtex_is_eGene": int(is_egene),
                "piggtex_lead_variant": eqtl["variant_id"] if eqtl else "",
                "piggtex_tss_distance": int(eqtl["tss_distance"]) if eqtl else "",
                "piggtex_eqtl_slope": float(eqtl["slope"]) if eqtl else math.nan,
                "piggtex_eqtl_pval_adj_BH": float(eqtl["pval_adj_BH"]) if eqtl else math.nan,
                "external_support_class": support,
            }
        )
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recurrence-dir", type=Path, required=True)
    parser.add_argument("--specificity-dir", type=Path, required=True)
    parser.add_argument("--piggtex-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        if not args.force:
            raise ValidationError(f"output directory is not empty: {args.out_dir}; use --force")
        shutil.rmtree(args.out_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    formal_path = args.specificity_dir / "primary_formal_tissue_enriched_with_sensitivity.tsv"
    ase_path = args.piggtex_dir / "6020.ASE.gene.sum.sort.txt.gz"
    gtf_path = args.piggtex_dir / "Sus_scrofa.Sscrofa11.1.100.gtf.gz"
    eqtl_root = args.piggtex_dir / "PigGTEx_v0.permutations_eQTL"
    for path in (formal_path, ase_path, gtf_path, eqtl_root):
        if not path.exists():
            raise ValidationError(f"missing input: {path}")

    print("[1/6] loading internal annotations and formal calls", flush=True)
    annotations = load_internal_annotations(args.recurrence_dir)
    formal = load_formal(formal_path)

    print("[2/6] loading PigGTEx eQTL tables and v100 annotation", flush=True)
    eqtl_tables, eqtl_qc = load_eqtl_tables(eqtl_root)
    external_gene_ids = set().union(*(set(table) for table in eqtl_tables.values()))
    gtf_gene_to_name, gtf_name_to_genes = parse_gtf_gene_names(gtf_path)
    mapping, crosswalk_rows = build_crosswalk(
        annotations,
        gtf_gene_to_name,
        gtf_name_to_genes,
        external_gene_ids,
    )
    mapping_status = {
        str(row["internal_gene_id"]): str(row["mapping_status"])
        for row in crosswalk_rows
    }

    formal_external_keys = {
        (
            mapping.get(strip_version(row["gene_id"]), ""),
            TISSUE_MAP[normalize_tissue(row["target_tissue"])],
        )
        for row in formal
    }
    print("[3/6] streaming significant-only PigGTEx ASE table", flush=True)
    ase_summary, ase_qc = load_ase_significant_summary(
        ase_path,
        set(mapping.values()),
        formal_external_keys,
    )

    all_universe = []
    all_strata = []
    all_enrichment = []
    source_files = [
        {"role": "formal_tissue_enrichment", "path": str(formal_path.resolve()), "sha256": sha256(formal_path)},
        {"role": "piggtex_ase_significant_only", "path": str(ase_path.resolve()), "sha256": sha256(ase_path)},
        {"role": "piggtex_gtf_v100", "path": str(gtf_path.resolve()), "sha256": sha256(gtf_path)},
    ]

    print("[4/6] tissue-stratified enrichment at thresholds 15/20/30", flush=True)
    for threshold in THRESHOLDS:
        recurrence_path = args.recurrence_dir / f"threshold_{threshold}_recurrence_all.tsv.gz"
        opportunity_path = args.recurrence_dir / f"threshold_{threshold}_eligible_gene_measurements.tsv.gz"
        source_files.extend(
            [
                {"role": f"threshold_{threshold}_recurrence", "path": str(recurrence_path.resolve()), "sha256": sha256(recurrence_path)},
                {"role": f"threshold_{threshold}_opportunity", "path": str(opportunity_path.resolve()), "sha256": sha256(opportunity_path)},
            ]
        )
        recurrence_rows = load_recurrence_rows(recurrence_path)
        opportunities = load_opportunity(opportunity_path)
        universe, strata, enrichment = analyze_mode(
            f"PRIMARY_THRESHOLD_{threshold}",
            threshold,
            recurrence_rows,
            opportunities,
            mapping,
            mapping_status,
            eqtl_tables,
            ase_summary,
        )
        all_universe.extend(universe)
        all_strata.extend(strata)
        all_enrichment.extend(enrichment)

        if threshold == 15:
            matched_path = args.recurrence_dir / "opportunity_matched_recurrence_all.tsv.gz"
            source_files.append(
                {"role": "opportunity_matched_recurrence_labels", "path": str(matched_path.resolve()), "sha256": sha256(matched_path)}
            )
            matched_labels = load_matched_labels(matched_path)
            universe, strata, enrichment = analyze_mode(
                "OPPORTUNITY_MATCHED_RECURRENCE_LABEL",
                threshold,
                recurrence_rows,
                opportunities,
                mapping,
                mapping_status,
                eqtl_tables,
                ase_summary,
                matched_labels=matched_labels,
            )
            all_universe.extend(universe)
            all_strata.extend(strata)
            all_enrichment.extend(enrichment)

    print("[5/6] writing candidate and cohort validation tables", flush=True)
    formal_rows = formal_external_rows(
        formal, mapping, mapping_status, eqtl_tables, ase_summary
    )
    write_tsv(
        args.out_dir / "tissue_mapping.tsv",
        ["internal_tissue", "piggtex_tissue", "mapping_role"],
        (
            {
                "internal_tissue": tissue,
                "piggtex_tissue": TISSUE_MAP[tissue],
                "mapping_role": "PRIMARY_PREDEFINED",
            }
            for tissue in TISSUE_ORDER
        ),
    )
    write_tsv(
        args.out_dir / "gene_crosswalk.tsv.gz",
        [
            "internal_gene_id",
            "internal_gene_name",
            "internal_gene_biotype",
            "piggtex_gene_id",
            "piggtex_gene_name_v100",
            "mapping_status",
        ],
        crosswalk_rows,
    )
    write_tsv(args.out_dir / "formal_8_external_validation.tsv", FORMAL_FIELDS, formal_rows)
    write_tsv(
        args.out_dir / "matched_double_testable_universe.tsv.gz",
        UNIVERSE_FIELDS,
        all_universe,
    )
    write_tsv(
        args.out_dir / "opportunity_strata_counts.tsv.gz",
        STRATA_FIELDS,
        all_strata,
    )
    write_tsv(
        args.out_dir / "tissue_eqtl_enrichment.tsv",
        ENRICHMENT_FIELDS,
        all_enrichment,
    )

    mapping_counts = Counter(row["mapping_status"] for row in crosswalk_rows)
    primary_15 = [
        row for row in all_enrichment if row["analysis_mode"] == "PRIMARY_THRESHOLD_15"
    ]
    both_pass = [
        row
        for row in primary_15
        if row["primary_fdr_pass"] and row["opportunity_stratified_fdr_pass"]
    ]
    candidate_counts = Counter(row["external_support_class"] for row in formal_rows)
    qc = {
        "status": "PASS",
        "version": VERSION,
        "definitions": {
            "piggtex_ase": "descriptive significant-only matched-tissue support; absence is not non-ASE",
            "formal_external_test": "one-sided enrichment of PigGTEx cis-eGenes among internal recurrent versus eligible non-recurrent gene-tissue hypotheses",
            "double_testable_universe": "internal recurrence hypothesis plus a tested phenotype row in the mapped PigGTEx tissue",
            "opportunity_adjustment": "exact stratified hypergeometric tail conditioning on internal eligible-animal, fragment, SNP, and block bins",
            "multiple_testing": "BH across the 11 predefined internal-tissue tests within each analysis mode",
        },
        "tissue_mapping": TISSUE_MAP,
        "gene_mapping_counts": dict(sorted(mapping_counts.items())),
        "n_internal_annotated_genes": len(annotations),
        "n_internal_genes_mapped": len(mapping),
        "piggtex_ase_qc": ase_qc,
        "piggtex_eqtl_qc": eqtl_qc,
        "n_formal_pairs": len(formal_rows),
        "formal_external_support_counts": dict(sorted(candidate_counts.items())),
        "primary_threshold15": {
            "n_tissues_primary_fdr_pass": sum(int(row["primary_fdr_pass"]) for row in primary_15),
            "n_tissues_opportunity_stratified_fdr_pass": sum(
                int(row["opportunity_stratified_fdr_pass"]) for row in primary_15
            ),
            "n_tissues_both_pass": len(both_pass),
            "both_pass_tissues": [row["internal_tissue"] for row in both_pass],
        },
        "row_counts": {
            "gene_crosswalk": len(crosswalk_rows),
            "formal_8_external_validation": len(formal_rows),
            "matched_double_testable_universe": len(all_universe),
            "opportunity_strata_counts": len(all_strata),
            "tissue_eqtl_enrichment": len(all_enrichment),
        },
        "source_files": source_files,
    }
    (args.out_dir / "analysis_qc.json").write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("[6/6] complete", flush=True)
    print(
        json.dumps(
            {
                "status": "PASS",
                "formal_support": dict(candidate_counts),
                "threshold15_both_pass_tissues": [
                    row["internal_tissue"] for row in both_pass
                ],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
