#!/usr/bin/env python3
"""Audit canonical recurrent ASE SNPs against independent GATK and PigGTEx evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "ase_canonical_snp_piggtex_evidence_v0.1.0"

PAIR_REQUIRED = {"gene_id", "gene_name", "gene_biotype", "tissue"}
TISSUE_PAIR_REQUIRED = {"gene_id", "gene_name", "gene_biotype", "target_tissue", "pair_level_primary"}
CANDIDATE_REQUIRED = {
    "internal_gene_id", "internal_tissue", "piggtex_gene_id", "piggtex_tissue",
    "chrom", "position", "allele_1", "allele_2", "internal_variant_unordered_key",
}
CANONICAL_SNP_REQUIRED = {
    "gene_id", "gene_name", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "n_animals_testable", "animals_testable", "n_animals_with_block_concordant_snp_ase",
    "animals_with_block_concordant_snp_ase", "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global", "is_statistically_recurrent_snp_ase",
}
BLOCK_REQUIRED = {
    "sample", "tissue", "gene_id", "gene_name", "variant_id", "contig", "position",
    "ref_allele", "alt_allele", "gatk_site_testable", "gatk_snp_level_ase",
    "gatk_block_direction_status", "gatk_q_bh_within_animal_tissue",
}
GATK_RECURRENCE_REQUIRED = {
    "tissue", "contig", "position", "ref_allele", "alt_allele", "n_animals_testable",
    "animals_testable", "n_animals_with_snp_level_ase", "animals_with_snp_level_ase",
    "recurrence_p_poisson_binomial", "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
}
GATK_SITE_REQUIRED = {
    "internal_tissue", "contig", "position", "ref_allele", "alt_allele",
    "n_internal_gatk_ase_measurements", "internal_gatk_ase_samples",
}
EQTL_REQUIRED = {
    "evidence_source", "piggtex_tissue", "piggtex_gene_id", "piggtex_variant_id",
    "internal_variant_unordered_key", "internal_gene_id", "internal_tissue",
    "pval_nominal", "slope", "slope_se", "fine_mapping_prob", "fine_mapping_cs",
}
PIGGTEX_SITE_REQUIRED = {"chr", "position", "refAllele", "altAllele", "FDR", "SampleID", "Tissue"}

SITE_HEADER = [
    "gene_id", "gene_name", "gene_biotype", "internal_tissue", "piggtex_tissue",
    "piggtex_gene_ids", "variant_unordered_key", "chrom", "position", "allele_1", "allele_2",
    "current_recurrent_gene_tissue_pair", "current_tissue_enriched_pair",
    "canonical_recurrent_gene_tissue_snp", "canonical_snp_n_testable_animals",
    "canonical_snp_animals_testable", "canonical_snp_n_positive_animals",
    "canonical_snp_positive_animals", "canonical_snp_p", "canonical_snp_q",
    "n_primary_block_contexts", "n_primary_block_testable_samples", "primary_block_testable_samples",
    "n_primary_block_gatk_ase_samples", "primary_block_gatk_ase_samples",
    "n_primary_block_direction_concordant_samples", "primary_block_direction_concordant_samples",
    "n_primary_block_direction_discordant_samples", "primary_block_direction_discordant_samples",
    "n_primary_block_gatk_testable_not_ase_samples", "primary_block_gatk_testable_not_ase_samples",
    "n_primary_block_not_gatk_testable_samples", "primary_block_not_gatk_testable_samples",
    "n_independent_gatk_ase_samples", "independent_gatk_ase_samples",
    "independent_gatk_recurrence_evaluable", "independent_gatk_n_testable_animals",
    "independent_gatk_n_positive_animals", "independent_gatk_p", "independent_gatk_q",
    "independent_gatk_statistically_recurrent",
    "piggtex_significant_site_rows_any_tissue", "piggtex_significant_site_rows_same_mapped_tissue",
    "piggtex_significant_site_samples_same_mapped_tissue", "piggtex_min_fdr_any_tissue",
    "piggtex_min_fdr_same_mapped_tissue", "piggtex_site_ase_support_status",
    "server_eqtl_audit_site_status", "piggtex_exact_significant_eqtl",
    "piggtex_exact_significant_eqtl_min_p", "piggtex_exact_significant_eqtl_slopes",
    "piggtex_exact_finemapped_eqtl", "piggtex_exact_finemapped_eqtl_max_pip",
    "piggtex_exact_finemapped_eqtl_credible_sets", "piggtex_exact_finemapped_eqtl_pip_ge_0_5",
    "piggtex_exact_finemapped_eqtl_pip_ge_0_9", "support_signature",
]

CANONICAL_AUDIT_HEADER = [
    "gene_id", "gene_name", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "n_animals_testable", "animals_testable", "n_animals_with_block_concordant_snp_ase",
    "animals_with_block_concordant_snp_ase", "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global", "pair_in_current_646_recurrence", "pair_in_current_39_tissue_enrichment",
    "n_primary_block_gatk_ase_samples", "primary_block_gatk_ase_samples",
    "n_primary_block_direction_concordant_samples", "primary_block_direction_concordant_samples",
    "n_primary_block_direction_discordant_samples", "primary_block_direction_discordant_samples",
    "n_independent_gatk_ase_samples", "independent_gatk_ase_samples",
    "independent_gatk_recurrence_evaluable", "independent_gatk_n_testable_animals",
    "independent_gatk_n_positive_animals", "independent_gatk_p", "independent_gatk_q",
    "independent_gatk_statistically_recurrent", "piggtex_mapped_tissue",
    "piggtex_significant_site_rows_same_mapped_tissue", "piggtex_min_fdr_same_mapped_tissue",
    "piggtex_exact_significant_eqtl", "piggtex_exact_finemapped_eqtl_max_pip",
    "server_eqtl_audit_site_status", "support_signature",
]

PAIR_SUMMARY_HEADER = [
    "gene_id", "gene_name", "target_tissue", "piggtex_tissue", "n_server_audit_candidate_sites",
    "n_canonical_recurrent_snp_sites", "n_sites_with_primary_block_gatk_ase",
    "n_sites_with_primary_block_direction_concordant_gatk_ase",
    "n_independent_gatk_recurrent_sites", "n_sites_with_piggtex_significant_ase_same_mapped_tissue",
    "n_sites_with_exact_significant_eqtl", "n_sites_with_exact_finemapped_eqtl",
    "n_sites_with_exact_finemapped_eqtl_pip_ge_0_5", "n_sites_with_exact_finemapped_eqtl_pip_ge_0_9",
    "n_full_site_resolved_candidates_pip_ge_0_5", "canonical_recurrent_snp_site_ids",
    "direction_concordant_gatk_site_ids", "same_tissue_piggtex_site_ase_ids",
    "exact_significant_eqtl_site_ids", "high_pip_full_site_resolved_candidate_ids",
    "best_supported_site_id", "pair_site_evidence_interpretation",
]

STAGE_HEADER = [
    "scope", "n_sites", "n_canonical_recurrent_snp", "n_with_primary_block_gatk_ase",
    "n_with_primary_block_direction_concordant_gatk_ase", "n_independent_gatk_recurrent",
    "n_with_piggtex_significant_site_ase_same_mapped_tissue", "n_with_exact_significant_eqtl",
    "n_with_exact_finemapped_eqtl", "n_with_exact_finemapped_eqtl_pip_ge_0_5",
    "n_with_exact_finemapped_eqtl_pip_ge_0_9", "n_full_site_resolved_candidates_pip_ge_0_5",
]


class AnalysisError(RuntimeError):
    """Raised when an input or cross-source identity invariant is violated."""


@contextmanager
def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def iter_tsv(path: Path, required: Set[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AnalysisError(f"{path}: missing columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value or "" for key, value in row.items()}


def normalize_tissue(value: str) -> str:
    stripped = value.strip()
    return {
        "blood": "Blood", "l-blood": "Blood", "tenderlo": "Tenderloin",
        "tenderloin": "Tenderloin",
    }.get(stripped.lower(), stripped)


def parse_bool(value: str, field: str, context: str) -> int:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return 1
    if normalized in {"0", "false", "no"}:
        return 0
    raise AnalysisError(f"{context}: invalid Boolean {field}={value!r}")


def unordered_site(contig: str, position: int, allele_a: str, allele_b: str) -> Tuple[str, int, str, str]:
    allele_1, allele_2 = sorted((allele_a.upper(), allele_b.upper()))
    return contig.removeprefix("chr"), int(position), allele_1, allele_2


def parse_variant_id(value: str) -> Tuple[str, int, str, str]:
    fields = value.rsplit("_", 3)
    if len(fields) != 4:
        raise AnalysisError(f"cannot parse variant ID: {value!r}")
    return unordered_site(fields[0], int(fields[1]), fields[2], fields[3])


def site_id(site: Tuple[str, int, str, str]) -> str:
    return f"{site[0]}_{site[1]}_{site[2]}_{site[3]}"


def fmt_float(value: object) -> str:
    if value is None or value == "":
        return ""
    number = float(value)
    return "" if not math.isfinite(number) else format(number, ".12g")


def sorted_csv(values: Iterable[str]) -> str:
    return ",".join(sorted(set(value for value in values if value)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_metadata(path: Path) -> Dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Mapping[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in header})
            count += 1
    return count


def new_block_record() -> Dict[str, object]:
    return {
        "contexts": set(), "testable": set(), "ase": set(), "concordant": set(),
        "discordant": set(), "testable_not_ase": set(), "not_testable": set(),
    }


def aggregate_block_rows(path: Path) -> Dict[Tuple[str, str, str, int, str, str], Dict[str, object]]:
    result: Dict[Tuple[str, str, str, int, str, str], Dict[str, object]] = {}
    seen_sample_keys: Dict[Tuple[str, str, str, str, int, str, str], str] = {}
    for row in iter_tsv(path, BLOCK_REQUIRED):
        sample = row["sample"].strip()
        tissue = normalize_tissue(row["tissue"])
        site = unordered_site(row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"])
        key = (row["gene_id"].strip(), tissue, *site)
        record = result.setdefault(key, new_block_record())
        context = (sample, row["variant_id"], row.get("haplotype_block_id", ""), row.get("measurement_id", ""))
        record["contexts"].add(context)
        sample_key = (sample, *key)
        status = row["gatk_block_direction_status"].strip()
        previous = seen_sample_keys.get(sample_key)
        if previous is not None and previous != status:
            raise AnalysisError(f"conflicting GATK block status for {sample_key}: {previous} versus {status}")
        seen_sample_keys[sample_key] = status
        if parse_bool(row["gatk_site_testable"], "gatk_site_testable", str(sample_key)):
            record["testable"].add(sample)
        else:
            record["not_testable"].add(sample)
        if parse_bool(row["gatk_snp_level_ase"], "gatk_snp_level_ase", str(sample_key)):
            record["ase"].add(sample)
        if status == "DIRECTION_CONCORDANT":
            record["concordant"].add(sample)
        elif status == "DIRECTION_DISCORDANT":
            record["discordant"].add(sample)
        elif status == "GATK_TESTABLE_NOT_ASE":
            record["testable_not_ase"].add(sample)
        elif status != "NOT_GATK_TESTABLE":
            raise AnalysisError(f"unexpected GATK block status: {status}")
    return result


def load_gatk_site_calls(path: Path) -> Dict[Tuple[str, str, int, str, str], Dict[str, str]]:
    result = {}
    for row in iter_tsv(path, GATK_SITE_REQUIRED):
        tissue = normalize_tissue(row["internal_tissue"])
        site = unordered_site(row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"])
        key = (tissue, *site)
        if key in result:
            raise AnalysisError(f"duplicate independent GATK site-call row: {key}")
        result[key] = row
    return result


def load_gatk_recurrence(path: Path) -> Dict[Tuple[str, str, int, str, str], Dict[str, str]]:
    result = {}
    for row in iter_tsv(path, GATK_RECURRENCE_REQUIRED):
        tissue = normalize_tissue(row["tissue"])
        site = unordered_site(row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"])
        key = (tissue, *site)
        if key in result:
            raise AnalysisError(f"duplicate independent GATK recurrence row: {key}")
        result[key] = row
    return result


def aggregate_eqtl(path: Path) -> Dict[Tuple[str, str, str, int, str, str], Dict[str, object]]:
    result = defaultdict(lambda: {
        "significant_rows": 0, "p_values": [], "slopes": set(), "finemapped_rows": 0,
        "pips": [], "credible_sets": set(), "piggtex_gene_ids": set(), "piggtex_tissues": set(),
    })
    for row in iter_tsv(path, EQTL_REQUIRED):
        tissue = normalize_tissue(row["internal_tissue"])
        site = parse_variant_id(row["internal_variant_unordered_key"])
        key = (row["internal_gene_id"].strip(), tissue, *site)
        record = result[key]
        record["piggtex_gene_ids"].add(row["piggtex_gene_id"].strip())
        record["piggtex_tissues"].add(row["piggtex_tissue"].strip())
        source = row["evidence_source"].strip()
        if source == "significant_eQTL":
            record["significant_rows"] += 1
            if row["pval_nominal"]:
                record["p_values"].append(float(row["pval_nominal"]))
            if row["slope"]:
                record["slopes"].add(fmt_float(row["slope"]))
        elif source == "finemapped_eQTL":
            record["finemapped_rows"] += 1
            if row["fine_mapping_prob"]:
                record["pips"].append(float(row["fine_mapping_prob"]))
            if row["fine_mapping_cs"]:
                record["credible_sets"].add(row["fine_mapping_cs"])
        else:
            raise AnalysisError(f"unexpected eQTL evidence source: {source}")
    return dict(result)


def audit_piggtex_sites(
    path: Path, target_sites: Set[Tuple[str, int, str, str]], tissue_map: Mapping[str, str]
) -> Tuple[Dict[Tuple[str, int, str, str], Dict[str, object]], Dict[str, object]]:
    result = defaultdict(lambda: {
        "rows_any": 0, "min_fdr_any": math.inf, "rows_by_tissue": Counter(),
        "min_fdr_by_tissue": {}, "samples_by_tissue": defaultdict(set),
    })
    total_rows = 0
    matched_rows = 0
    max_fdr = -math.inf
    tissues_seen = set()
    for row in iter_tsv(path, PIGGTEX_SITE_REQUIRED):
        total_rows += 1
        fdr = float(row["FDR"])
        max_fdr = max(max_fdr, fdr)
        external_tissue = row["Tissue"].strip()
        tissues_seen.add(external_tissue)
        site = unordered_site(row["chr"], int(row["position"]), row["refAllele"], row["altAllele"])
        if site not in target_sites:
            continue
        matched_rows += 1
        record = result[site]
        record["rows_any"] += 1
        record["min_fdr_any"] = min(record["min_fdr_any"], fdr)
        record["rows_by_tissue"][external_tissue] += 1
        record["samples_by_tissue"][external_tissue].add(row["SampleID"].strip())
        current = record["min_fdr_by_tissue"].get(external_tissue, math.inf)
        record["min_fdr_by_tissue"][external_tissue] = min(current, fdr)
    if max_fdr > 0.05 + 1e-12:
        raise AnalysisError(f"PigGTEx site table is not significant-only: maximum FDR={max_fdr}")
    audit = {
        "n_rows": total_rows, "n_target_matched_rows": matched_rows,
        "n_target_sites_with_any_support": len(result), "maximum_fdr": max_fdr,
        "n_tissues": len(tissues_seen), "external_boundary": (
            "The PigGTEx input is significant-only. Absence is not evidence of a tested negative."
        ),
        "mapped_internal_tissues": dict(sorted(tissue_map.items())),
    }
    return dict(result), audit


def support_signature(row: Mapping[str, object]) -> str:
    labels = ["CANONICAL_SNP_RECURRENCE"] if int(row.get("canonical_recurrent_gene_tissue_snp", 0)) else []
    if int(row.get("n_primary_block_direction_concordant_samples", 0)) > 0:
        labels.append("INDEPENDENT_GATK_DIRECTION_CONCORDANT")
    if int(row.get("independent_gatk_statistically_recurrent", 0)):
        labels.append("INDEPENDENT_GATK_RECURRENCE")
    if int(row.get("piggtex_significant_site_rows_same_mapped_tissue", 0)) > 0:
        labels.append("PIGGTEX_SAME_TISSUE_SITE_ASE")
    if int(row.get("piggtex_exact_significant_eqtl", 0)):
        labels.append("PIGGTEX_EXACT_SIGNIFICANT_EQTL")
    if int(row.get("piggtex_exact_finemapped_eqtl_pip_ge_0_5", 0)):
        labels.append("PIGGTEX_EXACT_FINEMAPPED_EQTL_PIP_GE_0.5")
    if int(row.get("piggtex_exact_finemapped_eqtl_pip_ge_0_9", 0)):
        labels.append("PIGGTEX_EXACT_FINEMAPPED_EQTL_PIP_GE_0.9")
    return ";".join(labels) if labels else "NO_LISTED_SUPPORT"


def is_full_high_pip(row: Mapping[str, object]) -> bool:
    return (
        int(row.get("canonical_recurrent_gene_tissue_snp", 0)) == 1
        and int(row.get("n_primary_block_direction_concordant_samples", 0)) > 0
        and int(row.get("piggtex_significant_site_rows_same_mapped_tissue", 0)) > 0
        and int(row.get("piggtex_exact_significant_eqtl", 0)) == 1
        and int(row.get("piggtex_exact_finemapped_eqtl_pip_ge_0_5", 0)) == 1
    )


def summarize_scope(name: str, rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    return {
        "scope": name,
        "n_sites": len(rows),
        "n_canonical_recurrent_snp": sum(int(row.get("canonical_recurrent_gene_tissue_snp", 0)) for row in rows),
        "n_with_primary_block_gatk_ase": sum(int(row.get("n_primary_block_gatk_ase_samples", 0)) > 0 for row in rows),
        "n_with_primary_block_direction_concordant_gatk_ase": sum(int(row.get("n_primary_block_direction_concordant_samples", 0)) > 0 for row in rows),
        "n_independent_gatk_recurrent": sum(int(row.get("independent_gatk_statistically_recurrent", 0)) for row in rows),
        "n_with_piggtex_significant_site_ase_same_mapped_tissue": sum(int(row.get("piggtex_significant_site_rows_same_mapped_tissue", 0)) > 0 for row in rows),
        "n_with_exact_significant_eqtl": sum(int(row.get("piggtex_exact_significant_eqtl", 0)) for row in rows),
        "n_with_exact_finemapped_eqtl": sum(int(row.get("piggtex_exact_finemapped_eqtl", 0)) for row in rows),
        "n_with_exact_finemapped_eqtl_pip_ge_0_5": sum(int(row.get("piggtex_exact_finemapped_eqtl_pip_ge_0_5", 0)) for row in rows),
        "n_with_exact_finemapped_eqtl_pip_ge_0_9": sum(int(row.get("piggtex_exact_finemapped_eqtl_pip_ge_0_9", 0)) for row in rows),
        "n_full_site_resolved_candidates_pip_ge_0_5": sum(is_full_high_pip(row) for row in rows),
    }


def build_site_row(
    key: Tuple[str, str, str, int, str, str], candidate: Mapping[str, object],
    pair_info: Mapping[Tuple[str, str], Mapping[str, str]], tissue_pairs: Set[Tuple[str, str]],
    canonical_snps: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, str]],
    blocks: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, object]],
    gatk_calls: Mapping[Tuple[str, str, int, str, str], Mapping[str, str]],
    gatk_recurrence: Mapping[Tuple[str, str, int, str, str], Mapping[str, str]],
    piggtex_sites: Mapping[Tuple[str, int, str, str], Mapping[str, object]],
    eqtl: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, object]],
) -> Dict[str, object]:
    gene_id, tissue, chrom, position, allele_1, allele_2 = key
    site = (chrom, position, allele_1, allele_2)
    physical_key = (tissue, *site)
    pair = (gene_id, tissue)
    info = pair_info[pair]
    canonical = canonical_snps.get(key, {})
    block = blocks.get(key, new_block_record())
    call = gatk_calls.get(physical_key, {})
    recurrence = gatk_recurrence.get(physical_key, {})
    external = piggtex_sites.get(site, {})
    external_tissue = str(candidate["piggtex_tissue"])
    same_rows = int(external.get("rows_by_tissue", {}).get(external_tissue, 0))
    same_samples = external.get("samples_by_tissue", {}).get(external_tissue, set())
    same_min = external.get("min_fdr_by_tissue", {}).get(external_tissue, math.inf)
    e = eqtl.get(key, {})
    pips = list(e.get("pips", []))
    max_pip = max(pips) if pips else None
    if key not in eqtl:
        eqtl_status = "IN_SERVER_AUDIT_SITE_UNIVERSE_NO_EXACT_EQTL_MATCH"
    elif int(e.get("finemapped_rows", 0)) > 0:
        eqtl_status = "EXACT_FINEMAPPED_EQTL_MATCH"
    else:
        eqtl_status = "EXACT_SIGNIFICANT_EQTL_MATCH"
    row = {
        "gene_id": gene_id, "gene_name": info["gene_name"], "gene_biotype": info["gene_biotype"],
        "internal_tissue": tissue, "piggtex_tissue": external_tissue,
        "piggtex_gene_ids": sorted_csv(candidate["piggtex_gene_ids"]),
        "variant_unordered_key": site_id(site), "chrom": chrom, "position": position,
        "allele_1": allele_1, "allele_2": allele_2,
        "current_recurrent_gene_tissue_pair": 1, "current_tissue_enriched_pair": int(pair in tissue_pairs),
        "canonical_recurrent_gene_tissue_snp": int(bool(canonical)),
        "canonical_snp_n_testable_animals": canonical.get("n_animals_testable", ""),
        "canonical_snp_animals_testable": canonical.get("animals_testable", ""),
        "canonical_snp_n_positive_animals": canonical.get("n_animals_with_block_concordant_snp_ase", ""),
        "canonical_snp_positive_animals": canonical.get("animals_with_block_concordant_snp_ase", ""),
        "canonical_snp_p": canonical.get("recurrence_p_poisson_binomial", ""),
        "canonical_snp_q": canonical.get("recurrence_q_bh_global", ""),
        "n_primary_block_contexts": len(block["contexts"]),
        "n_primary_block_testable_samples": len(block["testable"]),
        "primary_block_testable_samples": sorted_csv(block["testable"]),
        "n_primary_block_gatk_ase_samples": len(block["ase"]),
        "primary_block_gatk_ase_samples": sorted_csv(block["ase"]),
        "n_primary_block_direction_concordant_samples": len(block["concordant"]),
        "primary_block_direction_concordant_samples": sorted_csv(block["concordant"]),
        "n_primary_block_direction_discordant_samples": len(block["discordant"]),
        "primary_block_direction_discordant_samples": sorted_csv(block["discordant"]),
        "n_primary_block_gatk_testable_not_ase_samples": len(block["testable_not_ase"]),
        "primary_block_gatk_testable_not_ase_samples": sorted_csv(block["testable_not_ase"]),
        "n_primary_block_not_gatk_testable_samples": len(block["not_testable"]),
        "primary_block_not_gatk_testable_samples": sorted_csv(block["not_testable"]),
        "n_independent_gatk_ase_samples": call.get("n_internal_gatk_ase_measurements", "0"),
        "independent_gatk_ase_samples": call.get("internal_gatk_ase_samples", ""),
        "independent_gatk_recurrence_evaluable": int(bool(recurrence)),
        "independent_gatk_n_testable_animals": recurrence.get("n_animals_testable", ""),
        "independent_gatk_n_positive_animals": recurrence.get("n_animals_with_snp_level_ase", ""),
        "independent_gatk_p": recurrence.get("recurrence_p_poisson_binomial", ""),
        "independent_gatk_q": recurrence.get("recurrence_q_bh_global", ""),
        "independent_gatk_statistically_recurrent": parse_bool(
            recurrence.get("is_statistically_recurrent_snp_ase", "0"),
            "is_statistically_recurrent_snp_ase", site_id(site),
        ),
        "piggtex_significant_site_rows_any_tissue": int(external.get("rows_any", 0)),
        "piggtex_significant_site_rows_same_mapped_tissue": same_rows,
        "piggtex_significant_site_samples_same_mapped_tissue": len(same_samples),
        "piggtex_min_fdr_any_tissue": fmt_float(external.get("min_fdr_any", math.inf)),
        "piggtex_min_fdr_same_mapped_tissue": fmt_float(same_min),
        "piggtex_site_ase_support_status": (
            "SIGNIFICANT_ASE_SITE_SAME_MAPPED_TISSUE" if same_rows > 0 else
            "SIGNIFICANT_ASE_SITE_OTHER_TISSUE_ONLY" if int(external.get("rows_any", 0)) > 0 else
            "NOT_OBSERVED_IN_SIGNIFICANT_ONLY_SITE_TABLE"
        ),
        "server_eqtl_audit_site_status": eqtl_status,
        "piggtex_exact_significant_eqtl": int(int(e.get("significant_rows", 0)) > 0),
        "piggtex_exact_significant_eqtl_min_p": fmt_float(min(e.get("p_values", []))) if e.get("p_values") else "",
        "piggtex_exact_significant_eqtl_slopes": sorted_csv(e.get("slopes", set())),
        "piggtex_exact_finemapped_eqtl": int(int(e.get("finemapped_rows", 0)) > 0),
        "piggtex_exact_finemapped_eqtl_max_pip": fmt_float(max_pip),
        "piggtex_exact_finemapped_eqtl_credible_sets": sorted_csv(e.get("credible_sets", set())),
        "piggtex_exact_finemapped_eqtl_pip_ge_0_5": int(max_pip is not None and max_pip >= 0.5),
        "piggtex_exact_finemapped_eqtl_pip_ge_0_9": int(max_pip is not None and max_pip >= 0.9),
    }
    row["support_signature"] = support_signature(row)
    return row


def best_site(rows: Sequence[Mapping[str, object]]) -> str:
    def rank(row: Mapping[str, object]):
        return (
            int(is_full_high_pip(row)),
            int(row["canonical_recurrent_gene_tissue_snp"]),
            int(row["n_primary_block_direction_concordant_samples"]),
            int(row["piggtex_significant_site_rows_same_mapped_tissue"] > 0),
            int(row["piggtex_exact_significant_eqtl"]),
            float(row["piggtex_exact_finemapped_eqtl_max_pip"] or -1),
            -int(row["position"]),
        )
    return str(max(rows, key=rank)["variant_unordered_key"]) if rows else ""


def pair_summary(pair_rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    first = pair_rows[0]
    ids = lambda predicate: sorted_csv(str(row["variant_unordered_key"]) for row in pair_rows if predicate(row))
    canonical_ids = ids(lambda row: int(row["canonical_recurrent_gene_tissue_snp"]) == 1)
    concordant_ids = ids(lambda row: int(row["n_primary_block_direction_concordant_samples"]) > 0)
    same_site_ids = ids(lambda row: int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0)
    eqtl_ids = ids(lambda row: int(row["piggtex_exact_significant_eqtl"]) == 1)
    full_ids = ids(is_full_high_pip)
    if full_ids:
        interpretation = "FULL_SITE_RESOLVED_SUPPORT_WITH_EXACT_PIP_GE_0.5"
    elif canonical_ids and concordant_ids and same_site_ids:
        interpretation = "CANONICAL_GATK_AND_PIGGTEX_SITE_ASE_SUPPORT_WITHOUT_HIGH_PIP_FULL_MATCH"
    elif concordant_ids and same_site_ids:
        interpretation = "GATK_AND_PIGGTEX_SITE_ASE_SUPPORT"
    elif concordant_ids:
        interpretation = "INDEPENDENT_GATK_SUPPORT_ONLY_AT_AUDITED_SITES"
    else:
        interpretation = "NO_DIRECTION_CONCORDANT_GATK_SITE_SUPPORT_IN_PRIMARY_BLOCK_CONTEXTS"
    return {
        "gene_id": first["gene_id"], "gene_name": first["gene_name"],
        "target_tissue": first["internal_tissue"], "piggtex_tissue": first["piggtex_tissue"],
        "n_server_audit_candidate_sites": len(pair_rows),
        "n_canonical_recurrent_snp_sites": sum(int(row["canonical_recurrent_gene_tissue_snp"]) for row in pair_rows),
        "n_sites_with_primary_block_gatk_ase": sum(int(row["n_primary_block_gatk_ase_samples"]) > 0 for row in pair_rows),
        "n_sites_with_primary_block_direction_concordant_gatk_ase": sum(int(row["n_primary_block_direction_concordant_samples"]) > 0 for row in pair_rows),
        "n_independent_gatk_recurrent_sites": sum(int(row["independent_gatk_statistically_recurrent"]) for row in pair_rows),
        "n_sites_with_piggtex_significant_ase_same_mapped_tissue": sum(int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0 for row in pair_rows),
        "n_sites_with_exact_significant_eqtl": sum(int(row["piggtex_exact_significant_eqtl"]) for row in pair_rows),
        "n_sites_with_exact_finemapped_eqtl": sum(int(row["piggtex_exact_finemapped_eqtl"]) for row in pair_rows),
        "n_sites_with_exact_finemapped_eqtl_pip_ge_0_5": sum(int(row["piggtex_exact_finemapped_eqtl_pip_ge_0_5"]) for row in pair_rows),
        "n_sites_with_exact_finemapped_eqtl_pip_ge_0_9": sum(int(row["piggtex_exact_finemapped_eqtl_pip_ge_0_9"]) for row in pair_rows),
        "n_full_site_resolved_candidates_pip_ge_0_5": sum(is_full_high_pip(row) for row in pair_rows),
        "canonical_recurrent_snp_site_ids": canonical_ids,
        "direction_concordant_gatk_site_ids": concordant_ids,
        "same_tissue_piggtex_site_ase_ids": same_site_ids,
        "exact_significant_eqtl_site_ids": eqtl_ids,
        "high_pip_full_site_resolved_candidate_ids": full_ids,
        "best_supported_site_id": best_site(pair_rows),
        "pair_site_evidence_interpretation": interpretation,
    }


def canonical_audit_row(
    canonical: Mapping[str, str], current_pairs: Set[Tuple[str, str]], tissue_pairs: Set[Tuple[str, str]],
    site_rows: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, object]],
    blocks: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, object]],
    gatk_calls: Mapping[Tuple[str, str, int, str, str], Mapping[str, str]],
    gatk_recurrence: Mapping[Tuple[str, str, int, str, str], Mapping[str, str]],
    piggtex_sites: Mapping[Tuple[str, int, str, str], Mapping[str, object]],
    tissue_map: Mapping[str, str], eqtl: Mapping[Tuple[str, str, str, int, str, str], Mapping[str, object]],
) -> Dict[str, object]:
    tissue = normalize_tissue(canonical["tissue"])
    site = unordered_site(canonical["chrom"], int(canonical["position"]), canonical["ref"], canonical["alt"])
    key = (canonical["gene_id"], tissue, *site)
    if key in site_rows:
        source = site_rows[key]
        return {
            "gene_id": canonical["gene_id"], "gene_name": canonical["gene_name"], "tissue": tissue,
            "variant_id": canonical["variant_id"], "chrom": canonical["chrom"], "position": canonical["position"],
            "ref": canonical["ref"], "alt": canonical["alt"],
            "n_animals_testable": canonical["n_animals_testable"], "animals_testable": canonical["animals_testable"],
            "n_animals_with_block_concordant_snp_ase": canonical["n_animals_with_block_concordant_snp_ase"],
            "animals_with_block_concordant_snp_ase": canonical["animals_with_block_concordant_snp_ase"],
            "recurrence_p_poisson_binomial": canonical["recurrence_p_poisson_binomial"],
            "recurrence_q_bh_global": canonical["recurrence_q_bh_global"],
            "pair_in_current_646_recurrence": 1, "pair_in_current_39_tissue_enrichment": source["current_tissue_enriched_pair"],
            "n_primary_block_gatk_ase_samples": source["n_primary_block_gatk_ase_samples"],
            "primary_block_gatk_ase_samples": source["primary_block_gatk_ase_samples"],
            "n_primary_block_direction_concordant_samples": source["n_primary_block_direction_concordant_samples"],
            "primary_block_direction_concordant_samples": source["primary_block_direction_concordant_samples"],
            "n_primary_block_direction_discordant_samples": source["n_primary_block_direction_discordant_samples"],
            "primary_block_direction_discordant_samples": source["primary_block_direction_discordant_samples"],
            "n_independent_gatk_ase_samples": source["n_independent_gatk_ase_samples"],
            "independent_gatk_ase_samples": source["independent_gatk_ase_samples"],
            "independent_gatk_recurrence_evaluable": source["independent_gatk_recurrence_evaluable"],
            "independent_gatk_n_testable_animals": source["independent_gatk_n_testable_animals"],
            "independent_gatk_n_positive_animals": source["independent_gatk_n_positive_animals"],
            "independent_gatk_p": source["independent_gatk_p"], "independent_gatk_q": source["independent_gatk_q"],
            "independent_gatk_statistically_recurrent": source["independent_gatk_statistically_recurrent"],
            "piggtex_mapped_tissue": source["piggtex_tissue"],
            "piggtex_significant_site_rows_same_mapped_tissue": source["piggtex_significant_site_rows_same_mapped_tissue"],
            "piggtex_min_fdr_same_mapped_tissue": source["piggtex_min_fdr_same_mapped_tissue"],
            "piggtex_exact_significant_eqtl": source["piggtex_exact_significant_eqtl"],
            "piggtex_exact_finemapped_eqtl_max_pip": source["piggtex_exact_finemapped_eqtl_max_pip"],
            "server_eqtl_audit_site_status": source["server_eqtl_audit_site_status"],
            "support_signature": source["support_signature"],
        }
    block = blocks.get(key, new_block_record())
    physical_key = (tissue, *site)
    call = gatk_calls.get(physical_key, {})
    recurrence = gatk_recurrence.get(physical_key, {})
    external_tissue = tissue_map[tissue]
    external = piggtex_sites.get(site, {})
    same_rows = int(external.get("rows_by_tissue", {}).get(external_tissue, 0))
    same_min = external.get("min_fdr_by_tissue", {}).get(external_tissue, math.inf)
    pair_in_646 = int((canonical["gene_id"], tissue) in current_pairs)
    row_for_signature = {
        "canonical_recurrent_gene_tissue_snp": 1,
        "n_primary_block_direction_concordant_samples": len(block["concordant"]),
        "independent_gatk_statistically_recurrent": parse_bool(recurrence.get("is_statistically_recurrent_snp_ase", "0"), "is_statistically_recurrent_snp_ase", canonical["variant_id"]),
        "piggtex_significant_site_rows_same_mapped_tissue": same_rows,
        "piggtex_exact_significant_eqtl": int(int(eqtl.get(key, {}).get("significant_rows", 0)) > 0),
        "piggtex_exact_finemapped_eqtl_pip_ge_0_5": 0,
        "piggtex_exact_finemapped_eqtl_pip_ge_0_9": 0,
    }
    return {
        "gene_id": canonical["gene_id"], "gene_name": canonical["gene_name"], "tissue": tissue,
        "variant_id": canonical["variant_id"], "chrom": canonical["chrom"], "position": canonical["position"],
        "ref": canonical["ref"], "alt": canonical["alt"], "n_animals_testable": canonical["n_animals_testable"],
        "animals_testable": canonical["animals_testable"],
        "n_animals_with_block_concordant_snp_ase": canonical["n_animals_with_block_concordant_snp_ase"],
        "animals_with_block_concordant_snp_ase": canonical["animals_with_block_concordant_snp_ase"],
        "recurrence_p_poisson_binomial": canonical["recurrence_p_poisson_binomial"],
        "recurrence_q_bh_global": canonical["recurrence_q_bh_global"],
        "pair_in_current_646_recurrence": pair_in_646,
        "pair_in_current_39_tissue_enrichment": int((canonical["gene_id"], tissue) in tissue_pairs),
        "n_primary_block_gatk_ase_samples": len(block["ase"]), "primary_block_gatk_ase_samples": sorted_csv(block["ase"]),
        "n_primary_block_direction_concordant_samples": len(block["concordant"]),
        "primary_block_direction_concordant_samples": sorted_csv(block["concordant"]),
        "n_primary_block_direction_discordant_samples": len(block["discordant"]),
        "primary_block_direction_discordant_samples": sorted_csv(block["discordant"]),
        "n_independent_gatk_ase_samples": call.get("n_internal_gatk_ase_measurements", "0"),
        "independent_gatk_ase_samples": call.get("internal_gatk_ase_samples", ""),
        "independent_gatk_recurrence_evaluable": int(bool(recurrence)),
        "independent_gatk_n_testable_animals": recurrence.get("n_animals_testable", ""),
        "independent_gatk_n_positive_animals": recurrence.get("n_animals_with_snp_level_ase", ""),
        "independent_gatk_p": recurrence.get("recurrence_p_poisson_binomial", ""),
        "independent_gatk_q": recurrence.get("recurrence_q_bh_global", ""),
        "independent_gatk_statistically_recurrent": row_for_signature["independent_gatk_statistically_recurrent"],
        "piggtex_mapped_tissue": external_tissue,
        "piggtex_significant_site_rows_same_mapped_tissue": same_rows,
        "piggtex_min_fdr_same_mapped_tissue": fmt_float(same_min),
        "piggtex_exact_significant_eqtl": row_for_signature["piggtex_exact_significant_eqtl"],
        "piggtex_exact_finemapped_eqtl_max_pip": "",
        "server_eqtl_audit_site_status": (
            "PAIR_OUTSIDE_SERVER_AUDIT_CANDIDATE_UNIVERSE" if not pair_in_646
            else "SITE_OUTSIDE_SERVER_AUDIT_CANDIDATE_UNIVERSE"
        ),
        "support_signature": support_signature(row_for_signature),
    }


def run(args: argparse.Namespace) -> Dict[str, object]:
    paths = {name: Path(getattr(args, name)).resolve() for name in [
        "current_recurrent_pairs", "tissue_enriched_pairs", "canonical_recurrent_snps",
        "server_candidate_universe", "block_validation", "gatk_site_support",
        "gatk_recurrence_all", "piggtex_sites", "exact_eqtl_matches", "server_audit_metadata",
    ]}
    for name, path in paths.items():
        if not path.is_file():
            raise AnalysisError(f"missing input {name}: {path}")
    server_metadata = json.loads(paths["server_audit_metadata"].read_text(encoding="utf-8"))
    output = Path(args.out_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)

    pair_info = {}
    for row in iter_tsv(paths["current_recurrent_pairs"], PAIR_REQUIRED):
        key = (row["gene_id"].strip(), normalize_tissue(row["tissue"]))
        if key in pair_info:
            raise AnalysisError(f"duplicate current recurrent pair: {key}")
        pair_info[key] = row
    if len(pair_info) != 646:
        raise AnalysisError(f"expected 646 current recurrent pairs, observed {len(pair_info)}")

    tissue_pair_info = {}
    for row in iter_tsv(paths["tissue_enriched_pairs"], TISSUE_PAIR_REQUIRED):
        if parse_bool(row["pair_level_primary"], "pair_level_primary", row["gene_id"]) != 1:
            raise AnalysisError("current tissue-enriched pair file contains a non-primary row")
        key = (row["gene_id"].strip(), normalize_tissue(row["target_tissue"]))
        if key in tissue_pair_info:
            raise AnalysisError(f"duplicate tissue-enriched pair: {key}")
        tissue_pair_info[key] = row
    tissue_pairs = set(tissue_pair_info)
    if len(tissue_pairs) != 39 or not tissue_pairs <= set(pair_info):
        raise AnalysisError("current 39-pair set is not an exact subset of the current 646-pair recurrence set")

    candidates = {}
    candidate_pairs = set()
    tissue_map = {}
    for row in iter_tsv(paths["server_candidate_universe"], CANDIDATE_REQUIRED):
        tissue = normalize_tissue(row["internal_tissue"])
        site = unordered_site(row["chrom"], int(row["position"]), row["allele_1"], row["allele_2"])
        key = (row["internal_gene_id"].strip(), tissue, *site)
        pair = key[:2]
        candidate_pairs.add(pair)
        mapped = row["piggtex_tissue"].strip()
        if tissue in tissue_map and tissue_map[tissue] != mapped:
            raise AnalysisError(f"inconsistent tissue map for {tissue}")
        tissue_map[tissue] = mapped
        record = candidates.setdefault(key, {
            "piggtex_gene_ids": set(), "piggtex_tissue": mapped,
            "variant_unordered_keys": set(),
        })
        if record["piggtex_tissue"] != mapped:
            raise AnalysisError(f"inconsistent candidate mapped tissue for {key}")
        record["piggtex_gene_ids"].add(row["piggtex_gene_id"].strip())
        record["variant_unordered_keys"].add(row["internal_variant_unordered_key"].strip())
    if candidate_pairs != set(pair_info):
        raise AnalysisError(
            f"server candidate pair universe mismatch: missing={len(set(pair_info)-candidate_pairs)}, "
            f"extra={len(candidate_pairs-set(pair_info))}"
        )
    candidate_external_keys = {
        (piggtex_gene_id, str(record["piggtex_tissue"]), *key[2:])
        for key, record in candidates.items()
        for piggtex_gene_id in record["piggtex_gene_ids"]
    }
    expected_external_keys = int(server_metadata.get("n_candidate_exact_gene_tissue_variant_keys", -1))
    if len(candidate_external_keys) != expected_external_keys:
        raise AnalysisError(
            "local server-candidate external-key count does not reproduce server metadata: "
            f"local={len(candidate_external_keys)}, metadata={expected_external_keys}"
        )

    canonical_rows = []
    canonical_snps = {}
    for row in iter_tsv(paths["canonical_recurrent_snps"], CANONICAL_SNP_REQUIRED):
        if parse_bool(row["is_statistically_recurrent_snp_ase"], "is_statistically_recurrent_snp_ase", row["variant_id"]) != 1:
            raise AnalysisError("canonical recurrent SNP input contains a non-recurrent row")
        tissue = normalize_tissue(row["tissue"])
        site = unordered_site(row["chrom"], int(row["position"]), row["ref"], row["alt"])
        key = (row["gene_id"].strip(), tissue, *site)
        if key in canonical_snps:
            raise AnalysisError(f"duplicate canonical recurrent SNP hypothesis: {key}")
        canonical_snps[key] = row
        canonical_rows.append(row)
    if len(canonical_rows) != 1619:
        raise AnalysisError(f"expected 1,619 canonical recurrent SNP hypotheses, observed {len(canonical_rows)}")

    blocks = aggregate_block_rows(paths["block_validation"])
    gatk_calls = load_gatk_site_calls(paths["gatk_site_support"])
    gatk_recurrence = load_gatk_recurrence(paths["gatk_recurrence_all"])
    eqtl = aggregate_eqtl(paths["exact_eqtl_matches"])
    if any(key not in candidates for key in eqtl):
        raise AnalysisError("exact eQTL match output contains a key outside the server candidate universe")

    target_sites = {key[2:] for key in candidates} | {key[2:] for key in canonical_snps}
    piggtex_site_support, piggtex_audit = audit_piggtex_sites(paths["piggtex_sites"], target_sites, tissue_map)

    site_rows = []
    site_rows_by_key = {}
    for key in sorted(candidates, key=lambda value: (value[1], value[0], int(value[2]) if value[2].isdigit() else 999, value[3], value[4], value[5])):
        row = build_site_row(
            key, candidates[key], pair_info, tissue_pairs, canonical_snps, blocks, gatk_calls,
            gatk_recurrence, piggtex_site_support, eqtl,
        )
        site_rows.append(row)
        site_rows_by_key[key] = row

    current_39_rows = [row for row in site_rows if int(row["current_tissue_enriched_pair"]) == 1]
    pair_rows_by_pair = defaultdict(list)
    for row in current_39_rows:
        pair_rows_by_pair[(row["gene_id"], row["internal_tissue"])].append(row)
    if set(pair_rows_by_pair) != tissue_pairs:
        raise AnalysisError("not all current 39 tissue-enriched pairs have server-audit candidate sites")
    pair_summaries = [pair_summary(pair_rows_by_pair[pair]) for pair in sorted(tissue_pairs, key=lambda value: (value[1], value[0]))]

    canonical_audit = [
        canonical_audit_row(
            row, set(pair_info), tissue_pairs, site_rows_by_key, blocks, gatk_calls,
            gatk_recurrence, piggtex_site_support, tissue_map, eqtl,
        )
        for row in canonical_rows
    ]
    canonical_audit.sort(key=lambda row: (row["tissue"], row["gene_id"], int(row["chrom"]) if str(row["chrom"]).isdigit() else 999, int(row["position"])))

    high_candidates = [row for row in site_rows if is_full_high_pip(row)]
    high_candidates.sort(key=lambda row: (-float(row["piggtex_exact_finemapped_eqtl_max_pip"]), row["internal_tissue"], row["gene_id"], int(row["position"])))
    thyn1_rows = [row for row in site_rows if row["gene_name"] == "THYN1"]

    stage_rows = [
        summarize_scope("current_646_pair_server_audit_site_universe", site_rows),
        summarize_scope("current_39_tissue_enriched_pair_site_universe", current_39_rows),
        summarize_scope("current_646_pairs_canonical_recurrent_SNP_subset", [row for row in site_rows if int(row["canonical_recurrent_gene_tissue_snp"]) == 1]),
    ]
    canonical_stage_rows = [{
        "scope": "all_1619_canonical_recurrent_gene_tissue_SNP_hypotheses",
        "n_sites": len(canonical_audit),
        "n_canonical_recurrent_snp": len(canonical_audit),
        "n_with_primary_block_gatk_ase": sum(int(row["n_primary_block_gatk_ase_samples"]) > 0 for row in canonical_audit),
        "n_with_primary_block_direction_concordant_gatk_ase": sum(int(row["n_primary_block_direction_concordant_samples"]) > 0 for row in canonical_audit),
        "n_independent_gatk_recurrent": sum(int(row["independent_gatk_statistically_recurrent"]) for row in canonical_audit),
        "n_with_piggtex_significant_site_ase_same_mapped_tissue": sum(int(row["piggtex_significant_site_rows_same_mapped_tissue"]) > 0 for row in canonical_audit),
        "n_with_exact_significant_eqtl": sum(int(row["piggtex_exact_significant_eqtl"]) for row in canonical_audit),
        "n_with_exact_finemapped_eqtl": sum(bool(row["piggtex_exact_finemapped_eqtl_max_pip"]) for row in canonical_audit),
        "n_with_exact_finemapped_eqtl_pip_ge_0_5": sum(float(row["piggtex_exact_finemapped_eqtl_max_pip"] or -1) >= 0.5 for row in canonical_audit),
        "n_with_exact_finemapped_eqtl_pip_ge_0_9": sum(float(row["piggtex_exact_finemapped_eqtl_max_pip"] or -1) >= 0.9 for row in canonical_audit),
        "n_full_site_resolved_candidates_pip_ge_0_5": len(high_candidates),
    }]
    stage_rows.extend(canonical_stage_rows)

    outputs = {
        "current_646_pair_site_evidence_audit.tsv.gz": (SITE_HEADER, site_rows),
        "current_39_pair_site_evidence_audit.tsv": (SITE_HEADER, current_39_rows),
        "current_39_pair_site_evidence_summary.tsv": (PAIR_SUMMARY_HEADER, pair_summaries),
        "canonical_recurrent_snp_gatk_piggtex_audit.tsv": (CANONICAL_AUDIT_HEADER, canonical_audit),
        "high_pip_full_site_resolved_candidates.tsv": (SITE_HEADER, high_candidates),
        "thyn1_site_resolved_audit.tsv": (SITE_HEADER, thyn1_rows),
        "evidence_stage_summary.tsv": (STAGE_HEADER, stage_rows),
    }
    output_counts = {name: write_tsv(output / name, header, rows) for name, (header, rows) in outputs.items()}

    metadata = {
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS_LOCAL_REJOIN_PENDING_CANONICAL_REVIEW",
        "exact_site_key": "Sscrofa11.1 contig + 1-based position + unordered biallelic allele pair",
        "canonical_boundaries": {
            "gene_level_recurrence_pairs": 646,
            "tissue_enriched_pairs": 39,
            "canonical_recurrent_gene_tissue_snp_hypotheses": 1619,
            "independent_GATK_recurrence_is_validation_not_replacement": True,
            "retired_eight_pair_internal_statistics_used": False,
        },
        "external_provenance_boundary": {
            "PigGTEx_site_ASE": piggtex_audit["external_boundary"],
            "PigGTEx_gene_eQTL": (
                "Exact significant and fine-mapped eQTL rows were locally rejoined from the existing server audit. "
                "The raw server archives were not independently rescanned in this local run."
            ),
            "server_variant_matching": server_metadata.get("variant_matching", ""),
        },
        "counts": {
            "n_current_recurrent_pairs": len(pair_info),
            "n_current_tissue_enriched_pairs": len(tissue_pairs),
            "n_server_candidate_rows_before_collapse": sum(1 for _ in iter_tsv(paths["server_candidate_universe"], CANDIDATE_REQUIRED)),
            "n_server_candidate_unique_internal_sites": len(candidates),
            "n_server_candidate_unique_external_gene_tissue_site_keys": len(candidate_external_keys),
            "n_canonical_recurrent_snp_hypotheses": len(canonical_rows),
            "n_canonical_recurrent_snp_hypotheses_in_current_646_pairs": sum(key[:2] in pair_info for key in canonical_snps),
            "n_high_pip_full_site_resolved_candidates": len(high_candidates),
            "n_THYN1_site_rows_in_current_646_pair_audit": len(thyn1_rows),
        },
        "piggtex_site_input_audit": piggtex_audit,
        "inputs": {name: file_metadata(path) for name, path in paths.items()},
        "outputs": {},
    }
    metadata_path = output / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata["outputs"] = {
        name: {"rows": output_counts[name], **file_metadata(output / name)} for name in outputs
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-recurrent-pairs", required=True)
    parser.add_argument("--tissue-enriched-pairs", required=True)
    parser.add_argument("--canonical-recurrent-snps", required=True)
    parser.add_argument("--server-candidate-universe", required=True)
    parser.add_argument("--block-validation", required=True)
    parser.add_argument("--gatk-site-support", required=True)
    parser.add_argument("--gatk-recurrence-all", required=True)
    parser.add_argument("--piggtex-sites", required=True)
    parser.add_argument("--exact-eqtl-matches", required=True)
    parser.add_argument("--server-audit-metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    try:
        metadata = run(args)
    except (AnalysisError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    print(f"[{metadata['status']}] {Path(args.out_dir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
