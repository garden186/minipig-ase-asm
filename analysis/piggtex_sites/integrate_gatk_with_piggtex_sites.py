#!/usr/bin/env python3
"""Integrate independent minipig GATK SNP-level ASE with PigGTEx sites."""

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


VERSION = "ase_gatk_piggtex_site_integration_v0.1.0"

GATK_REQUIRED = {
    "sample", "tissue", "contig", "position", "ref_allele", "alt_allele",
    "major_allele", "q_bh_within_animal_tissue", "is_snp_level_ase",
}
RECURRENCE_REQUIRED = {
    "tissue", "contig", "position", "ref_allele", "alt_allele",
    "n_animals_testable", "n_animals_with_snp_level_ase",
    "recurrence_q_bh_global", "is_statistically_recurrent_snp_ase",
}
BLOCK_REQUIRED = {
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id",
    "variant_id", "contig", "position", "ref_allele", "alt_allele",
    "gatk_snp_level_ase", "gatk_block_direction_status",
}
TISSUE_PAIR_REQUIRED = {"gene_id", "gene_name", "target_tissue", "pair_level_primary"}
TISSUE_MAP_REQUIRED = {"internal_tissue", "piggtex_tissue"}
PIGGTEX_REQUIRED = {
    "chr", "position", "refAllele", "altAllele", "refCount", "altCount",
    "FDR", "SampleID", "Breed", "Tissue", "pCADD",
}

OUTPUT_HEADER = [
    "internal_tissue", "piggtex_mapped_tissue", "physical_snp_id", "contig",
    "position", "ref_allele", "alt_allele", "allele_1", "allele_2",
    "n_internal_gatk_ase_measurements", "internal_gatk_ase_samples",
    "n_internal_major_ref", "n_internal_major_alt", "internal_major_allele_pattern",
    "is_statistically_recurrent_gatk_snp_ase", "gatk_recurrence_n_testable_animals",
    "gatk_recurrence_n_ase_positive_animals", "gatk_recurrence_q_bh_global",
    "n_primary_block_gene_contexts", "n_primary_block_direction_concordant_contexts",
    "n_primary_block_direction_discordant_contexts", "primary_block_samples",
    "primary_block_gene_ids", "primary_block_gene_names", "primary_block_ids",
    "tissue_enriched_gene_tissue_pairs", "piggtex_significant_rows_any_tissue",
    "piggtex_significant_rows_same_mapped_tissue", "piggtex_samples_any_tissue",
    "piggtex_samples_same_mapped_tissue", "piggtex_breeds_any_tissue",
    "piggtex_breeds_same_mapped_tissue", "piggtex_tissues_any",
    "piggtex_min_fdr_any_tissue", "piggtex_min_fdr_same_mapped_tissue",
    "piggtex_max_pcadd", "piggtex_same_tissue_major_internal_ref_rows",
    "piggtex_same_tissue_major_internal_alt_rows", "piggtex_same_tissue_tie_rows",
    "piggtex_same_tissue_major_allele_pattern", "piggtex_site_support_status",
]

STAGE_HEADER = [
    "stage", "n_internal_tissue_snp_sites", "n_distinct_physical_snps",
    "n_with_PigGTEx_significant_site_ASE_any_tissue",
    "n_with_PigGTEx_significant_site_ASE_same_mapped_tissue",
    "n_not_observed_in_PigGTEx_significant_only_site_table",
]

PAIR_HEADER = [
    "gene_id", "gene_name", "target_tissue", "piggtex_mapped_tissue",
    "n_gatk_ase_sites_in_primary_blocks", "n_direction_concordant_gatk_ase_sites",
    "n_direction_discordant_gatk_ase_sites", "n_recurrent_direction_concordant_sites",
    "n_direction_concordant_sites_with_PigGTEx_significant_ASE_any_tissue",
    "n_direction_concordant_sites_with_PigGTEx_significant_ASE_same_mapped_tissue",
    "n_minipig_animals_with_direction_concordant_site",
    "direction_concordant_site_ids", "same_tissue_PigGTEx_supported_site_ids",
    "PigGTEx_site_evidence_interpretation",
]


class AnalysisError(RuntimeError):
    """Raised when an input or exact-site invariant is violated."""


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


def ordered_key(contig: str, position: int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    return contig.removeprefix("chr"), int(position), ref.upper(), alt.upper()


def unordered_key(contig: str, position: int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    allele_1, allele_2 = sorted((ref.upper(), alt.upper()))
    return contig.removeprefix("chr"), int(position), allele_1, allele_2


def physical_snp_id(contig: str, position: int, ref: str, alt: str) -> str:
    return f"{contig}_{position}_{ref}_{alt}"


def fmt_float(value: object) -> str:
    if value is None or value == "" or not math.isfinite(float(value)):
        return ""
    return format(float(value), ".12g")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_metadata(path: Path) -> Dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def new_record(tissue: str, mapped_tissue: str, contig: str, position: int, ref: str, alt: str) -> Dict[str, object]:
    key = unordered_key(contig, position, ref, alt)
    return {
        "internal_tissue": tissue,
        "mapped_tissue": mapped_tissue,
        "contig": contig,
        "position": position,
        "ref": ref,
        "alt": alt,
        "allele_1": key[2],
        "allele_2": key[3],
        "internal_samples": set(),
        "n_internal_major_ref": 0,
        "n_internal_major_alt": 0,
        "recurrence": None,
        "block_contexts": set(),
        "block_concordant_contexts": set(),
        "block_discordant_contexts": set(),
        "block_samples": set(),
        "block_gene_ids": set(),
        "block_gene_names": set(),
        "block_ids": set(),
        "tissue_enriched_pairs": set(),
        "piggtex_rows_any": 0,
        "piggtex_rows_same": 0,
        "piggtex_samples_any": set(),
        "piggtex_samples_same": set(),
        "piggtex_breeds_any": set(),
        "piggtex_breeds_same": set(),
        "piggtex_tissues_any": set(),
        "piggtex_min_fdr_any": math.inf,
        "piggtex_min_fdr_same": math.inf,
        "piggtex_pcadd": set(),
        "piggtex_same_major_ref": 0,
        "piggtex_same_major_alt": 0,
        "piggtex_same_tie": 0,
    }


def load_tissue_map(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for row in iter_tsv(path, TISSUE_MAP_REQUIRED):
        internal = normalize_tissue(row["internal_tissue"])
        if internal in mapping:
            raise AnalysisError(f"duplicate tissue mapping: {internal}")
        mapping[internal] = row["piggtex_tissue"].strip()
    return mapping


def load_gatk_calls(
    path: Path, tissue_map: Mapping[str, str]
) -> Dict[Tuple[str, str, int, str, str], Dict[str, object]]:
    records: Dict[Tuple[str, str, int, str, str], Dict[str, object]] = {}
    seen_measurements = set()
    for row in iter_tsv(path, GATK_REQUIRED):
        if parse_bool(row["is_snp_level_ase"], "is_snp_level_ase", row["sample"]) != 1:
            raise AnalysisError("GATK calls input contains a non-call row")
        sample = row["sample"].strip()
        tissue = normalize_tissue(row["tissue"])
        if tissue not in tissue_map:
            raise AnalysisError(f"unmapped internal tissue: {tissue}")
        contig, position, ref, alt = ordered_key(
            row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"]
        )
        measurement_key = (sample, tissue, contig, position, ref, alt)
        if measurement_key in seen_measurements:
            raise AnalysisError(f"duplicate GATK ASE measurement: {measurement_key}")
        seen_measurements.add(measurement_key)
        key = (tissue, *unordered_key(contig, position, ref, alt))
        record = records.setdefault(
            key, new_record(tissue, tissue_map[tissue], contig, position, ref, alt)
        )
        if (record["ref"], record["alt"]) != (ref, alt):
            raise AnalysisError(f"inconsistent internal REF/ALT orientation: {key}")
        record["internal_samples"].add(sample)
        major = row["major_allele"].strip().upper()
        if major == ref:
            record["n_internal_major_ref"] += 1
        elif major == alt:
            record["n_internal_major_alt"] += 1
        else:
            raise AnalysisError(f"internal major allele outside REF/ALT: {measurement_key}")
    return records


def attach_recurrence(
    path: Path, records: MutableMapping[Tuple[str, str, int, str, str], Dict[str, object]]
) -> int:
    n = 0
    for row in iter_tsv(path, RECURRENCE_REQUIRED):
        if parse_bool(
            row["is_statistically_recurrent_snp_ase"],
            "is_statistically_recurrent_snp_ase",
            row["tissue"],
        ) != 1:
            raise AnalysisError("recurrent input contains a non-recurrent row")
        tissue = normalize_tissue(row["tissue"])
        key = (tissue, *unordered_key(row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"]))
        if key not in records:
            raise AnalysisError(f"recurrent SNP absent from GATK call universe: {key}")
        if records[key]["recurrence"] is not None:
            raise AnalysisError(f"duplicate recurrent tissue-SNP: {key}")
        records[key]["recurrence"] = dict(row)
        n += 1
    return n


def load_tissue_enriched_pairs(path: Path) -> Dict[Tuple[str, str], str]:
    pairs: Dict[Tuple[str, str], str] = {}
    for row in iter_tsv(path, TISSUE_PAIR_REQUIRED):
        if parse_bool(row["pair_level_primary"], "pair_level_primary", row["gene_id"]) != 1:
            continue
        key = (row["gene_id"].strip(), normalize_tissue(row["target_tissue"]))
        pairs[key] = row["gene_name"].strip()
    return pairs


def attach_block_contexts(
    path: Path,
    records: MutableMapping[Tuple[str, str, int, str, str], Dict[str, object]],
    enriched_pairs: Mapping[Tuple[str, str], str],
) -> Dict[str, int]:
    counts = Counter()
    seen = set()
    for row in iter_tsv(path, BLOCK_REQUIRED):
        if parse_bool(row["gatk_snp_level_ase"], "gatk_snp_level_ase", row["variant_id"]) != 1:
            continue
        tissue = normalize_tissue(row["tissue"])
        key = (tissue, *unordered_key(row["contig"], int(row["position"]), row["ref_allele"], row["alt_allele"]))
        if key not in records:
            raise AnalysisError(f"block-linked GATK call absent from GATK call universe: {key}")
        context_key = (
            row["sample"], tissue, row["gene_id"], row["haplotype_block_id"], row["variant_id"]
        )
        if context_key in seen:
            raise AnalysisError(f"duplicate block context: {context_key}")
        seen.add(context_key)
        record = records[key]
        status = row["gatk_block_direction_status"]
        record["block_contexts"].add(context_key)
        if status == "DIRECTION_CONCORDANT":
            record["block_concordant_contexts"].add(context_key)
            counts["n_direction_concordant_contexts"] += 1
        elif status == "DIRECTION_DISCORDANT":
            record["block_discordant_contexts"].add(context_key)
            counts["n_direction_discordant_contexts"] += 1
        else:
            raise AnalysisError(f"unexpected block direction status for GATK call: {status}")
        record["block_samples"].add(row["sample"])
        record["block_gene_ids"].add(row["gene_id"])
        record["block_gene_names"].add(row["gene_name"])
        record["block_ids"].add(row["haplotype_block_id"])
        pair = (row["gene_id"], tissue)
        if pair in enriched_pairs:
            record["tissue_enriched_pairs"].add(pair)
        counts["n_primary_block_gatk_ase_contexts"] += 1
    return dict(counts)


def stream_piggtex(
    path: Path,
    records: MutableMapping[Tuple[str, str, int, str, str], Dict[str, object]],
) -> Dict[str, object]:
    physical_lookup: Dict[Tuple[str, int, str, str], List[Dict[str, object]]] = defaultdict(list)
    for record in records.values():
        key = unordered_key(
            str(record["contig"]), int(record["position"]), str(record["ref"]), str(record["alt"])
        )
        physical_lookup[key].append(record)
    stats: Dict[str, object] = {
        "n_rows": 0,
        "n_distinct_sites_sorted_key": 0,
        "n_malformed_rows": 0,
        "minimum_fdr": math.inf,
        "maximum_fdr": -math.inf,
        "samples": set(),
        "breeds": set(),
        "tissues": set(),
    }
    previous = None
    for row in iter_tsv(path, PIGGTEX_REQUIRED):
        stats["n_rows"] += 1
        try:
            chrom = row["chr"].removeprefix("chr")
            position = int(row["position"])
            pig_ref, pig_alt = row["refAllele"].upper(), row["altAllele"].upper()
            fdr = float(row["FDR"])
            ref_count, alt_count = int(row["refCount"]), int(row["altCount"])
        except (ValueError, KeyError):
            stats["n_malformed_rows"] += 1
            continue
        key = unordered_key(chrom, position, pig_ref, pig_alt)
        if key != previous:
            stats["n_distinct_sites_sorted_key"] += 1
            previous = key
        stats["minimum_fdr"] = min(stats["minimum_fdr"], fdr)
        stats["maximum_fdr"] = max(stats["maximum_fdr"], fdr)
        stats["samples"].add(row["SampleID"])
        stats["breeds"].add(row["Breed"])
        stats["tissues"].add(row["Tissue"])
        for record in physical_lookup.get(key, []):
            record["piggtex_rows_any"] += 1
            record["piggtex_samples_any"].add(row["SampleID"])
            record["piggtex_breeds_any"].add(row["Breed"])
            record["piggtex_tissues_any"].add(row["Tissue"])
            record["piggtex_min_fdr_any"] = min(record["piggtex_min_fdr_any"], fdr)
            try:
                record["piggtex_pcadd"].add(float(row["pCADD"]))
            except ValueError:
                pass
            if row["Tissue"] != record["mapped_tissue"]:
                continue
            record["piggtex_rows_same"] += 1
            record["piggtex_samples_same"].add(row["SampleID"])
            record["piggtex_breeds_same"].add(row["Breed"])
            record["piggtex_min_fdr_same"] = min(record["piggtex_min_fdr_same"], fdr)
            if ref_count == alt_count:
                record["piggtex_same_tie"] += 1
            else:
                pig_major = pig_ref if ref_count > alt_count else pig_alt
                if pig_major == record["ref"]:
                    record["piggtex_same_major_ref"] += 1
                elif pig_major == record["alt"]:
                    record["piggtex_same_major_alt"] += 1
                else:
                    raise AnalysisError(f"PigGTEx major allele outside internal alleles: {key}")
    if stats["n_malformed_rows"]:
        raise AnalysisError(f"PigGTEx site table contains malformed rows: {stats['n_malformed_rows']}")
    if float(stats["maximum_fdr"]) >= 0.05:
        raise AnalysisError(f"PigGTEx site input is not significant-only: max FDR={stats['maximum_fdr']}")
    return {
        "n_rows": stats["n_rows"],
        "n_distinct_sites_sorted_key": stats["n_distinct_sites_sorted_key"],
        "n_malformed_rows": stats["n_malformed_rows"],
        "minimum_fdr": stats["minimum_fdr"],
        "maximum_fdr": stats["maximum_fdr"],
        "n_samples": len(stats["samples"]),
        "n_breeds": len(stats["breeds"]),
        "n_tissues": len(stats["tissues"]),
    }


def allele_pattern(ref_count: int, alt_count: int) -> str:
    if ref_count > 0 and alt_count > 0:
        return "MIXED_REF_AND_ALT_MAJOR"
    if ref_count > 0:
        return "REF_MAJOR_ONLY"
    if alt_count > 0:
        return "ALT_MAJOR_ONLY"
    return "NONE"


def support_status(record: Mapping[str, object]) -> str:
    if int(record["piggtex_rows_same"]) > 0:
        return "SIGNIFICANT_ASE_SITE_SAME_MAPPED_TISSUE"
    if int(record["piggtex_rows_any"]) > 0:
        return "SIGNIFICANT_ASE_SITE_OTHER_TISSUE_ONLY"
    return "NOT_OBSERVED_IN_SIGNIFICANT_ONLY_SITE_TABLE"


def output_row(record: Mapping[str, object]) -> Sequence[object]:
    recurrence = record["recurrence"] or {}
    ref, alt = str(record["ref"]), str(record["alt"])
    return (
        record["internal_tissue"], record["mapped_tissue"],
        physical_snp_id(str(record["contig"]), int(record["position"]), ref, alt),
        record["contig"], record["position"], ref, alt, record["allele_1"], record["allele_2"],
        len(record["internal_samples"]), ",".join(sorted(record["internal_samples"])),
        record["n_internal_major_ref"], record["n_internal_major_alt"],
        allele_pattern(int(record["n_internal_major_ref"]), int(record["n_internal_major_alt"])),
        int(bool(recurrence)), recurrence.get("n_animals_testable", ""),
        recurrence.get("n_animals_with_snp_level_ase", ""),
        recurrence.get("recurrence_q_bh_global", ""), len(record["block_contexts"]),
        len(record["block_concordant_contexts"]), len(record["block_discordant_contexts"]),
        ",".join(sorted(record["block_samples"])), ",".join(sorted(record["block_gene_ids"])),
        ",".join(sorted(record["block_gene_names"])), ",".join(sorted(record["block_ids"])),
        ";".join(f"{gene}|{tissue}" for gene, tissue in sorted(record["tissue_enriched_pairs"])),
        record["piggtex_rows_any"], record["piggtex_rows_same"],
        ",".join(sorted(record["piggtex_samples_any"])),
        ",".join(sorted(record["piggtex_samples_same"])),
        ",".join(sorted(record["piggtex_breeds_any"])),
        ",".join(sorted(record["piggtex_breeds_same"])),
        ",".join(sorted(record["piggtex_tissues_any"])),
        fmt_float(record["piggtex_min_fdr_any"]), fmt_float(record["piggtex_min_fdr_same"]),
        fmt_float(max(record["piggtex_pcadd"]) if record["piggtex_pcadd"] else None),
        record["piggtex_same_major_ref"], record["piggtex_same_major_alt"],
        record["piggtex_same_tie"],
        allele_pattern(int(record["piggtex_same_major_ref"]), int(record["piggtex_same_major_alt"])),
        support_status(record),
    )


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with open_text(path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def stage_summary(records: Sequence[Mapping[str, object]]) -> List[Sequence[object]]:
    stages = [
        ("all_independent_GATK_SNP_ASE_tissue_sites", lambda record: True),
        ("statistically_recurrent_GATK_SNP_ASE", lambda record: bool(record["recurrence"])),
        ("primary_block_direction_concordant_GATK_SNP_ASE", lambda record: len(record["block_concordant_contexts"]) > 0),
        ("recurrent_and_primary_block_direction_concordant", lambda record: bool(record["recurrence"]) and len(record["block_concordant_contexts"]) > 0),
        ("tissue_enriched_pair_primary_block_direction_concordant", lambda record: len(record["tissue_enriched_pairs"]) > 0 and len(record["block_concordant_contexts"]) > 0),
    ]
    rows = []
    for name, predicate in stages:
        selected = [record for record in records if predicate(record)]
        rows.append(
            (
                name, len(selected),
                len({unordered_key(str(r["contig"]), int(r["position"]), str(r["ref"]), str(r["alt"])) for r in selected}),
                sum(int(r["piggtex_rows_any"]) > 0 for r in selected),
                sum(int(r["piggtex_rows_same"]) > 0 for r in selected),
                sum(int(r["piggtex_rows_any"]) == 0 for r in selected),
            )
        )
    return rows


def pair_summary(
    records: Sequence[Mapping[str, object]],
    pairs: Mapping[Tuple[str, str], str],
    tissue_map: Mapping[str, str],
) -> List[Sequence[object]]:
    rows = []
    for (gene_id, tissue), gene_name in sorted(pairs.items(), key=lambda item: (item[0][1], item[1], item[0][0])):
        selected = [
            record for record in records if (gene_id, tissue) in record["tissue_enriched_pairs"]
        ]
        concordant = [
            record for record in selected
            if any(context[2] == gene_id for context in record["block_concordant_contexts"])
        ]
        discordant = [
            record for record in selected
            if any(context[2] == gene_id for context in record["block_discordant_contexts"])
        ]
        concordant_ids = sorted(
            physical_snp_id(str(r["contig"]), int(r["position"]), str(r["ref"]), str(r["alt"]))
            for r in concordant
        )
        same_supported = sorted(
            physical_snp_id(str(r["contig"]), int(r["position"]), str(r["ref"]), str(r["alt"]))
            for r in concordant if int(r["piggtex_rows_same"]) > 0
        )
        animals = {
            context[0]
            for record in concordant for context in record["block_concordant_contexts"]
            if context[2] == gene_id
        }
        rows.append(
            (
                gene_id, gene_name, tissue, tissue_map[tissue], len(selected), len(concordant),
                len(discordant), sum(bool(r["recurrence"]) for r in concordant),
                sum(int(r["piggtex_rows_any"]) > 0 for r in concordant),
                len(same_supported), len(animals), ",".join(concordant_ids),
                ",".join(same_supported),
                "PigGTEx site table is significant-only; absence is not evidence of no ASE",
            )
        )
    return rows


def prepare_output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise AnalysisError(f"output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def run(args: argparse.Namespace) -> None:
    output = args.out_dir.resolve()
    prepare_output_directory(output)
    tissue_map = load_tissue_map(args.tissue_map.resolve())
    pairs = load_tissue_enriched_pairs(args.tissue_enriched_pairs.resolve())
    if len(pairs) != 39:
        raise AnalysisError(f"expected 39 primary tissue-enriched pairs, found {len(pairs)}")
    records_map = load_gatk_calls(args.gatk_calls.resolve(), tissue_map)
    n_recurrent = attach_recurrence(args.gatk_recurrent.resolve(), records_map)
    block_counts = attach_block_contexts(
        args.block_validation.resolve(), records_map, pairs
    )
    piggtex_stats = stream_piggtex(args.piggtex_sites.resolve(), records_map)
    records = sorted(
        records_map.values(),
        key=lambda record: (
            str(record["internal_tissue"]), str(record["contig"]), int(record["position"]),
            str(record["ref"]), str(record["alt"]),
        ),
    )

    all_path = output / "gatk_snp_piggtex_site_support.tsv.gz"
    recurrent_path = output / "recurrent_gatk_snp_piggtex_site_support.tsv"
    block_path = output / "primary_block_concordant_gatk_snp_piggtex_site_support.tsv"
    stage_path = output / "piggtex_site_support_stage_summary.tsv"
    pair_path = output / "tissue_enriched_pair_gatk_piggtex_site_summary.tsv"
    write_tsv(all_path, OUTPUT_HEADER, (output_row(record) for record in records))
    write_tsv(
        recurrent_path, OUTPUT_HEADER,
        (output_row(record) for record in records if record["recurrence"]),
    )
    write_tsv(
        block_path, OUTPUT_HEADER,
        (output_row(record) for record in records if record["block_concordant_contexts"]),
    )
    stage_rows = stage_summary(records)
    pair_rows = pair_summary(records, pairs, tissue_map)
    write_tsv(stage_path, STAGE_HEADER, stage_rows)
    write_tsv(pair_path, PAIR_HEADER, pair_rows)

    counts = {
        "n_gatk_ase_measurements": sum(len(record["internal_samples"]) for record in records),
        "n_gatk_ase_tissue_snp_sites": len(records),
        "n_recurrent_gatk_tissue_snp_sites": n_recurrent,
        "n_primary_tissue_enriched_pairs": len(pairs),
        **block_counts,
    }
    metadata = {
        "version": VERSION,
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "exact_site_key": "Sscrofa11.1 contig + 1-based position + unordered biallelic allele pair",
        "external_data_boundary": (
            "PigGTEx site ASE is a significant-only table with FDR<0.05; absence is recorded "
            "as not observed, never as a tested negative."
        ),
        "counts": counts,
        "piggtex_input_audit": piggtex_stats,
        "inputs": {
            "gatk_calls": file_metadata(args.gatk_calls.resolve()),
            "gatk_recurrent": file_metadata(args.gatk_recurrent.resolve()),
            "block_validation": file_metadata(args.block_validation.resolve()),
            "tissue_enriched_pairs": file_metadata(args.tissue_enriched_pairs.resolve()),
            "tissue_map": file_metadata(args.tissue_map.resolve()),
            "piggtex_sites": file_metadata(args.piggtex_sites.resolve()),
        },
        "outputs": {},
    }
    for path in (all_path, recurrent_path, block_path, stage_path, pair_path):
        metadata["outputs"][path.name] = file_metadata(path)
    with (output / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gatk-calls", type=Path, required=True)
    parser.add_argument("--gatk-recurrent", type=Path, required=True)
    parser.add_argument("--block-validation", type=Path, required=True)
    parser.add_argument("--tissue-enriched-pairs", type=Path, required=True)
    parser.add_argument("--tissue-map", type=Path, required=True)
    parser.add_argument("--piggtex-sites", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    try:
        run(parse_args(argv))
    except (AnalysisError, OSError, ValueError, KeyError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
