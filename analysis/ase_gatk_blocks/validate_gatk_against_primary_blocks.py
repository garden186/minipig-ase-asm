#!/usr/bin/env python3
"""Link independent GATK SNP-level ASE to current v0.2.6 primary blocks."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
import tarfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Tuple


VERSION = "ase_gatk_primary_block_validation_v0.1.0"

SITE_REQUIRED = {
    "sample", "tissue", "contig", "position", "variant_id", "ref_allele",
    "alt_allele", "ref_count", "alt_count", "total_count",
    "major_allele_fraction", "major_allele", "imbalance_direction",
    "p_exact_two_sided", "q_bh_within_animal_tissue", "is_snp_level_ase",
}

CORE_MEASUREMENT_REQUIRED = {
    "sample", "tissue", "gene_id", "gene_name", "top_haplotype_block_id",
    "core_haplotype_block_ids", "gene_ase_status_v026",
}

BLOCK_REQUIRED = {
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id",
    "measurement_id", "block_snps", "gene_variants", "haplotypeA",
    "haplotypeB", "v026_ase_direction", "v026_evidence_class",
    "v026_is_core_ase",
}

RECURRENCE_REQUIRED = {
    "tissue", "physical_snp_id", "n_animals_testable",
    "n_animals_with_snp_level_ase", "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
}

SITE_OUTPUT_HEADER = [
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id",
    "measurement_id", "is_representative_block", "v026_evidence_class",
    "v026_ase_direction", "variant_id", "contig", "position", "ref_allele",
    "alt_allele", "haplotype_a_allele", "haplotype_b_allele",
    "predicted_major_allele", "predicted_minor_allele", "gatk_site_testable",
    "gatk_ref_count", "gatk_alt_count", "gatk_total_count",
    "gatk_major_allele_fraction", "gatk_major_allele", "gatk_direction",
    "gatk_p_exact_two_sided", "gatk_q_bh_within_animal_tissue",
    "gatk_snp_level_ase", "gatk_block_direction_status",
    "gatk_site_is_statistically_recurrent", "gatk_recurrence_n_testable_animals",
    "gatk_recurrence_n_ase_positive_animals", "gatk_recurrence_q_bh_global",
]

BLOCK_SUMMARY_HEADER = [
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id",
    "measurement_id", "is_representative_block", "v026_evidence_class",
    "v026_ase_direction", "n_gene_variants", "n_gatk_testable_variants",
    "n_gatk_snp_level_ase_variants", "n_direction_concordant_variants",
    "n_direction_discordant_variants", "n_recurrent_gatk_ase_variants",
    "gatk_block_validation_status",
]

MEASUREMENT_SUMMARY_HEADER = [
    "sample", "tissue", "gene_id", "gene_name", "gene_ase_status_v026",
    "representative_block_id", "n_core_blocks", "n_gene_variant_contexts",
    "n_gatk_testable_variants", "n_gatk_snp_level_ase_variants",
    "n_direction_concordant_variants", "n_direction_discordant_variants",
    "n_recurrent_gatk_ase_variants", "gatk_measurement_validation_status",
]

TISSUE_SUMMARY_HEADER = [
    "tissue", "n_primary_gene_measurements", "n_core_block_gene_contexts",
    "n_core_gene_variant_contexts", "n_gatk_testable_variant_contexts",
    "n_gatk_snp_level_ase_variant_contexts", "n_direction_concordant_contexts",
    "n_direction_discordant_contexts", "n_recurrent_gatk_ase_contexts",
    "n_measurements_with_gatk_testable_variant", "n_measurements_with_gatk_ase",
    "n_measurements_with_direction_concordant_gatk_ase",
    "n_measurements_with_direction_discordant_gatk_ase",
]


class AnalysisError(RuntimeError):
    """Raised when an input or block-linkage invariant is violated."""


@contextmanager
def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def iter_tsv(path: Path, required: set[str]) -> Iterator[Dict[str, str]]:
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
        "blood": "Blood",
        "l-blood": "Blood",
        "tenderlo": "Tenderloin",
        "tenderloin": "Tenderloin",
    }.get(stripped.lower(), stripped)


def parse_variant_id(value: str) -> Tuple[str, int, str, str]:
    parts = value.strip().split("_")
    if len(parts) < 4:
        raise AnalysisError(f"invalid variant_id: {value!r}")
    contig = "_".join(parts[:-3])
    try:
        position = int(parts[-3])
    except ValueError as error:
        raise AnalysisError(f"invalid variant position: {value!r}") from error
    ref, alt = parts[-2].upper(), parts[-1].upper()
    return contig, position, ref, alt


def physical_snp_id(contig: str, position: int, ref: str, alt: str) -> str:
    return f"{contig}_{position}_{ref}_{alt}"


def parse_bool(value: str, field: str, context: str) -> int:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return 1
    if normalized in {"0", "false", "no"}:
        return 0
    raise AnalysisError(f"{context}: invalid Boolean {field}={value!r}")


def fmt_float(value: object) -> str:
    if value is None or value == "":
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


@dataclass
class Context:
    sample: str
    tissue: str
    gene_id: str
    gene_name: str
    block_id: str
    measurement_id: str
    is_representative: int
    evidence_class: str
    block_direction: str
    variant_id: str
    contig: str
    position: int
    ref: str
    alt: str
    hap_a: str
    hap_b: str
    predicted_major: str
    predicted_minor: str
    gatk: Dict[str, str] | None = None
    recurrence: Dict[str, str] | None = None


def load_measurement_registry(
    path: Path,
) -> Tuple[Dict[Tuple[str, str, str], Dict[str, str]], set[Tuple[str, str, str, str]]]:
    measurements: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    expected_blocks: set[Tuple[str, str, str, str]] = set()
    for row in iter_tsv(path, CORE_MEASUREMENT_REQUIRED):
        key = (row["sample"].strip(), normalize_tissue(row["tissue"]), row["gene_id"].strip())
        if key in measurements:
            raise AnalysisError(f"duplicate primary gene measurement: {key}")
        normalized = dict(row)
        normalized["tissue"] = key[1]
        measurements[key] = normalized
        for block_id in (item for item in row["core_haplotype_block_ids"].split(",") if item):
            expected_blocks.add((*key, block_id))
    return measurements, expected_blocks


def load_core_contexts(
    archive: Path,
    measurements: Mapping[Tuple[str, str, str], Mapping[str, str]],
    expected_blocks: set[Tuple[str, str, str, str]],
) -> Tuple[List[Context], Dict[Tuple[str, str, str], List[int]], Dict[str, int]]:
    contexts: List[Context] = []
    lookup: Dict[Tuple[str, str, str], List[int]] = defaultdict(list)
    observed_blocks: set[Tuple[str, str, str, str]] = set()
    counts = Counter()
    with tarfile.open(archive, "r:gz") as handle:
        members = sorted(
            (
                member for member in handle.getmembers()
                if member.isfile() and member.name.endswith("/ase_core_haplotype_blocks.tsv")
            ),
            key=lambda member: member.name,
        )
        if len(members) != 10:
            raise AnalysisError(f"expected 10 core block tables, found {len(members)}")
        for member in members:
            extracted = handle.extractfile(member)
            if extracted is None:
                raise AnalysisError(f"cannot read archive member: {member.name}")
            reader = csv.DictReader((line.decode("utf-8") for line in extracted), delimiter="\t")
            missing = sorted(BLOCK_REQUIRED - set(reader.fieldnames or []))
            if missing:
                raise AnalysisError(f"{member.name}: missing columns: {', '.join(missing)}")
            for row in reader:
                sample = row["sample"].strip()
                tissue = normalize_tissue(row["tissue"])
                gene_id = row["gene_id"].strip()
                block_id = row["haplotype_block_id"].strip()
                measurement_key = (sample, tissue, gene_id)
                block_key = (*measurement_key, block_id)
                if measurement_key not in measurements:
                    raise AnalysisError(f"core block lacks primary gene measurement: {block_key}")
                if parse_bool(row["v026_is_core_ase"], "v026_is_core_ase", str(block_key)) != 1:
                    raise AnalysisError(f"non-core row in core archive: {block_key}")
                observed_blocks.add(block_key)
                block_snps = [item for item in row["block_snps"].split(",") if item]
                gene_variants = [item for item in row["gene_variants"].split(",") if item]
                hap_a = [item.upper() for item in row["haplotypeA"].split(",") if item]
                hap_b = [item.upper() for item in row["haplotypeB"].split(",") if item]
                if not (len(block_snps) == len(hap_a) == len(hap_b)):
                    raise AnalysisError(f"block haplotype vector length mismatch: {block_key}")
                index = {variant_id: offset for offset, variant_id in enumerate(block_snps)}
                if len(index) != len(block_snps):
                    raise AnalysisError(f"duplicate variant within block: {block_key}")
                direction = row["v026_ase_direction"].strip()
                if direction not in {"HAP_A", "HAP_B"}:
                    raise AnalysisError(f"unresolved current core block direction: {block_key}")
                is_representative = int(
                    measurements[measurement_key]["top_haplotype_block_id"] == block_id
                )
                counts["n_core_block_gene_contexts"] += 1
                counts["n_representative_block_gene_contexts"] += is_representative
                for variant_id in gene_variants:
                    if variant_id not in index:
                        raise AnalysisError(f"gene variant absent from block SNPs: {block_key}|{variant_id}")
                    offset = index[variant_id]
                    contig, position, ref, alt = parse_variant_id(variant_id)
                    allele_a, allele_b = hap_a[offset], hap_b[offset]
                    if {allele_a, allele_b} != {ref, alt}:
                        raise AnalysisError(f"haplotype alleles do not match SNP: {block_key}|{variant_id}")
                    major = allele_a if direction == "HAP_A" else allele_b
                    minor = allele_b if direction == "HAP_A" else allele_a
                    context = Context(
                        sample, tissue, gene_id, row["gene_name"].strip(), block_id,
                        row["measurement_id"].strip(), is_representative,
                        row["v026_evidence_class"].strip(), direction, variant_id,
                        contig, position, ref, alt, allele_a, allele_b, major, minor,
                    )
                    context_index = len(contexts)
                    contexts.append(context)
                    lookup[(sample, tissue, variant_id)].append(context_index)
                    counts["n_core_gene_variant_contexts"] += 1
    if observed_blocks != expected_blocks:
        missing = sorted(expected_blocks - observed_blocks)[:10]
        extra = sorted(observed_blocks - expected_blocks)[:10]
        raise AnalysisError(f"core block registry mismatch; missing={missing}, extra={extra}")
    counts["n_primary_gene_measurements"] = len(measurements)
    counts["n_unique_core_block_registry_keys"] = len(observed_blocks)
    return contexts, dict(lookup), dict(counts)


def attach_gatk_sites(
    path: Path,
    contexts: List[Context],
    lookup: Mapping[Tuple[str, str, str], List[int]],
) -> Dict[str, int]:
    counts = Counter()
    matched_contexts: set[int] = set()
    for row in iter_tsv(path, SITE_REQUIRED):
        key = (row["sample"].strip(), normalize_tissue(row["tissue"]), row["variant_id"].strip())
        target_indices = lookup.get(key)
        if not target_indices:
            continue
        counts["n_gatk_site_rows_matching_any_core_context"] += 1
        for index in target_indices:
            if index in matched_contexts:
                raise AnalysisError(f"duplicate GATK site for core context: {key}")
            context = contexts[index]
            if (
                row["contig"].strip(), int(row["position"]), row["ref_allele"].upper(),
                row["alt_allele"].upper()
            ) != (context.contig, context.position, context.ref, context.alt):
                raise AnalysisError(f"GATK SNP identity mismatch: {key}")
            context.gatk = dict(row)
            matched_contexts.add(index)
            counts["n_gatk_testable_core_contexts"] += 1
            success = parse_bool(row["is_snp_level_ase"], "is_snp_level_ase", str(key))
            counts["n_gatk_ase_core_contexts"] += success
            if success:
                major = row["major_allele"].strip().upper()
                if major == context.predicted_major:
                    counts["n_direction_concordant_core_contexts"] += 1
                elif major == context.predicted_minor:
                    counts["n_direction_discordant_core_contexts"] += 1
                else:
                    raise AnalysisError(f"GATK major allele is outside block haplotypes: {key}")
    return dict(counts)


def load_recurrence(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    recurrence: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in iter_tsv(path, RECURRENCE_REQUIRED):
        if parse_bool(
            row["is_statistically_recurrent_snp_ase"],
            "is_statistically_recurrent_snp_ase",
            row["physical_snp_id"],
        ) != 1:
            continue
        key = (normalize_tissue(row["tissue"]), row["physical_snp_id"].strip())
        if key in recurrence:
            raise AnalysisError(f"duplicate recurrent tissue-SNP: {key}")
        recurrence[key] = dict(row)
    return recurrence


def attach_recurrence(
    contexts: List[Context], recurrence: Mapping[Tuple[str, str], Dict[str, str]]
) -> int:
    count = 0
    for context in contexts:
        key = (
            context.tissue,
            physical_snp_id(context.contig, context.position, context.ref, context.alt),
        )
        context.recurrence = recurrence.get(key)
        if context.recurrence is not None and context.gatk is not None:
            if parse_bool(context.gatk["is_snp_level_ase"], "is_snp_level_ase", str(key)):
                count += 1
    return count


def direction_status(context: Context) -> str:
    if context.gatk is None:
        return "NOT_GATK_TESTABLE"
    if parse_bool(context.gatk["is_snp_level_ase"], "is_snp_level_ase", context.variant_id) == 0:
        return "GATK_TESTABLE_NOT_ASE"
    major = context.gatk["major_allele"].strip().upper()
    if major == context.predicted_major:
        return "DIRECTION_CONCORDANT"
    if major == context.predicted_minor:
        return "DIRECTION_DISCORDANT"
    raise AnalysisError(f"unresolved direction status: {context.variant_id}")


def context_row(context: Context) -> Sequence[object]:
    gatk = context.gatk or {}
    recurrence = context.recurrence or {}
    gatk_success = (
        parse_bool(gatk["is_snp_level_ase"], "is_snp_level_ase", context.variant_id)
        if gatk else 0
    )
    recurrent_success = int(bool(recurrence) and gatk_success == 1)
    return (
        context.sample, context.tissue, context.gene_id, context.gene_name,
        context.block_id, context.measurement_id, context.is_representative,
        context.evidence_class, context.block_direction, context.variant_id,
        context.contig, context.position, context.ref, context.alt, context.hap_a,
        context.hap_b, context.predicted_major, context.predicted_minor,
        int(bool(gatk)), gatk.get("ref_count", ""), gatk.get("alt_count", ""),
        gatk.get("total_count", ""), gatk.get("major_allele_fraction", ""),
        gatk.get("major_allele", ""), gatk.get("imbalance_direction", ""),
        gatk.get("p_exact_two_sided", ""), gatk.get("q_bh_within_animal_tissue", ""),
        gatk_success, direction_status(context), recurrent_success,
        recurrence.get("n_animals_testable", ""),
        recurrence.get("n_animals_with_snp_level_ase", ""),
        recurrence.get("recurrence_q_bh_global", ""),
    )


def status_from_counts(testable: int, ase: int, concordant: int, discordant: int) -> str:
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
    raise AnalysisError("invalid validation-status counts")


def summarize_contexts(
    contexts: Sequence[Context],
    key_fields: str,
) -> List[Tuple[Tuple[str, ...], Dict[str, object]]]:
    summaries: MutableMapping[Tuple[str, ...], Dict[str, object]] = {}
    for context in contexts:
        if key_fields == "block":
            key = (context.sample, context.tissue, context.gene_id, context.block_id)
        elif key_fields == "measurement":
            key = (context.sample, context.tissue, context.gene_id)
        else:
            raise ValueError(key_fields)
        summary = summaries.setdefault(
            key,
            {
                "sample": context.sample,
                "tissue": context.tissue,
                "gene_id": context.gene_id,
                "gene_name": context.gene_name,
                "block_ids": set(),
                "measurement_ids": set(),
                "representative_block_ids": set(),
                "evidence_classes": set(),
                "block_directions": set(),
                "n_contexts": 0,
                "n_testable": 0,
                "n_ase": 0,
                "n_concordant": 0,
                "n_discordant": 0,
                "n_recurrent": 0,
            },
        )
        summary["block_ids"].add(context.block_id)
        summary["measurement_ids"].add(context.measurement_id)
        if context.is_representative:
            summary["representative_block_ids"].add(context.block_id)
        summary["evidence_classes"].add(context.evidence_class)
        summary["block_directions"].add(context.block_direction)
        summary["n_contexts"] += 1
        status = direction_status(context)
        if context.gatk is not None:
            summary["n_testable"] += 1
        if status in {"DIRECTION_CONCORDANT", "DIRECTION_DISCORDANT"}:
            summary["n_ase"] += 1
        if status == "DIRECTION_CONCORDANT":
            summary["n_concordant"] += 1
        if status == "DIRECTION_DISCORDANT":
            summary["n_discordant"] += 1
        if context.recurrence is not None and status in {"DIRECTION_CONCORDANT", "DIRECTION_DISCORDANT"}:
            summary["n_recurrent"] += 1
    return sorted(summaries.items())


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with open_text(path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def prepare_output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise AnalysisError(f"output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def run(args: argparse.Namespace) -> None:
    output = args.out_dir.resolve()
    prepare_output_directory(output)
    measurements, expected_blocks = load_measurement_registry(args.core_measurements.resolve())
    contexts, lookup, counts = load_core_contexts(
        args.block_archive.resolve(), measurements, expected_blocks
    )
    counts.update(attach_gatk_sites(args.site_tests.resolve(), contexts, lookup))
    recurrence = load_recurrence(args.recurrent_snps.resolve())
    counts["n_recurrent_gatk_ase_core_contexts"] = attach_recurrence(contexts, recurrence)
    contexts.sort(
        key=lambda item: (
            item.sample, item.tissue, item.gene_id, item.block_id,
            item.contig, item.position, item.ref, item.alt,
        )
    )

    block_summaries = summarize_contexts(contexts, "block")
    measurement_summaries = summarize_contexts(contexts, "measurement")
    if len(measurement_summaries) != len(measurements):
        raise AnalysisError("primary gene measurement summary count mismatch")

    site_path = output / "primary_block_gatk_snp_validation.tsv.gz"
    block_path = output / "primary_block_validation_summary.tsv"
    measurement_path = output / "primary_gene_measurement_validation_summary.tsv"
    tissue_path = output / "validation_summary_by_tissue.tsv"
    write_tsv(site_path, SITE_OUTPUT_HEADER, (context_row(context) for context in contexts))

    block_rows = []
    for _, summary in block_summaries:
        if len(summary["block_ids"]) != 1 or len(summary["measurement_ids"]) != 1:
            raise AnalysisError("block summary key is not unique")
        block_id = next(iter(summary["block_ids"]))
        block_rows.append(
            (
                summary["sample"], summary["tissue"], summary["gene_id"],
                summary["gene_name"], block_id, next(iter(summary["measurement_ids"])),
                int(block_id in summary["representative_block_ids"]),
                next(iter(summary["evidence_classes"])),
                next(iter(summary["block_directions"])), summary["n_contexts"],
                summary["n_testable"], summary["n_ase"], summary["n_concordant"],
                summary["n_discordant"], summary["n_recurrent"],
                status_from_counts(
                    summary["n_testable"], summary["n_ase"],
                    summary["n_concordant"], summary["n_discordant"],
                ),
            )
        )
    write_tsv(block_path, BLOCK_SUMMARY_HEADER, block_rows)

    measurement_rows = []
    measurement_status_by_tissue: Dict[str, List[Tuple[object, ...]]] = defaultdict(list)
    for key, summary in measurement_summaries:
        registry = measurements[key]
        representative = registry["top_haplotype_block_id"]
        if summary["representative_block_ids"] != {representative}:
            raise AnalysisError(f"representative block mismatch: {key}")
        row = (
            summary["sample"], summary["tissue"], summary["gene_id"], summary["gene_name"],
            registry["gene_ase_status_v026"], representative, len(summary["block_ids"]),
            summary["n_contexts"], summary["n_testable"], summary["n_ase"],
            summary["n_concordant"], summary["n_discordant"], summary["n_recurrent"],
            status_from_counts(
                summary["n_testable"], summary["n_ase"],
                summary["n_concordant"], summary["n_discordant"],
            ),
        )
        measurement_rows.append(row)
        measurement_status_by_tissue[str(summary["tissue"])].append(row)
    write_tsv(measurement_path, MEASUREMENT_SUMMARY_HEADER, measurement_rows)

    context_by_tissue: Dict[str, List[Context]] = defaultdict(list)
    blocks_by_tissue: Counter[str] = Counter()
    for context in contexts:
        context_by_tissue[context.tissue].append(context)
    for row in block_rows:
        blocks_by_tissue[str(row[1])] += 1
    tissue_rows = []
    for tissue in sorted(measurement_status_by_tissue):
        tissue_contexts = context_by_tissue[tissue]
        rows = measurement_status_by_tissue[tissue]
        statuses = [direction_status(context) for context in tissue_contexts]
        tissue_rows.append(
            (
                tissue, len(rows), blocks_by_tissue[tissue], len(tissue_contexts),
                sum(context.gatk is not None for context in tissue_contexts),
                sum(status in {"DIRECTION_CONCORDANT", "DIRECTION_DISCORDANT"} for status in statuses),
                statuses.count("DIRECTION_CONCORDANT"), statuses.count("DIRECTION_DISCORDANT"),
                sum(
                    context.recurrence is not None
                    and direction_status(context) in {"DIRECTION_CONCORDANT", "DIRECTION_DISCORDANT"}
                    for context in tissue_contexts
                ),
                sum(int(row[8]) > 0 for row in rows), sum(int(row[9]) > 0 for row in rows),
                sum(int(row[10]) > 0 for row in rows), sum(int(row[11]) > 0 for row in rows),
            )
        )
    write_tsv(tissue_path, TISSUE_SUMMARY_HEADER, tissue_rows)

    counts["n_core_block_summary_rows"] = len(block_rows)
    counts["n_primary_measurement_summary_rows"] = len(measurement_rows)
    counts["n_tissues"] = len(tissue_rows)
    counts["n_independent_recurrent_tissue_snps"] = len(recurrence)
    metadata = {
        "version": VERSION,
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent GATK SNP-level ASE support within current v0.2.6 primary local haplotype blocks",
        "interpretation_boundary": (
            "This validates site support and local block direction; it does not establish "
            "causal SNPs, parental origin, or independence of biological samples."
        ),
        "counts": counts,
        "inputs": {
            "site_tests": file_metadata(args.site_tests.resolve()),
            "core_measurements": file_metadata(args.core_measurements.resolve()),
            "block_archive": file_metadata(args.block_archive.resolve()),
            "recurrent_snps": file_metadata(args.recurrent_snps.resolve()),
        },
        "outputs": {},
    }
    for path in (site_path, block_path, measurement_path, tissue_path):
        metadata["outputs"][path.name] = file_metadata(path)
    with (output / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-tests", type=Path, required=True)
    parser.add_argument("--core-measurements", type=Path, required=True)
    parser.add_argument("--block-archive", type=Path, required=True)
    parser.add_argument("--recurrent-snps", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    try:
        run(parse_args(argv))
    except (AnalysisError, OSError, tarfile.TarError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
