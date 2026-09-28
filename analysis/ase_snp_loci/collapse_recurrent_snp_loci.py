#!/usr/bin/env python3
"""Collapse recurrent gene-tissue-SNP hypotheses into core-block-connected ASE loci.

Two recurrent SNPs are connected only when, in at least one successful animal,
both SNPs are supported by the same final v0.2.6 core ASE measurement for the
same gene and tissue. Connected components of this evidence graph are reported
as loci. Genomic distance is never used to create an edge.
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
import statistics
import sys
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "ase_snp_locus_collapse_v0.1.0"
GROUPING_METHOD = "SAME_GENE_TISSUE_CORE_BLOCK_CONNECTED_COMPONENT"
TISSUE_NORMALIZATION = {
    "Tenderlo": "Tenderloin",
    "blood": "Blood",
    "L-blood": "Blood",
}
VARIANT_RE = re.compile(r"^([^_]+)_(\d+)_([ACGTN]+)_([ACGTN]+)$", re.I)

RECURRENT_REQUIRED = {
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_animals_testable",
    "animals_testable",
    "n_animals_with_block_concordant_snp_ase",
    "animals_with_block_concordant_snp_ase",
    "recurrence_fraction",
    "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
    "original_ase_assignment_class",
    "canonical_transcript_features",
    "annotation_status",
}
SITE_REQUIRED = {
    "sample",
    "tissue",
    "variant_id",
    "ref_fragment_count",
    "alt_fragment_count",
    "total_informative_fragments",
    "major_allele",
    "site_level_ase_signal",
    "core_parent_measurement_ids",
}
BLOCK_REQUIRED = {
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "haplotype_block_id",
    "block_start",
    "block_end",
    "block_snps",
    "haplotypeA",
    "haplotypeB",
    "gene_variants_with_reads",
    "measurement_id",
    "analysis_set",
    "v026_technical_qc_status",
    "v026_is_core_ase",
    "v026_ase_direction",
}

LOCUS_FIELDS = [
    "locus_id",
    "locus_grouping_method",
    "gene_id",
    "gene_name",
    "tissue",
    "chrom",
    "locus_start",
    "locus_end",
    "locus_span_bp",
    "n_recurrent_snps",
    "recurrent_snp_ids",
    "lead_snp_id",
    "lead_snp_position",
    "lead_recurrence_q_bh_global",
    "lead_recurrence_p_poisson_binomial",
    "lead_recurrence_fraction",
    "lead_n_animals_with_block_concordant_snp_ase",
    "lead_median_success_total_informative_fragments",
    "lead_original_ase_assignment_class",
    "lead_canonical_transcript_features",
    "n_animals_with_any_recurrent_snp",
    "animals_with_any_recurrent_snp",
    "n_animals_with_all_recurrent_snps",
    "animals_with_all_recurrent_snps",
    "n_success_parent_block_links",
    "n_distinct_success_sample_blocks",
    "n_direct_core_block_edges",
    "n_animals_with_multi_snp_shared_core_block",
    "animals_with_multi_snp_shared_core_block",
    "original_ase_assignment_classes",
    "canonical_transcript_features",
    "locus_is_descriptive_collapse_of_significant_snps",
]

LINK_FIELDS = [
    "locus_id",
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "sample",
    "measurement_id",
    "haplotype_block_id",
    "block_start",
    "block_end",
    "v026_ase_direction",
    "parent_major_allele_at_snp",
    "site_major_allele",
    "ref_fragment_count",
    "alt_fragment_count",
    "total_informative_fragments",
]

EDGE_FIELDS = [
    "locus_id",
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id_1",
    "position_1",
    "variant_id_2",
    "position_2",
    "n_shared_success_sample_blocks",
    "n_shared_success_animals",
    "shared_success_animals",
    "shared_measurement_ids",
    "shared_haplotype_block_ids",
]


class AnalysisError(RuntimeError):
    pass


def normalize_sample(value: str) -> str:
    text = str(value).strip()
    return text.zfill(4) if text.isdigit() else text


def normalize_tissue(value: str) -> str:
    text = str(value).strip()
    return TISSUE_NORMALIZATION.get(text, text)


def normalize_gene_id(value: str) -> str:
    return str(value).strip().split(".")[0]


def split_values(value: str) -> List[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def parse_bool(value: str, field: str, context: str) -> bool:
    text = str(value).strip().upper()
    if text in {"1", "TRUE", "YES"}:
        return True
    if text in {"0", "FALSE", "NO"}:
        return False
    raise AnalysisError(f"{context}: invalid Boolean in {field}: {value!r}")


def parse_int(value: str, field: str, context: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{context}: invalid integer in {field}: {value!r}") from error


def parse_float(value: str, field: str, context: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{context}: invalid number in {field}: {value!r}") from error
    if not math.isfinite(result):
        raise AnalysisError(f"{context}: non-finite number in {field}: {value!r}")
    return result


def fmt_float(value: float | int) -> str:
    return f"{float(value):.12g}"


def parse_variant(value: str) -> Tuple[str, int, str, str]:
    match = VARIANT_RE.match(str(value).strip())
    if not match:
        raise AnalysisError(f"malformed SNP identifier: {value!r}")
    chrom, position, ref, alt = match.groups()
    ref, alt = ref.upper(), alt.upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        raise AnalysisError(f"non-biallelic-SNP identifier: {value!r}")
    return chrom.removeprefix("chr"), int(position), ref, alt


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def iter_tsv(path: Path, required: Set[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AnalysisError(f"{path}: missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value if value is not None else "" for key, value in row.items()}


def input_header(path: Path) -> List[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as error:
            raise AnalysisError(f"empty input: {path}") from error


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def fingerprint(path: Path) -> Dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    status = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": status.st_size,
        "modified_utc": datetime.fromtimestamp(status.st_mtime, timezone.utc).isoformat(),
        "sha256": digest.hexdigest(),
    }


@dataclass(frozen=True)
class Block:
    sample: str
    tissue: str
    gene_id: str
    gene_name: str
    block_id: str
    measurement_id: str
    block_start: int
    block_end: int
    direction: str
    variants_with_reads: frozenset[str]
    parent_major_alleles: Mapping[str, str]


@dataclass(frozen=True)
class Link:
    gene_id: str
    gene_name: str
    tissue: str
    variant_id: str
    chrom: str
    position: int
    sample: str
    measurement_id: str
    block_id: str
    block_start: int
    block_end: int
    direction: str
    parent_major_allele: str
    site_major_allele: str
    ref_count: int
    alt_count: int
    total_count: int


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            next_value = self.parent[value]
            self.parent[value] = root
            value = next_value
        return root

    def union(self, first: str, second: str) -> None:
        root_first, root_second = self.find(first), self.find(second)
        if root_first != root_second:
            low, high = sorted((root_first, root_second))
            self.parent[high] = low


def load_recurrent(
    path: Path,
) -> Tuple[List[str], Dict[Tuple[str, str, str], Dict[str, str]], Set[Tuple[str, str, str, str]]]:
    header = input_header(path)
    rows: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    success_keys: Set[Tuple[str, str, str, str]] = set()
    for row in iter_tsv(path, RECURRENT_REQUIRED):
        gene_id = normalize_gene_id(row["gene_id"])
        tissue = normalize_tissue(row["tissue"])
        variant = row["variant_id"].strip()
        context = f"{gene_id}|{tissue}|{variant}"
        if not parse_bool(
            row["is_statistically_recurrent_snp_ase"],
            "is_statistically_recurrent_snp_ase",
            context,
        ):
            raise AnalysisError(f"{path}: non-recurrent row in recurrent input: {context}")
        if row["annotation_status"] != "PASS":
            raise AnalysisError(f"{path}: non-PASS feature annotation: {context}")
        chrom, position, ref, alt = parse_variant(variant)
        if (
            chrom != row["chrom"].removeprefix("chr")
            or position != parse_int(row["position"], "position", context)
            or ref != row["ref"].upper()
            or alt != row["alt"].upper()
        ):
            raise AnalysisError(f"{path}: variant columns disagree for {context}")
        row = dict(row)
        row["gene_id"], row["tissue"] = gene_id, tissue
        key = (gene_id, tissue, variant)
        if key in rows:
            raise AnalysisError(f"{path}: duplicate recurrent hypothesis: {context}")
        successes = split_values(row["animals_with_block_concordant_snp_ase"])
        n_success = parse_int(
            row["n_animals_with_block_concordant_snp_ase"],
            "n_animals_with_block_concordant_snp_ase",
            context,
        )
        if n_success < 2 or len(successes) != n_success or len(successes) != len(set(successes)):
            raise AnalysisError(f"{path}: invalid success count/list for {context}")
        rows[key] = row
        for sample in successes:
            success_key = (normalize_sample(sample), tissue, gene_id, variant)
            if success_key in success_keys:
                raise AnalysisError(f"{path}: duplicate success key: {success_key}")
            success_keys.add(success_key)
    if not rows:
        raise AnalysisError(f"{path}: no recurrent hypotheses")
    return header, rows, success_keys


def load_core_blocks(path: Path) -> Tuple[Dict[str, Block], Dict[str, int]]:
    blocks: Dict[str, Block] = {}
    counts = Counter()
    with tarfile.open(path, "r:*") as archive:
        members = sorted(
            (
                member
                for member in archive.getmembers()
                if member.isfile() and member.name.endswith("/ase_core_haplotype_blocks.tsv")
            ),
            key=lambda member: member.name,
        )
        if not members:
            raise AnalysisError(f"{path}: no ase_core_haplotype_blocks.tsv members")
        counts["n_core_archive_members"] = len(members)
        for member in members:
            binary = archive.extractfile(member)
            if binary is None:
                raise AnalysisError(f"{path}: cannot read {member.name}")
            with io.TextIOWrapper(binary, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                missing = sorted(BLOCK_REQUIRED - set(reader.fieldnames or []))
                if missing:
                    raise AnalysisError(
                        f"{path}:{member.name}: missing columns: {', '.join(missing)}"
                    )
                for line_number, row in enumerate(reader, start=2):
                    counts["n_core_archive_rows"] += 1
                    context = f"{member.name}:{line_number}"
                    if row["analysis_set"] != "MAIN" or row["v026_technical_qc_status"] != "PASS":
                        raise AnalysisError(f"{context}: core archive row is not MAIN/technical-PASS")
                    if not parse_bool(row["v026_is_core_ase"], "v026_is_core_ase", context):
                        raise AnalysisError(f"{context}: non-core row in core block file")
                    variants = split_values(row["block_snps"])
                    hap_a = [value.upper() for value in split_values(row["haplotypeA"])]
                    hap_b = [value.upper() for value in split_values(row["haplotypeB"])]
                    if not (len(variants) == len(hap_a) == len(hap_b)):
                        raise AnalysisError(f"{context}: unaligned block SNP/haplotype columns")
                    allele_a = dict(zip(variants, hap_a))
                    allele_b = dict(zip(variants, hap_b))
                    direction = row["v026_ase_direction"].strip()
                    if direction == "HAP_A":
                        major = allele_a
                    elif direction == "HAP_B":
                        major = allele_b
                    else:
                        raise AnalysisError(f"{context}: invalid core ASE direction {direction!r}")
                    observed = frozenset(split_values(row["gene_variants_with_reads"]))
                    if not observed or not observed.issubset(major):
                        raise AnalysisError(f"{context}: invalid gene_variants_with_reads")
                    measurement_id = row["measurement_id"].strip()
                    if not measurement_id or measurement_id in blocks:
                        raise AnalysisError(f"{context}: blank or duplicate measurement_id")
                    sample = normalize_sample(row["sample"])
                    expected_sample = member.name.split("/", 1)[0]
                    if sample != expected_sample:
                        raise AnalysisError(f"{context}: sample does not match archive member")
                    blocks[measurement_id] = Block(
                        sample=sample,
                        tissue=normalize_tissue(row["tissue"]),
                        gene_id=normalize_gene_id(row["gene_id"]),
                        gene_name=row["gene_name"].strip(),
                        block_id=row["haplotype_block_id"].strip(),
                        measurement_id=measurement_id,
                        block_start=parse_int(row["block_start"], "block_start", context),
                        block_end=parse_int(row["block_end"], "block_end", context),
                        direction=direction,
                        variants_with_reads=observed,
                        parent_major_alleles={variant: major[variant] for variant in observed},
                    )
    return blocks, dict(counts)


def map_successes_to_blocks(
    site_path: Path,
    recurrent_rows: Mapping[Tuple[str, str, str], Mapping[str, str]],
    success_keys: Set[Tuple[str, str, str, str]],
    blocks: Mapping[str, Block],
) -> Tuple[List[Link], Dict[Tuple[str, str, str, str], int], Dict[str, int]]:
    targets: MutableMapping[Tuple[str, str, str], Set[str]] = defaultdict(set)
    for sample, tissue, gene_id, variant in success_keys:
        targets[(sample, tissue, variant)].add(gene_id)
    links: List[Link] = []
    coverage: Dict[Tuple[str, str, str, str], int] = {}
    mapped: Set[Tuple[str, str, str, str]] = set()
    seen_target_sites: Set[Tuple[str, str, str]] = set()
    counts = Counter()

    for row in iter_tsv(site_path, SITE_REQUIRED):
        counts["n_core_site_rows_scanned"] += 1
        site_key = (
            normalize_sample(row["sample"]),
            normalize_tissue(row["tissue"]),
            row["variant_id"].strip(),
        )
        genes = targets.get(site_key)
        if not genes:
            continue
        if site_key in seen_target_sites:
            raise AnalysisError(f"{site_path}: duplicate target physical site row: {site_key}")
        seen_target_sites.add(site_key)
        counts["n_target_success_physical_site_rows"] += 1
        if not parse_bool(row["site_level_ase_signal"], "site_level_ase_signal", str(site_key)):
            raise AnalysisError(f"{site_path}: recurrent success is not site-level ASE: {site_key}")
        ref_count = parse_int(row["ref_fragment_count"], "ref_fragment_count", str(site_key))
        alt_count = parse_int(row["alt_fragment_count"], "alt_fragment_count", str(site_key))
        total = parse_int(
            row["total_informative_fragments"], "total_informative_fragments", str(site_key)
        )
        if ref_count + alt_count != total:
            raise AnalysisError(f"{site_path}: REF+ALT != total for {site_key}")
        major_allele = row["major_allele"].strip().upper()
        parent_ids = split_values(row["core_parent_measurement_ids"])
        for gene_id in sorted(genes):
            success_key = (site_key[0], site_key[1], gene_id, site_key[2])
            recurrent = recurrent_rows[(gene_id, site_key[1], site_key[2])]
            matches: List[Block] = []
            for measurement_id in parent_ids:
                block = blocks.get(measurement_id)
                if block is None:
                    raise AnalysisError(
                        f"{site_path}: core parent measurement absent from archive: {measurement_id}"
                    )
                if (
                    block.sample == site_key[0]
                    and block.tissue == site_key[1]
                    and block.gene_id == gene_id
                    and site_key[2] in block.variants_with_reads
                    and block.parent_major_alleles.get(site_key[2], "") == major_allele
                ):
                    matches.append(block)
            if not matches:
                raise AnalysisError(
                    "successful gene-specific SNP cannot be linked to a direction-concordant "
                    f"v0.2.6 core parent block: {success_key}"
                )
            coverage[(gene_id, site_key[1], site_key[2], site_key[0])] = total
            mapped.add(success_key)
            chrom = recurrent["chrom"].removeprefix("chr")
            position = parse_int(recurrent["position"], "position", str(success_key))
            for block in sorted(matches, key=lambda value: value.measurement_id):
                links.append(
                    Link(
                        gene_id=gene_id,
                        gene_name=recurrent["gene_name"],
                        tissue=site_key[1],
                        variant_id=site_key[2],
                        chrom=chrom,
                        position=position,
                        sample=site_key[0],
                        measurement_id=block.measurement_id,
                        block_id=block.block_id,
                        block_start=block.block_start,
                        block_end=block.block_end,
                        direction=block.direction,
                        parent_major_allele=block.parent_major_alleles[site_key[2]],
                        site_major_allele=major_allele,
                        ref_count=ref_count,
                        alt_count=alt_count,
                        total_count=total,
                    )
                )
    missing = success_keys - mapped
    if missing:
        raise AnalysisError(
            f"{len(missing)} recurrent success keys were not found in core site input; "
            f"examples={sorted(missing)[:10]}"
        )
    counts["n_recurrent_success_keys"] = len(success_keys)
    counts["n_recurrent_success_parent_block_links"] = len(links)
    counts["n_recurrent_success_keys_with_parent_block"] = len(mapped)
    return links, coverage, dict(counts)


def make_locus_id(
    gene_id: str, tissue: str, chrom: str, start: int, end: int, variants: Sequence[str]
) -> str:
    payload = "|".join((gene_id, tissue, chrom, str(start), str(end), *variants))
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:10]
    safe_tissue = re.sub(r"[^A-Za-z0-9]+", "_", tissue).strip("_")
    return f"ASEL_{gene_id}_{safe_tissue}_{chrom}_{start}_{end}_{digest}"


def lead_sort_key(
    row: Mapping[str, str], median_coverage: float
) -> Tuple[float, float, float, int, float, int, str]:
    context = f"{row['gene_id']}|{row['tissue']}|{row['variant_id']}"
    return (
        parse_float(row["recurrence_q_bh_global"], "recurrence_q_bh_global", context),
        parse_float(
            row["recurrence_p_poisson_binomial"], "recurrence_p_poisson_binomial", context
        ),
        -parse_float(row["recurrence_fraction"], "recurrence_fraction", context),
        -parse_int(
            row["n_animals_with_block_concordant_snp_ase"],
            "n_animals_with_block_concordant_snp_ase",
            context,
        ),
        -median_coverage,
        parse_int(row["position"], "position", context),
        row["variant_id"],
    )


def collapse_loci(
    recurrent_rows: Mapping[Tuple[str, str, str], Dict[str, str]],
    links: Sequence[Link],
    coverage: Mapping[Tuple[str, str, str, str], int],
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    Dict[str, int],
]:
    rows_by_group: MutableMapping[Tuple[str, str], Dict[str, Dict[str, str]]] = defaultdict(dict)
    for (gene_id, tissue, variant), row in recurrent_rows.items():
        rows_by_group[(gene_id, tissue)][variant] = row

    links_by_sample_block: MutableMapping[Tuple[str, str, str, str], Set[str]] = defaultdict(set)
    links_by_hypothesis: MutableMapping[Tuple[str, str, str], List[Link]] = defaultdict(list)
    for link in links:
        links_by_sample_block[(link.gene_id, link.tissue, link.sample, link.measurement_id)].add(
            link.variant_id
        )
        links_by_hypothesis[(link.gene_id, link.tissue, link.variant_id)].append(link)

    edge_evidence: MutableMapping[
        Tuple[str, str, str, str], Set[Tuple[str, str, str]]
    ] = defaultdict(set)
    for (gene_id, tissue, sample, measurement_id), variants in links_by_sample_block.items():
        if len(variants) < 2:
            continue
        block_id = next(
            link.block_id
            for variant in variants
            for link in links_by_hypothesis[(gene_id, tissue, variant)]
            if link.sample == sample and link.measurement_id == measurement_id
        )
        for first, second in combinations(sorted(variants), 2):
            edge_evidence[(gene_id, tissue, first, second)].add(
                (sample, measurement_id, block_id)
            )

    locus_for_hypothesis: Dict[Tuple[str, str, str], str] = {}
    locus_members: Dict[str, List[str]] = {}
    loci: List[Dict[str, object]] = []
    counts = Counter()

    for (gene_id, tissue), grouped_rows in sorted(rows_by_group.items()):
        union_find = UnionFind(grouped_rows)
        for edge_gene, edge_tissue, first, second in edge_evidence:
            if edge_gene == gene_id and edge_tissue == tissue:
                union_find.union(first, second)
        components: MutableMapping[str, List[str]] = defaultdict(list)
        for variant in grouped_rows:
            components[union_find.find(variant)].append(variant)
        ordered_components = sorted(
            components.values(),
            key=lambda variants: min(parse_variant(variant)[1] for variant in variants),
        )
        for variants in ordered_components:
            variants = sorted(variants, key=lambda value: (parse_variant(value)[1], value))
            component_rows = [grouped_rows[variant] for variant in variants]
            chroms = {row["chrom"].removeprefix("chr") for row in component_rows}
            if len(chroms) != 1:
                raise AnalysisError(f"component crosses chromosomes: {gene_id}|{tissue}|{variants}")
            chrom = next(iter(chroms))
            positions = [parse_int(row["position"], "position", gene_id) for row in component_rows]
            start, end = min(positions), max(positions)
            locus_id = make_locus_id(gene_id, tissue, chrom, start, end, variants)
            if locus_id in locus_members:
                raise AnalysisError(f"locus ID collision: {locus_id}")
            locus_members[locus_id] = variants
            for variant in variants:
                locus_for_hypothesis[(gene_id, tissue, variant)] = locus_id

            medians: Dict[str, float] = {}
            success_sets: List[Set[str]] = []
            for row in component_rows:
                samples = {normalize_sample(value) for value in split_values(
                    row["animals_with_block_concordant_snp_ase"]
                )}
                success_sets.append(samples)
                values = [coverage[(gene_id, tissue, row["variant_id"], sample)] for sample in samples]
                medians[row["variant_id"]] = float(statistics.median(values))
            lead = min(
                component_rows,
                key=lambda row: lead_sort_key(row, medians[row["variant_id"]]),
            )
            union_samples = sorted(set().union(*success_sets))
            intersection_samples = sorted(set.intersection(*success_sets))
            component_links = [
                link
                for variant in variants
                for link in links_by_hypothesis[(gene_id, tissue, variant)]
            ]
            component_edges = {
                key: evidence
                for key, evidence in edge_evidence.items()
                if key[0] == gene_id
                and key[1] == tissue
                and key[2] in variants
                and key[3] in variants
            }
            multi_snp_samples = sorted(
                {
                    sample
                    for evidence in component_edges.values()
                    for sample, _, _ in evidence
                }
            )
            loci.append(
                {
                    "locus_id": locus_id,
                    "locus_grouping_method": GROUPING_METHOD,
                    "gene_id": gene_id,
                    "gene_name": lead["gene_name"],
                    "tissue": tissue,
                    "chrom": chrom,
                    "locus_start": start,
                    "locus_end": end,
                    "locus_span_bp": end - start + 1,
                    "n_recurrent_snps": len(variants),
                    "recurrent_snp_ids": ",".join(variants),
                    "lead_snp_id": lead["variant_id"],
                    "lead_snp_position": lead["position"],
                    "lead_recurrence_q_bh_global": lead["recurrence_q_bh_global"],
                    "lead_recurrence_p_poisson_binomial": lead[
                        "recurrence_p_poisson_binomial"
                    ],
                    "lead_recurrence_fraction": lead["recurrence_fraction"],
                    "lead_n_animals_with_block_concordant_snp_ase": lead[
                        "n_animals_with_block_concordant_snp_ase"
                    ],
                    "lead_median_success_total_informative_fragments": fmt_float(
                        medians[lead["variant_id"]]
                    ),
                    "lead_original_ase_assignment_class": lead[
                        "original_ase_assignment_class"
                    ],
                    "lead_canonical_transcript_features": lead[
                        "canonical_transcript_features"
                    ],
                    "n_animals_with_any_recurrent_snp": len(union_samples),
                    "animals_with_any_recurrent_snp": ",".join(union_samples),
                    "n_animals_with_all_recurrent_snps": len(intersection_samples),
                    "animals_with_all_recurrent_snps": ",".join(intersection_samples),
                    "n_success_parent_block_links": len(component_links),
                    "n_distinct_success_sample_blocks": len(
                        {(link.sample, link.measurement_id) for link in component_links}
                    ),
                    "n_direct_core_block_edges": len(component_edges),
                    "n_animals_with_multi_snp_shared_core_block": len(multi_snp_samples),
                    "animals_with_multi_snp_shared_core_block": ",".join(multi_snp_samples),
                    "original_ase_assignment_classes": ",".join(
                        sorted({row["original_ase_assignment_class"] for row in component_rows})
                    ),
                    "canonical_transcript_features": ",".join(
                        sorted(
                            {
                                value
                                for row in component_rows
                                for value in split_values(row["canonical_transcript_features"])
                            }
                        )
                    ),
                    "locus_is_descriptive_collapse_of_significant_snps": 1,
                }
            )
            counts["n_recurrent_ase_loci"] += 1
            counts["n_single_snp_loci"] += int(len(variants) == 1)
            counts["n_multi_snp_loci"] += int(len(variants) > 1)

    membership: List[Dict[str, object]] = []
    locus_by_id = {row["locus_id"]: row for row in loci}
    for key, row in recurrent_rows.items():
        gene_id, tissue, variant = key
        locus_id = locus_for_hypothesis[key]
        locus = locus_by_id[locus_id]
        sample_values = [normalize_sample(value) for value in split_values(
            row["animals_with_block_concordant_snp_ase"]
        )]
        counts_for_success = [coverage[(gene_id, tissue, variant, sample)] for sample in sample_values]
        hypothesis_links = links_by_hypothesis[key]
        augmented = dict(row)
        augmented.update(
            {
                "locus_id": locus_id,
                "locus_grouping_method": GROUPING_METHOD,
                "n_locus_recurrent_snps": locus["n_recurrent_snps"],
                "is_locus_lead_snp": int(variant == locus["lead_snp_id"]),
                "n_success_parent_block_links": len(hypothesis_links),
                "n_distinct_success_parent_blocks": len(
                    {(link.sample, link.measurement_id) for link in hypothesis_links}
                ),
                "success_total_informative_fragments": ",".join(
                    str(value) for value in counts_for_success
                ),
                "median_success_total_informative_fragments": fmt_float(
                    statistics.median(counts_for_success)
                ),
            }
        )
        membership.append(augmented)

    link_rows: List[Dict[str, object]] = []
    for link in links:
        locus_id = locus_for_hypothesis[(link.gene_id, link.tissue, link.variant_id)]
        link_rows.append(
            {
                "locus_id": locus_id,
                "gene_id": link.gene_id,
                "gene_name": link.gene_name,
                "tissue": link.tissue,
                "variant_id": link.variant_id,
                "chrom": link.chrom,
                "position": link.position,
                "sample": link.sample,
                "measurement_id": link.measurement_id,
                "haplotype_block_id": link.block_id,
                "block_start": link.block_start,
                "block_end": link.block_end,
                "v026_ase_direction": link.direction,
                "parent_major_allele_at_snp": link.parent_major_allele,
                "site_major_allele": link.site_major_allele,
                "ref_fragment_count": link.ref_count,
                "alt_fragment_count": link.alt_count,
                "total_informative_fragments": link.total_count,
            }
        )

    edge_rows: List[Dict[str, object]] = []
    for (gene_id, tissue, first, second), evidence in sorted(edge_evidence.items()):
        locus_id = locus_for_hypothesis[(gene_id, tissue, first)]
        if locus_id != locus_for_hypothesis[(gene_id, tissue, second)]:
            raise AnalysisError("direct block edge split across loci")
        first_row, second_row = recurrent_rows[(gene_id, tissue, first)], recurrent_rows[
            (gene_id, tissue, second)
        ]
        animals = sorted({sample for sample, _, _ in evidence})
        edge_rows.append(
            {
                "locus_id": locus_id,
                "gene_id": gene_id,
                "gene_name": first_row["gene_name"],
                "tissue": tissue,
                "variant_id_1": first,
                "position_1": first_row["position"],
                "variant_id_2": second,
                "position_2": second_row["position"],
                "n_shared_success_sample_blocks": len(evidence),
                "n_shared_success_animals": len(animals),
                "shared_success_animals": ",".join(animals),
                "shared_measurement_ids": ",".join(sorted({value[1] for value in evidence})),
                "shared_haplotype_block_ids": ",".join(sorted({value[2] for value in evidence})),
            }
        )

    tissue_rows: List[Dict[str, object]] = []
    tissues = sorted({row["tissue"] for row in membership})
    for tissue in tissues:
        tissue_members = [row for row in membership if row["tissue"] == tissue]
        tissue_loci = [row for row in loci if row["tissue"] == tissue]
        tissue_rows.append(
            {
                "tissue": tissue,
                "n_recurrent_gene_tissue_snp_hypotheses": len(tissue_members),
                "n_recurrent_ase_loci": len(tissue_loci),
                "n_single_snp_loci": sum(int(row["n_recurrent_snps"] == 1) for row in tissue_loci),
                "n_multi_snp_loci": sum(int(row["n_recurrent_snps"] > 1) for row in tissue_loci),
                "n_recurrent_genes": len({row["gene_id"] for row in tissue_members}),
            }
        )

    counts["n_recurrent_gene_tissue_snp_hypotheses"] = len(membership)
    counts["n_locus_lead_snps"] = sum(int(row["is_locus_lead_snp"]) for row in membership)
    counts["n_direct_core_block_edges"] = len(edge_rows)
    counts["n_recurrent_genes"] = len({row["gene_id"] for row in membership})
    counts["n_tissues"] = len(tissue_rows)
    return loci, membership, link_rows, edge_rows, tissue_rows, dict(counts)


def run_analysis(args: argparse.Namespace) -> Dict[str, object]:
    for path in (args.recurrent_annotation, args.core_site_measurements, args.block_archive):
        if not path.is_file() or path.stat().st_size == 0:
            raise AnalysisError(f"missing or empty input: {path}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    recurrent_header, recurrent_rows, success_keys = load_recurrent(args.recurrent_annotation)
    blocks, block_counts = load_core_blocks(args.block_archive)
    links, coverage, mapping_counts = map_successes_to_blocks(
        args.core_site_measurements, recurrent_rows, success_keys, blocks
    )
    loci, membership, link_rows, edge_rows, tissue_rows, collapse_counts = collapse_loci(
        recurrent_rows, links, coverage
    )

    loci.sort(key=lambda row: (row["tissue"], row["gene_id"], row["chrom"], row["locus_start"]))
    membership.sort(
        key=lambda row: (row["tissue"], row["gene_id"], row["chrom"], int(row["position"]))
    )
    link_rows.sort(
        key=lambda row: (
            row["tissue"], row["gene_id"], row["chrom"], int(row["position"]),
            row["sample"], row["measurement_id"],
        )
    )
    edge_rows.sort(
        key=lambda row: (
            row["tissue"], row["gene_id"], int(row["position_1"]), int(row["position_2"])
        )
    )

    membership_fields = [
        "locus_id",
        "locus_grouping_method",
        "n_locus_recurrent_snps",
        "is_locus_lead_snp",
        "n_success_parent_block_links",
        "n_distinct_success_parent_blocks",
        "success_total_informative_fragments",
        "median_success_total_informative_fragments",
        *recurrent_header,
    ]
    write_tsv(args.out_dir / "recurrent_ase_loci.tsv", LOCUS_FIELDS, loci)
    write_tsv(args.out_dir / "recurrent_snp_locus_membership.tsv", membership_fields, membership)
    write_tsv(
        args.out_dir / "recurrent_snp_parent_block_links.tsv.gz", LINK_FIELDS, link_rows
    )
    write_tsv(args.out_dir / "recurrent_snp_block_edges.tsv.gz", EDGE_FIELDS, edge_rows)
    write_tsv(
        args.out_dir / "locus_counts_by_tissue.tsv",
        [
            "tissue",
            "n_recurrent_gene_tissue_snp_hypotheses",
            "n_recurrent_ase_loci",
            "n_single_snp_loci",
            "n_multi_snp_loci",
            "n_recurrent_genes",
        ],
        tissue_rows,
    )

    counts = {**block_counts, **mapping_counts, **collapse_counts}
    if counts["n_locus_lead_snps"] != counts["n_recurrent_ase_loci"]:
        raise AnalysisError("each locus must have exactly one lead SNP")
    if counts["n_recurrent_gene_tissue_snp_hypotheses"] != len(recurrent_rows):
        raise AnalysisError("not all recurrent hypotheses were assigned to a locus")

    qc = {
        "status": "PASS",
        "checks": {
            "all_recurrent_hypotheses_are_significant_input_rows": True,
            "every_recurrent_success_maps_to_direction_concordant_v026_core_parent": True,
            "graph_edges_require_same_gene_tissue_sample_and_core_measurement": True,
            "genomic_distance_not_used_for_grouping": True,
            "every_recurrent_hypothesis_assigned_to_exactly_one_locus": True,
            "every_locus_has_exactly_one_deterministic_lead_snp": True,
            "locus_collapse_does_not_create_new_statistical_test": True,
        },
        "counts": counts,
    }
    metadata = {
        "version": VERSION,
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "unit_input": "statistically recurrent gene-tissue-SNP hypothesis",
        "unit_output": "same-gene-and-tissue core-block-connected recurrent ASE locus",
        "definitions": {
            "locus": (
                "a connected component of recurrent SNPs for one gene and tissue; an edge exists "
                "only if two SNPs occur together in the same direction-concordant final v0.2.6 "
                "core ASE parent measurement in at least one successful animal"
            ),
            "singleton_locus": (
                "a statistically recurrent gene-tissue-SNP with no shared successful core parent "
                "measurement linking it to another recurrent SNP for that gene and tissue"
            ),
            "lead_snp": (
                "deterministically selected within each locus by recurrence q, recurrence p, "
                "higher recurrence fraction, more successful animals, higher median successful-site "
                "coverage, genomic position, then variant ID"
            ),
            "statistical_scope": (
                "descriptive collapse of already significant recurrent SNP hypotheses; no locus-level "
                "p-value or new significance claim is generated"
            ),
        },
        "grouping_method": GROUPING_METHOD,
        "distance_rule": "NONE",
        "inputs": {
            "recurrent_annotation": fingerprint(args.recurrent_annotation),
            "core_site_measurements": fingerprint(args.core_site_measurements),
            "final_v026_block_archive": fingerprint(args.block_archive),
        },
        "outputs": {
            "loci": str((args.out_dir / "recurrent_ase_loci.tsv").resolve()),
            "membership": str((args.out_dir / "recurrent_snp_locus_membership.tsv").resolve()),
            "parent_block_links": str(
                (args.out_dir / "recurrent_snp_parent_block_links.tsv.gz").resolve()
            ),
            "block_edges": str((args.out_dir / "recurrent_snp_block_edges.tsv.gz").resolve()),
            "tissue_counts": str((args.out_dir / "locus_counts_by_tissue.tsv").resolve()),
        },
        "counts": counts,
        "software": {"python": sys.version},
    }
    for name, value in (("analysis_qc.json", qc), ("run_metadata.json", metadata)):
        with (args.out_dir / name).open("w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
    return metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recurrent-annotation", required=True, type=Path)
    parser.add_argument("--core-site-measurements", required=True, type=Path)
    parser.add_argument("--block-archive", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        metadata = run_analysis(args)
    except (AnalysisError, KeyError, OSError, tarfile.TarError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] {VERSION}: {metadata['counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
