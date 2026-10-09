#!/usr/bin/env python3
"""Annotate unique testable exact CpG-K windows without changing ASM calls."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import gzip
import os
import re
import sqlite3
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from asm_annotation_common import AnnotationError as CpGKError, VERSION, open_text, write_json


GENE_BIOTYPE = "protein_coding"
DEFAULT_CANONICAL_TAG = "ensembl_canonical"
PRIMARY_FEATURE_ORDER = (
    "PROMOTER",
    "FIVE_PRIME_UTR",
    "CDS",
    "THREE_PRIME_UTR",
    "INTRON",
    "OTHER_GENE_BODY",
)
DISPLAY_FEATURE_ORDER = (
    "PROMOTER",
    "FIVE_PRIME_UTR",
    "CDS",
    "THREE_PRIME_UTR",
    "EXON",
    "INTRON",
    "TRANSCRIPT_BODY",
)
DISPLAY_INDEX = {value: index for index, value in enumerate(DISPLAY_FEATURE_ORDER)}

MEASUREMENT_REQUIRED = {
    "sample", "tissue", "chrom", "k", "window_id", "first_cpg_index",
    "last_cpg_index", "cpg_positions_1based", "start0", "end0", "span_bp",
    "ps", "individual_ASM_call",
}

GENE_FIELDS = (
    "chrom", "k", "window_id", "first_cpg_index", "last_cpg_index",
    "cpg_positions_1based", "start0", "end0", "span_bp",
    "annotation_gene_count", "multi_gene_overlap", "gene_id", "gene_name",
    "gene_biotype", "gene_strand", "gene_relationship", "gene_body_overlap_bp",
    "canonical_status", "canonical_transcript_ids", "canonical_tss_1based",
    "primary_feature", "canonical_features", "all_transcript_features",
    "canonical_promoter_overlap", "any_transcript_promoter_overlap",
    "include_primary_gene_feature_summary", "has_direct_protein_coding_annotation",
)

TRANSCRIPT_FIELDS = (
    "chrom", "k", "window_id", "cpg_positions_1based", "start0", "end0",
    "gene_id", "gene_name", "gene_strand", "transcript_id", "transcript_name",
    "transcript_biotype", "is_ensembl_canonical", "transcript_tss_1based",
    "features", "promoter_overlap_bp", "transcript_body_overlap_bp",
    "five_prime_utr_overlap_bp", "cds_overlap_bp", "three_prime_utr_overlap_bp",
    "exon_overlap_bp", "intron_overlap_bp",
)


@dataclasses.dataclass(frozen=True)
class Interval:
    start0: int
    end0: int


@dataclasses.dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    transcript_name: str = ""
    transcript_biotype: str = ""
    start0: int | None = None
    end0: int | None = None
    tags: set[str] = dataclasses.field(default_factory=set)
    exons: list[Interval] = dataclasses.field(default_factory=list)
    cds: list[Interval] = dataclasses.field(default_factory=list)
    five_prime_utr: list[Interval] = dataclasses.field(default_factory=list)
    three_prime_utr: list[Interval] = dataclasses.field(default_factory=list)

    def absorb(self, start0: int, end0: int) -> None:
        self.start0 = start0 if self.start0 is None else min(self.start0, start0)
        self.end0 = end0 if self.end0 is None else max(self.end0, end0)

    @property
    def tss0(self) -> int:
        if self.start0 is None or self.end0 is None:
            raise CpGKError(f"transcript has no span: {self.transcript_id}")
        return self.start0 if self.strand == "+" else self.end0 - 1


@dataclasses.dataclass
class Gene:
    gene_id: str
    chrom: str
    strand: str
    gene_name: str = ""
    gene_biotype: str = ""
    start0: int | None = None
    end0: int | None = None
    relation_start0: int | None = None
    relation_end0: int | None = None
    gene_feature_seen: bool = False
    transcripts: dict[str, Transcript] = dataclasses.field(default_factory=dict)

    def absorb(self, start0: int, end0: int) -> None:
        self.start0 = start0 if self.start0 is None else min(self.start0, start0)
        self.end0 = end0 if self.end0 is None else max(self.end0, end0)


@dataclasses.dataclass(frozen=True)
class Window:
    chrom: str
    k: int
    window_id: str
    first_cpg_index: int
    last_cpg_index: int
    cpg_positions_1based: str
    start0: int
    end0: int
    span_bp: int


def normalize_chrom(value: str) -> str:
    return value.strip().removeprefix("chr")


def normalize_tag(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def parse_attributes(text: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for key, value in re.findall(r"([^\s;]+)\s+\"([^\"]*)\"", text):
        result.setdefault(key, []).append(value)
    return result


def first_attribute(attributes: Mapping[str, Sequence[str]], *keys: str) -> str:
    for key in keys:
        values = attributes.get(key)
        if values:
            return values[0]
    return ""


def interval_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2))


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    result: list[Interval] = []
    for current in sorted(set(intervals), key=lambda item: (item.start0, item.end0)):
        if not result or current.start0 > result[-1].end0:
            result.append(current)
        else:
            previous = result[-1]
            result[-1] = Interval(previous.start0, max(previous.end0, current.end0))
    return result


def union_overlap(window: Window, intervals: Iterable[Interval]) -> int:
    return sum(
        interval_overlap(window.start0, window.end0, item.start0, item.end0)
        for item in merge_intervals(intervals)
    )


def ordered_features(features: Iterable[str]) -> list[str]:
    return sorted(
        set(features),
        key=lambda value: (DISPLAY_INDEX.get(value, len(DISPLAY_INDEX)), value),
    )


def promoter_interval(transcript: Transcript, upstream: int, downstream: int) -> Interval:
    if transcript.strand == "+":
        return Interval(max(0, transcript.tss0 - upstream), transcript.tss0 + downstream + 1)
    return Interval(max(0, transcript.tss0 - downstream), transcript.tss0 + upstream + 1)


def reject_system_tmp(path: Path) -> None:
    normalized = str(path.resolve()).replace("\\", "/")
    if normalized == "/tmp" or normalized.startswith("/tmp/"):
        raise CpGKError(f"system /tmp is forbidden: {path}")


def cleanup_sqlite(path: Path) -> None:
    for candidate in (path, Path(str(path) + "-journal"), Path(str(path) + "-wal"), Path(str(path) + "-shm")):
        try:
            candidate.unlink()
        except FileNotFoundError:
            pass


@contextmanager
def atomic_tsv_writer(path: Path, fieldnames: Sequence[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(partial, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            yield writer
        os.replace(partial, path)
    finally:
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


def validate_measurement_header(fieldnames: Sequence[str] | None, path: Path) -> None:
    if fieldnames is None:
        raise CpGKError(f"empty measurement input: {path}")
    missing = sorted(MEASUREMENT_REQUIRED.difference(fieldnames))
    if missing:
        raise CpGKError(f"{path}: missing measurement fields: {','.join(missing)}")


def load_unique_windows(measurements: Path, connection: sqlite3.Connection) -> dict[str, object]:
    connection.execute(
        """
        CREATE TABLE windows (
          chrom TEXT NOT NULL, k INTEGER NOT NULL, window_id TEXT NOT NULL,
          first_cpg_index INTEGER NOT NULL, last_cpg_index INTEGER NOT NULL,
          cpg_positions_1based TEXT NOT NULL, start0 INTEGER NOT NULL,
          end0 INTEGER NOT NULL, span_bp INTEGER NOT NULL, conflict INTEGER NOT NULL DEFAULT 0,
          PRIMARY KEY (chrom, k, window_id)
        )
        """
    )
    sql = """
        INSERT INTO windows VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 0)
        ON CONFLICT(chrom, k, window_id) DO UPDATE SET conflict = CASE WHEN
          windows.first_cpg_index != excluded.first_cpg_index OR
          windows.last_cpg_index != excluded.last_cpg_index OR
          windows.cpg_positions_1based != excluded.cpg_positions_1based OR
          windows.start0 != excluded.start0 OR windows.end0 != excluded.end0 OR
          windows.span_bp != excluded.span_bp THEN 1 ELSE windows.conflict END
    """
    measurement_rows = 0
    asm_rows = 0
    units: set[tuple[str, str]] = set()
    chroms: set[str] = set()
    batch: list[tuple[object, ...]] = []
    with open_text(measurements) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_measurement_header(reader.fieldnames, measurements)
        for raw in reader:
            measurement_rows += 1
            sample = raw["sample"].strip()
            tissue = raw["tissue"].strip()
            if not sample or not tissue:
                raise CpGKError(f"blank sample/tissue at measurement row {measurement_rows}")
            chrom = normalize_chrom(raw["chrom"])
            k = int(raw["k"])
            first = int(raw["first_cpg_index"])
            last = int(raw["last_cpg_index"])
            positions = [int(value) for value in raw["cpg_positions_1based"].split(",")]
            start0, end0, span = int(raw["start0"]), int(raw["end0"]), int(raw["span_bp"])
            if k not in {3, 4, 5} or last - first + 1 != k or len(positions) != k:
                raise CpGKError(f"invalid exact CpG-K contract: {raw['window_id']}")
            if start0 != positions[0] - 1 or end0 != positions[-1] + 1 or span != end0 - start0:
                raise CpGKError(f"invalid exact-window coordinates: {raw['window_id']}")
            call = int(raw["individual_ASM_call"])
            if call not in {0, 1}:
                raise CpGKError(f"invalid individual_ASM_call: {raw['window_id']}")
            asm_rows += call
            units.add((sample, tissue))
            chroms.add(chrom)
            batch.append((chrom, k, raw["window_id"], first, last, raw["cpg_positions_1based"], start0, end0, span))
            if len(batch) >= 10000:
                connection.executemany(sql, batch)
                batch.clear()
        if batch:
            connection.executemany(sql, batch)
    connection.commit()
    conflict = connection.execute("SELECT window_id FROM windows WHERE conflict=1 LIMIT 1").fetchone()
    if conflict:
        raise CpGKError(f"same exact window has inconsistent coordinates: {conflict[0]}")
    unique_windows = int(connection.execute("SELECT COUNT(*) FROM windows").fetchone()[0])
    if measurement_rows == 0 or unique_windows == 0:
        raise CpGKError(f"no testable exact-window measurements in {measurements}")
    return {
        "measurement_rows": measurement_rows,
        "individual_ASM_rows": asm_rows,
        "unique_windows": unique_windows,
        "units": sorted(f"{sample}|{tissue}" for sample, tissue in units),
        "chromosomes": sorted(chroms, key=lambda value: int(value) if value.isdigit() else 10**9),
    }


def load_gtf(
    path: Path,
    requested_chroms: set[str],
    canonical_tag: str,
    promoter_upstream: int,
    promoter_downstream: int,
) -> tuple[dict[str, list[Gene]], dict[str, object]]:
    genes: dict[tuple[str, str], Gene] = {}
    malformed = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                malformed += 1
                continue
            chrom, _source, feature, start_text, end_text, _score, strand, _frame, attribute_text = fields
            chrom = normalize_chrom(chrom)
            if chrom not in requested_chroms or strand not in {"+", "-"}:
                continue
            attributes = parse_attributes(attribute_text)
            gene_id = first_attribute(attributes, "gene_id")
            if not gene_id:
                malformed += 1
                continue
            try:
                start0, end0 = int(start_text) - 1, int(end_text)
            except ValueError:
                malformed += 1
                continue
            if start0 < 0 or end0 <= start0:
                malformed += 1
                continue
            key = (chrom, gene_id)
            gene = genes.get(key)
            if gene is None:
                gene = Gene(gene_id=gene_id, chrom=chrom, strand=strand)
                genes[key] = gene
            elif gene.strand != strand:
                raise CpGKError(f"inconsistent strand for {gene_id} at GTF line {line_number}")
            gene.absorb(start0, end0)
            gene.gene_name = gene.gene_name or first_attribute(attributes, "gene_name")
            feature_lower = feature.lower()
            if feature_lower == "gene":
                gene.gene_feature_seen = True
                gene.gene_biotype = first_attribute(attributes, "gene_biotype", "gene_type")
            transcript_id = first_attribute(attributes, "transcript_id")
            if not transcript_id:
                continue
            transcript = gene.transcripts.get(transcript_id)
            if transcript is None:
                transcript = Transcript(transcript_id, gene_id, chrom, strand)
                gene.transcripts[transcript_id] = transcript
            transcript.absorb(start0, end0)
            transcript.transcript_name = transcript.transcript_name or first_attribute(attributes, "transcript_name")
            transcript.transcript_biotype = transcript.transcript_biotype or first_attribute(
                attributes, "transcript_biotype", "transcript_type"
            )
            transcript.tags.update(attributes.get("tag", ()))
            interval = Interval(start0, end0)
            if feature_lower == "exon":
                transcript.exons.append(interval)
            elif feature_lower == "cds":
                transcript.cds.append(interval)
            elif feature_lower in {"five_prime_utr", "5utr", "5_prime_utr"}:
                transcript.five_prime_utr.append(interval)
            elif feature_lower in {"three_prime_utr", "3utr", "3_prime_utr"}:
                transcript.three_prime_utr.append(interval)

    normalized_canonical = normalize_tag(canonical_tag)
    by_chrom: dict[str, list[Gene]] = {}
    canonical_transcripts = 0
    canonical_genes = 0
    transcript_count = 0
    for gene in genes.values():
        if not gene.gene_feature_seen or gene.gene_biotype != GENE_BIOTYPE:
            continue
        valid: dict[str, Transcript] = {}
        for transcript_id, transcript in gene.transcripts.items():
            if transcript.start0 is None or transcript.end0 is None:
                continue
            transcript.exons = merge_intervals(transcript.exons)
            transcript.cds = merge_intervals(transcript.cds)
            transcript.five_prime_utr = merge_intervals(transcript.five_prime_utr)
            transcript.three_prime_utr = merge_intervals(transcript.three_prime_utr)
            valid[transcript_id] = transcript
        gene.transcripts = valid
        transcript_count += len(valid)
        canonical = [
            transcript for transcript in valid.values()
            if normalized_canonical in {normalize_tag(tag) for tag in transcript.tags}
        ]
        canonical_transcripts += len(canonical)
        canonical_genes += int(bool(canonical))
        if gene.start0 is None or gene.end0 is None:
            continue
        starts, ends = [gene.start0], [gene.end0]
        for transcript in valid.values():
            promoter = promoter_interval(transcript, promoter_upstream, promoter_downstream)
            starts.append(promoter.start0)
            ends.append(promoter.end0)
        gene.relation_start0, gene.relation_end0 = min(starts), max(ends)
        by_chrom.setdefault(gene.chrom, []).append(gene)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda gene: (gene.relation_start0 or 0, gene.relation_end0 or 0, gene.gene_id))
    if not by_chrom:
        raise CpGKError(f"no {GENE_BIOTYPE} genes loaded from {path}")
    return by_chrom, {
        "malformed_rows": malformed,
        "genes_loaded": sum(len(values) for values in by_chrom.values()),
        "transcripts_loaded": transcript_count,
        "genes_with_ensembl_canonical": canonical_genes,
        "ensembl_canonical_transcripts": canonical_transcripts,
        "gene_biotype_filter": GENE_BIOTYPE,
        "requested_chromosomes": sorted(requested_chroms),
    }


def transcript_relation(
    window: Window,
    gene: Gene,
    transcript: Transcript,
    canonical_tag: str,
    promoter_upstream: int,
    promoter_downstream: int,
) -> dict[str, object] | None:
    if transcript.start0 is None or transcript.end0 is None:
        return None
    promoter = promoter_interval(transcript, promoter_upstream, promoter_downstream)
    promoter_bp = interval_overlap(window.start0, window.end0, promoter.start0, promoter.end0)
    body_bp = interval_overlap(window.start0, window.end0, transcript.start0, transcript.end0)
    if promoter_bp == 0 and body_bp == 0:
        return None
    exon_bp = union_overlap(window, transcript.exons)
    cds_bp = union_overlap(window, transcript.cds)
    five_bp = union_overlap(window, transcript.five_prime_utr)
    three_bp = union_overlap(window, transcript.three_prime_utr)
    intron_bp = max(0, body_bp - exon_bp) if transcript.exons else 0
    features: list[str] = []
    for feature, overlap in (
        ("PROMOTER", promoter_bp),
        ("FIVE_PRIME_UTR", five_bp),
        ("CDS", cds_bp),
        ("THREE_PRIME_UTR", three_bp),
        ("EXON", exon_bp),
        ("INTRON", intron_bp),
    ):
        if overlap > 0:
            features.append(feature)
    if body_bp > 0 and exon_bp == 0 and intron_bp == 0:
        features.append("TRANSCRIPT_BODY")
    normalized_canonical = normalize_tag(canonical_tag)
    is_canonical = normalized_canonical in {normalize_tag(tag) for tag in transcript.tags}
    return {
        "chrom": window.chrom,
        "k": window.k,
        "window_id": window.window_id,
        "cpg_positions_1based": window.cpg_positions_1based,
        "start0": window.start0,
        "end0": window.end0,
        "gene_id": gene.gene_id,
        "gene_name": gene.gene_name or "NA",
        "gene_strand": gene.strand,
        "transcript_id": transcript.transcript_id,
        "transcript_name": transcript.transcript_name or "NA",
        "transcript_biotype": transcript.transcript_biotype or "NA",
        "is_ensembl_canonical": int(is_canonical),
        "transcript_tss_1based": transcript.tss0 + 1,
        "features": ";".join(ordered_features(features)) or "NA",
        "promoter_overlap_bp": promoter_bp,
        "transcript_body_overlap_bp": body_bp,
        "five_prime_utr_overlap_bp": five_bp,
        "cds_overlap_bp": cds_bp,
        "three_prime_utr_overlap_bp": three_bp,
        "exon_overlap_bp": exon_bp,
        "intron_overlap_bp": intron_bp,
    }


def primary_feature(features: set[str], gene_body_overlap_bp: int) -> str:
    for feature in PRIMARY_FEATURE_ORDER[:-1]:
        if feature in features:
            return feature
    if gene_body_overlap_bp > 0 or features.intersection({"EXON", "TRANSCRIPT_BODY"}):
        return "OTHER_GENE_BODY"
    return ""


def annotate_gene(
    window: Window,
    gene: Gene,
    canonical_tag: str,
    promoter_upstream: int,
    promoter_downstream: int,
) -> tuple[dict[str, object] | None, list[dict[str, object]]]:
    tx_rows = [
        row for transcript in sorted(gene.transcripts.values(), key=lambda value: value.transcript_id)
        if (row := transcript_relation(
            window, gene, transcript, canonical_tag, promoter_upstream, promoter_downstream
        )) is not None
    ]
    gene_body_bp = interval_overlap(
        window.start0, window.end0, gene.start0 or 0, gene.end0 or 0
    )
    if not tx_rows and gene_body_bp == 0:
        return None, []
    normalized = normalize_tag(canonical_tag)
    canonical = sorted(
        [
            transcript for transcript in gene.transcripts.values()
            if normalized in {normalize_tag(tag) for tag in transcript.tags}
        ],
        key=lambda value: value.transcript_id,
    )
    canonical_ids = {transcript.transcript_id for transcript in canonical}
    canonical_rows = [row for row in tx_rows if str(row["transcript_id"]) in canonical_ids]
    alternative_rows = [row for row in tx_rows if str(row["transcript_id"]) not in canonical_ids]
    canonical_features = {
        feature for row in canonical_rows for feature in str(row["features"]).split(";")
        if feature and feature != "NA"
    }
    related_features = {
        feature for row in tx_rows for feature in str(row["features"]).split(";")
        if feature and feature != "NA"
    }
    chosen = primary_feature(canonical_features, gene_body_bp)
    if not canonical:
        canonical_status = "NO_ENSEMBL_CANONICAL"
        relationship = "NO_ENSEMBL_CANONICAL"
        include = 0
    else:
        canonical_status = "SINGLE_ENSEMBL_CANONICAL" if len(canonical) == 1 else "MULTIPLE_ENSEMBL_CANONICAL"
        if canonical_rows:
            relationship = "CANONICAL_TRANSCRIPT"
            include = int(bool(chosen))
        elif gene_body_bp > 0:
            relationship = "GENE_BODY_ONLY"
            chosen = "OTHER_GENE_BODY"
            include = 1
        elif alternative_rows:
            relationship = "ALTERNATIVE_ISOFORM_ONLY"
            include = 0
        else:
            return None, []
    if not chosen:
        chosen = "ALTERNATIVE_TRANSCRIPT_RELATION" if alternative_rows else "OTHER_GENE_BODY"
    row = {
        "chrom": window.chrom,
        "k": window.k,
        "window_id": window.window_id,
        "first_cpg_index": window.first_cpg_index,
        "last_cpg_index": window.last_cpg_index,
        "cpg_positions_1based": window.cpg_positions_1based,
        "start0": window.start0,
        "end0": window.end0,
        "span_bp": window.span_bp,
        "gene_id": gene.gene_id,
        "gene_name": gene.gene_name or "NA",
        "gene_biotype": gene.gene_biotype,
        "gene_strand": gene.strand,
        "gene_relationship": relationship,
        "gene_body_overlap_bp": gene_body_bp,
        "canonical_status": canonical_status,
        "canonical_transcript_ids": ";".join(transcript.transcript_id for transcript in canonical) or "NA",
        "canonical_tss_1based": ";".join(str(transcript.tss0 + 1) for transcript in canonical) or "NA",
        "primary_feature": chosen,
        "canonical_features": ";".join(ordered_features(canonical_features)) or "NA",
        "all_transcript_features": ";".join(ordered_features(related_features)) or "NA",
        "canonical_promoter_overlap": int("PROMOTER" in canonical_features),
        "any_transcript_promoter_overlap": int("PROMOTER" in related_features),
        "include_primary_gene_feature_summary": include,
        "has_direct_protein_coding_annotation": 1,
    }
    return row, tx_rows


def unannotated_row(window: Window) -> dict[str, object]:
    return {
        "chrom": window.chrom, "k": window.k, "window_id": window.window_id,
        "first_cpg_index": window.first_cpg_index, "last_cpg_index": window.last_cpg_index,
        "cpg_positions_1based": window.cpg_positions_1based,
        "start0": window.start0, "end0": window.end0, "span_bp": window.span_bp,
        "annotation_gene_count": 0, "multi_gene_overlap": 0,
        "gene_id": "NA", "gene_name": "NA", "gene_biotype": "NA", "gene_strand": "NA",
        "gene_relationship": "NO_DIRECT_PROTEIN_CODING_ANNOTATION",
        "gene_body_overlap_bp": 0, "canonical_status": "NA",
        "canonical_transcript_ids": "NA", "canonical_tss_1based": "NA",
        "primary_feature": "NO_DIRECT_PROTEIN_CODING_ANNOTATION",
        "canonical_features": "NA", "all_transcript_features": "NA",
        "canonical_promoter_overlap": 0, "any_transcript_promoter_overlap": 0,
        "include_primary_gene_feature_summary": 0,
        "has_direct_protein_coding_annotation": 0,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measurements", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--work-db", required=True)
    parser.add_argument("--canonical-tag", default=DEFAULT_CANONICAL_TAG)
    parser.add_argument("--promoter-upstream", type=int, default=2000)
    parser.add_argument("--promoter-downstream", type=int, default=500)
    parser.add_argument("--keep-work-db", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    measurements = Path(args.measurements)
    gtf = Path(args.gtf)
    output_root = Path(args.output_root)
    work_db = Path(args.work_db)
    if not measurements.is_file() or not gtf.is_file():
        raise CpGKError("measurements and GTF must exist")
    if args.promoter_upstream < 0 or args.promoter_downstream < 0:
        raise CpGKError("promoter boundaries must be non-negative")
    reject_system_tmp(work_db)
    work_db.parent.mkdir(parents=True, exist_ok=True)
    cleanup_sqlite(work_db)
    connection = sqlite3.connect(work_db)
    success = False
    try:
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        input_qc = load_unique_windows(measurements, connection)
        genes_by_chrom, gtf_qc = load_gtf(
            gtf, set(input_qc["chromosomes"]), args.canonical_tag,
            args.promoter_upstream, args.promoter_downstream,
        )
        gene_path = output_root / "exact_window_gene_annotations.tsv.gz"
        transcript_path = output_root / "exact_window_transcript_annotations.tsv.gz"
        qc = {
            "unique_windows": 0, "gene_annotation_rows": 0, "transcript_annotation_rows": 0,
            "windows_with_direct_protein_coding_annotation": 0,
            "windows_without_direct_protein_coding_annotation": 0,
            "windows_with_multiple_gene_annotations": 0,
            "primary_summary_gene_rows": 0,
        }
        with atomic_tsv_writer(gene_path, GENE_FIELDS) as gene_writer, atomic_tsv_writer(
            transcript_path, TRANSCRIPT_FIELDS
        ) as transcript_writer:
            for chrom in input_qc["chromosomes"]:
                genes = genes_by_chrom.get(chrom, [])
                active: list[Gene] = []
                gene_index = 0
                cursor = connection.execute(
                    """
                    SELECT chrom,k,window_id,first_cpg_index,last_cpg_index,
                           cpg_positions_1based,start0,end0,span_bp
                    FROM windows WHERE chrom=? ORDER BY start0,end0,k,window_id
                    """,
                    (chrom,),
                )
                for values in cursor:
                    window = Window(*values)
                    qc["unique_windows"] += 1
                    active = [gene for gene in active if (gene.relation_end0 or 0) > window.start0]
                    while gene_index < len(genes) and (genes[gene_index].relation_start0 or 0) < window.end0:
                        gene = genes[gene_index]
                        if (gene.relation_end0 or 0) > window.start0:
                            active.append(gene)
                        gene_index += 1
                    gene_rows: list[dict[str, object]] = []
                    transcript_rows: list[dict[str, object]] = []
                    for gene in sorted(active, key=lambda value: value.gene_id):
                        gene_row, current_transcripts = annotate_gene(
                            window, gene, args.canonical_tag,
                            args.promoter_upstream, args.promoter_downstream,
                        )
                        if gene_row is not None:
                            gene_rows.append(gene_row)
                            transcript_rows.extend(current_transcripts)
                    if not gene_rows:
                        gene_writer.writerow(unannotated_row(window))
                        qc["windows_without_direct_protein_coding_annotation"] += 1
                        qc["gene_annotation_rows"] += 1
                        continue
                    count = len(gene_rows)
                    qc["windows_with_direct_protein_coding_annotation"] += 1
                    qc["windows_with_multiple_gene_annotations"] += int(count > 1)
                    for row in gene_rows:
                        row["annotation_gene_count"] = count
                        row["multi_gene_overlap"] = int(count > 1)
                        gene_writer.writerow(row)
                        qc["gene_annotation_rows"] += 1
                        qc["primary_summary_gene_rows"] += int(row["include_primary_gene_feature_summary"])
                    for row in transcript_rows:
                        transcript_writer.writerow(row)
                        qc["transcript_annotation_rows"] += 1
        if qc["unique_windows"] != input_qc["unique_windows"]:
            raise CpGKError("annotation did not visit every unique exact window")
        metadata_path = output_root / "exact_window_annotation_metadata.json"
        write_json(
            metadata_path,
            {
                "program": "annotate_exact_cpgk_windows",
                "version": VERSION,
                "status": "PASS",
                "inputs": {
                    "measurements": str(measurements),
                    "measurements_size_bytes": measurements.stat().st_size,
                    "gtf": str(gtf),
                    "gtf_size_bytes": gtf.stat().st_size,
                },
                "annotation_contract": {
                    "unit": "unique testable reference-consecutive exact CpG-K window",
                    "gene_biotype_filter": GENE_BIOTYPE,
                    "gene_assignment": "direct gene-body or transcript-promoter overlap; no nearest-gene assignment",
                    "canonical_tag": args.canonical_tag,
                    "canonical_fallback": "none",
                    "promoter_relative_to_TSS_inclusive": [-args.promoter_upstream, args.promoter_downstream],
                    "primary_feature_hierarchy": list(PRIMARY_FEATURE_ORDER),
                    "all_overlapping_features_retained": True,
                    "alternative_transcripts_retained_in_relational_output": True,
                    "no_direct_annotation_label": "NO_DIRECT_PROTEIN_CODING_ANNOTATION",
                },
                "input_qc": input_qc,
                "gtf_qc": gtf_qc,
                "annotation_qc": qc,
                "outputs": {"gene": str(gene_path), "transcript": str(transcript_path)},
            },
        )
        success = True
        print(
            f"[PASS] exact-window annotation: windows={qc['unique_windows']:,}; "
            f"gene rows={qc['gene_annotation_rows']:,}; no-direct={qc['windows_without_direct_protein_coding_annotation']:,}",
            file=sys.stderr,
        )
        return 0
    finally:
        connection.close()
        if success and not args.keep_work_db:
            cleanup_sqlite(work_db)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CpGKError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(2)
