#!/usr/bin/env python3
"""Annotate ASE SNPs with the exact v0.2.5 block-level GTF rules.

This is deliberately not a consequence predictor.  It reproduces the original
ASE annotation model:

* select genes whose GTF gene_biotype is protein_coding;
* build each gene's exon union from protein-coding (or unspecified) transcripts;
* classify sites as UNIQUE_EXONIC, AMBIGUOUS_SHARED_EXON, INTRONIC_ONLY, or
  INTERGENIC across those selected genes;
* classify every transcript position with the original v0.2.5 precedence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import json
import math
import re
import sys
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "ase_snp_feature_annotation_v0.1.0"
SELECTED_GENE_BIOTYPE = "protein_coding"

INPUT_REQUIRED = {
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_animals_testable",
    "n_animals_with_block_concordant_snp_ase",
    "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
}

ANNOTATION_FIELDS = [
    "gtf_gene_name",
    "gtf_gene_biotype",
    "gtf_gene_contig",
    "gtf_gene_start",
    "gtf_gene_end",
    "gtf_gene_strand",
    "original_ase_assignment_class",
    "site_gene_assignment_requirement",
    "target_gene_in_protein_coding_exon_union",
    "n_exonic_protein_coding_genes_at_site",
    "exonic_protein_coding_gene_ids_at_site",
    "gene_span_ids_at_site",
    "n_transcripts_for_gene",
    "n_exon_union_transcripts_for_gene",
    "n_canonical_transcripts_for_gene",
    "canonical_transcript_ids",
    "canonical_transcript_features",
    "all_transcript_features",
    "exon_union_transcript_features",
    "n_transcripts_exonic_at_site",
    "n_exon_union_transcripts_exonic_at_site",
    "annotation_status",
]

TRANSCRIPT_HEADER = [
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_animals_testable",
    "n_animals_with_block_concordant_snp_ase",
    "recurrence_q_bh_global",
    "transcript_id",
    "transcript_name",
    "transcript_biotype",
    "canonical_transcript",
    "included_in_original_exon_union",
    "transcript_feature",
]

COUNT_HEADER = ["dataset", "dimension", "category", "n_gene_tissue_snp_hypotheses"]

EXONIC_TRANSCRIPT_FEATURES = {
    "CDS",
    "FIVE_PRIME_UTR",
    "THREE_PRIME_UTR",
    "START_CODON",
    "STOP_CODON",
    "EXON_NONCODING_OR_UNSPECIFIED",
}


class AnalysisError(RuntimeError):
    pass


@contextmanager
def open_text(path: Path, mode: str = "rt"):
    if "b" in mode:
        raise ValueError("open_text supports text mode only")
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def input_header(path: Path) -> List[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise AnalysisError(f"empty input: {path}") from error
    missing = sorted(INPUT_REQUIRED - set(header))
    if missing:
        raise AnalysisError(f"{path}: missing required columns: {', '.join(missing)}")
    return header


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(INPUT_REQUIRED - fields)
        if missing:
            raise AnalysisError(f"{path}: missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value if value is not None else "" for key, value in row.items()}


def normalize_gene_id(value: str) -> str:
    return value.strip().split(".")[0]


def parse_int(value: str, field: str, context: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{context}: invalid integer in {field}: {value!r}") from error


def parse_gtf_attributes(text: str) -> Dict[str, object]:
    """Exact parser used by the v0.2.5 block-level ASE caller."""
    attributes: Dict[str, object] = {}
    for key, quoted, bare in re.findall(
        r'(?:^|;\s*)([^\s;]+)\s+(?:"([^"]*)"|([^;\s]+))', text
    ):
        value = quoted if quoted != "" else bare
        if key in attributes:
            current = attributes[key]
            if not isinstance(current, list):
                current = [current]
                attributes[key] = current
            current.append(value)
        else:
            attributes[key] = value
    return attributes


def attribute_first(attributes: Mapping[str, object], *keys: str) -> str:
    for key in keys:
        value = attributes.get(key)
        if isinstance(value, list):
            return str(value[0]) if value else ""
        if value:
            return str(value)
    return ""


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge 1-based inclusive intervals exactly as in v0.2.5."""
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(interval) for interval in merged]


def load_gtf_annotation(
    path: Path,
) -> Tuple[Dict[str, Dict[str, object]], Dict[str, List[Dict[str, object]]], Dict[str, int]]:
    """Reproduce v0.2.5 GTF selection and exon-union construction."""
    genes: Dict[str, Dict[str, object]] = {}
    transcripts: Dict[str, Dict[str, object]] = {}
    counts = Counter()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                counts["n_malformed_gtf_rows"] += 1
                continue
            contig, source, feature = fields[0], fields[1], fields[2]
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise AnalysisError(f"{path}:{line_number}: invalid GTF coordinates") from error
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])
            gene_id = normalize_gene_id(attribute_first(attributes, "gene_id"))
            if not gene_id:
                continue
            gene_biotype = attribute_first(
                attributes, "gene_biotype", "gene_type", "biotype"
            )
            if gene_biotype != SELECTED_GENE_BIOTYPE:
                continue
            counts["n_selected_gtf_rows"] += 1

            if feature == "gene":
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": attribute_first(attributes, "gene_name") or gene_id,
                    "gene_biotype": gene_biotype,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source": source,
                }

            transcript_id = attribute_first(attributes, "transcript_id")
            if not transcript_id:
                continue
            transcript = transcripts.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "transcript_name": attribute_first(attributes, "transcript_name")
                    or transcript_id,
                    "transcript_biotype": attribute_first(
                        attributes, "transcript_biotype", "transcript_type"
                    ),
                    "gene_id": gene_id,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "canonical": False,
                    "tags": set(),
                    "features": defaultdict(list),
                },
            )
            transcript["start"] = min(int(transcript["start"]), start)
            transcript["end"] = max(int(transcript["end"]), end)
            tags = attributes.get("tag", [])
            if not isinstance(tags, list):
                tags = [tags] if tags else []
            transcript["tags"].update(tags)  # type: ignore[union-attr]
            if "Ensembl_canonical" in transcript["tags"]:  # type: ignore[operator]
                transcript["canonical"] = True
            if feature != "transcript":
                transcript["features"][feature].append((start, end))  # type: ignore[index]

    transcripts = {
        transcript_id: transcript
        for transcript_id, transcript in transcripts.items()
        if transcript["gene_id"] in genes
    }
    transcripts_by_gene: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for transcript in transcripts.values():
        transcripts_by_gene[str(transcript["gene_id"])].append(transcript)
    for gene_id in transcripts_by_gene:
        transcripts_by_gene[gene_id].sort(
            key=lambda item: (not bool(item["canonical"]), str(item["transcript_id"]))
        )
        exon_intervals: List[Tuple[int, int]] = []
        for transcript in transcripts_by_gene[gene_id]:
            transcript_biotype = str(transcript["transcript_biotype"])
            if transcript_biotype and transcript_biotype != SELECTED_GENE_BIOTYPE:
                continue
            exon_intervals.extend(transcript["features"].get("exon", []))  # type: ignore[union-attr]
        genes[gene_id]["exon_intervals"] = merge_intervals(exon_intervals)

    counts["n_selected_genes"] = len(genes)
    counts["n_selected_gene_transcripts"] = len(transcripts)
    counts["n_genes_with_exon_union"] = sum(bool(gene.get("exon_intervals")) for gene in genes.values())
    return genes, dict(transcripts_by_gene), dict(counts)


def classify_transcript_position(transcript: Mapping[str, object], position: int) -> str:
    """Exact transcript-feature precedence used by v0.2.5."""
    if position < int(transcript["start"]) or position > int(transcript["end"]):
        return "OUTSIDE_TRANSCRIPT"
    precedence = [
        ("CDS", "CDS"),
        ("five_prime_utr", "FIVE_PRIME_UTR"),
        ("three_prime_utr", "THREE_PRIME_UTR"),
        ("start_codon", "START_CODON"),
        ("stop_codon", "STOP_CODON"),
        ("exon", "EXON_NONCODING_OR_UNSPECIFIED"),
    ]
    features = transcript["features"]
    for feature_name, label in precedence:
        for start, end in features.get(feature_name, []):  # type: ignore[union-attr]
            if int(start) <= position <= int(end):
                return label
    return "INTRON"


def collect_variants(
    path: Path,
) -> Tuple[Dict[str, Tuple[str, int, str, str]], Dict[str, Set[str]], int]:
    variant_info: Dict[str, Tuple[str, int, str, str]] = {}
    target_genes: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    for row in iter_tsv(path):
        n_rows += 1
        variant_id = row["variant_id"]
        context = f"{row['gene_id']}|{row['tissue']}|{variant_id}"
        info = (
            row["chrom"],
            parse_int(row["position"], "position", context),
            row["ref"],
            row["alt"],
        )
        if variant_id in variant_info and variant_info[variant_id] != info:
            raise AnalysisError(f"inconsistent variant metadata for {variant_id}")
        variant_info[variant_id] = info
        target_genes[variant_id].add(normalize_gene_id(row["gene_id"]))
    return variant_info, dict(target_genes), n_rows


def map_variants_to_selected_genes(
    variant_info: Mapping[str, Tuple[str, int, str, str]],
    genes: Mapping[str, Mapping[str, object]],
) -> Dict[str, Dict[str, Tuple[str, ...]]]:
    """Sweep relevant positions across gene spans and v0.2.5 exon unions."""
    variants_by_contig_position: Dict[str, Dict[int, List[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for variant_id, (contig, position, _, _) in variant_info.items():
        variants_by_contig_position[contig][position].append(variant_id)

    genes_by_contig: Dict[str, List[Mapping[str, object]]] = defaultdict(list)
    exons_by_contig: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for gene in genes.values():
        contig = str(gene["contig"])
        genes_by_contig[contig].append(gene)
        for start, end in gene.get("exon_intervals", []):  # type: ignore[assignment]
            exons_by_contig[contig].append((int(start), int(end), str(gene["gene_id"])))
    for contig in genes_by_contig:
        genes_by_contig[contig].sort(key=lambda gene: (int(gene["start"]), int(gene["end"])))
    for contig in exons_by_contig:
        exons_by_contig[contig].sort()

    span: Dict[str, Tuple[str, ...]] = {}
    exonic: Dict[str, Tuple[str, ...]] = {}
    for contig, positions in variants_by_contig_position.items():
        contig_genes = genes_by_contig.get(contig, [])
        active_genes: Dict[str, Mapping[str, object]] = {}
        gene_end_heap: List[Tuple[int, str]] = []
        gene_index = 0
        for position in sorted(positions):
            while gene_index < len(contig_genes) and int(contig_genes[gene_index]["start"]) <= position:
                gene = contig_genes[gene_index]
                gene_id = str(gene["gene_id"])
                active_genes[gene_id] = gene
                heapq.heappush(gene_end_heap, (int(gene["end"]), gene_id))
                gene_index += 1
            while gene_end_heap and gene_end_heap[0][0] < position:
                end, gene_id = heapq.heappop(gene_end_heap)
                current = active_genes.get(gene_id)
                if current is not None and int(current["end"]) == end:
                    active_genes.pop(gene_id, None)
            overlapping = tuple(sorted(active_genes))
            for variant_id in positions[position]:
                span[variant_id] = overlapping

        intervals = exons_by_contig.get(contig, [])
        active_exons: Dict[Tuple[str, int, int], str] = {}
        exon_end_heap: List[Tuple[int, str, int]] = []
        interval_index = 0
        for position in sorted(positions):
            while interval_index < len(intervals) and intervals[interval_index][0] <= position:
                start, end, gene_id = intervals[interval_index]
                active_exons[(gene_id, start, end)] = gene_id
                heapq.heappush(exon_end_heap, (end, gene_id, start))
                interval_index += 1
            while exon_end_heap and exon_end_heap[0][0] < position:
                end, gene_id, start = heapq.heappop(exon_end_heap)
                active_exons.pop((gene_id, start, end), None)
            overlapping = tuple(sorted(set(active_exons.values())))
            for variant_id in positions[position]:
                exonic[variant_id] = overlapping

    result: Dict[str, Dict[str, Tuple[str, ...]]] = {}
    for variant_id in variant_info:
        result[variant_id] = {
            "exonic": exonic.get(variant_id, ()),
            "span": span.get(variant_id, ()),
        }
    return result


def is_exon_union_transcript(transcript: Mapping[str, object]) -> bool:
    biotype = str(transcript["transcript_biotype"])
    return not biotype or biotype == SELECTED_GENE_BIOTYPE


def gene_variant_annotation(
    gene_id: str,
    variant_id: str,
    position: int,
    genes: Mapping[str, Mapping[str, object]],
    transcripts_by_gene: Mapping[str, Sequence[Mapping[str, object]]],
    variant_gene_map: Mapping[str, Mapping[str, Tuple[str, ...]]],
) -> Dict[str, object]:
    gene = genes.get(gene_id)
    if gene is None:
        return {field: "" for field in ANNOTATION_FIELDS} | {"annotation_status": "TARGET_GENE_NOT_IN_SELECTED_GTF"}
    site_map = variant_gene_map[variant_id]
    exonic = site_map["exonic"]
    span = site_map["span"]
    if len(exonic) == 1:
        assignment_class = "UNIQUE_EXONIC"
        requirement = "DIRECT"
    elif len(exonic) > 1:
        assignment_class = "AMBIGUOUS_SHARED_EXON"
        requirement = "STRAND_RESOLUTION_REQUIRED_FOR_SITE"
    elif span:
        assignment_class = "INTRONIC_ONLY"
        requirement = "NOT_USED_FOR_ASE_GENE_ASSIGNMENT"
    else:
        assignment_class = "INTERGENIC"
        requirement = "NOT_USED_FOR_ASE_GENE_ASSIGNMENT"

    transcripts = list(transcripts_by_gene.get(gene_id, []))
    transcript_features: List[Tuple[Mapping[str, object], str]] = [
        (transcript, classify_transcript_position(transcript, position))
        for transcript in transcripts
    ]
    union_features = [
        feature for transcript, feature in transcript_features if is_exon_union_transcript(transcript)
    ]
    canonical = [
        (transcript, feature)
        for transcript, feature in transcript_features
        if bool(transcript["canonical"])
    ]
    target_exonic = gene_id in exonic
    annotation_status = "PASS" if target_exonic else "TARGET_GENE_NOT_IN_EXON_UNION"
    return {
        "gtf_gene_name": gene["gene_name"],
        "gtf_gene_biotype": gene["gene_biotype"],
        "gtf_gene_contig": gene["contig"],
        "gtf_gene_start": gene["start"],
        "gtf_gene_end": gene["end"],
        "gtf_gene_strand": gene["strand"],
        "original_ase_assignment_class": assignment_class,
        "site_gene_assignment_requirement": requirement,
        "target_gene_in_protein_coding_exon_union": int(target_exonic),
        "n_exonic_protein_coding_genes_at_site": len(exonic),
        "exonic_protein_coding_gene_ids_at_site": ",".join(exonic),
        "gene_span_ids_at_site": ",".join(span),
        "n_transcripts_for_gene": len(transcripts),
        "n_exon_union_transcripts_for_gene": sum(
            is_exon_union_transcript(transcript) for transcript in transcripts
        ),
        "n_canonical_transcripts_for_gene": len(canonical),
        "canonical_transcript_ids": ",".join(
            str(transcript["transcript_id"]) for transcript, _ in canonical
        ),
        "canonical_transcript_features": ",".join(sorted({feature for _, feature in canonical})),
        "all_transcript_features": ",".join(sorted({feature for _, feature in transcript_features})),
        "exon_union_transcript_features": ",".join(sorted(set(union_features))),
        "n_transcripts_exonic_at_site": sum(
            feature in EXONIC_TRANSCRIPT_FEATURES for _, feature in transcript_features
        ),
        "n_exon_union_transcripts_exonic_at_site": sum(
            feature in EXONIC_TRANSCRIPT_FEATURES for feature in union_features
        ),
        "annotation_status": annotation_status,
    }


def format_row(row: Mapping[str, object], header: Sequence[str]) -> Dict[str, object]:
    return {field: row.get(field, "") for field in header}


def annotate_file(
    input_path: Path,
    output_path: Path,
    input_fields: Sequence[str],
    cache: Mapping[Tuple[str, str], Mapping[str, object]],
    counts: MutableMapping[Tuple[str, str], Counter],
    dataset: str,
) -> int:
    output_fields = list(input_fields) + ANNOTATION_FIELDS
    with open_text(output_path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        n_rows = 0
        for row in iter_tsv(input_path):
            gene_id = normalize_gene_id(row["gene_id"])
            annotation = cache[(gene_id, row["variant_id"])]
            output = dict(row)
            output.update(annotation)
            writer.writerow(format_row(output, output_fields))
            n_rows += 1
            counts[(dataset, "original_ase_assignment_class")][
                str(annotation["original_ase_assignment_class"])
            ] += 1
            canonical = str(annotation["canonical_transcript_features"]) or "NO_CANONICAL_TRANSCRIPT_FEATURE"
            counts[(dataset, "canonical_transcript_features")][canonical] += 1
            counts[(dataset, "annotation_status")][str(annotation["annotation_status"])] += 1
    return n_rows


def write_recurrent_transcript_rows(
    input_path: Path,
    output_path: Path,
    transcripts_by_gene: Mapping[str, Sequence[Mapping[str, object]]],
) -> int:
    n_rows = 0
    with open_text(output_path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSCRIPT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in iter_tsv(input_path):
            gene_id = normalize_gene_id(row["gene_id"])
            position = parse_int(row["position"], "position", f"{gene_id}|{row['variant_id']}")
            for transcript in transcripts_by_gene.get(gene_id, []):
                output = {
                    **{field: row.get(field, "") for field in TRANSCRIPT_HEADER},
                    "transcript_id": transcript["transcript_id"],
                    "transcript_name": transcript["transcript_name"],
                    "transcript_biotype": transcript["transcript_biotype"],
                    "canonical_transcript": int(bool(transcript["canonical"])),
                    "included_in_original_exon_union": int(is_exon_union_transcript(transcript)),
                    "transcript_feature": classify_transcript_position(transcript, position),
                }
                writer.writerow(format_row(output, TRANSCRIPT_HEADER))
                n_rows += 1
    return n_rows


def prepare_output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise AnalysisError(f"output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def file_metadata(path: Path) -> Dict[str, object]:
    stat = path.stat()
    return {"path": str(path.resolve()), "size_bytes": stat.st_size}


def run_analysis(args: argparse.Namespace) -> Dict[str, object]:
    for path in (args.all_hypotheses, args.recurrent, args.gtf):
        if not path.is_file():
            raise AnalysisError(f"missing input: {path}")
    prepare_output_directory(args.out_dir)

    all_fields = input_header(args.all_hypotheses)
    recurrent_fields = input_header(args.recurrent)
    variant_info, targets_by_variant, n_input_hypotheses = collect_variants(args.all_hypotheses)
    genes, transcripts_by_gene, gtf_counts = load_gtf_annotation(args.gtf)
    variant_gene_map = map_variants_to_selected_genes(variant_info, genes)

    cache: Dict[Tuple[str, str], Dict[str, object]] = {}
    for variant_id, target_genes in targets_by_variant.items():
        position = variant_info[variant_id][1]
        for gene_id in target_genes:
            cache[(gene_id, variant_id)] = gene_variant_annotation(
                gene_id,
                variant_id,
                position,
                genes,
                transcripts_by_gene,
                variant_gene_map,
            )

    failed = [key for key, annotation in cache.items() if annotation["annotation_status"] != "PASS"]
    if failed:
        raise AnalysisError(
            f"{len(failed)} tested gene-SNP keys do not reproduce the original exon-union assignment; "
            f"examples={failed[:10]}"
        )

    output_paths = {
        "all_tested": args.out_dir / "all_tested_gene_tissue_snp_feature_annotations.tsv.gz",
        "recurrent": args.out_dir / "recurrent_gene_tissue_snp_feature_annotations.tsv",
        "recurrent_transcripts": args.out_dir / "recurrent_snp_transcript_annotations.tsv.gz",
        "category_counts": args.out_dir / "annotation_category_counts.tsv",
    }
    category_counts: Dict[Tuple[str, str], Counter] = defaultdict(Counter)
    n_all_output = annotate_file(
        args.all_hypotheses,
        output_paths["all_tested"],
        all_fields,
        cache,
        category_counts,
        "ALL_TESTED",
    )
    n_recurrent_output = annotate_file(
        args.recurrent,
        output_paths["recurrent"],
        recurrent_fields,
        cache,
        category_counts,
        "STATISTICALLY_RECURRENT",
    )
    n_transcript_rows = write_recurrent_transcript_rows(
        args.recurrent, output_paths["recurrent_transcripts"], transcripts_by_gene
    )

    with output_paths["category_counts"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COUNT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (dataset, dimension), counter in sorted(category_counts.items()):
            for category, count in sorted(counter.items()):
                writer.writerow(
                    {
                        "dataset": dataset,
                        "dimension": dimension,
                        "category": category,
                        "n_gene_tissue_snp_hypotheses": count,
                    }
                )

    recurrent_assignment = category_counts[("STATISTICALLY_RECURRENT", "original_ase_assignment_class")]
    counts = {
        **gtf_counts,
        "n_input_all_tested_hypotheses": n_input_hypotheses,
        "n_unique_physical_variants": len(variant_info),
        "n_unique_tested_gene_snp_keys": len(cache),
        "n_all_tested_annotation_rows": n_all_output,
        "n_recurrent_annotation_rows": n_recurrent_output,
        "n_recurrent_transcript_annotation_rows": n_transcript_rows,
        "n_recurrent_unique_exonic": recurrent_assignment.get("UNIQUE_EXONIC", 0),
        "n_recurrent_ambiguous_shared_exon": recurrent_assignment.get("AMBIGUOUS_SHARED_EXON", 0),
    }
    qc = {
        "version": VERSION,
        "status": "PASS",
        "counts": counts,
        "invariants": {
            "all_tested_target_genes_reproduced_in_original_exon_union": True,
            "only_original_v025_assignment_classes_used": True,
            "transcript_feature_precedence_matches_v025": True,
            "no_vep_or_snpeff_consequence_annotation": True,
            "no_new_promoter_or_splice_category": True,
        },
    }
    with (args.out_dir / "analysis_qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")

    metadata = {
        "version": VERSION,
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "annotation_model": {
            "source": "phaser_ase_standalone_v0.2.5/call_phaser_ase_candidates.py",
            "selected_gene_biotype": SELECTED_GENE_BIOTYPE,
            "gene_assignment": "union of exons from protein-coding or unspecified transcripts",
            "site_assignment_classes": [
                "UNIQUE_EXONIC",
                "AMBIGUOUS_SHARED_EXON",
                "INTRONIC_ONLY",
                "INTERGENIC",
            ],
            "transcript_feature_precedence": [
                "CDS",
                "FIVE_PRIME_UTR",
                "THREE_PRIME_UTR",
                "START_CODON",
                "STOP_CODON",
                "EXON_NONCODING_OR_UNSPECIFIED",
                "INTRON",
                "OUTSIDE_TRANSCRIPT",
            ],
            "excluded_scope": [
                "VEP",
                "SnpEff",
                "protein consequence",
                "new promoter category",
                "new splice category",
            ],
        },
        "inputs": {
            "all_hypotheses": file_metadata(args.all_hypotheses),
            "recurrent": file_metadata(args.recurrent),
            "gtf": file_metadata(args.gtf),
        },
        "outputs": {key: str(path.resolve()) for key, path in output_paths.items()},
        "counts": counts,
        "software": {"python": sys.version},
    }
    with (args.out_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    return metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-hypotheses", required=True, type=Path)
    parser.add_argument("--recurrent", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        metadata = run_analysis(args)
    except AnalysisError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] recurrent annotated rows: {metadata['counts']['n_recurrent_annotation_rows']}")
    print(f"[PASS] output: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
