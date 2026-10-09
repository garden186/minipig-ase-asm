#!/usr/bin/env python3
"""
Call WGS-filtered local ASE candidates from GTF, phASER outputs and RNA BAM.

Primary evidence
----------------
* Gene/transcript mapping: GTF intervals and phASER variant coordinates.
* ASE counts: strict BAM/QNAME fragment reconstruction.  A fragment must have
  one haplotype, one inferred transcript strand and one gene.
* Shared-exon strand resolution: a shared exon is assigned only when the
  fragment strand matches exactly one gene.  Such fragments are merged with
  direct gene-specific evidence; their provenance remains supplementary QC.
* Local haplotype-block QC: phASER haplotypes and variant connections.
* WGS SNP QC: phased WGS genotype plus phASER-WGS REF/ALT balance.

Used only for direction-resolved interpretation, not ASE existence
------------------------------------------------------------------
* mapping of local phASER Hap A/B to one WGS phase-set orientation

Explicitly not used
-------------------
* blockGWPhase, gwStat, phase_concordant
* phaser_gene_ae gene-level aggregation
* any precomputed custom ASE block master

The statistical and graph-QC code remains Python 3.6 compatible.  BAM
reconstruction additionally requires pysam.
"""

from __future__ import print_function

import argparse
import csv
import functools
import gzip
import hashlib
import heapq
import json
import math
import os
import re
import shutil
import statistics
import sys
from collections import defaultdict

import strand_bam_engine as bam_engine
import wgs_snp_qc

VERSION = "0.2.5"
MIN_POSITIVE_P = sys.float_info.min


def raise_csv_field_limit():
    """Allow phASER read-ID fields that can exceed csv's 128 KiB default."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = limit // 10


raise_csv_field_limit()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def split_csv(value):
    if not value:
        return []
    return [x for x in value.split(",") if x]


def normalize_bam(value):
    value = os.path.basename(value or "")
    for suffix in (".bam", ".cram"):
        if value.endswith(suffix):
            value = value[:-len(suffix)]
    return value


def pair_key(a, b):
    a = "" if a is None else str(a)
    b = "" if b is None else str(b)
    return (a, b) if a <= b else (b, a)


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return "{:.12g}".format(value)
    return str(value)


@functools.lru_cache(maxsize=100000)
def exact_binom_two_sided_half(a_count, b_count):
    """Exact two-sided binomial p-value for H0 p=0.5.

    At p=0.5 the PMF is symmetric, so the probability-ordering two-sided test
    equals 2 * P[X <= min(a,b)], capped at 1. The tail is summed in log space.
    """
    n = int(a_count) + int(b_count)
    if n <= 0:
        return None
    k = min(int(a_count), int(b_count))
    logs = []
    base = math.lgamma(n + 1) - n * math.log(2.0)
    for i in range(k + 1):
        logs.append(base - math.lgamma(i + 1) - math.lgamma(n - i + 1))
    maximum = max(logs)
    log_cdf = maximum + math.log(sum(math.exp(x - maximum) for x in logs))
    pvalue = min(1.0, 2.0 * math.exp(log_cdf))
    # For very deep and extremely imbalanced measurements exp(log_cdf) can
    # underflow to 0. A mathematical p-value is never exactly zero. Preserve a
    # finite, rankable value for BH adjustment and -log10(p) plotting while
    # keeping the result on the conservative side of machine precision.
    return max(MIN_POSITIVE_P, pvalue)


def bh_adjust(indexed_pvalues):
    """Return {row_index: BH q-value}; indexed_pvalues is [(index, p), ...]."""
    clean = [(idx, p) for idx, p in indexed_pvalues if p is not None]
    clean.sort(key=lambda x: x[1])
    m = len(clean)
    result = {}
    running = 1.0
    for rank_from_end in range(m - 1, -1, -1):
        idx, pvalue = clean[rank_from_end]
        rank = rank_from_end + 1
        adjusted = min(1.0, pvalue * m / float(rank))
        running = min(running, adjusted)
        result[idx] = running
    return result


def read_master(path, sample_filter=None):
    rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty block master: {}".format(path))
        required = {
            "sample", "tissue", "bam", "gene_id", "component_id",
            "component_variants", "gene_variants", "local_a_count",
            "local_b_count", "local_total_count"
        }
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError("Component master is missing: {}".format(",".join(missing)))
        for row in reader:
            if sample_filter and row.get("sample") != sample_filter:
                continue
            rows.append(row)
    if not rows:
        raise ValueError("No component rows selected from {}".format(path))
    samples = sorted(set(row["sample"] for row in rows))
    if len(samples) != 1:
        raise ValueError(
            "variant_connections.txt is per individual; select one sample. "
            "Found: {}".format(",".join(samples))
        )
    return rows, reader.fieldnames


def component_catalog(master_rows):
    components = {}
    variant_to_components = defaultdict(set)
    for row in master_rows:
        component_id = row["component_id"]
        variants = tuple(split_csv(row["component_variants"]))
        if not variants:
            continue
        entry = components.setdefault(component_id, {
            "component_id": component_id,
            "variants": variants,
            "variant_set": set(variants),
            "genes": set(),
            "rows": 0,
            "start": row.get("component_start", ""),
            "end": row.get("component_end", ""),
            "contig": row.get("gene_contig", row.get("contig", "")),
            "phaser_pi": row.get("phaser_pi", "")
        })
        if entry["variants"] != variants:
            raise ValueError("Inconsistent variants for {}".format(component_id))
        entry["genes"].add(row["gene_id"])
        entry["rows"] += 1
        for variant in variants:
            variant_to_components[variant].add(component_id)
    return components, variant_to_components


def parse_gtf_attributes(text):
    attributes = {}
    for key, quoted, bare in re.findall(
        r'(?:^|;\s*)([^\s;]+)\s+(?:"([^"]*)"|([^;\s]+))', text
    ):
        value = quoted if quoted != "" else bare
        if key in attributes:
            if not isinstance(attributes[key], list):
                attributes[key] = [attributes[key]]
            attributes[key].append(value)
        else:
            attributes[key] = value
    return attributes


def attribute_first(attributes, *keys):
    for key in keys:
        value = attributes.get(key)
        if isinstance(value, list):
            return value[0] if value else ""
        if value:
            return value
    return ""


def load_gtf_annotation(path, selected_biotypes):
    """Load gene and transcript metadata needed by the phASER-only analysis."""
    genes = {}
    transcripts = {}
    selected = set(selected_biotypes)
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            contig, source, feature = fields[0], fields[1], fields[2]
            start, end, strand = parse_int(fields[3]), parse_int(fields[4]), fields[6]
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attribute_first(attributes, "gene_id")
            if not gene_id:
                continue
            gene_biotype = attribute_first(
                attributes, "gene_biotype", "gene_type", "biotype"
            )
            if selected and gene_biotype not in selected:
                continue

            if feature == "gene":
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": attribute_first(attributes, "gene_name") or gene_id,
                    "gene_biotype": gene_biotype,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source": source
                }

            transcript_id = attribute_first(attributes, "transcript_id")
            if not transcript_id:
                continue
            transcript = transcripts.setdefault(transcript_id, {
                "transcript_id": transcript_id,
                "transcript_name": (
                    attribute_first(attributes, "transcript_name") or transcript_id
                ),
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
                "features": defaultdict(list)
            })
            transcript["start"] = min(transcript["start"], start)
            transcript["end"] = max(transcript["end"], end)
            tags = attributes.get("tag", [])
            if not isinstance(tags, list):
                tags = [tags] if tags else []
            transcript["tags"].update(tags)
            if "Ensembl_canonical" in transcript["tags"]:
                transcript["canonical"] = True
            if feature != "transcript":
                transcript["features"][feature].append((start, end))

    transcripts = {
        transcript_id: transcript
        for transcript_id, transcript in transcripts.items()
        if transcript["gene_id"] in genes
    }
    transcripts_by_gene = defaultdict(list)
    for transcript in transcripts.values():
        transcripts_by_gene[transcript["gene_id"]].append(transcript)
    for gene_id in transcripts_by_gene:
        transcripts_by_gene[gene_id].sort(
            key=lambda x: (not x["canonical"], x["transcript_id"])
        )
        exon_intervals = []
        for transcript in transcripts_by_gene[gene_id]:
            transcript_biotype = transcript["transcript_biotype"]
            if transcript_biotype and transcript_biotype != "protein_coding":
                continue
            exon_intervals.extend(transcript["features"].get("exon", []))
        genes[gene_id]["exon_intervals"] = merge_intervals(exon_intervals)
    return genes, transcripts_by_gene


def merge_intervals(intervals):
    """Merge 1-based inclusive intervals."""
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(interval) for interval in merged]


def component_identity_from_row(row):
    variants = tuple(split_csv(row.get("variants", "")))
    return (
        row.get("contig", ""),
        parse_int(row.get("start")),
        parse_int(row.get("stop")),
        variants
    )


def make_component_id(sample, identity):
    contig, start, stop, variants = identity
    digest = hashlib.sha1(",".join(variants).encode("utf-8")).hexdigest()[:10]
    return "{}:{}:{}-{}:C{}".format(sample, contig, start, stop, digest)


def scan_phaser_components(path, sample):
    """First pass over haplotypic_counts: collect definitions, not read IDs."""
    components = {}
    identity_to_component = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "contig", "start", "stop", "variants", "variantCount",
            "haplotypeA", "haplotypeB", "aCount", "bCount", "bam"
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "Invalid haplotypic_counts; missing: {}".format(
                    ",".join(sorted(required.difference(reader.fieldnames or [])))
                )
            )
        for row in reader:
            identity = component_identity_from_row(row)
            variants = identity[3]
            if not variants:
                continue
            component_id = identity_to_component.get(identity)
            if component_id is None:
                component_id = make_component_id(sample, identity)
                identity_to_component[identity] = component_id
                components[component_id] = {
                    "component_id": component_id,
                    "identity": identity,
                    "contig": identity[0],
                    "start": identity[1],
                    "end": identity[2],
                    "phaser_pi": "",
                    "variants": variants,
                    "variant_set": set(variants),
                    "haplotypeA": row.get("haplotypeA", ""),
                    "haplotypeB": row.get("haplotypeB", ""),
                    "genes": set(),
                    "gene_variants": {},
                    "n_bam_rows": 0
                }
            components[component_id]["n_bam_rows"] += 1
    if not components:
        raise ValueError("No phASER components found in {}".format(path))
    return components, identity_to_component


def parse_variant_id(variant_id):
    parts = variant_id.rsplit("_", 3)
    if len(parts) != 4:
        return None
    contig, position, ref, alt = parts
    try:
        position = int(position)
    except ValueError:
        return None
    return contig, position, ref, alt


def map_components_to_genes(components, genes):
    """Map SNPs by selected transcript-exon unions.

    v0.2.4 preserves both uniquely exonic and shared-exon SNPs.  Shared-exon
    SNPs are not counted from annotation alone; their candidate genes are
    carried into the BAM strand-resolution step.  Consequently a component
    containing only shared-exon SNPs is retained instead of being discarded.
    """
    variants_by_contig_position = defaultdict(lambda: defaultdict(list))
    invalid_variants = []
    for component in components.values():
        for variant in component["variants"]:
            parsed = parse_variant_id(variant)
            if parsed is None:
                invalid_variants.append(variant)
                continue
            contig, position = parsed[0], parsed[1]
            variants_by_contig_position[contig][position].append(variant)

    genes_by_contig = defaultdict(list)
    exons_by_contig = defaultdict(list)
    for gene in genes.values():
        genes_by_contig[gene["contig"]].append(gene)
        for start, end in gene.get("exon_intervals", []):
            exons_by_contig[gene["contig"]].append(
                (start, end, gene["gene_id"])
            )
    for contig in genes_by_contig:
        genes_by_contig[contig].sort(key=lambda x: (x["start"], x["end"]))
    for contig in exons_by_contig:
        exons_by_contig[contig].sort()

    variant_to_span_genes = defaultdict(list)
    for contig, positions in variants_by_contig_position.items():
        contig_genes = genes_by_contig.get(contig, [])
        active = {}
        end_heap = []
        gene_index = 0
        for position in sorted(positions):
            while (
                gene_index < len(contig_genes) and
                contig_genes[gene_index]["start"] <= position
            ):
                gene = contig_genes[gene_index]
                active[gene["gene_id"]] = gene
                heapq.heappush(end_heap, (gene["end"], gene["gene_id"]))
                gene_index += 1
            while end_heap and end_heap[0][0] < position:
                end, gene_id = heapq.heappop(end_heap)
                current = active.get(gene_id)
                if current is not None and current["end"] == end:
                    active.pop(gene_id, None)
            overlapping = sorted(active)
            for variant in positions[position]:
                variant_to_span_genes[variant].extend(overlapping)

    variant_to_exonic_genes = defaultdict(list)
    for contig, positions in variants_by_contig_position.items():
        intervals = exons_by_contig.get(contig, [])
        active = {}
        end_heap = []
        interval_index = 0
        for position in sorted(positions):
            while (
                interval_index < len(intervals) and
                intervals[interval_index][0] <= position
            ):
                start, end, gene_id = intervals[interval_index]
                active[(gene_id, start, end)] = gene_id
                heapq.heappush(end_heap, (end, gene_id, start))
                interval_index += 1
            while end_heap and end_heap[0][0] < position:
                end, gene_id, start = heapq.heappop(end_heap)
                active.pop((gene_id, start, end), None)
            overlapping = sorted(set(active.values()))
            for variant in positions[position]:
                variant_to_exonic_genes[variant] = overlapping

    relevant_components = {}
    for component_id, component in components.items():
        gene_variants = defaultdict(list)
        candidate_gene_variants = defaultdict(list)
        ambiguous_gene_variants = defaultdict(list)
        ambiguous_variants = {}
        intronic_variants = {}
        for variant in component["variants"]:
            exonic_genes = variant_to_exonic_genes.get(variant, [])
            if len(exonic_genes) == 1:
                gene_variants[exonic_genes[0]].append(variant)
                candidate_gene_variants[exonic_genes[0]].append(variant)
            elif len(exonic_genes) > 1:
                ambiguous_variants[variant] = list(exonic_genes)
                for gene_id in exonic_genes:
                    candidate_gene_variants[gene_id].append(variant)
                    ambiguous_gene_variants[gene_id].append(variant)
            else:
                span_genes = variant_to_span_genes.get(variant, [])
                if span_genes:
                    intronic_variants[variant] = list(span_genes)
        if not candidate_gene_variants:
            continue
        component["gene_variants"] = dict(gene_variants)
        component["candidate_gene_variants"] = dict(
            candidate_gene_variants
        )
        component["ambiguous_gene_variants"] = dict(
            ambiguous_gene_variants
        )
        component["genes"] = set(candidate_gene_variants)
        component["ambiguous_variants"] = ambiguous_variants
        component["intronic_variants"] = intronic_variants
        relevant_components[component_id] = component
    assignment = {
        "exonic": dict(variant_to_exonic_genes),
        "span": dict(variant_to_span_genes)
    }
    return relevant_components, assignment, sorted(set(invalid_variants))


def catalog_from_components(components):
    variant_to_components = defaultdict(set)
    for component_id, component in components.items():
        for variant in component["variants"]:
            variant_to_components[variant].add(component_id)
    return variant_to_components


def prepare_component_haplotype_alleles(components):
    """Attach parsed SNP/HapA/HapB records required by the BAM engine."""
    for component in components.values():
        variants = component["variants"]
        hap_a = tuple(split_csv(component.get("haplotypeA", "")))
        hap_b = tuple(split_csv(component.get("haplotypeB", "")))
        if not (len(variants) == len(hap_a) == len(hap_b)):
            raise ValueError(
                "Variant/haplotype length mismatch for {}".format(
                    component["component_id"]
                )
            )
        parsed_variants = []
        for variant_id, allele_a, allele_b in zip(
            variants, hap_a, hap_b
        ):
            parsed = parse_variant_id(variant_id)
            if parsed is None:
                raise ValueError("Cannot parse variant {}".format(variant_id))
            parsed_variants.append(
                {
                    "variant_id": variant_id,
                    "contig": parsed[0],
                    "position": parsed[1],
                    "ref": parsed[2].upper(),
                    "alt": parsed[3].upper(),
                    "hapA": allele_a.upper(),
                    "hapB": allele_b.upper()
                }
            )
        component["parsed_variants"] = parsed_variants
    return components


def build_variant_assignment_index(components, variant_assignment):
    """Build the annotation evidence consumed for each BAM SNP observation."""
    result = {}
    variants = {
        variant
        for component in components.values()
        for variant in component["variants"]
    }
    for variant in variants:
        exonic = list(variant_assignment["exonic"].get(variant, []))
        if len(exonic) == 1:
            assignment_class = "UNIQUE_EXONIC"
        elif len(exonic) > 1:
            assignment_class = "AMBIGUOUS_SHARED_EXON"
        elif variant_assignment["span"].get(variant):
            assignment_class = "INTRONIC_ONLY"
        else:
            assignment_class = "INTERGENIC"
        result[variant] = {
            "assignment_class": assignment_class,
            "gene_ids": exonic
        }
    return result


def load_bam_tissue_map(path):
    mapping = {}
    if not path:
        return mapping
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"bam", "tissue"}.issubset(reader.fieldnames):
            raise ValueError("--bam-tissue-map requires bam and tissue columns")
        for row in reader:
            mapping[normalize_bam(row["bam"])] = row["tissue"]
    return mapping


def infer_tissue(bam, sample, bam_tissue_map):
    normalized = normalize_bam(bam)
    if normalized in bam_tissue_map:
        return bam_tissue_map[normalized]
    stem = normalized
    if stem.endswith(".phaseready"):
        stem = stem[:-len(".phaseready")]
    prefix = sample + "_"
    if stem.startswith(prefix):
        return stem[len(prefix):]
    return stem


def calculate_gene_read_metrics(row, component, gene_variants, args):
    variants = component["variants"]
    variant_index = {variant: index for index, variant in enumerate(variants)}
    selected = [variant_index[v] for v in gene_variants if v in variant_index]
    a_segments = row.get("aReads", "").split(";")
    b_segments = row.get("bReads", "").split(";")
    while len(a_segments) < len(variants):
        a_segments.append("")
    while len(b_segments) < len(variants):
        b_segments.append("")
    a_sets = [set(split_csv(x)) for x in a_segments[:len(variants)]]
    b_sets = [set(split_csv(x)) for x in b_segments[:len(variants)]]
    ids_available = any(a_sets) or any(b_sets)
    reported_a = parse_int(row.get("aCount"))
    reported_b = parse_int(row.get("bCount"))
    all_component_variants_in_gene = set(gene_variants) == set(variants)

    if not ids_available:
        if all_component_variants_in_gene:
            return {
                "count_available": True,
                "a": reported_a,
                "b": reported_b,
                "total": reported_a + reported_b,
                "count_source": "PHASER_REPORTED_COMPONENT_COUNTS",
                "count_validation": "READ_IDS_UNAVAILABLE_FOR_ROW",
                "n_gene_variants_with_reads": None,
                "gene_variants_with_reads": [],
                "support_class": (
                    "SINGLE_SNP_ONLY"
                    if len(gene_variants) == 1
                    else "READ_IDS_UNAVAILABLE_COMPONENT_TOTAL"
                ),
                "top_variant": "",
                "top_variant_unique_reads": None,
                "top_variant_unique_fraction": None,
                "top_loo": None,
                "loo_all_robust": None,
                "loo_direction_consistent": None,
                "loo_min_total": None,
                "loo_min_major_fraction": None,
                "loo_max_p": None
            }
        return {
            "count_available": False,
            "a": None,
            "b": None,
            "total": None,
            "count_source": "UNAVAILABLE",
            "count_validation": "PARTIAL_COMPONENT_WITHOUT_READ_IDS",
            "n_gene_variants_with_reads": None,
            "gene_variants_with_reads": [],
            "support_class": "UNRESOLVED_LOCAL_COUNT",
            "top_variant": "",
            "top_variant_unique_reads": None,
            "top_variant_unique_fraction": None,
            "top_loo": None,
            "loo_all_robust": None,
            "loo_direction_consistent": None,
            "loo_min_total": None,
            "loo_min_major_fraction": None,
            "loo_max_p": None
        }

    per_variant = {}
    a_union = set()
    b_union = set()
    for index in selected:
        variant = variants[index]
        per_variant[variant] = {
            "a": set(a_sets[index]),
            "b": set(b_sets[index])
        }
        a_union.update(a_sets[index])
        b_union.update(b_sets[index])
    contributing = [
        variant for variant in gene_variants
        if variant in per_variant and
        (per_variant[variant]["a"] or per_variant[variant]["b"])
    ]

    top_variant = ""
    top_unique_reads = 0
    top_reduction = -1
    top_loo = None
    all_loo = []
    for variant in contributing:
        other_a = set()
        other_b = set()
        for other, read_sets in per_variant.items():
            if other == variant:
                continue
            other_a.update(read_sets["a"])
            other_b.update(read_sets["b"])
        unique_reads = len(per_variant[variant]["a"] - other_a)
        unique_reads += len(per_variant[variant]["b"] - other_b)
        reduction = (
            len(a_union) + len(b_union) - len(other_a) - len(other_b)
        )
        loo_total = len(other_a) + len(other_b)
        loo_mhf = (
            max(len(other_a), len(other_b)) / float(loo_total)
            if loo_total else None
        )
        loo_p = (
            exact_binom_two_sided_half(len(other_a), len(other_b))
            if loo_total else None
        )
        loo = {
            "variant": variant,
            "a": len(other_a),
            "b": len(other_b),
            "total": loo_total,
            "mhf": loo_mhf,
            "p": loo_p,
            "unique_reads": unique_reads,
            "reduction": reduction
        }
        all_loo.append(loo)
        if top_loo is None or (reduction, unique_reads) > (
            top_reduction, top_unique_reads
        ):
            top_variant = variant
            top_unique_reads = unique_reads
            top_reduction = reduction
            top_loo = loo

    total = len(a_union) + len(b_union)
    unique_fraction = top_unique_reads / float(total) if total else None
    loo_effect_persists = False
    if top_loo is not None and top_loo["total"] >= args.min_total:
        loo_effect_persists = (
            top_loo["mhf"] is not None and
            top_loo["mhf"] >= args.min_major_fraction and
            top_loo["p"] is not None and
            top_loo["p"] <= 0.05
        )
    full_direction = 1 if len(a_union) > len(b_union) else -1
    usable_loo = [loo for loo in all_loo if loo["total"] >= args.min_total]
    loo_direction_consistent = bool(all_loo)
    for loo in all_loo:
        loo_direction = (
            1 if loo["a"] > loo["b"] else
            -1 if loo["b"] > loo["a"] else 0
        )
        if loo_direction != full_direction:
            loo_direction_consistent = False
            break
    loo_all_robust = (
        len(contributing) > 1 and
        len(usable_loo) == len(all_loo) and
        loo_direction_consistent and
        all(
            loo["mhf"] is not None and
            loo["mhf"] >= args.min_major_fraction and
            loo["p"] is not None and
            loo["p"] <= 0.05
            for loo in all_loo
        )
    )
    if len(contributing) <= 1:
        support_class = "SINGLE_SNP_ONLY"
    elif (
        unique_fraction is not None and
        unique_fraction >= args.single_variant_driver_fraction and
        not loo_effect_persists
    ):
        support_class = "SINGLE_SNP_DOMINANT"
    elif loo_all_robust:
        support_class = "MULTI_SNP_LOO_ROBUST"
    else:
        support_class = "MULTI_SNP_SUPPORTED"

    validation = "READ_ID_UNION_COMPUTED"
    if all_component_variants_in_gene:
        validation = (
            "MATCH"
            if len(a_union) == reported_a and len(b_union) == reported_b
            else "MISMATCH"
        )
    return {
        "count_available": True,
        "a": len(a_union),
        "b": len(b_union),
        "total": total,
        "count_source": "PHASER_READ_INDEX_UNION",
        "count_validation": validation,
        "n_gene_variants_with_reads": len(contributing),
        "gene_variants_with_reads": contributing,
        "support_class": support_class,
        "top_variant": top_variant,
        "top_variant_unique_reads": top_unique_reads,
        "top_variant_unique_fraction": unique_fraction,
        "top_loo": top_loo,
        "loo_all_robust": loo_all_robust,
        "loo_direction_consistent": loo_direction_consistent,
        "loo_min_total": min(
            (loo["total"] for loo in all_loo), default=None
        ),
        "loo_min_major_fraction": min(
            (loo["mhf"] for loo in all_loo if loo["mhf"] is not None),
            default=None
        ),
        "loo_max_p": max(
            (loo["p"] for loo in all_loo if loo["p"] is not None),
            default=None
        )
    }


def classify_transcript_position(transcript, position):
    if position < transcript["start"] or position > transcript["end"]:
        return "OUTSIDE_TRANSCRIPT"
    precedence = [
        ("CDS", "CDS"),
        ("five_prime_utr", "FIVE_PRIME_UTR"),
        ("three_prime_utr", "THREE_PRIME_UTR"),
        ("start_codon", "START_CODON"),
        ("stop_codon", "STOP_CODON"),
        ("exon", "EXON_NONCODING_OR_UNSPECIFIED")
    ]
    for feature_name, label in precedence:
        for start, end in transcript["features"].get(feature_name, []):
            if start <= position <= end:
                return label
    return "INTRON"


def build_transcript_annotation_rows(components, genes, transcripts_by_gene):
    seen = set()
    rows = []
    for component in components.values():
        for gene_id, variants in component[
            "candidate_gene_variants"
        ].items():
            for variant in variants:
                parsed = parse_variant_id(variant)
                if parsed is None:
                    continue
                position = parsed[1]
                for transcript in transcripts_by_gene.get(gene_id, []):
                    key = (gene_id, variant, transcript["transcript_id"])
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append({
                        "gene_id": gene_id,
                        "gene_name": genes[gene_id]["gene_name"],
                        "variant_id": variant,
                        "contig": parsed[0],
                        "position": position,
                        "transcript_id": transcript["transcript_id"],
                        "transcript_name": transcript["transcript_name"],
                        "transcript_biotype": transcript["transcript_biotype"],
                        "canonical_transcript": transcript["canonical"],
                        "transcript_feature": classify_transcript_position(
                            transcript, position
                        )
                    })
    return rows


def write_transcript_annotation_stream(
    path, components, genes, transcripts_by_gene, fieldnames
):
    temp_path = path + ".tmp"
    count = 0
    seen = set()
    with open(temp_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for component in components.values():
            for gene_id, variants in component[
                "candidate_gene_variants"
            ].items():
                for variant in variants:
                    parsed = parse_variant_id(variant)
                    if parsed is None:
                        continue
                    for transcript in transcripts_by_gene.get(gene_id, []):
                        key = (gene_id, variant, transcript["transcript_id"])
                        if key in seen:
                            continue
                        seen.add(key)
                        row = {
                            "gene_id": gene_id,
                            "gene_name": genes[gene_id]["gene_name"],
                            "variant_id": variant,
                            "contig": parsed[0],
                            "position": parsed[1],
                            "transcript_id": transcript["transcript_id"],
                            "transcript_name": transcript["transcript_name"],
                            "transcript_biotype": (
                                transcript["transcript_biotype"]
                            ),
                            "canonical_transcript": transcript["canonical"],
                            "transcript_feature": classify_transcript_position(
                                transcript, parsed[1]
                            )
                        }
                        writer.writerow({
                            field: fmt(row.get(field)) for field in fieldnames
                        })
                        count += 1
    os.replace(temp_path, path)
    return count


def load_haplotypes(path, components):
    by_variant_tuple = defaultdict(list)
    for component_id, component in components.items():
        by_variant_tuple[tuple(component["variants"])].append(component_id)
        by_variant_tuple[tuple(sorted(component["variants"]))].append(component_id)

    found = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            variants = tuple(split_csv(row.get("variant_ids", "")))
            component_ids = set(by_variant_tuple.get(variants, []))
            component_ids.update(by_variant_tuple.get(tuple(sorted(variants)), []))
            for component_id in component_ids:
                found[component_id] = row
    return found


def load_haplotypic_read_ids(path, components):
    """Load only relevant phASER rows, keyed by normalized BAM and variant tuple."""
    relevant_sets = {}
    for component in components.values():
        relevant_sets[tuple(component["variants"])] = True
        relevant_sets[tuple(sorted(component["variants"]))] = True

    result = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variants", "bam", "aReads", "bReads"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "haplotypic_counts must contain variants,bam,aReads,bReads. "
                "The phASER run requires --output_read_ids 1."
            )
        for row in reader:
            variants = tuple(split_csv(row.get("variants", "")))
            if variants not in relevant_sets and tuple(sorted(variants)) not in relevant_sets:
                continue
            a_segments = row.get("aReads", "").split(";")
            b_segments = row.get("bReads", "").split(";")
            while len(a_segments) < len(variants):
                a_segments.append("")
            while len(b_segments) < len(variants):
                b_segments.append("")
            a_sets = [set(split_csv(segment)) for segment in a_segments[:len(variants)]]
            b_sets = [set(split_csv(segment)) for segment in b_segments[:len(variants)]]
            key = (normalize_bam(row.get("bam", "")), variants)
            result[key] = {
                "variants": variants,
                "a_sets": a_sets,
                "b_sets": b_sets,
                "reported_a": parse_int(row.get("aCount")),
                "reported_b": parse_int(row.get("bCount"))
            }
            result[(key[0], tuple(sorted(variants)))] = result[key]
    return result


def read_id_metrics(master_row, read_id_index, min_total, min_major_fraction,
                    driver_fraction):
    variants = tuple(split_csv(master_row["component_variants"]))
    bam = normalize_bam(master_row["bam"])
    record = read_id_index.get((bam, variants))
    if record is None:
        record = read_id_index.get((bam, tuple(sorted(variants))))
    if record is None:
        return None

    has_any_read_ids = any(record["a_sets"]) or any(record["b_sets"])
    if (
        not has_any_read_ids and
        record["reported_a"] + record["reported_b"] > 0
    ):
        # phASER may emit singleton haplotypic-count rows without read IDs even
        # when --output_read_ids 1 was used. Those rows remain usable through
        # their counts, but cannot be reconstructed from aReads/bReads.
        return {
            "ids_available": False,
            "raw_component_reported_a": record["reported_a"],
            "raw_component_reported_b": record["reported_b"]
        }

    index = {variant: i for i, variant in enumerate(record["variants"])}
    gene_variants = [v for v in split_csv(master_row["gene_variants"]) if v in index]
    selected = [index[v] for v in gene_variants]
    if not selected:
        return None

    a_union = set()
    b_union = set()
    for i in selected:
        a_union.update(record["a_sets"][i])
        b_union.update(record["b_sets"][i])

    per_variant = {}
    for i in selected:
        variant = record["variants"][i]
        per_variant[variant] = {
            "a": set(record["a_sets"][i]),
            "b": set(record["b_sets"][i])
        }

    contributing = [
        variant for variant, sets in per_variant.items()
        if sets["a"] or sets["b"]
    ]

    top_variant = ""
    top_unique_reads = 0
    top_reduction = -1
    top_loo = None
    all_loo = []
    for variant in contributing:
        other_a = set()
        other_b = set()
        for other, sets in per_variant.items():
            if other == variant:
                continue
            other_a.update(sets["a"])
            other_b.update(sets["b"])
        unique_reads = len(per_variant[variant]["a"] - other_a)
        unique_reads += len(per_variant[variant]["b"] - other_b)
        reduction = (len(a_union) + len(b_union)) - (len(other_a) + len(other_b))
        loo_total = len(other_a) + len(other_b)
        loo_mhf = (
            max(len(other_a), len(other_b)) / float(loo_total)
            if loo_total else None
        )
        loo_p = exact_binom_two_sided_half(len(other_a), len(other_b)) if loo_total else None
        loo = {
            "variant": variant,
            "a": len(other_a),
            "b": len(other_b),
            "total": loo_total,
            "mhf": loo_mhf,
            "p": loo_p,
            "unique_reads": unique_reads,
            "reduction": reduction
        }
        all_loo.append(loo)
        score = (reduction, unique_reads)
        best_score = (
            top_reduction,
            top_unique_reads
        )
        if top_loo is None or score > best_score:
            top_variant = variant
            top_unique_reads = unique_reads
            top_reduction = reduction
            top_loo = loo

    total = len(a_union) + len(b_union)
    unique_fraction = top_unique_reads / float(total) if total else None
    loo_effect_persists = False
    if top_loo is not None and top_loo["total"] >= min_total:
        loo_effect_persists = (
            top_loo["mhf"] is not None and
            top_loo["mhf"] >= min_major_fraction and
            top_loo["p"] is not None and
            top_loo["p"] <= 0.05
        )

    if len(contributing) <= 1:
        support_class = "SINGLE_SNP_ONLY"
    elif (
        unique_fraction is not None and
        unique_fraction >= driver_fraction and
        not loo_effect_persists
    ):
        support_class = "SINGLE_SNP_DOMINANT"
    elif loo_effect_persists:
        support_class = "MULTI_SNP_LOO_ROBUST"
    else:
        support_class = "MULTI_SNP_SUPPORTED"

    return {
        "ids_available": True,
        "a": len(a_union),
        "b": len(b_union),
        "total": total,
        "gene_variants_with_reads": contributing,
        "n_gene_variants_with_reads": len(contributing),
        "top_variant": top_variant,
        "top_variant_unique_reads": top_unique_reads,
        "top_variant_unique_fraction": unique_fraction,
        "top_loo": top_loo,
        "support_class": support_class,
        "raw_component_reported_a": record["reported_a"],
        "raw_component_reported_b": record["reported_b"]
    }


def load_connections(path, variant_to_components, components, cc_threshold):
    internal = defaultdict(list)
    cross = []
    relevant_pairs = {}

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "variant_a", "variant_b", "supporting_connections",
            "total_connections", "conflicting_configuration_p"
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("Invalid variant_connections file: {}".format(path))
        for row in reader:
            a = row["variant_a"]
            b = row["variant_b"]
            comps_a = variant_to_components.get(a)
            comps_b = variant_to_components.get(b)
            if not comps_a or not comps_b:
                continue
            support = parse_int(row["supporting_connections"])
            total = parse_int(row["total_connections"])
            pvalue = parse_float(row["conflicting_configuration_p"], 0.0)
            edge = {
                "variant_a": a,
                "variant_b": b,
                "supporting_connections": support,
                "total_connections": total,
                "conflicting_reads": max(0, total - support),
                "support_fraction": support / float(total) if total else None,
                "conflicting_configuration_p": pvalue,
                "retained": pvalue >= cc_threshold,
                "pair": pair_key(a, b)
            }
            relevant_pairs[edge["pair"]] = edge
            shared_components = comps_a.intersection(comps_b)
            if shared_components:
                for component_id in shared_components:
                    internal[component_id].append(edge)
            else:
                for component_a in sorted(comps_a):
                    for component_b in sorted(comps_b):
                        shared_genes = (
                            components[component_a]["genes"] &
                            components[component_b]["genes"]
                        )
                        if shared_genes:
                            item = dict(edge)
                            item["component_a"] = component_a
                            item["component_b"] = component_b
                            item["shared_genes"] = sorted(shared_genes)
                            cross.append(item)
    return internal, cross, relevant_pairs


def load_relevant_connection_pairs(path, relevant_variants, cc_threshold):
    pairs = {}
    if not path:
        return pairs
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            a = row.get("variant_a", "")
            b = row.get("variant_b", "")
            if a not in relevant_variants or b not in relevant_variants:
                continue
            support = parse_int(row.get("supporting_connections"))
            total = parse_int(row.get("total_connections"))
            pvalue = parse_float(row.get("conflicting_configuration_p"), 0.0)
            pairs[pair_key(a, b)] = {
                "supporting_connections": support,
                "total_connections": total,
                "conflicting_configuration_p": pvalue,
                "retained": pvalue >= cc_threshold
            }
    return pairs


def load_allele_config(path, relevant_pairs):
    result = {}
    qc = {
        "path": path or "",
        "rows_total": 0,
        "rows_relevant": 0,
        "rows_malformed_skipped": 0
    }
    if not path:
        return result, qc
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            qc["rows_total"] += 1
            variant_a = row.get("variant_a")
            variant_b = row.get("variant_b")
            configuration = row.get("configuration")
            if not variant_a or not variant_b or not configuration:
                qc["rows_malformed_skipped"] += 1
                continue
            key = pair_key(variant_a, variant_b)
            if key in relevant_pairs:
                result[key] = configuration
                qc["rows_relevant"] += 1
    if qc["rows_malformed_skipped"]:
        print(
            "[WARN] skipped {} malformed allele_config rows: {}".format(
                qc["rows_malformed_skipped"], path
            ),
            file=sys.stderr
        )
    return result, qc


def connected_components(nodes, retained_edges):
    adjacency = {node: set() for node in nodes}
    for edge in retained_edges:
        a, b = edge["variant_a"], edge["variant_b"]
        if a in adjacency and b in adjacency:
            adjacency[a].add(b)
            adjacency[b].add(a)
    seen = set()
    groups = []
    for start in nodes:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        group = []
        while stack:
            node = stack.pop()
            group.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        groups.append(group)
    return groups, adjacency


def bridge_pairs(nodes, retained_edges):
    adjacency = {node: [] for node in nodes}
    for edge in retained_edges:
        a, b = edge["variant_a"], edge["variant_b"]
        if a in adjacency and b in adjacency:
            adjacency[a].append(b)
            adjacency[b].append(a)

    sys.setrecursionlimit(max(10000, len(nodes) * 3 + 100))
    discovery = {}
    low = {}
    parent = {}
    bridges = set()
    counter = [0]

    def dfs(node):
        counter[0] += 1
        discovery[node] = counter[0]
        low[node] = counter[0]
        for neighbor in adjacency[node]:
            if neighbor not in discovery:
                parent[neighbor] = node
                dfs(neighbor)
                low[node] = min(low[node], low[neighbor])
                if low[neighbor] > discovery[node]:
                    bridges.add(pair_key(node, neighbor))
            elif parent.get(node) != neighbor:
                low[node] = min(low[node], discovery[neighbor])

    for node in nodes:
        if node not in discovery:
            dfs(node)
    return bridges


def evidence_class(full_edge, wgs_edge, full_config, wgs_config):
    if wgs_edge is None:
        return "RNA_ENABLED_OR_WGS_NOT_OBSERVED"
    if full_config and wgs_config and full_config != wgs_config:
        return "EVIDENCE_DISCORDANT"
    if full_edge["retained"] and not wgs_edge["retained"]:
        return "EVIDENCE_DISCORDANT"
    if full_edge["retained"] and wgs_edge["retained"]:
        if full_edge["total_connections"] > wgs_edge["total_connections"]:
            return "WGS_RNA_SUPPORTED"
        return "WGS_BACKED"
    if not full_edge["retained"] and wgs_edge["retained"]:
        return "FULL_REJECTED_WGS_RETAINED"
    return "REJECTED_BOTH"


def threshold_class(pvalue, cc_threshold, reference_threshold):
    if pvalue >= cc_threshold:
        return "PASSES_CONFLICT_FILTER"
    lower = min(cc_threshold, reference_threshold)
    upper = max(cc_threshold, reference_threshold)
    if lower <= pvalue < upper:
        return "THRESHOLD_SENSITIVE"
    return "STRONG_CONFLICT"


def build_graph_qc(components, internal_edges, haplotypes, wgs_pairs,
                   full_configs, wgs_configs, cc_threshold,
                   reference_threshold, min_bridge_support,
                   max_conflict_fraction):
    graph_rows = []
    edge_rows = []

    for component_id in sorted(components):
        component = components[component_id]
        nodes = list(component["variants"])
        observed = internal_edges.get(component_id, [])
        retained = [edge for edge in observed if edge["retained"]]
        rejected = [edge for edge in observed if not edge["retained"]]
        groups, adjacency = connected_components(nodes, retained)
        bridges = bridge_pairs(nodes, retained)
        retained_by_pair = {edge["pair"]: edge for edge in retained}
        bridge_supports = [
            retained_by_pair[pair]["supporting_connections"]
            for pair in bridges if pair in retained_by_pair
        ]
        total_connections = sum(edge["total_connections"] for edge in retained)
        total_conflicts = sum(edge["conflicting_reads"] for edge in retained)
        weighted_conflict_fraction = (
            total_conflicts / float(total_connections)
            if total_connections else None
        )
        support_values = [edge["supporting_connections"] for edge in retained]
        support_fractions = [
            edge["support_fraction"] for edge in retained
            if edge["support_fraction"] is not None
        ]
        hap = haplotypes.get(component_id)
        reported_edges_total = (
            parse_int(hap.get("edges_total")) if hap is not None else None
        )
        reported_edges_supporting = (
            parse_float(hap.get("edges_supporting")) if hap is not None else None
        )
        if len(nodes) == 1:
            reconstruction = "SINGLETON"
            phase_qc = "SINGLETON"
            phase_reason = "ONE_VARIANT"
        else:
            connected = len(groups) == 1
            edge_count_match = (
                reported_edges_total is not None and
                reported_edges_total == len(retained)
            )
            edge_configuration_match = (
                reported_edges_supporting is not None and
                reported_edges_total is not None and
                abs(
                    reported_edges_supporting - reported_edges_total
                ) < 1e-9
            )
            reconstruction = "MATCH" if connected and edge_count_match else "MISMATCH"
            if hap is None:
                phase_qc = "PHASE_LC"
                phase_reason = "HAPLOTYPE_ROW_MISSING"
            elif reconstruction != "MATCH":
                phase_qc = "PHASE_LC"
                phase_reason = "GRAPH_RECONSTRUCTION_MISMATCH"
            elif not edge_configuration_match:
                phase_qc = "PHASE_LC"
                phase_reason = "HAPLOTYPE_EDGE_CONFIGURATION_MISMATCH"
            elif bridge_supports and min(bridge_supports) < min_bridge_support:
                phase_qc = "PHASE_MC"
                phase_reason = "WEAK_BRIDGE"
            elif (
                weighted_conflict_fraction is not None and
                weighted_conflict_fraction > max_conflict_fraction
            ):
                phase_qc = "PHASE_MC"
                phase_reason = "HIGH_RETAINED_EDGE_CONFLICT"
            else:
                phase_qc = "PHASE_HC"
                phase_reason = "CONNECTED_GRAPH_WITH_SUPPORTED_BRIDGES"

        graph_rows.append({
            "component_id": component_id,
            "contig": component["contig"],
            "component_start": component["start"],
            "component_end": component["end"],
            "phaser_pi": component["phaser_pi"],
            "n_variants": len(nodes),
            "n_edges_observed": len(observed),
            "n_edges_retained": len(retained),
            "n_edges_rejected": len(rejected),
            "n_retained_connected_groups": len(groups),
            "retained_graph_connected": len(groups) == 1,
            "cycle_rank": len(retained) - len(nodes) + len(groups),
            "n_bridges": len(bridges),
            "min_bridge_support": min(bridge_supports) if bridge_supports else None,
            "n_single_read_bridges": sum(x <= 1 for x in bridge_supports),
            "min_retained_edge_support": min(support_values) if support_values else None,
            "median_retained_edge_support": (
                statistics.median(support_values) if support_values else None
            ),
            "min_retained_support_fraction": (
                min(support_fractions) if support_fractions else None
            ),
            "weighted_retained_conflict_fraction": weighted_conflict_fraction,
            "haplotypes_edges_supporting": reported_edges_supporting,
            "haplotypes_edges_total": reported_edges_total,
            "haplotype_edge_configuration_match": (
                len(nodes) == 1 or (
                    reported_edges_supporting is not None
                    and reported_edges_total is not None
                    and abs(
                        reported_edges_supporting - reported_edges_total
                    ) < 1e-9
                )
            ),
            "graph_reconstruction": reconstruction,
            "phase_qc": phase_qc,
            "phase_qc_reason": phase_reason,
            "component_variants": ",".join(nodes)
        })

        for edge in observed:
            wgs = wgs_pairs.get(edge["pair"])
            edge_rows.append({
                "component_id": component_id,
                "variant_a": edge["variant_a"],
                "variant_b": edge["variant_b"],
                "supporting_connections": edge["supporting_connections"],
                "total_connections": edge["total_connections"],
                "conflicting_reads": edge["conflicting_reads"],
                "support_fraction": edge["support_fraction"],
                "conflicting_configuration_p": edge["conflicting_configuration_p"],
                "cc_threshold": cc_threshold,
                "edge_status": threshold_class(
                    edge["conflicting_configuration_p"],
                    cc_threshold,
                    reference_threshold
                ),
                "is_bridge": edge["pair"] in bridges and edge["retained"],
                "full_configuration": full_configs.get(edge["pair"], ""),
                "wgs_supporting_connections": (
                    wgs["supporting_connections"] if wgs else None
                ),
                "wgs_total_connections": wgs["total_connections"] if wgs else None,
                "wgs_conflicting_configuration_p": (
                    wgs["conflicting_configuration_p"] if wgs else None
                ),
                "wgs_configuration": wgs_configs.get(edge["pair"], ""),
                "evidence_class": evidence_class(
                    edge,
                    wgs,
                    full_configs.get(edge["pair"], ""),
                    wgs_configs.get(edge["pair"], "")
                )
            })
    return graph_rows, edge_rows


def assess_local_phase_for_snps(
    component_id,
    supported_snps,
    internal_edges,
    haplotypes,
    min_bridge_support,
    max_conflict_fraction,
):
    """Assess the read graph induced by WGS-balanced RNA-supported SNPs."""
    nodes = sorted(set(supported_snps))
    if len(nodes) <= 1:
        return {
            "local_phase_status": "SINGLE_SNP",
            "local_phase_reason": "ONE_RNA_SUPPORTED_SNP",
            "local_phase_n_snps": len(nodes),
            "local_phase_n_edges": 0,
            "local_phase_min_bridge_support": None,
            "local_phase_conflict_fraction": None,
        }
    node_set = set(nodes)
    observed = [
        edge for edge in internal_edges.get(component_id, [])
        if edge["retained"]
        and edge["variant_a"] in node_set
        and edge["variant_b"] in node_set
    ]
    groups, _ = connected_components(nodes, observed)
    if len(groups) != 1:
        return {
            "local_phase_status": "FAIL",
            "local_phase_reason": "WGS_FILTERED_SNP_GRAPH_DISCONNECTED",
            "local_phase_n_snps": len(nodes),
            "local_phase_n_edges": len(observed),
            "local_phase_min_bridge_support": None,
            "local_phase_conflict_fraction": None,
        }
    hap = haplotypes.get(component_id)
    reported_total = (
        parse_int(hap.get("edges_total")) if hap is not None else None
    )
    reported_supporting = (
        parse_float(hap.get("edges_supporting"))
        if hap is not None else None
    )
    if not (
        reported_total is not None
        and reported_supporting is not None
        and abs(reported_supporting - reported_total) < 1e-9
    ):
        return {
            "local_phase_status": "FAIL",
            "local_phase_reason": "HAPLOTYPE_EDGE_CONFIGURATION_MISMATCH",
            "local_phase_n_snps": len(nodes),
            "local_phase_n_edges": len(observed),
            "local_phase_min_bridge_support": None,
            "local_phase_conflict_fraction": None,
        }
    bridges = bridge_pairs(nodes, observed)
    retained_by_pair = {edge["pair"]: edge for edge in observed}
    bridge_supports = [
        retained_by_pair[pair]["supporting_connections"]
        for pair in bridges if pair in retained_by_pair
    ]
    total_connections = sum(
        edge["total_connections"] for edge in observed
    )
    total_conflicts = sum(edge["conflicting_reads"] for edge in observed)
    conflict_fraction = (
        total_conflicts / float(total_connections)
        if total_connections else None
    )
    minimum_bridge = min(bridge_supports) if bridge_supports else None
    if bridge_supports and minimum_bridge < min_bridge_support:
        status = "REVIEW"
        reason = "WEAK_BRIDGE"
    elif (
        conflict_fraction is not None
        and conflict_fraction > max_conflict_fraction
    ):
        status = "REVIEW"
        reason = "HIGH_RETAINED_EDGE_CONFLICT"
    else:
        status = "PASS"
        reason = "CONNECTED_WGS_FILTERED_SNP_GRAPH"
    return {
        "local_phase_status": status,
        "local_phase_reason": reason,
        "local_phase_n_snps": len(nodes),
        "local_phase_n_edges": len(observed),
        "local_phase_min_bridge_support": minimum_bridge,
        "local_phase_conflict_fraction": conflict_fraction,
    }


def classify_candidate(row, graph, thresholds):
    """Return one compact biological class plus explicit QC reasons."""
    reasons = []
    if row["analysis_set"] == "EXCLUDED_TISSUE":
        reasons.append("EXCLUDED_TISSUE")
    if row["count_validation"] == "MISMATCH":
        reasons.append("COUNT_MISMATCH")
    if row["ase_test_status"] == "LOW_COVERAGE":
        reasons.append("LOW_COVERAGE")
    elif row["ase_test_status"] != "TESTED":
        reasons.append("NOT_TESTABLE")
    if row["ase_test_status"] == "TESTED" and (
        row["ase_q_bh"] is None
        or row["ase_q_bh"] > thresholds["fdr"]
        or row["major_haplotype_fraction"] < thresholds["min_major_fraction"]
    ):
        reasons.append("NO_SIGNIFICANT_ASE")
    if (
        row.get("fragment_hap_conflict_fraction") is not None
        and row["fragment_hap_conflict_fraction"]
        > thresholds.get("max_fragment_hap_conflict_fraction", 1.0)
    ):
        reasons.append("FRAGMENT_CONFLICT")
    if row.get("minor_support_class") == "EXTREME_MINOR":
        reasons.append("EXTREME_MINOR")
    elif row.get("minor_support_class") == "LOW_MINOR":
        reasons.append("LOW_MINOR")
    if row.get(
        "n_wgs_balanced_gene_snps",
        row.get("n_gene_variants_with_reads", 0),
    ) == 0:
        reasons.append("NO_WGS_BALANCED_SNP")

    support = row.get(
        "snp_support_status",
        row.get("count_support_class", ""),
    )
    single_snp = support in {
        "SINGLE_SNP_ONLY", "SINGLE_SNP_DOMINANT"
    }
    if not single_snp and row.get("local_phase_status") != "PASS":
        reasons.append("LOCAL_PHASE_REVIEW")

    hard_failures = {
        "EXCLUDED_TISSUE", "COUNT_MISMATCH", "LOW_COVERAGE",
        "NOT_TESTABLE", "NO_SIGNIFICANT_ASE", "FRAGMENT_CONFLICT",
        "EXTREME_MINOR", "LOW_MINOR", "NO_WGS_BALANCED_SNP",
    }
    if hard_failures.intersection(reasons):
        final_class = "NOT_PRIMARY"
    elif single_snp:
        final_class = "SECONDARY_SINGLE_SNP"
        reasons.append(support)
    elif "LOCAL_PHASE_REVIEW" in reasons:
        final_class = "NOT_PRIMARY"
    elif support == "MULTI_SNP_LOO_ROBUST":
        final_class = "PRIMARY_ROBUST"
    elif support == "MULTI_SNP_SUPPORTED":
        final_class = "PRIMARY_SUPPORTED"
    else:
        final_class = "NOT_PRIMARY"
        reasons.append("SNP_SUPPORT_UNRESOLVED")
    return final_class, ",".join(dict.fromkeys(reasons))


def build_candidate_rows(master_rows, graph_rows, read_id_index, args):
    graph_by_component = {row["component_id"]: row for row in graph_rows}
    rows = []
    excluded = set(split_csv(args.exclude_tissues))
    thresholds = {
        "fdr": args.fdr,
        "min_major_fraction": args.min_major_fraction
    }

    for master in master_rows:
        row = dict(master)
        graph = graph_by_component.get(master["component_id"])
        raw = (
            read_id_metrics(
                master,
                read_id_index,
                args.min_total,
                args.min_major_fraction,
                args.single_variant_driver_fraction
            )
            if read_id_index is not None else None
        )
        usable_raw = raw is not None and raw.get("ids_available", True)
        master_a = parse_int(master.get("local_a_count"))
        master_b = parse_int(master.get("local_b_count"))
        if usable_raw:
            a_count = raw["a"]
            b_count = raw["b"]
            count_source = "RAW_READ_ID_UNION"
            count_validation = (
                "MATCH" if a_count == master_a and b_count == master_b
                else "MISMATCH"
            )
            n_contributing = raw["n_gene_variants_with_reads"]
            contributing_variants = ",".join(raw["gene_variants_with_reads"])
            support_class = raw["support_class"]
            top_loo = raw["top_loo"]
        else:
            a_count = master_a
            b_count = master_b
            count_source = "MASTER_LOCAL_COUNTS"
            count_validation = (
                "READ_IDS_UNAVAILABLE_FOR_ROW"
                if raw is not None else "NOT_CHECKED"
            )
            n_contributing = None
            contributing_variants = ""
            support_class = "READ_IDS_NOT_AVAILABLE"
            top_loo = None

        total = a_count + b_count
        minor = min(a_count, b_count)
        major_fraction = max(a_count, b_count) / float(total) if total else None
        log2_ratio = math.log((a_count + 0.5) / float(b_count + 0.5), 2)
        allelic_imbalance = abs(a_count - b_count) / float(total) if total else None
        n_gene_variants = parse_int(master.get("n_gene_variants"))
        if n_contributing is not None:
            ase_unit = (
                "SINGLE_SNP_ASE"
                if n_contributing <= 1
                else "LOCAL_HAPLOTYPE_BLOCK_ASE"
            )
        else:
            ase_unit = (
                "SINGLE_SNP_ASE"
                if n_gene_variants <= 1
                else "LOCAL_HAPLOTYPE_BLOCK_ASE"
            )

        analysis_set = (
            "EXCLUDED_TISSUE" if master["tissue"] in excluded else "MAIN"
        )
        if total < args.min_total:
            test_status = "LOW_COVERAGE"
            pvalue = None
        elif minor <= args.extreme_minor_max:
            test_status = "EXTREME"
            pvalue = exact_binom_two_sided_half(a_count, b_count)
        elif minor < args.min_each_haplotype:
            test_status = "INSUFFICIENT_MINOR_HAP"
            pvalue = exact_binom_two_sided_half(a_count, b_count)
        else:
            test_status = "TESTABLE"
            pvalue = exact_binom_two_sided_half(a_count, b_count)

        row.update({
            "analysis_set": analysis_set,
            "ase_unit_type_local": ase_unit,
            "count_source": count_source,
            "count_validation": count_validation,
            "ase_a_count": a_count,
            "ase_b_count": b_count,
            "ase_total_count": total,
            "major_haplotype_fraction": major_fraction,
            "log2_a_over_b": log2_ratio,
            "absolute_allelic_imbalance": allelic_imbalance,
            "ase_test_status": test_status,
            "ase_p_exact": pvalue,
            "ase_q_bh": None,
            "n_gene_variants_with_reads": n_contributing,
            "gene_variants_with_reads": contributing_variants,
            "count_support_class": support_class,
            "top_count_driver_variant": raw["top_variant"] if usable_raw else "",
            "top_driver_unique_reads": (
                raw["top_variant_unique_reads"] if usable_raw else None
            ),
            "top_driver_unique_fraction": (
                raw["top_variant_unique_fraction"] if usable_raw else None
            ),
            "loo_a_count": top_loo["a"] if top_loo else None,
            "loo_b_count": top_loo["b"] if top_loo else None,
            "loo_total_count": top_loo["total"] if top_loo else None,
            "loo_major_haplotype_fraction": top_loo["mhf"] if top_loo else None,
            "loo_p_exact": top_loo["p"] if top_loo else None,
            "phase_qc": graph["phase_qc"] if graph else "PHASE_LC",
            "phase_qc_reason": (
                graph["phase_qc_reason"] if graph else "GRAPH_QC_MISSING"
            ),
            "graph_reconstruction": (
                graph["graph_reconstruction"] if graph else "MISSING"
            )
        })
        rows.append(row)

    by_family = defaultdict(list)
    for i, row in enumerate(rows):
        if row["analysis_set"] == "MAIN" and row["ase_test_status"] == "TESTABLE":
            key = (row["sample"], row["tissue"])
            by_family[key].append((i, row["ase_p_exact"]))
    for indexed in by_family.values():
        adjusted = bh_adjust(indexed)
        for index, qvalue in adjusted.items():
            rows[index]["ase_q_bh"] = qvalue

    for row in rows:
        graph = graph_by_component.get(row["component_id"])
        (
            row["candidate_class"],
            row["exclusion_reasons"],
        ) = classify_candidate(row, graph, thresholds)
        row["is_primary_candidate"] = row["candidate_class"] in {
            "PRIMARY_ROBUST", "PRIMARY_SUPPORTED"
        }
    return rows


def apply_unique_measurement_bh(rows):
    """Apply BH once per unique measurement within sample x tissue."""
    measurement_rows = defaultdict(list)
    for index, row in enumerate(rows):
        measurement_rows[row["measurement_id"]].append(index)
    for indices in measurement_rows.values():
        for index in indices:
            rows[index]["measurement_annotation_multiplicity"] = len(indices)

    by_family = defaultdict(list)
    measurement_index = {}
    for measurement_id, indices in measurement_rows.items():
        row = rows[indices[0]]
        signatures = {
            (
                rows[index]["ase_a_count"],
                rows[index]["ase_b_count"],
                rows[index]["ase_p_exact"]
            )
            for index in indices
        }
        if len(signatures) != 1:
            raise ValueError(
                "Inconsistent duplicate measurement: {}".format(
                    measurement_id
                )
            )
        if (
            row["analysis_set"] == "MAIN" and
            row["ase_test_status"] == "TESTED"
        ):
            family = (row["sample"], row["tissue"])
            local_index = len(by_family[family])
            by_family[family].append(
                (local_index, row["ase_p_exact"])
            )
            measurement_index[(family, local_index)] = measurement_id

    measurement_q = {}
    for family, indexed_pvalues in by_family.items():
        for local_index, qvalue in bh_adjust(indexed_pvalues).items():
            measurement_q[
                measurement_index[(family, local_index)]
            ] = qvalue
    for measurement_id, indices in measurement_rows.items():
        qvalue = measurement_q.get(measurement_id)
        for index in indices:
            rows[index]["ase_q_bh"] = qvalue
    return rows


def build_standalone_candidate_rows(
    haplotypic_counts_path,
    sample,
    components,
    identity_to_component,
    genes,
    transcripts_by_gene,
    graph_rows,
    bam_tissue_map,
    args
):
    """Second streaming pass: create component × gene × BAM ASE rows."""
    graph_by_component = {row["component_id"]: row for row in graph_rows}
    excluded = set(split_csv(args.exclude_tissues))
    rows = []

    with open_text(haplotypic_counts_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for phaser_row in reader:
            identity = component_identity_from_row(phaser_row)
            component_id = identity_to_component.get(identity)
            if component_id is None or component_id not in components:
                continue
            component = components[component_id]
            graph = graph_by_component.get(component_id)
            bam = phaser_row.get("bam", "")
            tissue = infer_tissue(bam, sample, bam_tissue_map)
            analysis_set = (
                "EXCLUDED_TISSUE" if tissue in excluded else "MAIN"
            )

            for gene_id in sorted(component["gene_variants"]):
                gene = genes[gene_id]
                gene_variants = component["gene_variants"][gene_id]
                metrics = calculate_gene_read_metrics(
                    phaser_row, component, gene_variants, args
                )
                a_count = metrics["a"]
                b_count = metrics["b"]
                total = metrics["total"]
                if metrics["count_available"]:
                    minor = min(a_count, b_count)
                    major_fraction = (
                        max(a_count, b_count) / float(total) if total else None
                    )
                    log2_ratio = math.log(
                        (a_count + 0.5) / float(b_count + 0.5), 2
                    )
                    allelic_imbalance = (
                        abs(a_count - b_count) / float(total) if total else None
                    )
                    if total < args.min_total:
                        test_status = "LOW_COVERAGE"
                        pvalue = None
                    else:
                        test_status = "TESTED"
                        pvalue = exact_binom_two_sided_half(a_count, b_count)
                    if minor <= args.extreme_minor_max:
                        minor_support_class = "EXTREME_MINOR"
                    elif minor < args.min_each_haplotype:
                        minor_support_class = "LOW_MINOR"
                    else:
                        minor_support_class = "ADEQUATE_MINOR"
                else:
                    major_fraction = None
                    log2_ratio = None
                    allelic_imbalance = None
                    test_status = "UNRESOLVED_LOCAL_COUNT"
                    pvalue = None
                    minor_support_class = "UNRESOLVED"

                n_contributing = metrics["n_gene_variants_with_reads"]
                if n_contributing is not None:
                    ase_unit = (
                        "SINGLE_SNP_ASE"
                        if n_contributing <= 1
                        else "LOCAL_HAPLOTYPE_BLOCK_ASE"
                    )
                else:
                    ase_unit = (
                        "SINGLE_SNP_ASE"
                        if len(gene_variants) <= 1
                        else "LOCAL_HAPLOTYPE_BLOCK_ASE"
                    )
                relation = (
                    "ALL_COMPONENT_VARIANTS_UNIQUELY_EXONIC_FOR_GENE"
                    if set(gene_variants) == set(component["variants"])
                    else "GENE_SPECIFIC_EXONIC_SUBSET"
                )
                transcripts = transcripts_by_gene.get(gene_id, [])
                canonical = [
                    transcript["transcript_id"]
                    for transcript in transcripts if transcript["canonical"]
                ]
                top_loo = metrics["top_loo"]
                measurement_key = "|".join([
                    sample,
                    tissue,
                    normalize_bam(bam),
                    component_id,
                    ",".join(sorted(gene_variants))
                ])
                measurement_id = "M" + hashlib.sha1(
                    measurement_key.encode("utf-8")
                ).hexdigest()[:16]
                row = {
                    "sample": sample,
                    "tissue": tissue,
                    "bam": bam,
                    "gene_id": gene_id,
                    "gene_name": gene["gene_name"],
                    "gene_biotype": gene["gene_biotype"],
                    "gene_contig": gene["contig"],
                    "gene_start": gene["start"],
                    "gene_end": gene["end"],
                    "strand": gene["strand"],
                    "component_id": component_id,
                    "component_start": component["start"],
                    "component_end": component["end"],
                    "component_gene_relation": relation,
                    "n_component_variants": len(component["variants"]),
                    "component_variants": ",".join(component["variants"]),
                    "n_gene_variants": len(gene_variants),
                    "gene_variants": ",".join(gene_variants),
                    "gene_assignment_policy": (
                        "UNIQUE_PROTEIN_CODING_EXON_UNION"
                    ),
                    "measurement_id": measurement_id,
                    "haplotypeA": phaser_row.get(
                        "haplotypeA", component["haplotypeA"]
                    ),
                    "haplotypeB": phaser_row.get(
                        "haplotypeB", component["haplotypeB"]
                    ),
                    "phaser_component_a_count": parse_int(
                        phaser_row.get("aCount")
                    ),
                    "phaser_component_b_count": parse_int(
                        phaser_row.get("bCount")
                    ),
                    "phaser_component_total_count": parse_int(
                        phaser_row.get("totalCount")
                    ),
                    "transcript_count": len(transcripts),
                    "transcript_ids": ",".join(
                        transcript["transcript_id"] for transcript in transcripts
                    ),
                    "canonical_transcript_ids": ",".join(canonical),
                    "analysis_set": analysis_set,
                    "ase_unit_type_local": ase_unit,
                    "count_source": metrics["count_source"],
                    "count_validation": metrics["count_validation"],
                    "ase_a_count": a_count,
                    "ase_b_count": b_count,
                    "ase_total_count": total,
                    "major_haplotype_fraction": major_fraction,
                    "log2_a_over_b": log2_ratio,
                    "absolute_allelic_imbalance": allelic_imbalance,
                    "ase_test_status": test_status,
                    "minor_support_class": minor_support_class,
                    "ase_p_exact": pvalue,
                    "ase_q_bh": None,
                    "n_gene_variants_with_reads": n_contributing,
                    "gene_variants_with_reads": ",".join(
                        metrics["gene_variants_with_reads"]
                    ),
                    "count_support_class": metrics["support_class"],
                    "top_count_driver_variant": metrics["top_variant"],
                    "top_driver_unique_reads": (
                        metrics["top_variant_unique_reads"]
                    ),
                    "top_driver_unique_fraction": (
                        metrics["top_variant_unique_fraction"]
                    ),
                    "loo_a_count": top_loo["a"] if top_loo else None,
                    "loo_b_count": top_loo["b"] if top_loo else None,
                    "loo_total_count": top_loo["total"] if top_loo else None,
                    "loo_major_haplotype_fraction": (
                        top_loo["mhf"] if top_loo else None
                    ),
                    "loo_p_exact": top_loo["p"] if top_loo else None,
                    "loo_all_robust": metrics["loo_all_robust"],
                    "loo_direction_consistent": (
                        metrics["loo_direction_consistent"]
                    ),
                    "loo_min_total": metrics["loo_min_total"],
                    "loo_min_major_fraction": (
                        metrics["loo_min_major_fraction"]
                    ),
                    "loo_max_p": metrics["loo_max_p"],
                    "phase_qc": graph["phase_qc"] if graph else "PHASE_LC",
                    "phase_qc_reason": (
                        graph["phase_qc_reason"]
                        if graph else "GRAPH_QC_MISSING"
                    ),
                    "graph_reconstruction": (
                        graph["graph_reconstruction"] if graph else "MISSING"
                    )
                }
                rows.append(row)

    apply_unique_measurement_bh(rows)

    thresholds = {
        "fdr": args.fdr,
        "min_major_fraction": args.min_major_fraction
    }
    for row in rows:
        graph = graph_by_component.get(row["component_id"])
        (
            row["candidate_class"],
            row["exclusion_reasons"],
        ) = classify_candidate(row, graph, thresholds)
        row["is_primary_candidate"] = row["candidate_class"] in {
            "PRIMARY_ROBUST", "PRIMARY_SUPPORTED"
        }
        row["candidate_tier"] = row["candidate_class"]
    return rows


def build_gene_summary(candidate_rows, fdr=0.05, min_major_fraction=0.65):
    groups = defaultdict(list)
    for row in candidate_rows:
        if row["analysis_set"] != "MAIN":
            continue
        key = (
            row["sample"], row["tissue"], row["gene_id"],
            row.get("gene_name", ""), row.get("gene_biotype", "")
        )
        groups[key].append(row)

    summaries = []
    for key in sorted(groups):
        rows = groups[key]
        informative = [
            row for row in rows
            if row["ase_total_count"] is not None and row["ase_total_count"] > 0
        ]
        tested = [row for row in rows if row["ase_test_status"] == "TESTED"]
        q_effect_pass = [
            row for row in rows
            if row["ase_test_status"] == "TESTED"
            and row["ase_q_bh"] is not None
            and row["ase_q_bh"] <= fdr
            and row["major_haplotype_fraction"] is not None
            and row["major_haplotype_fraction"] >= min_major_fraction
        ]
        primary = [row for row in rows if row["is_primary_candidate"]]
        secondary = [
            row for row in rows
            if row["candidate_class"] == "SECONDARY_SINGLE_SNP"
        ]
        blocks = sorted(
            set(row["haplotype_block_id"] for row in informative)
        )
        top_pool = primary or secondary or q_effect_pass or informative
        top = None
        if top_pool:
            top = max(
                top_pool,
                key=lambda row: (
                    row["major_haplotype_fraction"] or -1,
                    row["ase_total_count"]
                )
            )
        phase_scope = (
            "SINGLE_BLOCK" if len(blocks) <= 1
            else "MULTI_BLOCK_UNRESOLVED"
        )
        if any(
            row["candidate_class"] == "PRIMARY_ROBUST" for row in rows
        ):
            gene_status = "ASE_GENE_ROBUST"
        elif any(
            row["candidate_class"] == "PRIMARY_SUPPORTED" for row in rows
        ):
            gene_status = "ASE_GENE_SUPPORTED"
        elif secondary:
            gene_status = "SECONDARY_ONLY"
        else:
            gene_status = "NOT_ASE_GENE"

        summaries.append({
            "sample": key[0],
            "tissue": key[1],
            "gene_id": key[2],
            "gene_name": key[3],
            "gene_biotype": key[4],
            "n_informative_haplotype_blocks": len(blocks),
            "n_testable_haplotype_blocks": len(tested),
            "n_q_effect_pass_blocks": len(q_effect_pass),
            "n_primary_ase_blocks": len(primary),
            "n_primary_robust_blocks": sum(
                row["candidate_class"] == "PRIMARY_ROBUST"
                for row in rows
            ),
            "n_primary_supported_blocks": sum(
                row["candidate_class"] == "PRIMARY_SUPPORTED"
                for row in rows
            ),
            "n_secondary_single_snp_blocks": sum(
                row["candidate_class"] == "SECONDARY_SINGLE_SNP"
                for row in rows
            ),
            "n_extreme_minor_blocks": sum(
                row.get("minor_support_class") == "EXTREME_MINOR"
                for row in rows
            ),
            "gene_phase_scope": phase_scope,
            "gene_ase_status": gene_status,
            "top_haplotype_block_id": (
                top["haplotype_block_id"] if top else ""
            ),
            "top_candidate_class": top["candidate_class"] if top else "",
            "top_major_haplotype_fraction": (
                top["major_haplotype_fraction"] if top else None
            ),
            "top_total_count": top["ase_total_count"] if top else None,
            "primary_haplotype_block_ids": ",".join(
                sorted(
                    row["haplotype_block_id"] for row in primary
                )
            )
        })
    return summaries


def cross_edge_rows(cross_edges, wgs_pairs, full_configs, wgs_configs,
                    cc_threshold, reference_threshold):
    rows = []
    for edge in cross_edges:
        wgs = wgs_pairs.get(edge["pair"])
        rows.append({
            "component_a": edge["component_a"],
            "component_b": edge["component_b"],
            "shared_gene_ids": ",".join(edge["shared_genes"]),
            "variant_a": edge["variant_a"],
            "variant_b": edge["variant_b"],
            "supporting_connections": edge["supporting_connections"],
            "total_connections": edge["total_connections"],
            "conflicting_reads": edge["conflicting_reads"],
            "support_fraction": edge["support_fraction"],
            "conflicting_configuration_p": edge["conflicting_configuration_p"],
            "edge_status": threshold_class(
                edge["conflicting_configuration_p"],
                cc_threshold,
                reference_threshold
            ),
            "full_configuration": full_configs.get(edge["pair"], ""),
            "wgs_configuration": wgs_configs.get(edge["pair"], ""),
            "evidence_class": evidence_class(
                edge, wgs,
                full_configs.get(edge["pair"], ""),
                wgs_configs.get(edge["pair"], "")
            )
        })
    return rows


def write_tsv(path, rows, fieldnames):
    temp_path = path + ".tmp"
    with open(temp_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})
    os.replace(temp_path, path)


def build_variant_assignment_rows(all_components, assignment, genes):
    variant_components = defaultdict(set)
    for component_id, component in all_components.items():
        for variant in component["variants"]:
            variant_components[variant].add(component_id)
    rows = []
    for variant in sorted(variant_components):
        parsed = parse_variant_id(variant)
        exonic = sorted(set(assignment["exonic"].get(variant, [])))
        span = sorted(set(assignment["span"].get(variant, [])))
        if len(exonic) == 1:
            assignment_class = "UNIQUE_EXONIC"
        elif len(exonic) > 1:
            assignment_class = "AMBIGUOUS_SHARED_EXON"
        elif span:
            assignment_class = "INTRONIC_ONLY"
        else:
            assignment_class = "INTERGENIC"
        rows.append({
            "variant_id": variant,
            "contig": parsed[0] if parsed else "",
            "position": parsed[1] if parsed else None,
            "assignment_class": assignment_class,
            "exonic_gene_ids": ",".join(exonic),
            "exonic_gene_names": ",".join(
                genes[gene_id]["gene_name"] for gene_id in exonic
            ),
            "gene_span_ids": ",".join(span),
            "component_ids": ",".join(sorted(variant_components[variant]))
        })
    return rows


def build_measurement_rows(candidate_rows):
    grouped = defaultdict(list)
    for row in candidate_rows:
        grouped[row["measurement_id"]].append(row)
    result = []
    for measurement_id in sorted(grouped):
        rows = grouped[measurement_id]
        first = rows[0]
        result.append({
            "measurement_id": measurement_id,
            "sample": first["sample"],
            "tissue": first["tissue"],
            "bam": first["bam"],
            "haplotype_block_id": first["haplotype_block_id"],
            "gene_ids": ",".join(sorted(set(row["gene_id"] for row in rows))),
            "n_gene_annotations": len(rows),
            "gene_variants": first["gene_variants"],
            "strict_a_count": first.get("strict_a_count"),
            "strict_b_count": first.get("strict_b_count"),
            "strict_total_count": first.get("strict_total_count"),
            "direct_total_count": first.get("direct_total_count"),
            "strand_resolved_total": first.get("strand_resolved_total"),
            "rescue_fraction": first.get("rescue_fraction"),
            "rescue_dependency": first.get("rescue_dependency"),
            "ase_a_count": first["ase_a_count"],
            "ase_b_count": first["ase_b_count"],
            "ase_total_count": first["ase_total_count"],
            "ase_p_exact": first["ase_p_exact"],
            "ase_q_bh": first["ase_q_bh"],
            "minor_support_class": first["minor_support_class"],
            "candidate_class": first["candidate_class"],
            "is_primary_candidate": first["is_primary_candidate"]
        })
    return result


def prepare_output_directory(path, force):
    absolute = os.path.abspath(path)
    protected = {
        os.path.abspath(os.sep),
        os.path.abspath(os.path.expanduser("~")),
        os.path.abspath(os.getcwd())
    }
    if absolute in protected or len([x for x in absolute.split(os.sep) if x]) < 3:
        raise ValueError("Refusing unsafe output directory: {}".format(absolute))
    if os.path.exists(absolute):
        if not force:
            raise SystemExit(
                "Output directory exists; use --force: {}".format(absolute)
            )
        shutil.rmtree(absolute)
    os.makedirs(absolute)
    return absolute


def load_phaser_measurements(
    haplotypic_counts_path,
    sample,
    components,
    identity_to_component,
    bam_tissue_map,
):
    """Load one compact record per component x BAM row."""
    measurements = []
    seen = set()
    with open_text(haplotypic_counts_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            identity = component_identity_from_row(raw)
            component_id = identity_to_component.get(identity)
            if component_id is None or component_id not in components:
                continue
            bam = normalize_bam(raw.get("bam", ""))
            tissue = infer_tissue(bam, sample, bam_tissue_map)
            key = (component_id, bam, tissue)
            if key in seen:
                raise ValueError(
                    "Duplicate component/BAM measurement: {}".format(key)
                )
            seen.add(key)
            measurements.append(
                {
                    "sample": sample,
                    "tissue": tissue,
                    "bam": bam,
                    "component_id": component_id,
                    "reported_a": parse_int(raw.get("aCount")),
                    "reported_b": parse_int(raw.get("bCount")),
                    "reported_total": parse_int(raw.get("totalCount")),
                    "haplotypeA": raw.get(
                        "haplotypeA",
                        components[component_id].get("haplotypeA", "")
                    ),
                    "haplotypeB": raw.get(
                        "haplotypeB",
                        components[component_id].get("haplotypeB", "")
                    ),
                }
            )
    return measurements


def empty_strict_gene_count():
    return {
        "strict": {"A": set(), "B": set()},
        "direct": {"A": set(), "B": set()},
        "strand_resolved": {"A": set(), "B": set()},
        "variant_qnames": {
            "A": defaultdict(set), "B": defaultdict(set)
        },
        "strand_resolved_variant_qnames": {
            "A": defaultdict(set), "B": defaultdict(set)
        }
    }


def calculate_strict_variant_metrics(gene_count, args):
    """Calculate full-signal SNP driver and leave-one-variant-out metrics."""
    a_union = set(gene_count["strict"]["A"])
    b_union = set(gene_count["strict"]["B"])
    variant_ids = sorted(
        set(gene_count["variant_qnames"]["A"]).union(
            gene_count["variant_qnames"]["B"]
        )
    )
    per_variant = {}
    for variant in variant_ids:
        per_variant[variant] = {
            "a": set(gene_count["variant_qnames"]["A"].get(variant, set())),
            "b": set(gene_count["variant_qnames"]["B"].get(variant, set()))
        }
    contributing = [
        variant for variant in variant_ids
        if per_variant[variant]["a"] or per_variant[variant]["b"]
    ]

    top_variant = ""
    top_unique_reads = 0
    top_reduction = -1
    top_loo = None
    all_loo = []
    for variant in contributing:
        other_a = set()
        other_b = set()
        for other, read_sets in per_variant.items():
            if other == variant:
                continue
            other_a.update(read_sets["a"])
            other_b.update(read_sets["b"])
        unique_reads = len(per_variant[variant]["a"] - other_a)
        unique_reads += len(per_variant[variant]["b"] - other_b)
        reduction = (
            len(a_union) + len(b_union) - len(other_a) - len(other_b)
        )
        loo_total = len(other_a) + len(other_b)
        loo_mhf = (
            max(len(other_a), len(other_b)) / float(loo_total)
            if loo_total else None
        )
        loo_p = (
            exact_binom_two_sided_half(len(other_a), len(other_b))
            if loo_total else None
        )
        loo = {
            "variant": variant,
            "a": len(other_a),
            "b": len(other_b),
            "total": loo_total,
            "mhf": loo_mhf,
            "p": loo_p,
            "unique_reads": unique_reads,
            "reduction": reduction
        }
        all_loo.append(loo)
        if top_loo is None or (reduction, unique_reads, variant) > (
            top_reduction, top_unique_reads, top_variant
        ):
            top_variant = variant
            top_unique_reads = unique_reads
            top_reduction = reduction
            top_loo = loo

    total = len(a_union) + len(b_union)
    unique_fraction = top_unique_reads / float(total) if total else None
    loo_effect_persists = False
    if top_loo is not None and top_loo["total"] >= args.min_total:
        loo_effect_persists = (
            top_loo["mhf"] is not None and
            top_loo["mhf"] >= args.min_major_fraction and
            top_loo["p"] is not None and
            top_loo["p"] <= 0.05
        )
    if len(a_union) > len(b_union):
        full_direction = 1
    elif len(b_union) > len(a_union):
        full_direction = -1
    else:
        full_direction = 0
    usable_loo = [
        loo for loo in all_loo if loo["total"] >= args.min_total
    ]
    loo_direction_consistent = bool(all_loo) and full_direction != 0
    for loo in all_loo:
        loo_direction = (
            1 if loo["a"] > loo["b"] else
            -1 if loo["b"] > loo["a"] else 0
        )
        if loo_direction != full_direction:
            loo_direction_consistent = False
            break
    loo_all_robust = (
        len(contributing) > 1 and
        len(usable_loo) == len(all_loo) and
        loo_direction_consistent and
        all(
            loo["mhf"] is not None and
            loo["mhf"] >= args.min_major_fraction and
            loo["p"] is not None and
            loo["p"] <= 0.05
            for loo in all_loo
        )
    )
    if len(contributing) <= 1:
        support_class = "SINGLE_SNP_ONLY"
    elif (
        unique_fraction is not None and
        unique_fraction >= args.single_variant_driver_fraction and
        not loo_effect_persists
    ):
        support_class = "SINGLE_SNP_DOMINANT"
    elif loo_all_robust:
        support_class = "MULTI_SNP_LOO_ROBUST"
    else:
        support_class = "MULTI_SNP_SUPPORTED"

    return {
        "a": len(a_union),
        "b": len(b_union),
        "total": total,
        "n_gene_variants_with_reads": len(contributing),
        "gene_variants_with_reads": contributing,
        "support_class": support_class,
        "top_variant": top_variant,
        "top_variant_unique_reads": top_unique_reads,
        "top_variant_unique_fraction": unique_fraction,
        "top_loo": top_loo,
        "loo_all_robust": loo_all_robust,
        "loo_direction_consistent": loo_direction_consistent,
        "loo_min_total": min(
            (loo["total"] for loo in all_loo), default=None
        ),
        "loo_min_major_fraction": min(
            (loo["mhf"] for loo in all_loo if loo["mhf"] is not None),
            default=None
        ),
        "loo_max_p": max(
            (loo["p"] for loo in all_loo if loo["p"] is not None),
            default=None
        )
    }


def calculate_rescue_provenance_metrics(gene_count, args):
    """Summarize strand-required fragments for supplementary QC only."""
    rescue_a = set(gene_count["strand_resolved"]["A"])
    rescue_b = set(gene_count["strand_resolved"]["B"])
    rescue_union = rescue_a.union(rescue_b)
    direct_a = set(gene_count["direct"]["A"])
    direct_b = set(gene_count["direct"]["B"])
    direct_total = len(direct_a) + len(direct_b)
    direct_mhf = (
        max(len(direct_a), len(direct_b)) / float(direct_total)
        if direct_total else None
    )
    direct_p = (
        exact_binom_two_sided_half(len(direct_a), len(direct_b))
        if direct_total else None
    )
    full_a = set(gene_count["strict"]["A"])
    full_b = set(gene_count["strict"]["B"])
    full_direction = (
        1 if len(full_a) > len(full_b)
        else -1 if len(full_b) > len(full_a)
        else 0
    )
    direct_direction = (
        1 if len(direct_a) > len(direct_b)
        else -1 if len(direct_b) > len(direct_a)
        else 0
    )
    exclusion_robust = bool(
        rescue_union
        and direct_total >= args.min_total
        and direct_mhf is not None
        and direct_mhf >= args.min_major_fraction
        and direct_p is not None
        and direct_p <= 0.05
        and direct_direction == full_direction
        and full_direction != 0
    )
    if not rescue_union:
        dependency = "NO_STRAND_RESOLVED_FRAGMENTS"
    elif exclusion_robust:
        dependency = "NOT_DEPENDENT"
    else:
        dependency = "RESCUE_DEPENDENT_REVIEW"
    return {
        "strand_resolved_a": len(rescue_a),
        "strand_resolved_b": len(rescue_b),
        "strand_resolved_total": len(rescue_union),
        "direct_a": len(direct_a),
        "direct_b": len(direct_b),
        "direct_total": direct_total,
        "direct_major_haplotype_fraction": direct_mhf,
        "direct_p_exact": direct_p,
        "rescue_exclusion_robust": exclusion_robust,
        "rescue_dependency": dependency,
    }


def _flatten_filter_counts(values):
    return ",".join(
        "{}={}".format(key, value)
        for key, value in sorted(values.items())
    )


def build_bam_strict_candidate_rows(
    haplotypic_counts_path,
    sample,
    components,
    identity_to_component,
    genes,
    transcripts_by_gene,
    variant_assignment,
    graph_rows,
    internal_edges,
    haplotypes,
    bam_tissue_map,
    wgs_rows,
    phased_wgs_rows,
    wgs_balanced_variants,
    args,
):
    """Build v0.2.5 rows from WGS-filtered RNA QNAME fragments."""
    measurements = load_phaser_measurements(
        haplotypic_counts_path,
        sample,
        components,
        identity_to_component,
        bam_tissue_map,
    )
    if not measurements:
        raise ValueError("No component/BAM measurements were found")
    assignments = build_variant_assignment_index(
        components, variant_assignment
    )
    graph_by_component = {row["component_id"]: row for row in graph_rows}
    wgs_by_variant = {row["variant_id"]: row for row in wgs_rows}
    excluded = set(split_csv(args.exclude_tissues))
    manifest = bam_engine.load_bam_manifest(args.bam_manifest)
    as_cutoffs = bam_engine.load_phaser_as_cutoffs(args.phaser_log)
    grouped = defaultdict(list)
    for measurement in measurements:
        grouped[(measurement["bam"], measurement["tissue"])].append(
            measurement
        )

    fragment_handle = None
    fragment_writer = None
    audit_dir = os.path.join(args.out_dir, "audit")
    os.makedirs(audit_dir, exist_ok=True)
    fragment_tmp = os.path.join(
        audit_dir, "ase_fragment_assignment.tsv.tmp.gz"
    )
    fragment_path = os.path.join(
        audit_dir, "ase_fragment_assignment.tsv.gz"
    )
    if args.write_fragments:
        fragment_handle = gzip.open(fragment_tmp, "wt", newline="")
        fragment_writer = csv.DictWriter(
            fragment_handle,
            fieldnames=bam_engine.FRAGMENT_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        fragment_writer.writeheader()

    candidate_rows = []
    count_rows = []
    validation_rows = []
    fragment_status_totals = defaultdict(int)
    processed = 0
    fragment_success = False
    try:
        for group_key in sorted(grouped):
            local_measurements = grouped[group_key]
            bam_path, resolved_tissue, min_as = bam_engine.resolve_bam(
                local_measurements[0],
                sample,
                args.bam_root,
                manifest,
                as_cutoffs,
                args.min_as,
            )
            if resolved_tissue != local_measurements[0]["tissue"]:
                for measurement in local_measurements:
                    measurement["tissue"] = resolved_tissue
            sys.stderr.write(
                "[BAM] {} tissue={} components={}\n".format(
                    bam_path, resolved_tissue, len(local_measurements)
                )
            )
            for result in bam_engine.stream_bam_fragment_evidence(
                bam_path,
                local_measurements,
                components,
                assignments,
                genes,
                args,
                min_as=min_as,
                fragment_writer=fragment_writer,
                eligible_variants=wgs_balanced_variants,
            ):
                processed += 1
                if processed % 1000 == 0:
                    sys.stderr.write(
                        "[PROGRESS] reconstructed component/BAM "
                        "measurements={}\n".format(processed)
                    )
                measurement = result["measurement"]
                component = result["component"]
                comparison = result["comparison"]
                comparison["filter_counts"] = _flatten_filter_counts(
                    comparison.get("filter_counts", {})
                )
                validation_rows.append(comparison)
                graph = graph_by_component.get(component["component_id"])
                analysis_set = (
                    "EXCLUDED_TISSUE"
                    if measurement["tissue"] in excluded else "MAIN"
                )

                for gene_id in sorted(
                    component["candidate_gene_variants"]
                ):
                    gene = genes[gene_id]
                    gene_count = result["gene_counts"].get(
                        gene_id, empty_strict_gene_count()
                    )
                    metrics = calculate_strict_variant_metrics(
                        gene_count, args
                    )
                    rescue = calculate_rescue_provenance_metrics(
                        gene_count, args
                    )
                    a_count = metrics["a"]
                    b_count = metrics["b"]
                    total = metrics["total"]
                    minor = min(a_count, b_count)
                    major_fraction = (
                        max(a_count, b_count) / float(total)
                        if total else None
                    )
                    log2_ratio = math.log(
                        (a_count + 0.5) / float(b_count + 0.5), 2
                    )
                    allelic_imbalance = (
                        abs(a_count - b_count) / float(total)
                        if total else None
                    )
                    if total < args.min_total:
                        test_status = "LOW_COVERAGE"
                        pvalue = None
                    else:
                        test_status = "TESTED"
                        pvalue = exact_binom_two_sided_half(
                            a_count, b_count
                        )
                    if minor <= args.extreme_minor_max:
                        minor_support_class = "EXTREME_MINOR"
                    elif minor < args.min_each_haplotype:
                        minor_support_class = "LOW_MINOR"
                    else:
                        minor_support_class = "ADEQUATE_MINOR"

                    ase_unit = (
                        "SINGLE_SNP_ASE"
                        if metrics["n_gene_variants_with_reads"] <= 1
                        else "LOCAL_HAPLOTYPE_BLOCK_ASE"
                    )
                    gene_variants = component[
                        "candidate_gene_variants"
                    ][gene_id]
                    unique_variants = component["gene_variants"].get(
                        gene_id, []
                    )
                    ambiguous_variants = component[
                        "ambiguous_gene_variants"
                    ].get(gene_id, [])
                    relation = (
                        "ALL_COMPONENT_VARIANTS_EXONIC_FOR_GENE"
                        if set(gene_variants) == set(component["variants"])
                        else "GENE_SPECIFIC_EXONIC_SUBSET"
                    )
                    transcripts = transcripts_by_gene.get(gene_id, [])
                    canonical = [
                        transcript["transcript_id"]
                        for transcript in transcripts
                        if transcript["canonical"]
                    ]
                    measurement_key = "|".join(
                        [
                            sample,
                            measurement["tissue"],
                            normalize_bam(measurement["bam"]),
                            component["component_id"],
                            gene_id,
                        ]
                    )
                    measurement_id = "M" + hashlib.sha1(
                        measurement_key.encode("utf-8")
                    ).hexdigest()[:16]
                    existing_a = len(
                        gene_count["direct"]["A"]
                    )
                    existing_b = len(
                        gene_count["direct"]["B"]
                    )
                    rescue_fraction = (
                        rescue["strand_resolved_total"] / float(total)
                        if total else None
                    )
                    balanced_gene_variants = sorted(
                        set(gene_variants).intersection(
                            wgs_balanced_variants
                        )
                    )
                    wgs_failed_gene_variants = sorted(
                        set(gene_variants).difference(
                            wgs_balanced_variants
                        )
                    )
                    phase_agreement = wgs_snp_qc.phase_agreement(
                        component,
                        metrics["gene_variants_with_reads"],
                        phased_wgs_rows,
                    )
                    positions = sorted(
                        row["position"]
                        for row in component["parsed_variants"]
                    )
                    adjacent_distances = [
                        right - left
                        for left, right in zip(
                            positions[:-1], positions[1:]
                        )
                    ]
                    local_phase = assess_local_phase_for_snps(
                        component["component_id"],
                        metrics["gene_variants_with_reads"],
                        internal_edges,
                        haplotypes,
                        args.min_bridge_support,
                        args.max_conflict_fraction,
                    )
                    top_loo = metrics["top_loo"]
                    row = {
                        "sample": sample,
                        "tissue": measurement["tissue"],
                        "bam": measurement["bam"],
                        "gene_id": gene_id,
                        "gene_name": gene["gene_name"],
                        "gene_biotype": gene["gene_biotype"],
                        "gene_contig": gene["contig"],
                        "gene_start": gene["start"],
                        "gene_end": gene["end"],
                        "strand": gene["strand"],
                        "component_id": component["component_id"],
                        "haplotype_block_id": component["component_id"],
                        "component_start": component["start"],
                        "block_start": component["start"],
                        "component_end": component["end"],
                        "block_end": component["end"],
                        "block_span_bp": (
                            component["end"] - component["start"]
                        ),
                        "max_adjacent_snp_distance": (
                            max(adjacent_distances)
                            if adjacent_distances else None
                        ),
                        "component_gene_relation": relation,
                        "n_component_variants": len(
                            component["variants"]
                        ),
                        "n_block_snps": len(component["variants"]),
                        "component_variants": ",".join(
                            component["variants"]
                        ),
                        "block_snps": ",".join(component["variants"]),
                        "n_gene_variants": len(gene_variants),
                        "gene_variants": ",".join(gene_variants),
                        "n_unique_exonic_variants": len(unique_variants),
                        "unique_exonic_variants": ",".join(unique_variants),
                        "n_ambiguous_exonic_variants": len(
                            ambiguous_variants
                        ),
                        "ambiguous_exonic_variants": ",".join(
                            ambiguous_variants
                        ),
                        "gene_assignment_policy": (
                            "DIRECT_OR_STRAND_RESOLVED_EXON_UNION"
                        ),
                        "measurement_id": measurement_id,
                        "haplotypeA": measurement["haplotypeA"],
                        "haplotypeB": measurement["haplotypeB"],
                        "phaser_component_a_count": measurement[
                            "reported_a"
                        ],
                        "phaser_component_b_count": measurement[
                            "reported_b"
                        ],
                        "phaser_component_total_count": measurement[
                            "reported_total"
                        ],
                        "reconstructed_phaserlike_a": comparison[
                            "reconstructed_phaserlike_a"
                        ],
                        "reconstructed_phaserlike_b": comparison[
                            "reconstructed_phaserlike_b"
                        ],
                        "transcript_count": len(transcripts),
                        "transcript_ids": ",".join(
                            transcript["transcript_id"]
                            for transcript in transcripts
                        ),
                        "canonical_transcript_ids": ",".join(canonical),
                        "analysis_set": analysis_set,
                        "ase_unit_type_local": ase_unit,
                        "count_source": (
                            "WGS_BALANCED_SNP_RNA_FRAGMENT_UNION"
                        ),
                        "count_validation": (
                            "MATCH" if comparison["count_match_exact"]
                            else "MISMATCH"
                        ),
                        "existing_unique_a": existing_a,
                        "existing_unique_b": existing_b,
                        "existing_unique_total": existing_a + existing_b,
                        "strand_resolved_a": rescue["strand_resolved_a"],
                        "strand_resolved_b": rescue["strand_resolved_b"],
                        "strand_resolved_total": rescue[
                            "strand_resolved_total"
                        ],
                        "strict_a_count": a_count,
                        "strict_b_count": b_count,
                        "strict_total_count": total,
                        "rescue_fraction": rescue_fraction,
                        "rescue_dependency": rescue["rescue_dependency"],
                        "rescue_exclusion_robust": rescue[
                            "rescue_exclusion_robust"
                        ],
                        "direct_a_count": rescue["direct_a"],
                        "direct_b_count": rescue["direct_b"],
                        "direct_total_count": rescue["direct_total"],
                        "fragment_hap_conflict_fraction": comparison[
                            "analysis_fragment_hap_conflict_fraction"
                        ],
                        "n_block_haplotype_conflicts": comparison[
                            "analysis_haplotype_conflict"
                        ],
                        "n_wgs_balanced_gene_snps": len(
                            balanced_gene_variants
                        ),
                        "wgs_balanced_gene_snps": ",".join(
                            balanced_gene_variants
                        ),
                        "n_wgs_failed_gene_snps": len(
                            wgs_failed_gene_variants
                        ),
                        "wgs_failed_gene_snps": ",".join(
                            wgs_failed_gene_variants
                        ),
                        "ase_a_count": a_count,
                        "ase_b_count": b_count,
                        "ase_total_count": total,
                        "major_haplotype_fraction": major_fraction,
                        "log2_a_over_b": log2_ratio,
                        "absolute_allelic_imbalance": allelic_imbalance,
                        "ase_test_status": test_status,
                        "minor_support_class": minor_support_class,
                        "ase_p_exact": pvalue,
                        "ase_q_bh": None,
                        "n_gene_variants_with_reads": metrics[
                            "n_gene_variants_with_reads"
                        ],
                        "gene_variants_with_reads": ",".join(
                            metrics["gene_variants_with_reads"]
                        ),
                        "snp_support_status": metrics["support_class"],
                        "top_count_driver_variant": metrics[
                            "top_variant"
                        ],
                        "top_driver_unique_reads": metrics[
                            "top_variant_unique_reads"
                        ],
                        "top_driver_unique_fraction": metrics[
                            "top_variant_unique_fraction"
                        ],
                        "loo_a_count": (
                            top_loo["a"] if top_loo else None
                        ),
                        "loo_b_count": (
                            top_loo["b"] if top_loo else None
                        ),
                        "loo_total_count": (
                            top_loo["total"] if top_loo else None
                        ),
                        "loo_major_haplotype_fraction": (
                            top_loo["mhf"] if top_loo else None
                        ),
                        "loo_p_exact": (
                            top_loo["p"] if top_loo else None
                        ),
                        "loo_all_robust": metrics["loo_all_robust"],
                        "loo_direction_consistent": metrics[
                            "loo_direction_consistent"
                        ],
                        "loo_min_total": metrics["loo_min_total"],
                        "loo_min_major_fraction": metrics[
                            "loo_min_major_fraction"
                        ],
                        "loo_max_p": metrics["loo_max_p"],
                        "phase_qc": (
                            graph["phase_qc"]
                            if graph else "PHASE_LC"
                        ),
                        "phase_qc_reason": (
                            graph["phase_qc_reason"]
                            if graph else "GRAPH_QC_MISSING"
                        ),
                        "graph_reconstruction": (
                            graph["graph_reconstruction"]
                            if graph else "MISSING"
                        ),
                        **local_phase,
                        "min_bridge_support": (
                            graph["min_bridge_support"]
                            if graph else None
                        ),
                        **phase_agreement
                    }
                    candidate_rows.append(row)
                    count_rows.append(
                        {
                            key: row[key]
                            for key in (
                                "sample", "tissue", "bam", "component_id",
                                "haplotype_block_id",
                                "gene_id", "gene_name", "strand",
                                "direct_a_count", "direct_b_count",
                                "direct_total_count",
                                "strand_resolved_a", "strand_resolved_b",
                                "strand_resolved_total",
                                "strict_a_count", "strict_b_count",
                                "strict_total_count", "rescue_fraction",
                                "rescue_dependency",
                                "rescue_exclusion_robust"
                            )
                        }
                    )
        fragment_success = True
    finally:
        if fragment_handle is not None:
            fragment_handle.close()
            if fragment_success:
                os.replace(fragment_tmp, fragment_path)
            elif os.path.exists(fragment_tmp):
                os.remove(fragment_tmp)

    apply_unique_measurement_bh(candidate_rows)
    thresholds = {
        "fdr": args.fdr,
        "min_major_fraction": args.min_major_fraction,
        "max_fragment_hap_conflict_fraction": (
            args.max_fragment_hap_conflict_fraction
        )
    }
    for row in candidate_rows:
        graph = graph_by_component.get(row["component_id"])
        (
            row["candidate_class"],
            row["exclusion_reasons"],
        ) = classify_candidate(row, graph, thresholds)
        row["candidate_tier"] = row["candidate_class"]
        row["is_primary_candidate"] = row["candidate_class"] in {
            "PRIMARY_ROBUST", "PRIMARY_SUPPORTED"
        }
        fragment_status_totals[row["candidate_class"]] += 1
    return (
        candidate_rows,
        count_rows,
        validation_rows,
        measurements,
        processed,
    )


def parser():
    p = argparse.ArgumentParser(
        description=(
            "WGS-filtered read-linked SNP-block ASE caller: GTF + phASER "
            "outputs + indexed RNA BAM + phased WGS SNP QC"
        )
    )
    p.add_argument("--sample", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--variant-connections", required=True)
    p.add_argument("--haplotypes", required=True)
    p.add_argument("--haplotypic-counts", required=True)
    p.add_argument("--phased-wgs-vcf", required=True)
    p.add_argument("--phaser-wgs-allelic-counts", required=True)
    p.add_argument("--wgs-variant-connections")
    p.add_argument("--full-allele-config")
    p.add_argument("--wgs-allele-config")
    p.add_argument("--bam-tissue-map")
    p.add_argument("--bam-root")
    p.add_argument("--bam-manifest")
    p.add_argument("--phaser-log")
    p.add_argument(
        "--library-strand",
        choices=("reverse", "forward"),
        default="reverse"
    )
    p.add_argument("--min-mapq", type=int, default=255)
    p.add_argument("--min-baseq", type=int, default=10)
    p.add_argument("--max-isize", type=int, default=1000000)
    p.add_argument("--min-as", type=float)
    p.add_argument(
        "--require-proper-pair", type=int, choices=(0, 1), default=1
    )
    p.add_argument(
        "--exclude-duplicates", type=int, choices=(0, 1), default=1
    )
    p.add_argument(
        "--exclude-qcfail", type=int, choices=(0, 1), default=1
    )
    p.add_argument(
        "--allow-stale-index", type=int, choices=(0, 1), default=0
    )
    p.add_argument("--bam-threads", type=int, default=2)
    p.add_argument(
        "--merge-locus-gap",
        type=int,
        default=1000,
        help=(
            "Merge component fetch intervals separated by at most this "
            "many bp; affects speed only, not target SNPs"
        )
    )
    p.add_argument(
        "--write-fragments", type=int, choices=(0, 1), default=0
    )
    p.add_argument("--gene-biotypes", default="protein_coding")
    p.add_argument(
        "--write-transcript-annotation",
        type=int,
        choices=(0, 1),
        default=0
    )
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cc-threshold", type=float, default=0.05)
    p.add_argument("--reference-cc-threshold", type=float, default=0.01)
    p.add_argument("--min-total", type=int, default=15)
    p.add_argument("--min-each-haplotype", type=int, default=3)
    p.add_argument("--extreme-minor-max", type=int, default=1)
    p.add_argument("--min-major-fraction", type=float, default=0.65)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min-bridge-support", type=int, default=3)
    p.add_argument("--max-conflict-fraction", type=float, default=0.10)
    p.add_argument("--single-variant-driver-fraction", type=float, default=0.80)
    p.add_argument(
        "--max-fragment-hap-conflict-fraction",
        type=float,
        default=0.10
    )
    p.add_argument("--exclude-tissues", default="Cranial")
    p.add_argument("--min-wgs-gq", type=float, default=20.0)
    p.add_argument("--min-wgs-total", type=int, default=15)
    p.add_argument("--min-wgs-each", type=int, default=5)
    p.add_argument("--min-wgs-alt-fraction", type=float, default=0.30)
    p.add_argument("--max-wgs-alt-fraction", type=float, default=0.70)
    p.add_argument("--force", action="store_true")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    args.require_proper_pair = bool(args.require_proper_pair)
    args.exclude_duplicates = bool(args.exclude_duplicates)
    args.exclude_qcfail = bool(args.exclude_qcfail)
    args.allow_stale_index = bool(args.allow_stale_index)
    args.write_fragments = bool(args.write_fragments)
    if not (
        0.0 <= args.min_wgs_alt_fraction
        < args.max_wgs_alt_fraction <= 1.0
    ):
        raise SystemExit(
            "[ERROR] WGS ALT-fraction bounds must satisfy 0 <= min < max <= 1"
        )
    if not args.bam_root and not args.bam_manifest:
        raise SystemExit(
            "[ERROR] v0.2.5 requires --bam-root or --bam-manifest"
        )
    args.out_dir = prepare_output_directory(args.out_dir, args.force)

    selected_biotypes = split_csv(args.gene_biotypes)
    genes, transcripts_by_gene = load_gtf_annotation(
        args.gtf, selected_biotypes
    )
    if not genes:
        raise ValueError(
            "No selected genes found in GTF for biotypes: {}".format(
                ",".join(selected_biotypes)
            )
        )
    all_components, identity_to_component = scan_phaser_components(
        args.haplotypic_counts, args.sample
    )
    components, variant_assignment, invalid_variants = map_components_to_genes(
        all_components, genes
    )
    if not components:
        raise ValueError(
            "No phASER variants overlap selected GTF genes. Check contig naming."
        )
    identity_to_component = {
        component["identity"]: component_id
        for component_id, component in components.items()
    }
    prepare_component_haplotype_alleles(components)
    variant_to_components = catalog_from_components(components)
    requested_wgs_variants = {}
    for variant_id in variant_to_components:
        parsed = wgs_snp_qc.parse_variant_id(variant_id)
        if parsed is None:
            raise ValueError(
                "Cannot parse WGS target SNP {}".format(variant_id)
            )
        requested_wgs_variants[variant_id] = parsed
    phased_wgs_rows = wgs_snp_qc.load_phased_wgs_vcf(
        args.phased_wgs_vcf,
        args.sample,
        requested_wgs_variants,
    )
    phaser_wgs_counts = wgs_snp_qc.load_phaser_wgs_counts(
        args.phaser_wgs_allelic_counts,
        requested_wgs_variants,
    )
    wgs_rows, wgs_balanced_variants = wgs_snp_qc.evaluate_wgs_snps(
        requested_wgs_variants,
        phased_wgs_rows,
        phaser_wgs_counts,
        args.min_wgs_gq,
        args.min_wgs_total,
        args.min_wgs_each,
        args.min_wgs_alt_fraction,
        args.max_wgs_alt_fraction,
    )
    haplotypes = load_haplotypes(args.haplotypes, components)
    internal_edges, cross_edges, full_pairs = load_connections(
        args.variant_connections,
        variant_to_components,
        components,
        args.cc_threshold
    )
    relevant_variants = set(variant_to_components)
    wgs_pairs = load_relevant_connection_pairs(
        args.wgs_variant_connections,
        relevant_variants,
        args.cc_threshold
    )
    full_configs, full_config_qc = load_allele_config(
        args.full_allele_config, full_pairs
    )
    wgs_configs, wgs_config_qc = load_allele_config(
        args.wgs_allele_config, full_pairs
    )
    graph_rows, edge_rows = build_graph_qc(
        components,
        internal_edges,
        haplotypes,
        wgs_pairs,
        full_configs,
        wgs_configs,
        args.cc_threshold,
        args.reference_cc_threshold,
        args.min_bridge_support,
        args.max_conflict_fraction
    )
    bam_tissue_map = load_bam_tissue_map(args.bam_tissue_map)
    (
        candidate_rows,
        gene_component_count_rows,
        component_count_validation_rows,
        phaser_measurements,
        n_bam_measurements_processed,
    ) = build_bam_strict_candidate_rows(
        args.haplotypic_counts,
        args.sample,
        components,
        identity_to_component,
        genes,
        transcripts_by_gene,
        variant_assignment,
        graph_rows,
        internal_edges,
        haplotypes,
        bam_tissue_map,
        wgs_rows,
        phased_wgs_rows,
        wgs_balanced_variants,
        args
    )
    measurement_rows = build_measurement_rows(candidate_rows)
    assignment_rows = build_variant_assignment_rows(
        all_components, variant_assignment, genes
    )
    gene_rows = build_gene_summary(
        candidate_rows, args.fdr, args.min_major_fraction
    )
    cross_rows = cross_edge_rows(
        cross_edges,
        wgs_pairs,
        full_configs,
        wgs_configs,
        args.cc_threshold,
        args.reference_cc_threshold
    )
    graph_fields = [
        "component_id", "contig", "component_start", "component_end",
        "phaser_pi", "n_variants", "n_edges_observed", "n_edges_retained",
        "n_edges_rejected", "n_retained_connected_groups",
        "retained_graph_connected", "cycle_rank", "n_bridges",
        "min_bridge_support", "n_single_read_bridges",
        "min_retained_edge_support", "median_retained_edge_support",
        "min_retained_support_fraction",
        "weighted_retained_conflict_fraction",
        "haplotypes_edges_supporting", "haplotypes_edges_total",
        "haplotype_edge_configuration_match",
        "graph_reconstruction", "phase_qc", "phase_qc_reason",
        "component_variants"
    ]
    edge_fields = [
        "component_id", "variant_a", "variant_b",
        "supporting_connections", "total_connections", "conflicting_reads",
        "support_fraction", "conflicting_configuration_p", "cc_threshold",
        "edge_status", "is_bridge", "full_configuration",
        "wgs_supporting_connections", "wgs_total_connections",
        "wgs_conflicting_configuration_p", "wgs_configuration",
        "evidence_class"
    ]
    candidate_fields = [
        "sample", "tissue", "bam", "gene_id", "gene_name",
        "gene_biotype", "gene_contig", "gene_start", "gene_end", "strand",
        "haplotype_block_id", "block_start", "block_end",
        "block_span_bp", "max_adjacent_snp_distance",
        "n_block_snps", "block_snps", "n_gene_variants",
        "gene_variants", "measurement_id", "haplotypeA", "haplotypeB",
        "analysis_set", "ase_unit_type_local", "count_source",
        "count_validation", "n_wgs_balanced_gene_snps",
        "wgs_balanced_gene_snps", "n_wgs_failed_gene_snps",
        "wgs_failed_gene_snps",
        "fragment_hap_conflict_fraction",
        "n_block_haplotype_conflicts",
        "ase_a_count", "ase_b_count", "ase_total_count",
        "major_haplotype_fraction", "log2_a_over_b",
        "absolute_allelic_imbalance", "ase_test_status",
        "minor_support_class", "ase_p_exact",
        "ase_q_bh", "n_gene_variants_with_reads",
        "gene_variants_with_reads", "snp_support_status",
        "top_count_driver_variant", "top_driver_unique_reads",
        "top_driver_unique_fraction", "loo_a_count", "loo_b_count",
        "loo_total_count", "loo_major_haplotype_fraction", "loo_p_exact",
        "loo_all_robust", "loo_direction_consistent", "loo_min_total",
        "loo_min_major_fraction", "loo_max_p",
        "local_phase_status", "local_phase_reason", "local_phase_n_snps",
        "local_phase_n_edges", "local_phase_min_bridge_support",
        "local_phase_conflict_fraction", "wgs_phase_status",
        "wgs_phase_orientation", "wgs_phase_set",
        "n_wgs_phase_compared_snps", "n_wgs_phase_direct_snps",
        "n_wgs_phase_flipped_snps", "n_wgs_phase_mismatch_snps",
        "candidate_class", "is_primary_candidate", "exclusion_reasons"
    ]
    gene_fields = [
        "sample", "tissue", "gene_id", "gene_name", "gene_biotype",
        "n_informative_haplotype_blocks", "n_testable_haplotype_blocks",
        "n_q_effect_pass_blocks", "n_primary_ase_blocks",
        "n_primary_robust_blocks", "n_primary_supported_blocks",
        "n_secondary_single_snp_blocks", "n_extreme_minor_blocks",
        "gene_phase_scope", "gene_ase_status", "top_haplotype_block_id",
        "top_candidate_class", "top_major_haplotype_fraction",
        "top_total_count", "primary_haplotype_block_ids"
    ]
    cross_fields = [
        "component_a", "component_b", "shared_gene_ids",
        "variant_a", "variant_b", "supporting_connections",
        "total_connections", "conflicting_reads", "support_fraction",
        "conflicting_configuration_p", "edge_status",
        "full_configuration", "wgs_configuration", "evidence_class"
    ]
    component_gene_fields = [
        "component_id", "contig", "component_start", "component_end",
        "n_component_variants", "component_variants", "gene_id",
        "gene_name", "gene_biotype", "gene_start", "gene_end", "strand",
        "component_gene_relation", "n_gene_variants", "gene_variants",
        "n_unique_exonic_variants", "unique_exonic_variants",
        "n_ambiguous_exonic_variants", "ambiguous_exonic_variants",
        "gene_assignment_policy", "n_bam_rows"
    ]
    transcript_fields = [
        "gene_id", "gene_name", "variant_id", "contig", "position",
        "transcript_id", "transcript_name", "transcript_biotype",
        "canonical_transcript", "transcript_feature"
    ]
    measurement_fields = [
        "measurement_id", "sample", "tissue", "bam",
        "haplotype_block_id",
        "gene_ids", "n_gene_annotations", "gene_variants",
        "strict_a_count", "strict_b_count", "strict_total_count",
        "direct_total_count", "strand_resolved_total",
        "rescue_fraction", "rescue_dependency",
        "ase_a_count", "ase_b_count", "ase_total_count", "ase_p_exact",
        "ase_q_bh", "minor_support_class", "candidate_class",
        "is_primary_candidate"
    ]
    count_fields = [
        "sample", "tissue", "bam", "haplotype_block_id", "gene_id",
        "gene_name", "strand", "direct_a_count", "direct_b_count",
        "direct_total_count", "strand_resolved_a",
        "strand_resolved_b", "strand_resolved_total", "strict_a_count",
        "strict_b_count", "strict_total_count", "rescue_fraction",
        "rescue_dependency", "rescue_exclusion_robust"
    ]
    validation_fields = [
        "sample", "tissue", "bam", "bam_path", "component_id",
        "contig", "start", "stop", "n_variants", "reported_a",
        "reported_b", "reported_total", "reconstructed_phaserlike_a",
        "reconstructed_phaserlike_b", "reconstructed_phaserlike_total",
        "delta_a", "delta_b", "count_match_exact",
        "n_qnames_with_variant_support", "n_qnames_assigned",
        "n_qnames_with_strand_resolution_evidence",
        "n_qnames_strand_resolved", "n_haplotype_conflict",
        "fragment_hap_conflict_fraction", "n_gene_conflict",
        "analysis_qnames_with_variant_support",
        "analysis_qnames_assigned", "analysis_haplotype_conflict",
        "analysis_fragment_hap_conflict_fraction",
        "n_strand_conflict", "n_annotated_exon_strand_mismatch",
        "n_no_gene", "n_unresolved",
        "n_alignments_examined_in_locus",
        "n_allele_observations_in_locus", "min_as_applied",
        "filter_counts"
    ]
    wgs_snp_fields = [
        "variant_id", "contig", "position", "ref", "alt",
        "vcf_filter", "vcf_gt", "vcf_gq", "vcf_phase_set",
        "wgs_ref_count", "wgs_alt_count", "wgs_total_count",
        "wgs_alt_fraction", "wgs_snp_status", "wgs_snp_reason"
    ]
    assignment_fields = [
        "variant_id", "contig", "position", "assignment_class",
        "exonic_gene_ids", "exonic_gene_names", "gene_span_ids",
        "component_ids"
    ]
    component_gene_rows = []
    for component_id in sorted(components):
        component = components[component_id]
        for gene_id in sorted(component["candidate_gene_variants"]):
            gene = genes[gene_id]
            gene_variants = component["candidate_gene_variants"][gene_id]
            unique_variants = component["gene_variants"].get(gene_id, [])
            ambiguous_variants = component[
                "ambiguous_gene_variants"
            ].get(gene_id, [])
            component_gene_rows.append({
                "component_id": component_id,
                "contig": component["contig"],
                "component_start": component["start"],
                "component_end": component["end"],
                "n_component_variants": len(component["variants"]),
                "component_variants": ",".join(component["variants"]),
                "gene_id": gene_id,
                "gene_name": gene["gene_name"],
                "gene_biotype": gene["gene_biotype"],
                "gene_start": gene["start"],
                "gene_end": gene["end"],
                "strand": gene["strand"],
                "component_gene_relation": (
                    "ALL_COMPONENT_VARIANTS_EXONIC_FOR_GENE"
                    if set(gene_variants) == set(component["variants"])
                    else "GENE_SPECIFIC_EXONIC_SUBSET"
                ),
                "n_gene_variants": len(gene_variants),
                "gene_variants": ",".join(gene_variants),
                "n_unique_exonic_variants": len(unique_variants),
                "unique_exonic_variants": ",".join(unique_variants),
                "n_ambiguous_exonic_variants": len(ambiguous_variants),
                "ambiguous_exonic_variants": ",".join(
                    ambiguous_variants
                ),
                "gene_assignment_policy": (
                    "STRAND_RESOLVED_PROTEIN_CODING_EXON_UNION"
                ),
                "n_bam_rows": component["n_bam_rows"]
            })

    audit_dir = os.path.join(args.out_dir, "audit")
    os.makedirs(audit_dir, exist_ok=True)
    primary_rows = [
        row for row in candidate_rows if row["is_primary_candidate"]
    ]
    primary_gene_rows = [
        row for row in gene_rows
        if row["gene_ase_status"] in {
            "ASE_GENE_ROBUST", "ASE_GENE_SUPPORTED"
        }
    ]

    # Four manuscript-facing tables.
    write_tsv(
        os.path.join(args.out_dir, "ase_haplotype_block_results.tsv"),
        candidate_rows, candidate_fields
    )
    write_tsv(
        os.path.join(args.out_dir, "ase_primary_haplotype_blocks.tsv"),
        primary_rows, candidate_fields
    )
    write_tsv(
        os.path.join(args.out_dir, "ase_gene_summary.tsv"),
        gene_rows, gene_fields
    )
    write_tsv(
        os.path.join(args.out_dir, "ase_primary_genes.tsv"),
        primary_gene_rows, gene_fields
    )

    # Detailed provenance and reconstruction checks remain supplementary.
    write_tsv(
        os.path.join(audit_dir, "ase_local_phase_graph_qc.tsv"),
        graph_rows, graph_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_local_phase_edge_qc.tsv"),
        edge_rows, edge_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_cross_block_edges.tsv"),
        cross_rows, cross_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_strand_assignment_qc.tsv"),
        gene_component_count_rows, count_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_count_validation.tsv"),
        component_count_validation_rows, validation_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_measurement_unique.tsv"),
        measurement_rows, measurement_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_snp_qc.tsv"),
        wgs_rows, wgs_snp_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_variant_gene_assignment.tsv"),
        assignment_rows, assignment_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_ambiguous_gene_assignment.tsv"),
        [
            row for row in assignment_rows
            if row["assignment_class"] == "AMBIGUOUS_SHARED_EXON"
        ],
        assignment_fields
    )
    write_tsv(
        os.path.join(audit_dir, "ase_intronic_support.tsv"),
        [
            row for row in assignment_rows
            if row["assignment_class"] == "INTRONIC_ONLY"
        ],
        assignment_fields
    )
    write_tsv(
        os.path.join(audit_dir, "phaser_haplotype_block_gene_map.tsv"),
        component_gene_rows, component_gene_fields
    )
    if args.write_transcript_annotation:
        transcript_row_count = write_transcript_annotation_stream(
            os.path.join(
                audit_dir, "ase_variant_transcript_annotation.tsv"
            ),
            components,
            genes,
            transcripts_by_gene,
            transcript_fields
        )
    else:
        transcript_row_count = 0

    class_counts = defaultdict(int)
    for row in candidate_rows:
        class_counts[row["candidate_class"]] += 1
    phase_counts = defaultdict(int)
    for row in graph_rows:
        phase_counts[row["phase_qc"]] += 1
    summary = {
        "version": VERSION,
        "sample": args.sample,
        "parameters": vars(args),
        "n_gtf_selected_genes": len(genes),
        "n_gtf_selected_transcripts": sum(
            len(items) for items in transcripts_by_gene.values()
        ),
        "n_phaser_haplotype_blocks_all": len(all_components),
        "n_gene_overlapping_haplotype_blocks": len(components),
        "n_haplotype_block_gene_rows": len(component_gene_rows),
        "n_haplotype_block_result_rows": len(candidate_rows),
        "n_primary_haplotype_block_rows": len(primary_rows),
        "n_unique_measurements": len(measurement_rows),
        "n_phaser_block_bam_measurements": len(phaser_measurements),
        "n_bam_measurements_processed": n_bam_measurements_processed,
        "n_count_validation_exact_match": sum(
            row["count_match_exact"]
            for row in component_count_validation_rows
        ),
        "count_validation_exact_match_fraction": (
            sum(
                row["count_match_exact"]
                for row in component_count_validation_rows
            ) / float(len(component_count_validation_rows))
            if component_count_validation_rows else None
        ),
        "n_strand_resolved_fragments_across_gene_rows": sum(
            row["strand_resolved_total"]
            for row in gene_component_count_rows
        ),
        "n_analysis_fragment_haplotype_conflict": sum(
            row["analysis_haplotype_conflict"]
            for row in component_count_validation_rows
        ),
        "n_fragment_gene_conflict": sum(
            row["n_gene_conflict"]
            for row in component_count_validation_rows
        ),
        "n_fragment_strand_conflict": sum(
            row["n_strand_conflict"]
            for row in component_count_validation_rows
        ),
        "n_fragments_with_strand_resolution_evidence": sum(
            row["n_qnames_with_strand_resolution_evidence"]
            for row in component_count_validation_rows
        ),
        "n_annotated_exon_strand_mismatch_fragments": sum(
            row["n_annotated_exon_strand_mismatch"]
            for row in component_count_validation_rows
        ),
        "n_ambiguous_shared_exon_variants": sum(
            row["assignment_class"] == "AMBIGUOUS_SHARED_EXON"
            for row in assignment_rows
        ),
        "n_intronic_only_variants": sum(
            row["assignment_class"] == "INTRONIC_ONLY"
            for row in assignment_rows
        ),
        "n_graph_rows": len(graph_rows),
        "n_edge_rows": len(edge_rows),
        "n_cross_component_edges": len(cross_rows),
        "n_gene_summary_rows": len(gene_rows),
        "n_primary_gene_rows": len(primary_gene_rows),
        "n_wgs_snp_qc_rows": len(wgs_rows),
        "n_wgs_balanced_snps": len(wgs_balanced_variants),
        "n_transcript_annotation_rows": transcript_row_count,
        "n_invalid_variant_ids": len(invalid_variants),
        "invalid_variant_ids_first_20": invalid_variants[:20],
        "allele_config_qc": {
            "full": full_config_qc,
            "wgs": wgs_config_qc
        },
        "candidate_class_counts": dict(sorted(class_counts.items())),
        "phase_qc_counts": dict(sorted(phase_counts.items())),
        "primary_calling_uses": [
            "protein-coding transcript exon-union and gene strand",
            "indexed RNA BAM QNAME fragment reconstruction",
            "WGS-balanced heterozygous SNPs only",
            "direct exonic or uniquely strand-resolved gene assignment",
            "one HapA/HapB count per RNA fragment",
            "connected local phASER variant graph",
            "two-sided exact binomial test",
            "BH FDR within each sample-by-tissue family",
            "adequate support for the minor haplotype",
            "multi-SNP support"
        ],
        "primary_calling_does_not_use": [
            "WGS phase orientation as a significance gate",
            "single-SNP-only observations",
            "extreme or low minor-haplotype support",
            "gene_ae.txt"
        ],
        "wgs_phase_interpretation": [
            "WGS allelic balance is an ASE eligibility filter",
            "WGS phase agreement is reported for direction interpretation",
            "WGS phase agreement is not a primary-candidate gate"
        ]
    }
    with open(os.path.join(args.out_dir, "ase_v025_qc.json"), "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print("[OK] {}".format(args.out_dir))
    print("[OK] genes={} components={} rows={} primary_candidates={}".format(
        len(genes),
        len(components),
        len(candidate_rows),
        sum(1 for row in candidate_rows if row["is_primary_candidate"])
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
