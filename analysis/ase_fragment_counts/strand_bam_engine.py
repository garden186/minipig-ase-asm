#!/usr/bin/env python3
"""BAM/QNAME engine for WGS-filtered phASER local ASE.

The engine is deliberately independent of candidate statistics.  It performs
only the evidence-preserving operations that must precede an ASE test:

* reconstruct HapA/HapB support from RNA-seq BAM alignments;
* combine R1/R2 and all SNP observations sharing a QNAME into one fragment;
* require one haplotype, one inferred transcript strand and one gene;
* resolve opposite-strand shared exons when exactly one gene strand matches;
* distinguish direct gene assignment from assignment that requires strand;
* optionally reconstruct a second fragment set restricted to WGS-balanced SNPs.

Read-linked SNP blocks are grouped into non-overlapping genomic loci for each
BAM. Each locus is fetched once, avoiding a BAM fetch for every block.
"""

from __future__ import print_function

import csv
import gzip
import os
import re
from collections import defaultdict

try:
    import pysam
except ImportError:
    pysam = None


FRAGMENT_FIELDS = [
    "sample",
    "tissue",
    "bam",
    "haplotype_block_id",
    "qname",
    "transcript_strand",
    "gene_id",
    "gene_name",
    "strand_mismatch_gene_ids",
    "strand_mismatch_gene_names",
    "haplotype",
    "gene_assignment_method",
    "strand_resolution_required",
    "final_status",
    "supported_variants",
    "gene_supported_variants",
    "ambiguous_supported_variants",
    "observations",
]


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def split_csv(value):
    if value is None:
        return []
    return [
        item for item in str(value).split(",")
        if item not in ("", ".", "NA")
    ]


def normalize_bam(value):
    name = os.path.basename(str(value or ""))
    for suffix in (".bam", ".cram"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def load_bam_manifest(path):
    """Return BAM-label -> path/tissue/min-AS metadata."""
    manifest = {}
    if not path:
        return manifest
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"bam", "path"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("--bam-manifest requires bam and path columns")
        for row in reader:
            label = normalize_bam(row.get("bam"))
            manifest[label] = {
                "path": row.get("path", ""),
                "tissue": row.get("tissue", ""),
                "min_as": (
                    None
                    if row.get("min_as", "") in ("", ".", "NA")
                    else float(row["min_as"])
                ),
            }
    return manifest


def load_phaser_as_cutoffs(path):
    """Parse per-BAM alignment-score cutoffs from a phASER log."""
    cutoffs = {}
    if not path:
        return cutoffs
    current_bam = None
    file_pattern = re.compile(r"^\s*file:\s+(.+?)\s*$")
    cutoff_pattern = re.compile(
        r"using alignment score cutoff of\s+(-?[0-9]+(?:\.[0-9]+)?)"
    )
    with open_text(path) as handle:
        for line in handle:
            file_match = file_pattern.search(line)
            if file_match:
                current_bam = normalize_bam(file_match.group(1))
                continue
            cutoff_match = cutoff_pattern.search(line)
            if cutoff_match and current_bam:
                cutoffs[current_bam] = float(cutoff_match.group(1))
    return cutoffs


def resolve_bam(
    measurement,
    sample,
    bam_root,
    manifest,
    phaser_as_cutoffs,
    global_min_as,
):
    """Resolve the physical BAM path without guessing outside the project."""
    label = normalize_bam(measurement["bam"])
    entry = manifest.get(label)
    if entry:
        path = entry["path"]
        tissue = entry.get("tissue") or measurement["tissue"]
        if entry.get("min_as") is not None:
            min_as = entry["min_as"]
        else:
            min_as = phaser_as_cutoffs.get(label, global_min_as)
        return path, tissue, min_as
    if not bam_root:
        raise ValueError(
            "BAM {} is not in --bam-manifest and --bam-root was not given"
            .format(label)
        )
    tissue = measurement["tissue"]
    path = os.path.join(bam_root, tissue, sample, label + ".bam")
    return path, tissue, phaser_as_cutoffs.get(label, global_min_as)


def resolve_bam_contig(bam, contig):
    refs = set(bam.references)
    if contig in refs:
        return contig
    if contig.startswith("chr") and contig[3:] in refs:
        return contig[3:]
    candidate = "chr" + contig
    if candidate in refs:
        return candidate
    return None


def inferred_transcript_strand(read, library_strand):
    """Infer RNA transcript strand from a paired-end alignment."""
    if not read.is_paired:
        return None
    alignment_strand = "-" if read.is_reverse else "+"
    if library_strand == "reverse":
        if read.is_read1:
            return "+" if alignment_strand == "-" else "-"
        if read.is_read2:
            return alignment_strand
    elif library_strand == "forward":
        if read.is_read1:
            return alignment_strand
        if read.is_read2:
            return "+" if alignment_strand == "-" else "-"
    return None


def read_passes_filters(read, options, min_as):
    if read.is_unmapped:
        return False, "UNMAPPED"
    if read.is_secondary:
        return False, "SECONDARY"
    if read.is_supplementary:
        return False, "SUPPLEMENTARY"
    if options.exclude_qcfail and read.is_qcfail:
        return False, "QCFAIL"
    if options.exclude_duplicates and read.is_duplicate:
        return False, "DUPLICATE"
    if options.require_proper_pair and not read.is_proper_pair:
        return False, "NOT_PROPER_PAIR"
    if read.mapping_quality < options.min_mapq:
        return False, "LOW_MAPQ"
    if (
        options.max_isize > 0
        and abs(int(read.template_length)) > options.max_isize
    ):
        return False, "ISIZE"
    if min_as is not None:
        try:
            alignment_score = float(read.get_tag("AS"))
        except (KeyError, ValueError):
            alignment_score = 0.0
        if alignment_score < min_as:
            return False, "LOW_AS"
    return True, "PASS"


def new_fragment_record():
    return {
        "hap_support": set(),
        "gene_evidence": set(),
        "unique_gene_evidence": set(),
        "strand_resolved_gene_evidence": set(),
        "exon_strand_mismatch_gene_evidence": set(),
        "gene_variant_evidence": defaultdict(set),
        "supported_variants": set(),
        "ambiguous_supported_variants": set(),
        "observation_strings": [],
        "transcript_strands": set(),
    }


def add_gene_evidence(
    fragment,
    variant_id,
    assignment,
    transcript_strand,
    genes,
):
    """Add strand-concordant gene evidence for one SNP observation."""
    assignment_class = assignment.get("assignment_class", "")
    gene_ids = assignment.get("gene_ids", [])
    if assignment_class == "UNIQUE_EXONIC" and len(gene_ids) == 1:
        gene_id = gene_ids[0]
        gene_strand = genes.get(gene_id, {}).get("strand", "")
        if gene_strand not in ("+", "-"):
            return "", "ANNOTATED_GENE_STRAND_UNAVAILABLE"
        if transcript_strand != gene_strand:
            fragment["exon_strand_mismatch_gene_evidence"].add(gene_id)
            return "", "ANNOTATED_EXON_STRAND_MISMATCH"
        fragment["gene_evidence"].add(gene_id)
        fragment["unique_gene_evidence"].add(gene_id)
        fragment["gene_variant_evidence"][gene_id].add(variant_id)
        return gene_id, "UNIQUE_EXONIC"

    if assignment_class == "AMBIGUOUS_SHARED_EXON":
        fragment["ambiguous_supported_variants"].add(variant_id)
        matching = [
            gene_id
            for gene_id in gene_ids
            if genes.get(gene_id, {}).get("strand") == transcript_strand
        ]
        if len(matching) == 1:
            gene_id = matching[0]
            fragment["gene_evidence"].add(gene_id)
            fragment["strand_resolved_gene_evidence"].add(gene_id)
            fragment["gene_variant_evidence"][gene_id].add(variant_id)
            return gene_id, "STRAND_RESOLVED"
        if len(matching) > 1:
            return "", "STRAND_NONUNIQUE"
        return "", "STRAND_NO_MATCH"
    return "", "NO_EXONIC_GENE"


def _empty_gene_count():
    return {
        "strict": {"A": set(), "B": set()},
        "direct": {"A": set(), "B": set()},
        "strand_resolved": {"A": set(), "B": set()},
        "variant_qnames": {
            "A": defaultdict(set),
            "B": defaultdict(set),
        },
        "strand_resolved_variant_qnames": {
            "A": defaultdict(set),
            "B": defaultdict(set),
        },
    }


def build_target_interval_index(measurements, components, merge_gap=0):
    """Merge overlapping component spans into independently fetchable loci."""
    by_contig = defaultdict(list)
    for measurement in measurements:
        component_id = measurement["component_id"]
        component = components[component_id]
        by_contig[component["contig"]].append(component)

    loci = []
    for contig, items in sorted(by_contig.items()):
        unique = {
            component["component_id"]: component for component in items
        }
        ordered = sorted(
            unique.values(),
            key=lambda x: (x["start"], x["end"], x["component_id"]),
        )
        current = None
        for component in ordered:
            if (
                current is None
                or component["start"] > current["end"] + max(0, merge_gap)
            ):
                if current is not None:
                    loci.append(current)
                current = {
                    "contig": contig,
                    "start": component["start"],
                    "end": component["end"],
                    "component_ids": [component["component_id"]],
                }
            else:
                current["end"] = max(current["end"], component["end"])
                current["component_ids"].append(component["component_id"])
        if current is not None:
            loci.append(current)
    return loci


def _check_bam_and_index(bam_path, allow_stale_index):
    if not os.path.isfile(bam_path):
        raise FileNotFoundError("BAM not found: {}".format(bam_path))
    candidates = [bam_path + ".bai"]
    if bam_path.endswith(".bam"):
        candidates.append(bam_path[:-4] + ".bai")
    indexes = [path for path in candidates if os.path.isfile(path)]
    if not indexes:
        raise FileNotFoundError("BAM index not found for {}".format(bam_path))
    index_path = indexes[0]
    if (
        not allow_stale_index
        and os.path.getmtime(index_path) < os.path.getmtime(bam_path)
    ):
        raise RuntimeError(
            "BAM index is older than BAM; reindex before analysis: {}"
            .format(bam_path)
        )
    return index_path


def _position_targets(component_ids, components):
    targets = defaultdict(list)
    for component_id in component_ids:
        component = components[component_id]
        for variant in component["parsed_variants"]:
            targets[variant["position"]].append((component_id, variant))
    return targets


def finalize_fragment_assignment(
    qname,
    fragment,
    measurement,
    component,
    genes,
):
    """Finalize one QNAME and return its auditable assignment."""
    if len(fragment["hap_support"]) == 0:
        hap_assignment = ""
        hap_status = "NO_HAPLOTYPE_SUPPORT"
    elif len(fragment["hap_support"]) > 1:
        hap_assignment = ""
        hap_status = "CONFLICTING_HAPLOTYPE"
    else:
        hap_assignment = next(iter(fragment["hap_support"]))
        hap_status = "STRICT_HAPLOTYPE"

    if len(fragment["transcript_strands"]) == 1:
        transcript_strand = next(iter(fragment["transcript_strands"]))
        strand_status = "CONSISTENT"
    elif len(fragment["transcript_strands"]) > 1:
        transcript_strand = ""
        strand_status = "CONFLICTING_TRANSCRIPT_STRAND"
    else:
        transcript_strand = ""
        strand_status = "NO_TRANSCRIPT_STRAND"

    if len(fragment["gene_evidence"]) == 1 and strand_status == "CONSISTENT":
        gene_id = next(iter(fragment["gene_evidence"]))
        gene_status = "UNIQUE_GENE"
    elif len(fragment["gene_evidence"]) > 1:
        gene_id = ""
        gene_status = "CONFLICTING_GENE"
    else:
        gene_id = ""
        gene_status = "NO_GENE_EVIDENCE"

    mismatch_gene_ids = sorted(
        fragment["exon_strand_mismatch_gene_evidence"]
    )
    assignment_method = ""
    if gene_id:
        has_unique = gene_id in fragment["unique_gene_evidence"]
        has_resolved = gene_id in fragment["strand_resolved_gene_evidence"]
        if has_unique:
            # A gene-specific exonic observation is sufficient for assignment.
            # Any compatible shared-exon/strand observation is corroborating
            # evidence only and does not create a separate result class.
            assignment_method = "DIRECT"
        elif has_resolved:
            assignment_method = "STRAND_RESOLVED"
        else:
            assignment_method = "OTHER"
    elif mismatch_gene_ids:
        assignment_method = "ANNOTATED_EXON_STRAND_MISMATCH"

    strict_assigned = bool(
        hap_assignment and gene_id and strand_status == "CONSISTENT"
    )
    if strict_assigned:
        final_status = "ASSIGNED"
    elif hap_status == "CONFLICTING_HAPLOTYPE":
        final_status = "HAPLOTYPE_CONFLICT"
    elif gene_status == "CONFLICTING_GENE":
        final_status = "GENE_CONFLICT"
    elif strand_status == "CONFLICTING_TRANSCRIPT_STRAND":
        final_status = "STRAND_CONFLICT"
    elif mismatch_gene_ids:
        final_status = "ANNOTATED_EXON_STRAND_MISMATCH"
    elif gene_status == "NO_GENE_EVIDENCE":
        final_status = "NO_GENE"
    else:
        final_status = "UNRESOLVED"

    strand_resolution_required = bool(
        final_status == "ASSIGNED"
        and assignment_method == "STRAND_RESOLVED"
    )
    gene_supported = (
        sorted(fragment["gene_variant_evidence"].get(gene_id, set()))
        if gene_id else []
    )
    row = {
        "sample": measurement["sample"],
        "tissue": measurement["tissue"],
        "bam": measurement["bam"],
        "haplotype_block_id": component["component_id"],
        "qname": qname,
        "transcript_strand": transcript_strand,
        "gene_id": gene_id,
        "gene_name": genes.get(gene_id, {}).get("gene_name", gene_id),
        "strand_mismatch_gene_ids": ",".join(mismatch_gene_ids),
        "strand_mismatch_gene_names": ",".join(
            genes.get(item, {}).get("gene_name", item)
            for item in mismatch_gene_ids
        ),
        "haplotype": hap_assignment,
        "gene_assignment_method": assignment_method,
        "strand_resolution_required": int(strand_resolution_required),
        "final_status": final_status,
        "supported_variants": ",".join(
            sorted(fragment["supported_variants"])
        ),
        "gene_supported_variants": ",".join(gene_supported),
        "ambiguous_supported_variants": ",".join(
            sorted(fragment["ambiguous_supported_variants"])
        ),
        "observations": ";".join(fragment["observation_strings"]),
    }
    return row, gene_supported


def aggregate_strict_gene_component_counts(
    measurement,
    component,
    fragments,
    genes,
    fragment_writer=None,
):
    """Partition accepted QNAMEs and calculate component-level BAM QC."""
    phaserlike_a = set()
    phaserlike_b = set()
    gene_counts = defaultdict(_empty_gene_count)
    status_counts = defaultdict(int)
    n_strand_resolution_evidence = 0

    for qname, fragment in fragments.items():
        if "A" in fragment["hap_support"]:
            phaserlike_a.add(qname)
        if "B" in fragment["hap_support"]:
            phaserlike_b.add(qname)
        row, gene_supported = finalize_fragment_assignment(
            qname, fragment, measurement, component, genes
        )
        status_counts[row["final_status"]] += 1
        if fragment["strand_resolved_gene_evidence"]:
            n_strand_resolution_evidence += 1
        if fragment_writer is not None:
            fragment_writer.writerow(row)

        if row["final_status"] != "ASSIGNED":
            continue
        gene_id = row["gene_id"]
        hap = row["haplotype"]
        counts = gene_counts[gene_id]
        counts["strict"][hap].add(qname)
        if row["strand_resolution_required"]:
            counts["strand_resolved"][hap].add(qname)
        else:
            counts["direct"][hap].add(qname)
        for variant_id in gene_supported:
            counts["variant_qnames"][hap][variant_id].add(qname)
            if row["strand_resolution_required"]:
                counts["strand_resolved_variant_qnames"][hap][
                    variant_id
                ].add(qname)

    n_supported = len(fragments)
    n_hap_conflict = status_counts["HAPLOTYPE_CONFLICT"]
    conflict_fraction = (
        n_hap_conflict / float(n_supported) if n_supported else 0.0
    )
    reported_a = int(measurement.get("reported_a", 0))
    reported_b = int(measurement.get("reported_b", 0))
    comparison = {
        "sample": measurement["sample"],
        "tissue": measurement["tissue"],
        "bam": measurement["bam"],
        "component_id": component["component_id"],
        "contig": component["contig"],
        "start": component["start"],
        "stop": component["end"],
        "n_variants": len(component["variants"]),
        "reported_a": reported_a,
        "reported_b": reported_b,
        "reported_total": reported_a + reported_b,
        "reconstructed_phaserlike_a": len(phaserlike_a),
        "reconstructed_phaserlike_b": len(phaserlike_b),
        "reconstructed_phaserlike_total": (
            len(phaserlike_a) + len(phaserlike_b)
        ),
        "delta_a": len(phaserlike_a) - reported_a,
        "delta_b": len(phaserlike_b) - reported_b,
        "count_match_exact": (
            len(phaserlike_a) == reported_a
            and len(phaserlike_b) == reported_b
        ),
        "n_qnames_with_variant_support": n_supported,
        "n_qnames_assigned": status_counts["ASSIGNED"],
        "n_qnames_with_strand_resolution_evidence": (
            n_strand_resolution_evidence
        ),
        "n_qnames_strand_resolved": sum(
            len(rec["strand_resolved"]["A"])
            + len(rec["strand_resolved"]["B"])
            for rec in gene_counts.values()
        ),
        "n_haplotype_conflict": n_hap_conflict,
        "fragment_hap_conflict_fraction": conflict_fraction,
        "n_gene_conflict": status_counts["GENE_CONFLICT"],
        "n_strand_conflict": status_counts["STRAND_CONFLICT"],
        "n_annotated_exon_strand_mismatch": status_counts[
            "ANNOTATED_EXON_STRAND_MISMATCH"
        ],
        "n_no_gene": status_counts["NO_GENE"],
        "n_unresolved": status_counts["UNRESOLVED"],
    }
    return gene_counts, comparison


def stream_bam_fragment_evidence(
    bam_path,
    measurements,
    components,
    assignments,
    genes,
    options,
    min_as=None,
    fragment_writer=None,
    eligible_variants=None,
):
    """Yield all-SNP audit and WGS-eligible analysis counts per measurement."""
    if pysam is None:
        raise ImportError(
            "pysam is required for v0.2.5 BAM fragment reconstruction"
        )
    _check_bam_and_index(bam_path, options.allow_stale_index)
    by_component = {
        measurement["component_id"]: measurement
        for measurement in measurements
    }
    loci = build_target_interval_index(
        measurements,
        components,
        merge_gap=int(getattr(options, "merge_locus_gap", 0)),
    )
    yielded = set()

    with pysam.AlignmentFile(
        bam_path,
        "rb",
        threads=max(0, int(getattr(options, "bam_threads", 0))),
    ) as bam:
        for locus in loci:
            targets = _position_targets(locus["component_ids"], components)
            all_fragments = {
                component_id: defaultdict(new_fragment_record)
                for component_id in locus["component_ids"]
            }
            analysis_fragments = {
                component_id: defaultdict(new_fragment_record)
                for component_id in locus["component_ids"]
            }
            filter_counts = defaultdict(int)
            n_alignments_examined = 0
            n_allele_observations = 0
            bam_contig = resolve_bam_contig(bam, locus["contig"])
            if bam_contig is None:
                raise ValueError(
                    "Contig {} absent from {}".format(
                        locus["contig"], bam_path
                    )
                )
            for read in bam.fetch(
                bam_contig,
                max(0, locus["start"] - 1),
                locus["end"],
            ):
                n_alignments_examined += 1
                passed, reason = read_passes_filters(read, options, min_as)
                if not passed:
                    filter_counts[reason] += 1
                    continue
                sequence = read.query_sequence
                if sequence is None:
                    filter_counts["NO_SEQUENCE"] += 1
                    continue
                qualities = read.query_qualities
                transcript_strand = inferred_transcript_strand(
                    read, options.library_strand
                )
                if transcript_strand is None:
                    filter_counts["NO_TRANSCRIPT_STRAND"] += 1
                    continue

                for query_pos, reference_pos in read.get_aligned_pairs(
                    matches_only=True
                ):
                    if query_pos is None or reference_pos is None:
                        continue
                    position = reference_pos + 1
                    local_targets = targets.get(position)
                    if not local_targets:
                        continue
                    if (
                        qualities is not None
                        and qualities[query_pos] < options.min_baseq
                    ):
                        filter_counts["LOW_BASEQ_OBSERVATION"] += 1
                        continue
                    base = sequence[query_pos].upper()
                    for component_id, variant in local_targets:
                        if (
                            len(variant["hapA"]) != 1
                            or len(variant["hapB"]) != 1
                        ):
                            filter_counts["NON_SNP_VARIANT"] += 1
                            continue
                        if base == variant["hapA"] == variant["hapB"]:
                            hap = "AB"
                        elif base == variant["hapA"]:
                            hap = "A"
                        elif base == variant["hapB"]:
                            hap = "B"
                        else:
                            filter_counts["OTHER_ALLELE_OBSERVATION"] += 1
                            continue

                        n_allele_observations += 1
                        assignment = assignments.get(
                            variant["variant_id"],
                            {
                                "assignment_class": "UNAVAILABLE",
                                "gene_ids": [],
                            },
                        )
                        targets_to_update = [
                            all_fragments[component_id][read.query_name]
                        ]
                        if (
                            eligible_variants is None
                            or variant["variant_id"] in eligible_variants
                        ):
                            targets_to_update.append(
                                analysis_fragments[component_id][
                                    read.query_name
                                ]
                            )
                        for fragment in targets_to_update:
                            fragment["supported_variants"].add(
                                variant["variant_id"]
                            )
                            fragment["transcript_strands"].add(
                                transcript_strand
                            )
                            if hap in ("A", "B"):
                                fragment["hap_support"].add(hap)
                            gene_id, source = add_gene_evidence(
                                fragment,
                                variant["variant_id"],
                                assignment,
                                transcript_strand,
                                genes,
                            )
                            fragment["observation_strings"].append(
                                "{}:{}:{}:{}:{}".format(
                                    variant["variant_id"],
                                    base,
                                    hap,
                                    transcript_strand,
                                    gene_id or source,
                                )
                            )

            for component_id in locus["component_ids"]:
                measurement = by_component[component_id]
                component = components[component_id]
                _, comparison = (
                    aggregate_strict_gene_component_counts(
                        measurement,
                        component,
                        all_fragments[component_id],
                        genes,
                        fragment_writer=None,
                    )
                )
                gene_counts, analysis_comparison = (
                    aggregate_strict_gene_component_counts(
                        measurement,
                        component,
                        analysis_fragments[component_id],
                        genes,
                        fragment_writer=fragment_writer,
                    )
                )
                comparison.update(
                    {
                        "bam_path": bam_path,
                        "n_alignments_examined_in_locus": (
                            n_alignments_examined
                        ),
                        "n_allele_observations_in_locus": (
                            n_allele_observations
                        ),
                        "min_as_applied": min_as,
                        "filter_counts": dict(filter_counts),
                        "analysis_qnames_with_variant_support": (
                            analysis_comparison[
                                "n_qnames_with_variant_support"
                            ]
                        ),
                        "analysis_qnames_assigned": analysis_comparison[
                            "n_qnames_assigned"
                        ],
                        "analysis_haplotype_conflict": analysis_comparison[
                            "n_haplotype_conflict"
                        ],
                        "analysis_fragment_hap_conflict_fraction": (
                            analysis_comparison[
                                "fragment_hap_conflict_fraction"
                            ]
                        ),
                    }
                )
                yielded.add(component_id)
                yield {
                    "measurement": measurement,
                    "component": component,
                    "gene_counts": gene_counts,
                    "comparison": comparison,
                }

    missing = set(by_component).difference(yielded)
    if missing:
        raise RuntimeError(
            "BAM scan did not yield components: {}".format(
                ",".join(sorted(missing))
            )
        )
