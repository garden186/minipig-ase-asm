#!/usr/bin/env python3
"""Validate recurrent ASE SNP locus-collapse outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterator, Sequence


class ValidationError(RuntimeError):
    pass


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open("r", encoding="utf-8", newline="")
    with handle:
        yield from csv.DictReader(handle, delimiter="\t")


def validate(out_dir: Path) -> Dict[str, object]:
    required = [
        "recurrent_ase_loci.tsv",
        "recurrent_snp_locus_membership.tsv",
        "recurrent_snp_parent_block_links.tsv.gz",
        "recurrent_snp_block_edges.tsv.gz",
        "locus_counts_by_tissue.tsv",
        "analysis_qc.json",
        "run_metadata.json",
    ]
    for name in required:
        path = out_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValidationError(f"missing or empty output: {path}")
    with (out_dir / "run_metadata.json").open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    with (out_dir / "analysis_qc.json").open(encoding="utf-8") as handle:
        qc = json.load(handle)
    if metadata.get("status") != "PASS" or qc.get("status") != "PASS":
        raise ValidationError("metadata or QC status is not PASS")
    if metadata.get("distance_rule") != "NONE":
        raise ValidationError("distance rule must be NONE")

    loci: Dict[str, Dict[str, str]] = {}
    for row in iter_tsv(out_dir / "recurrent_ase_loci.tsv"):
        locus_id = row["locus_id"]
        if locus_id in loci:
            raise ValidationError(f"duplicate locus ID: {locus_id}")
        if row["locus_is_descriptive_collapse_of_significant_snps"] != "1":
            raise ValidationError(f"locus is not marked descriptive: {locus_id}")
        if int(row["n_recurrent_snps"]) != len(row["recurrent_snp_ids"].split(",")):
            raise ValidationError(f"locus member count mismatch: {locus_id}")
        loci[locus_id] = row

    membership_keys = set()
    membership_locus = {}
    members_by_locus = defaultdict(list)
    leads = Counter()
    expected_success_keys = set()
    for row in iter_tsv(out_dir / "recurrent_snp_locus_membership.tsv"):
        key = (row["gene_id"], row["tissue"], row["variant_id"])
        if key in membership_keys:
            raise ValidationError(f"duplicate membership hypothesis: {key}")
        membership_keys.add(key)
        locus_id = row["locus_id"]
        if locus_id not in loci:
            raise ValidationError(f"membership has unknown locus: {locus_id}")
        members_by_locus[locus_id].append(row["variant_id"])
        membership_locus[key] = locus_id
        leads[locus_id] += int(row["is_locus_lead_snp"])
        samples = [value for value in row["animals_with_block_concordant_snp_ase"].split(",") if value]
        if len(samples) != int(row["n_animals_with_block_concordant_snp_ase"]):
            raise ValidationError(f"membership success list mismatch: {key}")
        expected_success_keys.update((sample, row["tissue"], row["gene_id"], row["variant_id"]) for sample in samples)

    for locus_id, locus in loci.items():
        members = sorted(members_by_locus[locus_id])
        if len(members) != int(locus["n_recurrent_snps"]):
            raise ValidationError(f"membership count mismatch for {locus_id}")
        if sorted(locus["recurrent_snp_ids"].split(",")) != members:
            raise ValidationError(f"membership variants mismatch for {locus_id}")
        if leads[locus_id] != 1:
            raise ValidationError(f"locus does not have exactly one lead SNP: {locus_id}")

    observed_success_keys = set()
    variants_by_sample_block = defaultdict(set)
    link_count = 0
    for row in iter_tsv(out_dir / "recurrent_snp_parent_block_links.tsv.gz"):
        link_count += 1
        key = (row["sample"], row["tissue"], row["gene_id"], row["variant_id"])
        if (row["gene_id"], row["tissue"], row["variant_id"]) not in membership_keys:
            raise ValidationError(f"parent link has unknown hypothesis: {key}")
        hypothesis = (row["gene_id"], row["tissue"], row["variant_id"])
        if membership_locus[hypothesis] != row["locus_id"]:
            raise ValidationError(f"parent link locus mismatch: {key}")
        if row["parent_major_allele_at_snp"] != row["site_major_allele"]:
            raise ValidationError(f"direction-discordant parent link: {key}")
        if int(row["ref_fragment_count"]) + int(row["alt_fragment_count"]) != int(
            row["total_informative_fragments"]
        ):
            raise ValidationError(f"invalid fragment total in link: {key}")
        observed_success_keys.add(key)
        variants_by_sample_block[
            (row["gene_id"], row["tissue"], row["sample"], row["measurement_id"])
        ].add(row["variant_id"])
    if expected_success_keys != observed_success_keys:
        raise ValidationError(
            f"success-parent mapping mismatch: missing={len(expected_success_keys-observed_success_keys)}, "
            f"extra={len(observed_success_keys-expected_success_keys)}"
        )

    expected_edge_keys = set()
    for (gene_id, tissue, _, _), variants in variants_by_sample_block.items():
        for first, second in combinations(sorted(variants), 2):
            expected_edge_keys.add((gene_id, tissue, first, second))

    edge_count = 0
    observed_edge_keys = set()
    for row in iter_tsv(out_dir / "recurrent_snp_block_edges.tsv.gz"):
        edge_count += 1
        first = (row["gene_id"], row["tissue"], row["variant_id_1"])
        second = (row["gene_id"], row["tissue"], row["variant_id_2"])
        if first not in membership_keys or second not in membership_keys:
            raise ValidationError(f"edge has unknown hypothesis: {first}, {second}")
        edge_key = (row["gene_id"], row["tissue"], row["variant_id_1"], row["variant_id_2"])
        if edge_key in observed_edge_keys:
            raise ValidationError(f"duplicate graph edge: {edge_key}")
        observed_edge_keys.add(edge_key)
        if membership_locus[first] != row["locus_id"] or membership_locus[second] != row["locus_id"]:
            raise ValidationError(f"edge locus mismatch: {edge_key}")
        if int(row["n_shared_success_sample_blocks"]) < 1:
            raise ValidationError(f"edge has no shared block evidence: {row}")
    if observed_edge_keys != expected_edge_keys:
        raise ValidationError(
            f"edge set cannot be reproduced from parent links: "
            f"missing={len(expected_edge_keys-observed_edge_keys)}, "
            f"extra={len(observed_edge_keys-expected_edge_keys)}"
        )

    graph = {key: key for key in membership_keys}

    def find(value):
        while graph[value] != value:
            graph[value] = graph[graph[value]]
            value = graph[value]
        return value

    def union(first, second):
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            graph[root_second] = root_first

    for gene_id, tissue, first, second in observed_edge_keys:
        union((gene_id, tissue, first), (gene_id, tissue, second))
    component_loci = defaultdict(set)
    locus_components = defaultdict(set)
    for key in membership_keys:
        component_loci[find(key)].add(membership_locus[key])
        locus_components[membership_locus[key]].add(find(key))
    if any(len(values) != 1 for values in component_loci.values()):
        raise ValidationError("a graph component was split across output loci")
    if any(len(values) != 1 for values in locus_components.values()):
        raise ValidationError("an output locus combines disconnected graph components")

    observed = {
        "n_recurrent_gene_tissue_snp_hypotheses": len(membership_keys),
        "n_recurrent_ase_loci": len(loci),
        "n_locus_lead_snps": sum(leads.values()),
        "n_recurrent_success_keys": len(observed_success_keys),
        "n_recurrent_success_parent_block_links": link_count,
        "n_direct_core_block_edges": edge_count,
        "n_single_snp_loci": sum(int(row["n_recurrent_snps"] == "1") for row in loci.values()),
        "n_multi_snp_loci": sum(int(row["n_recurrent_snps"] != "1") for row in loci.values()),
    }
    expected = metadata["counts"]
    for key, value in observed.items():
        if int(expected[key]) != value:
            raise ValidationError(
                f"metadata count mismatch for {key}: expected={expected[key]}, observed={value}"
            )
    result = {
        "status": "PASS",
        "checks": {
            "all_required_outputs_present": True,
            "one_membership_row_per_recurrent_hypothesis": True,
            "exactly_one_lead_snp_per_locus": True,
            "all_successes_have_direction_concordant_v026_core_parent_links": True,
            "all_edges_have_shared_core_block_evidence": True,
            "edge_set_reproduced_from_parent_links": True,
            "loci_reproduced_as_graph_connected_components": True,
            "no_distance_grouping_rule": True,
            "metadata_counts_reproduced": True,
        },
        "observed_counts": observed,
    }
    with (out_dir / "validation.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = validate(args.out_dir)
    except (ValidationError, KeyError, OSError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] validation: {result['observed_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
