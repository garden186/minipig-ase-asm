"""Reassign saved fragment evidence using an explicit GQ cutoff; never read BAMs."""

from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import sys
from common import (
    CALLS,
    STATES,
    PHASE_FIELDS,
    PACKAGE,
    rows,
    sha,
    read_json,
    write_json,
    write_table,
    stamp,
    verify_package,
)
from profile import SITE_FIELDS
from evidence import reconstruct, numeric_gq

ASM_FIELDS = [
    "sample",
    "tissue",
    "chrom",
    "cpg_pos1",
    "ps",
    "genotype_disrupted",
    "H1_methylated",
    "H1_unmethylated",
    "H2_methylated",
    "H2_unmethylated",
    "H1_depth",
    "H2_depth",
    "H1_methylation_fraction",
    "H2_methylation_fraction",
    "delta_m",
    "assignment_GQ_cutoff",
    "observation_status",
]


def write_asm_inputs(directory, min_gq=None):
    """Prepare same-CpG fragment contingency tables, without significance tests."""
    masks = {
        (r["sample"], r["tissue"], r["chrom"], r["cpg_pos1"]): r["genotype_disrupted"]
        for r in rows(directory / "cpg_assignment.tsv.gz")
    }

    def values():
        for r in rows(directory / "cpg_haplotype_by_ps.tsv.gz"):
            key = tuple(r[k] for k in ("sample", "tissue", "chrom", "cpg_pos1"))
            h = [int(r[f"FRAGMENT_{hap}_{call}"]) for hap in ("H1", "H2") for call in CALLS[:2]]
            depths = [h[0] + h[1], h[2] + h[3]]
            if not sum(depths):
                continue
            fractions = [
                h[0] / depths[0] if depths[0] else None,
                h[2] / depths[1] if depths[1] else None,
            ]
            delta = fractions[0] - fractions[1] if all(depths) else None
            yield list(key) + [r["ps"], masks[key]] + h + depths + [
                v if v is not None else "." for v in fractions
            ] + [
                delta if delta is not None else ".",
                min_gq if min_gq is not None else "NONE",
                "BOTH_HAPLOTYPES_OBSERVED" if all(depths) else "ONE_HAPLOTYPE_UNOBSERVED",
            ]

    write_table(directory / "cpg_asm_inputs.tsv.gz", ASM_FIELDS, values())


def source_gate(directory):
    gate = read_json(directory / "validation.json")
    if gate["status"] != "PASS_DESCRIPTIVE_CPG_ASSIGNMENT" or gate["errors"]:
        raise ValueError("A completed evidence region is required")
    for n, h in gate["sha256"].items():
        if sha(directory / n) != h:
            raise ValueError("Source region checksum differs: " + n)
    evidence = read_json(directory / "evidence_validation.json")
    if evidence["status"] != "PASS_UNFILTERED_EVIDENCE_REPLAY" or evidence["errors"]:
        raise ValueError("No unfiltered evidence replay gate")
    return gate


def region(source, output, min_gq):
    source_gate(source)
    output.mkdir(parents=True, exist_ok=False)
    task = read_json(source / "task.json")
    sites, phases, transitions, nfragments, nlinks, rescued = reconstruct(
        source, min_gq, verify_baseline=min_gq is None
    )
    metadata, template_totals = {}, {}
    for r in rows(source / "cpg_assignment.tsv.gz"):
        p = int(r["cpg_pos1"]) - 1
        metadata[p] = [r[k] for k in SITE_FIELDS[:5]]
        for u in ("READ", "FRAGMENT"):
            for c in CALLS:
                before = sum(int(r[f"{u}_{s}_{c}"]) for s in STATES)
                after = sum(sites[p, u, s, c] for s in STATES)
                if before != after:
                    raise ValueError("GQ filtering changed the total CpG methylation observations")

    def site_rows():
        for p, meta in sorted(metadata.items()):
            yield meta + [
                sites[p, u, s, c] for u in ("READ", "FRAGMENT") for s in STATES for c in CALLS
            ] + [rescued[p]]

    write_table(output / "cpg_assignment.tsv.gz", SITE_FIELDS, site_rows())
    keys = sorted({(p, ps) for p, ps, u, h, c in phases})
    write_table(
        output / "cpg_haplotype_by_ps.tsv.gz",
        PHASE_FIELDS,
        (
            metadata[p][:4]
            + [ps]
            + [
                phases[p, ps, u, h, c]
                for u in ("READ", "FRAGMENT")
                for h in ("H1", "H2")
                for c in CALLS
            ]
            for p, ps in keys
        ),
    )
    write_asm_inputs(output, min_gq)
    from check_asm_inputs import check

    check(output, min_gq)
    write_table(
        output / "fragment_assignment_transitions.tsv",
        ["baseline_state", "filtered_state", "fragments"],
        (list(k) + [v] for k, v in sorted(transitions.items())),
    )
    catalog = read_json(source / "link_catalog.json")["links"]
    distribution = Counter(
        (v["gq_source"], "MISSING" if v["gq"] is None else str(v["gq"])) for v in catalog.values()
    )
    write_table(
        output / "link_GQ_distribution.tsv",
        ["gq_source", "GQ", "eligible_SNPs"],
        (list(k) + [v] for k, v in sorted(distribution.items())),
    )
    gate = dict(
        status="PASS_GQ_REASSIGNMENT_COUNTS",
        errors=[],
        sample=task["sample"],
        tissue=task["tissue"],
        chrom=task["chrom"],
        region_start0=task["start0"],
        region_end0=task["end0"],
        source=str(source.resolve()),
        source_gate_sha256=sha(source / "validation.json"),
        min_gq=min_gq,
        missing_gq="RETAINED_WHEN_NO_CUTOFF_OTHERWISE_EXCLUDED",
        retained_links=nlinks,
        baseline_links=len(catalog),
        fragment_rows=nfragments,
        total_CpG_call_counts_preserved=True,
        cpg_mask="Unchanged genotype-disruption flag; GQ filtering does not certify uncertain CpG genotypes.",
        limits="A positive GQ cutoff affects link SNPs only. No rephasing, new mapping/base-quality filter, ASM test or bias correction.",
        sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()},
    )
    write_json(output / "validation.json", gate)
    return gate


def parse_gq(value):
    if value.upper() == "NONE":
        return None
    number = numeric_gq(value)
    if number is None:
        raise argparse.ArgumentTypeError("Use NONE or a finite nonnegative GQ cutoff")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Completed fragment-count tile")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-gq", type=parse_gq, required=True)
    args = parser.parse_args()
    region(args.source.resolve(), args.output.resolve(), args.min_gq)


if __name__ == "__main__":
    main()
