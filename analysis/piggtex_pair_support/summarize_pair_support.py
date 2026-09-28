#!/usr/bin/env python3
"""Intersect 39 haplotype ASE tissue-enriched pairs with independent SNP evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


VERSION = "ase_tissue_enriched_gatk_snp_support_v0.1.0"
PAIR_REQUIRED = {
    "gene_id",
    "gene_name",
    "target_tissue",
    "tissue_enrichment_q_global_BH",
    "recurrence_q_bh_leave_one_gene_out",
    "pair_level_primary",
}
GENE_REQUIRED = {
    "tissue",
    "physical_snp_id",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "site_gene_relation",
    "best_transcript_feature",
    "n_animals_testable",
    "n_animals_with_snp_level_ase",
    "recurrence_fraction",
    "recurrence_q_bh_global",
}
DIRECTION_REQUIRED = {
    "internal_tissue",
    "physical_snp_id",
    "protein_coding_assignment_class",
    "primary_coordinate_gene_id",
    "primary_coordinate_gene_name",
    "piggtex_significant_rows_same_mapped_tissue",
    "piggtex_min_fdr_same_mapped_tissue",
    "piggtex_max_pcadd",
    "major_allele_direction_comparison",
}

PAIR_FIELDS = [
    "gene_id",
    "gene_name",
    "target_tissue",
    "tissue_enrichment_q_global_BH",
    "recurrence_q_bh_leave_one_gene_out",
    "n_coordinate_assigned_independent_recurrent_snps",
    "n_exonic_independent_recurrent_snps",
    "n_unique_protein_coding_exonic_independent_recurrent_snps",
    "n_unique_exonic_snps_with_piggtex_same_tissue_exact_support",
    "n_unique_exonic_snps_with_direction_comparable_piggtex_support",
    "n_unique_exonic_snps_with_direction_concordant_piggtex_support",
    "n_unique_exonic_snps_with_direction_discordant_piggtex_support",
    "max_independent_snp_ase_positive_animals",
    "min_independent_snp_recurrence_q_bh_global",
    "independent_snp_ids",
    "unique_exonic_independent_snp_ids",
    "piggtex_exact_supported_unique_exonic_snp_ids",
    "direction_concordant_unique_exonic_snp_ids",
    "best_available_snp_support_class",
]

SITE_FIELDS = [
    "target_gene_id",
    "target_gene_name",
    "target_tissue",
    "physical_snp_id",
    "site_gene_relation",
    "best_transcript_feature",
    "n_animals_testable",
    "n_animals_with_snp_level_ase",
    "recurrence_fraction",
    "recurrence_q_bh_global",
    "protein_coding_assignment_class",
    "is_unique_protein_coding_exonic_assignment_to_target",
    "piggtex_significant_rows_same_mapped_tissue",
    "piggtex_min_fdr_same_mapped_tissue",
    "piggtex_max_pcadd",
    "major_allele_direction_comparison",
]


class AnalysisError(RuntimeError):
    pass


def normalize_gene_id(value: str) -> str:
    return value.strip().split(".")[0]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def load_tsv(path: Path, required: set[str]) -> List[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AnalysisError(f"{path}: missing columns: {', '.join(missing)}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            rows.append({key: value if value is not None else "" for key, value in row.items()})
    return rows


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path: Path):
    return {"size_bytes": path.stat().st_size, "sha256": sha256(path)}


def best_support_class(
    has_any: bool,
    has_unique_exonic: bool,
    has_exact_support: bool,
    has_concordant_direction: bool,
) -> str:
    if has_concordant_direction:
        return "DIRECTION_CONCORDANT_PIGGTEX_EXACT_SITE"
    if has_exact_support:
        return "PIGGTEX_EXACT_SITE_WITHOUT_CONCORDANT_DIRECTION"
    if has_unique_exonic:
        return "INDEPENDENT_RECURRENT_UNIQUE_EXONIC_SNP_ONLY"
    if has_any:
        return "OTHER_COORDINATE_ASSIGNED_INDEPENDENT_RECURRENT_SNP"
    return "NO_COORDINATE_ASSIGNED_INDEPENDENT_RECURRENT_SNP"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissue-enriched-pairs", type=Path, required=True)
    parser.add_argument("--gene-assignments", type=Path, required=True)
    parser.add_argument("--direction-comparison", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    for path in (args.tissue_enriched_pairs, args.gene_assignments, args.direction_comparison):
        if not path.is_file():
            raise AnalysisError(f"missing input: {path}")
    if args.output_dir.exists():
        if not args.force:
            raise AnalysisError(f"output directory already exists: {args.output_dir}")
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True)

    pair_rows = load_tsv(args.tissue_enriched_pairs, PAIR_REQUIRED)
    pairs = {}
    for row in pair_rows:
        if row["pair_level_primary"] != "1":
            raise AnalysisError("non-primary pair present in the 39-pair input")
        key = (normalize_gene_id(row["gene_id"]), row["target_tissue"])
        if key in pairs:
            raise AnalysisError(f"duplicate pair: {key}")
        pairs[key] = row

    direction_rows = load_tsv(args.direction_comparison, DIRECTION_REQUIRED)
    direction = {}
    for row in direction_rows:
        key = (row["internal_tissue"], row["physical_snp_id"])
        if key in direction:
            raise AnalysisError(f"duplicate direction key: {key}")
        direction[key] = row

    assignments: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    seen_contexts = set()
    for row in load_tsv(args.gene_assignments, GENE_REQUIRED):
        key = (normalize_gene_id(row["gene_id"]), row["tissue"])
        if key not in pairs:
            continue
        context = key + (row["physical_snp_id"],)
        if context in seen_contexts:
            raise AnalysisError(f"duplicate gene assignment context: {context}")
        seen_contexts.add(context)
        direction_key = (row["tissue"], row["physical_snp_id"])
        if direction_key not in direction:
            raise AnalysisError(f"direction record missing: {direction_key}")
        assignments[key].append(row)

    pair_output = []
    site_output = []
    for key, pair in pairs.items():
        gene_id, tissue = key
        rows = sorted(
            assignments.get(key, []),
            key=lambda row: (
                int(row["physical_snp_id"].split("_")[0]),
                int(row["physical_snp_id"].split("_")[1]),
                row["physical_snp_id"],
            ),
        )
        all_sites = []
        exonic_sites = []
        unique_exonic_sites = []
        exact_sites = []
        comparable_sites = []
        concordant_sites = []
        discordant_sites = []
        for row in rows:
            external = direction[(row["tissue"], row["physical_snp_id"])]
            site = row["physical_snp_id"]
            all_sites.append(site)
            if row["site_gene_relation"] == "EXONIC":
                exonic_sites.append(site)
            is_unique_exonic = (
                external["protein_coding_assignment_class"] == "UNIQUE_EXONIC_GENE"
                and normalize_gene_id(external["primary_coordinate_gene_id"]) == gene_id
            )
            if is_unique_exonic:
                unique_exonic_sites.append(site)
                if int(external["piggtex_significant_rows_same_mapped_tissue"]) > 0:
                    exact_sites.append(site)
                if external["major_allele_direction_comparison"] in {"CONCORDANT", "DISCORDANT"}:
                    comparable_sites.append(site)
                if external["major_allele_direction_comparison"] == "CONCORDANT":
                    concordant_sites.append(site)
                if external["major_allele_direction_comparison"] == "DISCORDANT":
                    discordant_sites.append(site)
            site_output.append(
                {
                    "target_gene_id": gene_id,
                    "target_gene_name": pair["gene_name"],
                    "target_tissue": tissue,
                    "physical_snp_id": site,
                    "site_gene_relation": row["site_gene_relation"],
                    "best_transcript_feature": row["best_transcript_feature"],
                    "n_animals_testable": row["n_animals_testable"],
                    "n_animals_with_snp_level_ase": row["n_animals_with_snp_level_ase"],
                    "recurrence_fraction": row["recurrence_fraction"],
                    "recurrence_q_bh_global": row["recurrence_q_bh_global"],
                    "protein_coding_assignment_class": external[
                        "protein_coding_assignment_class"
                    ],
                    "is_unique_protein_coding_exonic_assignment_to_target": int(
                        is_unique_exonic
                    ),
                    "piggtex_significant_rows_same_mapped_tissue": external[
                        "piggtex_significant_rows_same_mapped_tissue"
                    ],
                    "piggtex_min_fdr_same_mapped_tissue": external[
                        "piggtex_min_fdr_same_mapped_tissue"
                    ],
                    "piggtex_max_pcadd": external["piggtex_max_pcadd"],
                    "major_allele_direction_comparison": external[
                        "major_allele_direction_comparison"
                    ],
                }
            )
        support_class = best_support_class(
            bool(all_sites),
            bool(unique_exonic_sites),
            bool(exact_sites),
            bool(concordant_sites),
        )
        pair_output.append(
            {
                "gene_id": gene_id,
                "gene_name": pair["gene_name"],
                "target_tissue": tissue,
                "tissue_enrichment_q_global_BH": pair["tissue_enrichment_q_global_BH"],
                "recurrence_q_bh_leave_one_gene_out": pair[
                    "recurrence_q_bh_leave_one_gene_out"
                ],
                "n_coordinate_assigned_independent_recurrent_snps": len(all_sites),
                "n_exonic_independent_recurrent_snps": len(exonic_sites),
                "n_unique_protein_coding_exonic_independent_recurrent_snps": len(
                    unique_exonic_sites
                ),
                "n_unique_exonic_snps_with_piggtex_same_tissue_exact_support": len(
                    exact_sites
                ),
                "n_unique_exonic_snps_with_direction_comparable_piggtex_support": len(
                    comparable_sites
                ),
                "n_unique_exonic_snps_with_direction_concordant_piggtex_support": len(
                    concordant_sites
                ),
                "n_unique_exonic_snps_with_direction_discordant_piggtex_support": len(
                    discordant_sites
                ),
                "max_independent_snp_ase_positive_animals": max(
                    (int(row["n_animals_with_snp_level_ase"]) for row in rows), default=0
                ),
                "min_independent_snp_recurrence_q_bh_global": f"{min((float(row['recurrence_q_bh_global']) for row in rows), default=1.0):.12g}",
                "independent_snp_ids": ",".join(all_sites),
                "unique_exonic_independent_snp_ids": ",".join(unique_exonic_sites),
                "piggtex_exact_supported_unique_exonic_snp_ids": ",".join(exact_sites),
                "direction_concordant_unique_exonic_snp_ids": ",".join(concordant_sites),
                "best_available_snp_support_class": support_class,
            }
        )
    pair_output.sort(
        key=lambda row: (
            -int(row["n_unique_exonic_snps_with_direction_concordant_piggtex_support"]),
            -int(row["n_unique_exonic_snps_with_piggtex_same_tissue_exact_support"]),
            -int(row["n_unique_protein_coding_exonic_independent_recurrent_snps"]),
            float(row["tissue_enrichment_q_global_BH"]),
            row["gene_name"],
            row["target_tissue"],
        )
    )
    site_output.sort(
        key=lambda row: (
            row["target_gene_name"],
            row["target_tissue"],
            int(row["physical_snp_id"].split("_")[0]),
            int(row["physical_snp_id"].split("_")[1]),
        )
    )

    stage_definitions = [
        ("PRIMARY_TISSUE_ENRICHED_PAIRS", lambda row: True),
        (
            "WITH_ANY_COORDINATE_ASSIGNED_INDEPENDENT_RECURRENT_SNP",
            lambda row: int(row["n_coordinate_assigned_independent_recurrent_snps"]) > 0,
        ),
        (
            "WITH_ANY_EXONIC_INDEPENDENT_RECURRENT_SNP",
            lambda row: int(row["n_exonic_independent_recurrent_snps"]) > 0,
        ),
        (
            "WITH_UNIQUE_PROTEIN_CODING_EXONIC_INDEPENDENT_RECURRENT_SNP",
            lambda row: int(
                row["n_unique_protein_coding_exonic_independent_recurrent_snps"]
            )
            > 0,
        ),
        (
            "WITH_PIGGTEX_SAME_TISSUE_EXACT_SITE_SUPPORT",
            lambda row: int(
                row["n_unique_exonic_snps_with_piggtex_same_tissue_exact_support"]
            )
            > 0,
        ),
        (
            "WITH_DIRECTION_COMPARABLE_PIGGTEX_EXACT_SITE",
            lambda row: int(
                row["n_unique_exonic_snps_with_direction_comparable_piggtex_support"]
            )
            > 0,
        ),
        (
            "WITH_DIRECTION_CONCORDANT_PIGGTEX_EXACT_SITE",
            lambda row: int(
                row["n_unique_exonic_snps_with_direction_concordant_piggtex_support"]
            )
            > 0,
        ),
    ]
    stage_rows = [
        {
            "stage_order": index,
            "stage": name,
            "n_of_39_pairs": sum(predicate(row) for row in pair_output),
            "fraction_of_39_pairs": f"{sum(predicate(row) for row in pair_output) / len(pair_output):.12g}",
        }
        for index, (name, predicate) in enumerate(stage_definitions, start=1)
    ]
    class_counts = Counter(row["best_available_snp_support_class"] for row in pair_output)
    class_rows = [
        {
            "best_available_snp_support_class": category,
            "n_of_39_pairs": count,
            "fraction_of_39_pairs": f"{count / len(pair_output):.12g}",
        }
        for category, count in sorted(class_counts.items())
    ]

    pair_path = args.output_dir / "tissue_enriched_pair_independent_snp_support.tsv"
    site_path = args.output_dir / "supporting_independent_recurrent_snps.tsv"
    stage_path = args.output_dir / "support_stage_summary.tsv"
    class_path = args.output_dir / "best_support_class_summary.tsv"
    write_tsv(pair_path, PAIR_FIELDS, pair_output)
    write_tsv(site_path, SITE_FIELDS, site_output)
    write_tsv(stage_path, list(stage_rows[0]), stage_rows)
    write_tsv(class_path, list(class_rows[0]), class_rows)

    counts = {
        "primary_tissue_enriched_pairs": len(pair_output),
        "pairs_with_any_coordinate_assigned_independent_recurrent_snp": sum(
            int(row["n_coordinate_assigned_independent_recurrent_snps"]) > 0
            for row in pair_output
        ),
        "pairs_with_unique_protein_coding_exonic_independent_recurrent_snp": sum(
            int(row["n_unique_protein_coding_exonic_independent_recurrent_snps"]) > 0
            for row in pair_output
        ),
        "pairs_with_piggtex_same_tissue_exact_site_support": sum(
            int(row["n_unique_exonic_snps_with_piggtex_same_tissue_exact_support"]) > 0
            for row in pair_output
        ),
        "pairs_with_direction_comparable_piggtex_exact_site": sum(
            int(row["n_unique_exonic_snps_with_direction_comparable_piggtex_support"]) > 0
            for row in pair_output
        ),
        "pairs_with_direction_concordant_piggtex_exact_site": sum(
            int(row["n_unique_exonic_snps_with_direction_concordant_piggtex_support"]) > 0
            for row in pair_output
        ),
        "coordinate_assigned_independent_recurrent_snp_rows": len(site_output),
    }
    metadata = {
        "version": VERSION,
        "status": "COMPLETED_PENDING_INDEPENDENT_VALIDATION",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_contract": {
            "primary_pair_universe": "39 pair-level tissue-enriched recurrent haplotype ASE gene--tissue pairs",
            "validation_order": "independent GATK SNP recurrence was completed and frozen before this intersection",
            "strong_coordinate_assignment": "unique protein-coding exonic overlap assigned to the target Ensembl 115 gene",
            "data_independence_caution": "the SNP and haplotype analyses use the same RNA-seq data and are computationally independent, not independent biological datasets",
            "piggtex_absence": "not observed in a significant-only table, not a tested negative",
        },
        "inputs": {
            "tissue_enriched_pairs": str(args.tissue_enriched_pairs.resolve()),
            "gene_assignments": str(args.gene_assignments.resolve()),
            "direction_comparison": str(args.direction_comparison.resolve()),
        },
        "input_fingerprints": {
            "tissue_enriched_pairs": fingerprint(args.tissue_enriched_pairs),
            "gene_assignments": fingerprint(args.gene_assignments),
            "direction_comparison": fingerprint(args.direction_comparison),
        },
        "counts": counts,
        "best_support_class_counts": dict(sorted(class_counts.items())),
    }
    metadata_path = args.output_dir / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"counts": counts, "classes": dict(sorted(class_counts.items()))}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(2)
