#!/usr/bin/env python3
"""Prepare all independently recurrent GATK SNPs for PigGTEx gene-eQTL lookup."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple


VERSION = "ase_gatk_recurrent_snp_piggtex_eqtl_full_universe_v0.1.0"
GENE_REQUIRED = {
    "tissue",
    "physical_snp_id",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "n_animals_testable",
    "n_animals_with_snp_level_ase",
    "recurrence_fraction",
    "recurrence_q_bh_global",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "site_gene_relation",
    "best_transcript_feature",
}
SITE_REQUIRED = {
    "tissue",
    "physical_snp_id",
    "all_gene_assignment_class",
    "n_gene_spans",
}
CROSSWALK_REQUIRED = {
    "internal_gene_id",
    "piggtex_gene_id",
    "piggtex_gene_name_v100",
    "mapping_status",
}
TISSUE_REQUIRED = {"internal_tissue", "piggtex_tissue"}
PAIR_REQUIRED = {"gene_id", "target_tissue", "pair_level_primary"}
DIRECTION_REQUIRED = {
    "internal_tissue",
    "physical_snp_id",
    "piggtex_significant_rows_same_mapped_tissue",
    "major_allele_direction_comparison",
}
EQTL_REQUIRED = {
    "phenotype_id",
    "variant_id",
    "qval",
    "pval_adj_BH",
    "pval_nominal",
    "slope",
    "slope_se",
    "is_eGene",
}

CANDIDATE_FIELDS = [
    "internal_gene_id",
    "internal_gene_name",
    "internal_gene_biotype",
    "internal_tissue",
    "piggtex_gene_id",
    "piggtex_gene_name_v100",
    "piggtex_tissue",
    "chrom",
    "position",
    "allele_1",
    "allele_2",
    "internal_variant_unordered_key",
    "site_gene_relation",
    "best_transcript_feature",
    "n_animals_testable",
    "n_animals_with_snp_level_ase",
    "recurrence_fraction",
    "recurrence_q_bh_global",
    "is_primary_tissue_enriched_pair",
    "is_formally_tissue_enriched_ASE",
    "piggtex_same_tissue_significant_site_ase_observed",
    "major_allele_direction_comparison",
]

MAPPING_FIELDS = CANDIDATE_FIELDS[:4] + [
    "piggtex_gene_id",
    "piggtex_gene_name_v100",
    "gene_crosswalk_status",
] + CANDIDATE_FIELDS[6:]

LEAD_FIELDS = CANDIDATE_FIELDS + [
    "piggtex_gene_eqtl_record_present",
    "piggtex_is_eGene",
    "piggtex_lead_variant_id",
    "piggtex_lead_variant_unordered_key",
    "internal_snp_is_exact_piggtex_lead_variant",
    "piggtex_eqtl_qval",
    "piggtex_eqtl_pval_adj_BH",
    "piggtex_eqtl_pval_nominal",
    "piggtex_eqtl_slope",
    "piggtex_eqtl_slope_se",
]


class AnalysisError(RuntimeError):
    pass


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def iter_tsv(path: Path, required: set[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AnalysisError(f"{path}: missing columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value if value is not None else "" for key, value in row.items()}


def clean_gene(value: str) -> str:
    return value.strip().split(".")[0]


def variant_key(contig: str, position: int, allele_1: str, allele_2: str) -> str:
    first, second = sorted((allele_1.upper(), allele_2.upper()))
    return f"{contig.removeprefix('chr')}_{position}_{first}_{second}"


def variant_key_from_id(value: str) -> str:
    fields = value.strip().split("_")
    if len(fields) != 4:
        return ""
    try:
        position = int(fields[1])
    except ValueError:
        return ""
    return variant_key(fields[0], position, fields[2], fields[3])


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> int:
    count = 0
    with open_text(path, "wt") as handle:
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


def load_lead_eqtl(eqtl_root: Path, piggtex_tissues: set[str]):
    result: Dict[Tuple[str, str], Dict[str, str]] = {}
    audit = []
    for tissue in sorted(piggtex_tissues):
        path = eqtl_root / f"{tissue}.cis_qtl_fdr0.05.txt.gz"
        if not path.is_file():
            audit.append({"piggtex_tissue": tissue, "file_present": False, "n_rows": 0})
            continue
        n_rows = 0
        for row in iter_tsv(path, EQTL_REQUIRED):
            n_rows += 1
            gene = clean_gene(row["phenotype_id"])
            key = (tissue, gene)
            if key in result:
                raise AnalysisError(f"duplicate permutation eQTL gene row: {key}")
            result[key] = row
        audit.append({"piggtex_tissue": tissue, "file_present": True, "n_rows": n_rows})
    return result, audit


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-assignments", type=Path, required=True)
    parser.add_argument("--site-annotations", type=Path, required=True)
    parser.add_argument("--gene-crosswalk", type=Path, required=True)
    parser.add_argument("--tissue-map", type=Path, required=True)
    parser.add_argument("--tissue-enriched-pairs", type=Path, required=True)
    parser.add_argument("--direction-comparison", type=Path, required=True)
    parser.add_argument("--eqtl-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    input_paths = [
        args.gene_assignments,
        args.site_annotations,
        args.gene_crosswalk,
        args.tissue_map,
        args.tissue_enriched_pairs,
        args.direction_comparison,
    ]
    for path in input_paths:
        if not path.is_file():
            raise AnalysisError(f"missing input: {path}")
    if not args.eqtl_root.is_dir():
        raise AnalysisError(f"missing eQTL directory: {args.eqtl_root}")
    if args.output_dir.exists():
        if not args.force:
            raise AnalysisError(f"output directory already exists: {args.output_dir}")
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True)

    crosswalk = {}
    for row in iter_tsv(args.gene_crosswalk, CROSSWALK_REQUIRED):
        gene = clean_gene(row["internal_gene_id"])
        if gene in crosswalk:
            raise AnalysisError(f"duplicate crosswalk gene: {gene}")
        crosswalk[gene] = row
    tissue_map = {}
    for row in iter_tsv(args.tissue_map, TISSUE_REQUIRED):
        tissue = row["internal_tissue"]
        if tissue in tissue_map:
            raise AnalysisError(f"duplicate tissue mapping: {tissue}")
        tissue_map[tissue] = row["piggtex_tissue"]
    formal_pairs = {
        (clean_gene(row["gene_id"]), row["target_tissue"])
        for row in iter_tsv(args.tissue_enriched_pairs, PAIR_REQUIRED)
        if row["pair_level_primary"] == "1"
    }
    direction = {}
    for row in iter_tsv(args.direction_comparison, DIRECTION_REQUIRED):
        key = (row["internal_tissue"], row["physical_snp_id"])
        if key in direction:
            raise AnalysisError(f"duplicate direction key: {key}")
        direction[key] = row

    site_keys = set()
    n_site_rows = 0
    n_intergenic = 0
    for row in iter_tsv(args.site_annotations, SITE_REQUIRED):
        n_site_rows += 1
        key = (row["tissue"], row["physical_snp_id"])
        if key in site_keys:
            raise AnalysisError(f"duplicate site annotation key: {key}")
        site_keys.add(key)
        n_intergenic += int(int(row["n_gene_spans"]) == 0)

    mapping_rows = []
    candidate_rows = []
    seen_assignments = set()
    assigned_site_keys = set()
    mapped_site_keys = set()
    for row in iter_tsv(args.gene_assignments, GENE_REQUIRED):
        tissue = row["tissue"]
        site_key = (tissue, row["physical_snp_id"])
        if site_key not in site_keys:
            raise AnalysisError(f"gene assignment outside site universe: {site_key}")
        assigned_site_keys.add(site_key)
        gene = clean_gene(row["gene_id"])
        context = (gene, tissue, row["physical_snp_id"])
        if context in seen_assignments:
            raise AnalysisError(f"duplicate gene--tissue--site assignment: {context}")
        seen_assignments.add(context)
        if tissue not in tissue_map:
            raise AnalysisError(f"unmapped internal tissue: {tissue}")
        if site_key not in direction:
            raise AnalysisError(f"direction row missing: {site_key}")
        external = crosswalk.get(gene, {})
        pig_gene = clean_gene(external.get("piggtex_gene_id", ""))
        position = int(row["position"])
        alleles = sorted((row["ref_allele"].upper(), row["alt_allele"].upper()))
        direction_row = direction[site_key]
        base = {
            "internal_gene_id": gene,
            "internal_gene_name": row["gene_name"],
            "internal_gene_biotype": row["gene_biotype"],
            "internal_tissue": tissue,
            "piggtex_gene_id": pig_gene,
            "piggtex_gene_name_v100": external.get("piggtex_gene_name_v100", ""),
            "piggtex_tissue": tissue_map[tissue],
            "chrom": row["contig"],
            "position": position,
            "allele_1": alleles[0],
            "allele_2": alleles[1],
            "internal_variant_unordered_key": variant_key(
                row["contig"], position, row["ref_allele"], row["alt_allele"]
            ),
            "site_gene_relation": row["site_gene_relation"],
            "best_transcript_feature": row["best_transcript_feature"],
            "n_animals_testable": row["n_animals_testable"],
            "n_animals_with_snp_level_ase": row["n_animals_with_snp_level_ase"],
            "recurrence_fraction": row["recurrence_fraction"],
            "recurrence_q_bh_global": row["recurrence_q_bh_global"],
            "is_primary_tissue_enriched_pair": int((gene, tissue) in formal_pairs),
            "is_formally_tissue_enriched_ASE": str((gene, tissue) in formal_pairs),
            "piggtex_same_tissue_significant_site_ase_observed": int(
                int(direction_row["piggtex_significant_rows_same_mapped_tissue"]) > 0
            ),
            "major_allele_direction_comparison": direction_row[
                "major_allele_direction_comparison"
            ],
        }
        mapping_rows.append(
            {
                **base,
                "gene_crosswalk_status": external.get("mapping_status", "UNMAPPED"),
            }
        )
        if pig_gene:
            candidate_rows.append(base)
            mapped_site_keys.add(site_key)

    candidate_rows.sort(
        key=lambda row: (
            row["piggtex_tissue"],
            row["piggtex_gene_id"],
            int(row["chrom"]),
            int(row["position"]),
            row["allele_1"],
            row["allele_2"],
        )
    )
    mapping_rows.sort(
        key=lambda row: (
            row["internal_tissue"],
            row["internal_gene_id"],
            int(row["chrom"]),
            int(row["position"]),
        )
    )

    lead_eqtl, eqtl_audit = load_lead_eqtl(
        args.eqtl_root, {row["piggtex_tissue"] for row in candidate_rows}
    )
    lead_rows = []
    for row in candidate_rows:
        lead = lead_eqtl.get((row["piggtex_tissue"], row["piggtex_gene_id"]))
        lead_key = variant_key_from_id(lead["variant_id"]) if lead else ""
        lead_rows.append(
            {
                **row,
                "piggtex_gene_eqtl_record_present": int(lead is not None),
                "piggtex_is_eGene": lead["is_eGene"] if lead else "",
                "piggtex_lead_variant_id": lead["variant_id"] if lead else "",
                "piggtex_lead_variant_unordered_key": lead_key,
                "internal_snp_is_exact_piggtex_lead_variant": int(
                    bool(lead_key) and lead_key == row["internal_variant_unordered_key"]
                ),
                "piggtex_eqtl_qval": lead["qval"] if lead else "",
                "piggtex_eqtl_pval_adj_BH": lead["pval_adj_BH"] if lead else "",
                "piggtex_eqtl_pval_nominal": lead["pval_nominal"] if lead else "",
                "piggtex_eqtl_slope": lead["slope"] if lead else "",
                "piggtex_eqtl_slope_se": lead["slope_se"] if lead else "",
            }
        )

    grouped: Dict[Tuple[str, str], List[Mapping[str, object]]] = defaultdict(list)
    for row in lead_rows:
        grouped[(str(row["internal_gene_id"]), str(row["internal_tissue"]))].append(row)
    pair_rows = []
    for (gene, tissue), rows in grouped.items():
        first = rows[0]
        pair_rows.append(
            {
                "internal_gene_id": gene,
                "internal_gene_name": first["internal_gene_name"],
                "internal_tissue": tissue,
                "piggtex_gene_id": first["piggtex_gene_id"],
                "piggtex_tissue": first["piggtex_tissue"],
                "is_primary_tissue_enriched_pair": first[
                    "is_primary_tissue_enriched_pair"
                ],
                "n_independent_recurrent_snps": len(rows),
                "piggtex_gene_eqtl_record_present": first[
                    "piggtex_gene_eqtl_record_present"
                ],
                "piggtex_is_eGene": first["piggtex_is_eGene"],
                "piggtex_lead_variant_id": first["piggtex_lead_variant_id"],
                "n_internal_snps_matching_piggtex_lead_variant": sum(
                    int(row["internal_snp_is_exact_piggtex_lead_variant"]) for row in rows
                ),
                "matching_internal_snp_ids": ",".join(
                    str(row["internal_variant_unordered_key"])
                    for row in rows
                    if int(row["internal_snp_is_exact_piggtex_lead_variant"])
                ),
                "piggtex_eqtl_qval": first["piggtex_eqtl_qval"],
                "piggtex_eqtl_pval_adj_BH": first["piggtex_eqtl_pval_adj_BH"],
                "piggtex_eqtl_pval_nominal": first["piggtex_eqtl_pval_nominal"],
                "piggtex_eqtl_slope": first["piggtex_eqtl_slope"],
            }
        )
    pair_rows.sort(
        key=lambda row: (
            -int(row["n_internal_snps_matching_piggtex_lead_variant"]),
            -int(str(row["piggtex_is_eGene"]).upper() == "TRUE"),
            row["internal_gene_name"],
            row["internal_tissue"],
        )
    )

    candidate_path = args.output_dir / "full_independent_recurrent_snp_gene_eqtl_input.tsv.gz"
    mapping_path = args.output_dir / "all_gene_assignment_eqtl_mapping_status.tsv.gz"
    lead_path = args.output_dir / "local_permutation_lead_eqtl_lookup.tsv.gz"
    pair_path = args.output_dir / "gene_tissue_lead_eqtl_summary.tsv"
    eqtl_audit_path = args.output_dir / "local_eqtl_file_audit.tsv"
    write_tsv(candidate_path, CANDIDATE_FIELDS, candidate_rows)
    write_tsv(mapping_path, MAPPING_FIELDS, mapping_rows)
    write_tsv(lead_path, LEAD_FIELDS, lead_rows)
    pair_fields = list(pair_rows[0]) if pair_rows else []
    write_tsv(pair_path, pair_fields, pair_rows)
    write_tsv(eqtl_audit_path, list(eqtl_audit[0]), eqtl_audit)

    exact_lead_rows = [
        row for row in lead_rows if int(row["internal_snp_is_exact_piggtex_lead_variant"])
    ]
    counts = {
        "independent_recurrent_tissue_snp_hypotheses": n_site_rows,
        "intergenic_tissue_snp_hypotheses": n_intergenic,
        "tissue_snp_hypotheses_with_any_gene_assignment": len(assigned_site_keys),
        "tissue_snp_hypotheses_with_any_v100_mapped_gene_assignment": len(mapped_site_keys),
        "gene_tissue_snp_assignment_rows": len(mapping_rows),
        "v100_mapped_gene_tissue_snp_candidate_rows": len(candidate_rows),
        "candidate_gene_tissue_pairs": len(pair_rows),
        "candidate_pairs_with_piggtex_gene_eqtl_record": sum(
            int(row["piggtex_gene_eqtl_record_present"]) for row in pair_rows
        ),
        "candidate_pairs_that_are_piggtex_eGenes": sum(
            str(row["piggtex_is_eGene"]).upper() == "TRUE" for row in pair_rows
        ),
        "candidate_gene_tissue_snp_rows_matching_piggtex_lead_variant": len(
            exact_lead_rows
        ),
        "candidate_gene_tissue_pairs_with_exact_lead_variant_match": sum(
            int(row["n_internal_snps_matching_piggtex_lead_variant"]) > 0 for row in pair_rows
        ),
        "primary_39_pairs_represented_in_mapped_candidates": len(
            {
                (row["internal_gene_id"], row["internal_tissue"])
                for row in candidate_rows
                if int(row["is_primary_tissue_enriched_pair"])
            }
        ),
    }
    metadata = {
        "version": VERSION,
        "status": "PASS_LOCAL_PREFLIGHT_AND_LEAD_EQTL_LOOKUP",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_contract": {
            "candidate_universe": "all coordinate-assigned gene contexts from 13,919 independently recurrent GATK tissue--SNP hypotheses",
            "gene_mapping": "Ensembl release 115 internal gene to PigGTEx Ensembl release 100 crosswalk",
            "local_eqtl_scope": "PigGTEx permutation-mode top association per protein-coding gene; not the full significant variant set",
            "pip_status": "not evaluated locally because the available finemapped_enQTL archive is enhancer-QTL, not gene eQTL",
            "required_server_archives": [
                "PigGTEx_v0.significant_eQTL.tar",
                "PigGTEx_v0.finemapped_eQTL.tar.gz",
            ],
        },
        "inputs": {name: str(path.resolve()) for name, path in zip(
            [
                "gene_assignments",
                "site_annotations",
                "gene_crosswalk",
                "tissue_map",
                "tissue_enriched_pairs",
                "direction_comparison",
            ],
            input_paths,
        )} | {"eqtl_root": str(args.eqtl_root.resolve())},
        "input_fingerprints": {
            name: fingerprint(path)
            for name, path in zip(
                [
                    "gene_assignments",
                    "site_annotations",
                    "gene_crosswalk",
                    "tissue_map",
                    "tissue_enriched_pairs",
                    "direction_comparison",
                ],
                input_paths,
            )
        },
        "counts": counts,
    }
    metadata_path = args.output_dir / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(counts, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(2)
