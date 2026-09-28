#!/usr/bin/env python3
"""Build multilayer PigGTEx evidence for the current 39 tissue-enriched pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


VERSION = "0.1.0"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def values(text: str) -> set[str]:
    return {value for value in text.split(",") if value}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="JSON mapping primary_39, tissue_map, piggtex_gene_profiles, independent_snp_support to TSV paths")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite-owned-outputs", action="store_true")
    args = parser.parse_args()

    inputs = {k: Path(v).expanduser().resolve() for k,v in json.loads(args.inputs.read_text()).items()}
    missing = [str(path) for path in inputs.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))

    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "tissue_enriched_pair_multilayer_evidence.tsv"
    summary_path = output_dir / "evidence_layer_summary.tsv"
    metadata_path = output_dir / "run_metadata.json"
    owned = [evidence_path, summary_path, metadata_path]
    existing = [str(path) for path in owned if path.exists()]
    if existing and not args.overwrite_owned_outputs:
        raise FileExistsError("Owned outputs already exist; use --overwrite-owned-outputs:\n" + "\n".join(existing))

    pairs = read_tsv(inputs["primary_39"])
    if len(pairs) != 39 or len({(row["gene_id"], row["target_tissue"]) for row in pairs}) != 39:
        raise ValueError("The primary pair input is not exactly 39 unique gene-tissue pairs")

    tissue_map = {row["internal_tissue"]: row["piggtex_tissue"] for row in read_tsv(inputs["tissue_map"])}
    profiles = {row["internal_gene_id"]: row for row in read_tsv(inputs["piggtex_gene_profiles"])}
    snp_rows = {
        (row["gene_id"], row["target_tissue"]): row
        for row in read_tsv(inputs["independent_snp_support"])
    }
    if len(snp_rows) != 39:
        raise ValueError("The independent SNP support input is not exactly 39 unique pairs")

    rows: list[dict[str, object]] = []
    for pair in pairs:
        key = (pair["gene_id"], pair["target_tissue"])
        if pair["target_tissue"] not in tissue_map:
            raise ValueError(f"No PigGTEx tissue mapping for {pair['target_tissue']}")
        if pair["gene_id"] not in profiles:
            raise ValueError(f"No PigGTEx gene profile for {pair['gene_id']}")
        if key not in snp_rows:
            raise ValueError(f"No independent SNP support row for {key}")

        profile = profiles[pair["gene_id"]]
        snp = snp_rows[key]
        target = tissue_map[pair["target_tissue"]]
        if profile["gene_name"] != pair["gene_name"] or snp["gene_name"] != pair["gene_name"]:
            raise ValueError(f"Gene-name mismatch for {key}")

        gene_ase = int(target in values(profile["PigGTEx_significant_ASE_tissues_among_10_comparable"]))
        eqtl_tested = int(target in values(profile["PigGTEx_cis_eQTL_tested_tissues_among_10_comparable"]))
        cis_egene = int(target in values(profile["PigGTEx_cis_eGene_tissues_among_10_comparable"]))
        exact_count = int(snp["n_unique_exonic_snps_with_piggtex_same_tissue_exact_support"])
        comparable_count = int(snp["n_unique_exonic_snps_with_direction_comparable_piggtex_support"])
        concordant_count = int(snp["n_unique_exonic_snps_with_direction_concordant_piggtex_support"])
        discordant_count = int(snp["n_unique_exonic_snps_with_direction_discordant_piggtex_support"])
        if cis_egene and not eqtl_tested:
            raise ValueError(f"cis-eGene without eQTL testability for {key}")
        if comparable_count > exact_count or concordant_count + discordant_count > comparable_count:
            raise ValueError(f"Invalid exact-site evidence hierarchy for {key}")

        pattern = []
        if gene_ase:
            pattern.append("PIGGTEX_GENE_ASE_OBSERVED")
        if cis_egene:
            pattern.append("PIGGTEX_CIS_EGENE")
        if exact_count:
            pattern.append("PIGGTEX_EXACT_SITE_OBSERVED")
        if concordant_count:
            pattern.append("DIRECTION_CONCORDANT_EXACT_SITE")
        if not pattern:
            pattern.append("NO_POSITIVE_PIGGTEX_EVIDENCE_OBSERVED")

        rows.append(
            {
                "gene_id": pair["gene_id"],
                "gene_name": pair["gene_name"],
                "target_tissue": pair["target_tissue"],
                "piggtex_mapped_tissue": target,
                "n_target_eligible_animals": pair["n_target_eligible_animals_total"],
                "n_target_ase_positive_animals": pair["n_target_core_animals_total"],
                "tissue_enrichment_q_global_bh": pair["tissue_enrichment_q_global_BH"],
                "recurrence_q_bh_leave_one_gene_out": pair["recurrence_q_bh_leave_one_gene_out"],
                "piggtex_target_gene_ase_observed_significant_only": gene_ase,
                "piggtex_target_eqtl_tested": eqtl_tested,
                "piggtex_target_is_cis_egene": cis_egene,
                "n_coordinate_assigned_independent_recurrent_snps": snp["n_coordinate_assigned_independent_recurrent_snps"],
                "n_unique_exonic_independent_recurrent_snps": snp["n_unique_protein_coding_exonic_independent_recurrent_snps"],
                "n_exact_same_tissue_piggtex_site_matches": exact_count,
                "n_direction_comparable_exact_site_matches": comparable_count,
                "n_direction_concordant_exact_site_matches": concordant_count,
                "n_direction_discordant_exact_site_matches": discordant_count,
                "has_exact_same_tissue_piggtex_site_match": int(exact_count > 0),
                "has_direction_comparable_exact_site_match": int(comparable_count > 0),
                "has_direction_concordant_exact_site_match": int(concordant_count > 0),
                "best_available_snp_support_class": snp["best_available_snp_support_class"],
                "positive_external_evidence_pattern": "+".join(pattern),
                "piggtex_ase_absence_interpretation": "not observed in significant-only catalog, not a tested negative",
            }
        )

    rows.sort(key=lambda row: (str(row["target_tissue"]), str(row["gene_name"]), str(row["gene_id"])))
    fieldnames = list(rows[0])
    write_tsv(evidence_path, rows, fieldnames)

    layers = [
        (
            "matched_tissue_gene_ase_observed",
            "piggtex_target_gene_ase_observed_significant_only",
            len(rows),
            "PigGTEx significant-only gene-ASE table; absence is not a tested negative",
        ),
        (
            "matched_tissue_cis_egene",
            "piggtex_target_is_cis_egene",
            sum(int(row["piggtex_target_eqtl_tested"]) for row in rows),
            "PigGTEx cis-eQTL tested background",
        ),
        (
            "exact_same_tissue_recurrent_snp_match",
            "has_exact_same_tissue_piggtex_site_match",
            sum(int(row["n_unique_exonic_independent_recurrent_snps"]) > 0 for row in rows),
            "Pairs with at least one unique exonic independent recurrent GATK SNP",
        ),
        (
            "direction_concordant_exact_snp_match",
            "has_direction_concordant_exact_site_match",
            sum(int(row["has_direction_comparable_exact_site_match"]) for row in rows),
            "Pairs with at least one allele-harmonized direction-comparable exact site",
        ),
    ]
    summary_rows = []
    for layer, field, evaluable, scope in layers:
        positive = sum(int(row[field]) for row in rows)
        summary_rows.append(
            {
                "evidence_layer": layer,
                "n_pairs_positive": positive,
                "n_pairs_evaluable": evaluable,
                "n_pairs_total_universe": len(rows),
                "fraction_positive_among_evaluable": f"{positive / evaluable:.12g}",
                "evidence_scope": scope,
            }
        )
    write_tsv(
        summary_path,
        summary_rows,
        [
            "evidence_layer",
            "n_pairs_positive",
            "n_pairs_evaluable",
            "n_pairs_total_universe",
            "fraction_positive_among_evaluable",
            "evidence_scope",
        ],
    )

    pattern_counts = Counter(str(row["positive_external_evidence_pattern"]) for row in rows)
    metadata = {
        "version": VERSION,
        "status": "COMPLETED_PENDING_VALIDATION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Descriptive multilayer PigGTEx annotation of the current 39 primary tissue-enriched haplotype-ASE pairs",
        "interpretation_boundaries": [
            "PigGTEx gene- and site-level ASE tables contain significant observations only; absence is not a tested negative.",
            "Exact-site support is descriptive and is not a global enrichment test.",
            "The independent GATK SNP and haplotype analyses use the same RNA-seq data and are not independent biological replication.",
        ],
        "counts": {
            "primary_pairs": len(rows),
            "matched_tissue_gene_ase_observed": sum(int(row["piggtex_target_gene_ase_observed_significant_only"]) for row in rows),
            "matched_tissue_eqtl_tested": sum(int(row["piggtex_target_eqtl_tested"]) for row in rows),
            "matched_tissue_cis_egene": sum(int(row["piggtex_target_is_cis_egene"]) for row in rows),
            "pairs_with_unique_exonic_independent_recurrent_snp": sum(
                int(row["n_unique_exonic_independent_recurrent_snps"]) > 0 for row in rows
            ),
            "pairs_with_exact_same_tissue_recurrent_snp_match": sum(int(row["has_exact_same_tissue_piggtex_site_match"]) for row in rows),
            "pairs_with_direction_comparable_exact_site_match": sum(
                int(row["has_direction_comparable_exact_site_match"]) for row in rows
            ),
            "pairs_with_direction_concordant_exact_site_match": sum(int(row["has_direction_concordant_exact_site_match"]) for row in rows),
            "positive_evidence_patterns": dict(sorted(pattern_counts.items())),
        },
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
            for name, path in inputs.items()
        },
        "outputs": {
            evidence_path.name: {"path": str(evidence_path), "sha256": sha256(evidence_path), "size_bytes": evidence_path.stat().st_size},
            summary_path.name: {"path": str(summary_path), "sha256": sha256(summary_path), "size_bytes": summary_path.stat().st_size},
        },
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": metadata["status"], "counts": metadata["counts"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
