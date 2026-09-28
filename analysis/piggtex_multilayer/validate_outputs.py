#!/usr/bin/env python3
"""Independently validate the multilayer PigGTEx evidence package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[4])
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = root / "analysis/results/ase/ase_tissue_enriched_piggtex_multilayer_v0.1.0/full_cohort"
    evidence_path = output_dir / "tissue_enriched_pair_multilayer_evidence.tsv"
    summary_path = output_dir / "evidence_layer_summary.tsv"
    metadata_path = output_dir / "run_metadata.json"
    validation_path = output_dir / "validation.json"

    errors: list[str] = []
    for path in [evidence_path, summary_path, metadata_path]:
        if not path.exists():
            errors.append(f"Missing output: {path}")
    if errors:
        validation_path.parent.mkdir(parents=True, exist_ok=True)
        validation_path.write_text(json.dumps({"status": "FAIL", "errors": errors}, indent=2) + "\n", encoding="utf-8")
        raise SystemExit(1)

    rows = read_tsv(evidence_path)
    summaries = {row["evidence_layer"]: row for row in read_tsv(summary_path)}
    with metadata_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    checks = {
        "exactly_39_unique_pairs": len(rows) == 39 and len({(row["gene_id"], row["target_tissue"]) for row in rows}) == 39,
        "exactly_38_pairs_eqtl_tested_in_target_tissue": sum(
            int(row["piggtex_target_eqtl_tested"]) for row in rows
        )
        == 38,
        "cis_egene_implies_tested": all(
            int(row["piggtex_target_is_cis_egene"]) <= int(row["piggtex_target_eqtl_tested"]) for row in rows
        ),
        "site_hierarchy_valid": all(
            int(row["n_direction_comparable_exact_site_matches"])
            <= int(row["n_exact_same_tissue_piggtex_site_matches"])
            and int(row["n_direction_concordant_exact_site_matches"])
            + int(row["n_direction_discordant_exact_site_matches"])
            <= int(row["n_direction_comparable_exact_site_matches"])
            for row in rows
        ),
        "matched_tissue_gene_ase_count_34": sum(
            int(row["piggtex_target_gene_ase_observed_significant_only"]) for row in rows
        )
        == 34,
        "matched_tissue_cis_egene_count_19": sum(int(row["piggtex_target_is_cis_egene"]) for row in rows) == 19,
        "exact_same_tissue_recurrent_snp_pair_count_22": sum(
            int(row["has_exact_same_tissue_piggtex_site_match"]) for row in rows
        )
        == 22,
        "direction_concordant_exact_site_pair_count_12": sum(
            int(row["has_direction_concordant_exact_site_match"]) for row in rows
        )
        == 12,
        "unique_exonic_recurrent_snp_opportunity_count_27": sum(
            int(row["n_unique_exonic_independent_recurrent_snps"]) > 0 for row in rows
        )
        == 27,
        "direction_comparable_exact_site_pair_count_16": sum(
            int(row["has_direction_comparable_exact_site_match"]) for row in rows
        )
        == 16,
        "summary_has_four_layers": len(summaries) == 4,
        "metadata_counts_match": metadata["counts"]["primary_pairs"] == 39
        and metadata["counts"]["matched_tissue_gene_ase_observed"] == 34
        and metadata["counts"]["matched_tissue_eqtl_tested"] == 38
        and metadata["counts"]["matched_tissue_cis_egene"] == 19
        and metadata["counts"]["pairs_with_unique_exonic_independent_recurrent_snp"] == 27
        and metadata["counts"]["pairs_with_exact_same_tissue_recurrent_snp_match"] == 22
        and metadata["counts"]["pairs_with_direction_comparable_exact_site_match"] == 16
        and metadata["counts"]["pairs_with_direction_concordant_exact_site_match"] == 12,
        "input_fingerprints_reproduced": all(
            Path(item["path"]).exists() and sha256(Path(item["path"])) == item["sha256"]
            for item in metadata["inputs"].values()
        ),
        "output_fingerprints_reproduced": sha256(evidence_path)
        == metadata["outputs"][evidence_path.name]["sha256"]
        and sha256(summary_path) == metadata["outputs"][summary_path.name]["sha256"],
    }
    errors.extend(name for name, passed in checks.items() if not passed)
    validation = {
        "status": "PASS" if not errors else "FAIL",
        "checks": checks,
        "errors": errors,
        "observed": {
            "primary_pairs": len(rows),
            "matched_tissue_gene_ase_observed": sum(
                int(row["piggtex_target_gene_ase_observed_significant_only"]) for row in rows
            ),
            "matched_tissue_eqtl_tested": sum(int(row["piggtex_target_eqtl_tested"]) for row in rows),
            "matched_tissue_cis_egene": sum(int(row["piggtex_target_is_cis_egene"]) for row in rows),
            "pairs_with_unique_exonic_independent_recurrent_snp": sum(
                int(row["n_unique_exonic_independent_recurrent_snps"]) > 0 for row in rows
            ),
            "pairs_with_exact_same_tissue_recurrent_snp_match": sum(
                int(row["has_exact_same_tissue_piggtex_site_match"]) for row in rows
            ),
            "pairs_with_direction_comparable_exact_site_match": sum(
                int(row["has_direction_comparable_exact_site_match"]) for row in rows
            ),
            "pairs_with_direction_concordant_exact_site_match": sum(
                int(row["has_direction_concordant_exact_site_match"]) for row in rows
            ),
        },
        "output_fingerprints": {
            evidence_path.name: {"sha256": sha256(evidence_path), "size_bytes": evidence_path.stat().st_size},
            summary_path.name: {"sha256": sha256(summary_path), "size_bytes": summary_path.stat().st_size},
            metadata_path.name: {"sha256": sha256(metadata_path), "size_bytes": metadata_path.stat().st_size},
        },
    }
    validation_path.write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": validation["status"], "observed": validation["observed"]}, ensure_ascii=False))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
