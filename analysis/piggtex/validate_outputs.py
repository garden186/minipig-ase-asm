#!/usr/bin/env python3
"""Independent structural and statistical validator for PigGTEx validation."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

from piggtex_validation import (
    ENRICHMENT_FIELDS,
    FORMAL_EXPECTED,
    TISSUE_ORDER,
    bh_values,
    fisher_upper_tail,
    stratified_exact_tail,
)


class OutputError(RuntimeError):
    pass


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(left: float, right: float, tolerance: float = 1e-10) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def validate(out_dir: Path) -> Dict[str, object]:
    errors = []
    required = [
        "tissue_mapping.tsv",
        "gene_crosswalk.tsv.gz",
        "formal_8_external_validation.tsv",
        "matched_double_testable_universe.tsv.gz",
        "opportunity_strata_counts.tsv.gz",
        "tissue_eqtl_enrichment.tsv",
        "analysis_qc.json",
    ]
    for name in required:
        if not (out_dir / name).is_file():
            errors.append(f"missing output: {name}")
    if errors:
        return {"status": "FAIL", "errors": errors}

    mapping = read_tsv(out_dir / "tissue_mapping.tsv")
    formal = read_tsv(out_dir / "formal_8_external_validation.tsv")
    enrichment = read_tsv(out_dir / "tissue_eqtl_enrichment.tsv")
    strata_rows = read_tsv(out_dir / "opportunity_strata_counts.tsv.gz")
    qc = json.loads((out_dir / "analysis_qc.json").read_text(encoding="utf-8"))

    if len(mapping) != len(TISSUE_ORDER):
        errors.append(f"tissue mapping rows: {len(mapping)}")
    if [row["internal_tissue"] for row in mapping] != list(TISSUE_ORDER):
        errors.append("tissue mapping order/content mismatch")
    if len(formal) != FORMAL_EXPECTED:
        errors.append(f"formal rows: {len(formal)}")
    if any(row["piggtex_ase_observed_significant_only"] != "1" for row in formal):
        errors.append("not all formal pairs have matched-tissue PigGTEx ASE support")
    if sum(int(row["piggtex_is_eGene"]) for row in formal) != 4:
        errors.append("formal eGene count is not 4")
    if len(enrichment) != 44:
        errors.append(f"enrichment rows: {len(enrichment)} expected 44")
    if set(enrichment[0]) != set(ENRICHMENT_FIELDS):
        errors.append("enrichment header mismatch")

    strata_by_key: Dict[Tuple[str, str], List[Tuple[int, int, int, int]]] = defaultdict(list)
    for row in strata_rows:
        key = (row["analysis_mode"], row["internal_tissue"])
        strata_by_key[key].append(
            tuple(
                int(row[field])
                for field in (
                    "recurrent_eGene",
                    "recurrent_non_eGene",
                    "nonrecurrent_eGene",
                    "nonrecurrent_non_eGene",
                )
            )
        )

    modes: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in enrichment:
        modes[row["analysis_mode"]].append(row)
        table = tuple(
            int(row[field])
            for field in (
                "recurrent_eGene",
                "recurrent_non_eGene",
                "nonrecurrent_eGene",
                "nonrecurrent_non_eGene",
            )
        )
        if sum(table) != int(row["n_double_testable"]):
            errors.append(f"table total mismatch: {row['analysis_mode']}|{row['internal_tissue']}")
        fisher = fisher_upper_tail(*table)
        if not close(fisher, float(row["fisher_p_one_sided"])):
            errors.append(f"Fisher mismatch: {row['analysis_mode']}|{row['internal_tissue']}")
        strata = strata_by_key[(row["analysis_mode"], row["internal_tissue"])]
        if tuple(map(sum, zip(*strata))) != table:
            errors.append(f"strata total mismatch: {row['analysis_mode']}|{row['internal_tissue']}")
        stratified = stratified_exact_tail(strata)
        if not close(stratified, float(row["opportunity_stratified_exact_p_one_sided"])):
            errors.append(f"stratified P mismatch: {row['analysis_mode']}|{row['internal_tissue']}")

    for mode, rows in modes.items():
        if len(rows) != len(TISSUE_ORDER):
            errors.append(f"mode {mode} has {len(rows)} tissue rows")
            continue
        fisher_q = bh_values([float(row["fisher_p_one_sided"]) for row in rows])
        stratified_q = bh_values(
            [float(row["opportunity_stratified_exact_p_one_sided"]) for row in rows]
        )
        for row, q1, q2 in zip(rows, fisher_q, stratified_q):
            if not close(q1, float(row["fisher_q_BH_across_tissues"])):
                errors.append(f"Fisher BH mismatch: {mode}|{row['internal_tissue']}")
            if not close(q2, float(row["opportunity_stratified_q_BH_across_tissues"])):
                errors.append(f"stratified BH mismatch: {mode}|{row['internal_tissue']}")

    if qc.get("status") != "PASS":
        errors.append("analysis_qc status is not PASS")
    return {
        "status": "PASS" if not errors else "FAIL",
        "version": "0.1.0",
        "errors": errors,
        "checks": {
            "n_tissue_mapping_rows": len(mapping),
            "n_formal_pairs": len(formal),
            "n_formal_with_matched_ASE": sum(
                int(row["piggtex_ase_observed_significant_only"]) for row in formal
            ),
            "n_formal_cis_eGenes": sum(int(row["piggtex_is_eGene"]) for row in formal),
            "n_enrichment_modes": len(modes),
            "n_enrichment_rows": len(enrichment),
            "n_strata_rows": len(strata_rows),
        },
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = validate(args.out_dir)
    (args.out_dir / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
