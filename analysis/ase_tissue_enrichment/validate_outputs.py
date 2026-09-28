#!/usr/bin/env python3
"""Independent cross-file validation for formal ASE tissue-enrichment outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple


VERSION = "0.1.0"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def bh_values(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 1.0
    for reverse_index in range(len(order) - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, min(1.0, values[index] * len(values) / rank))
        adjusted[index] = running
    return adjusted


def close(left: float, right: float, tolerance: float = 1e-9) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--recurrence-dir", type=Path, required=True)
    args = parser.parse_args()
    errors: List[str] = []
    checks: Dict[str, object] = {}
    metadata = json.loads((args.out_dir / "analysis_metadata.json").read_text(encoding="utf-8"))
    if metadata.get("status") != "PASS" or metadata.get("version") != VERSION:
        errors.append("analysis metadata status/version mismatch")
    parameters = metadata.get("parameters", {})
    fdr = float(parameters.get("fdr", 0.05))
    min_animals = int(parameters.get("min_informative_animals", 3))
    min_target_core = int(parameters.get("min_informative_target_core", 2))
    summary_rows = {
        int(row["fragment_threshold"]): row
        for row in rows(args.out_dir / "threshold_summary.tsv")
    }

    for threshold in (15, 20, 30):
        prefix = f"threshold_{threshold}"
        contrast_path = args.out_dir / f"{prefix}_tissue_contrasts_all.tsv.gz"
        omnibus_path = args.out_dir / f"{prefix}_gene_omnibus.tsv.gz"
        formal_path = args.out_dir / f"{prefix}_formal_tissue_enriched_ase.tsv"
        status_path = args.out_dir / f"{prefix}_recurrent_specificity_status.tsv"
        recurrence_path = args.recurrence_dir / f"{prefix}_statistically_recurrent.tsv"
        summary = summary_rows[threshold]
        contrast_rows = list(rows(contrast_path))
        contrast_keys = [(row["gene_id"], row["target_tissue"]) for row in contrast_rows]
        if len(contrast_keys) != len(set(contrast_keys)):
            errors.append(f"{prefix}: duplicate contrast keys")
        evaluable = [row for row in contrast_rows if row["is_contrast_evaluable"] == "1"]
        if any(int(row["n_informative_animals"]) < min_animals for row in evaluable):
            errors.append(f"{prefix}: evaluable contrast below minimum animals")
        raw_p = [float(row["tissue_enrichment_exact_p_one_sided"]) for row in evaluable]
        expected_q = bh_values(raw_p)
        for row, q_value in zip(evaluable, expected_q):
            if not close(float(row["tissue_enrichment_q_global_BH"]), q_value):
                errors.append(f"{prefix}: contrast BH mismatch")
                break
        non_evaluable = [row for row in contrast_rows if row["is_contrast_evaluable"] == "0"]
        if any(row["tissue_enrichment_q_global_BH"] != "NA" for row in non_evaluable):
            errors.append(f"{prefix}: non-evaluable contrast has q value")

        omnibus_rows = list(rows(omnibus_path))
        omnibus_keys = [row["gene_id"] for row in omnibus_rows]
        if len(omnibus_keys) != len(set(omnibus_keys)):
            errors.append(f"{prefix}: duplicate omnibus genes")
        omnibus_p = [float(row["gene_omnibus_p_bonferroni"]) for row in omnibus_rows]
        omnibus_q = bh_values(omnibus_p)
        for row, q_value in zip(omnibus_rows, omnibus_q):
            if not close(float(row["gene_omnibus_q_global_BH"]), q_value):
                errors.append(f"{prefix}: omnibus BH mismatch")
                break

        formal_rows = list(rows(formal_path))
        formal_keys = {(row["gene_id"], row["target_tissue"]) for row in formal_rows}
        flagged_keys = {
            (row["gene_id"], row["target_tissue"])
            for row in contrast_rows
            if row["is_formally_tissue_enriched_ase"] == "1"
        }
        if formal_keys != flagged_keys:
            errors.append(f"{prefix}: formal subset differs from contrast flags")
        for row in formal_rows:
            if row["is_primary_recurrent"] != "1":
                errors.append(f"{prefix}: formal call is not recurrent")
                break
            if float(row["tissue_enrichment_q_global_BH"]) > fdr:
                errors.append(f"{prefix}: formal target contrast fails FDR")
                break
            if float(row["gene_omnibus_q_global_BH"]) > fdr:
                errors.append(f"{prefix}: formal gene fails omnibus FDR")
                break
            if int(row["n_informative_animals"]) < min_animals:
                errors.append(f"{prefix}: formal call below animal gate")
                break
            if int(row["n_target_core_informative"]) < min_target_core:
                errors.append(f"{prefix}: formal call below core gate")
                break
            if float(row["conditional_core_excess"]) <= 0.0:
                errors.append(f"{prefix}: formal call lacks positive excess")
                break

        recurrence_count = sum(1 for _ in rows(recurrence_path))
        status_rows = list(rows(status_path))
        if len(status_rows) != recurrence_count:
            errors.append(f"{prefix}: recurrence status row count mismatch")
        status_formal = {
            (row["gene_id"], row["tissue"])
            for row in status_rows
            if row["is_formally_tissue_enriched_ase"] == "1"
        }
        if status_formal != formal_keys:
            errors.append(f"{prefix}: recurrence status formal keys mismatch")

        observed = {
            "n_cross_tissue_candidate_contrasts": len(contrast_rows),
            "n_tissue_contrasts_tested": len(evaluable),
            "n_genes_omnibus_tested": len(omnibus_rows),
            "n_genes_omnibus_fdr_pass": sum(
                row["gene_omnibus_fdr_pass"] == "1" for row in omnibus_rows
            ),
            "n_tissue_contrasts_global_fdr_pass": sum(
                float(row["tissue_enrichment_q_global_BH"]) <= fdr for row in evaluable
            ),
            "n_primary_recurrent_pairs": recurrence_count,
            "n_formal_tissue_enriched_recurrent_pairs": len(formal_rows),
            "n_formal_tissue_enriched_genes": len({row["gene_id"] for row in formal_rows}),
            "n_formal_tissue_restricted_recurrence_pairs": sum(
                row["is_formally_tissue_restricted_recurrence"] == "1"
                for row in formal_rows
            ),
        }
        for field, value in observed.items():
            if int(summary[field]) != value:
                errors.append(f"{prefix}: summary mismatch for {field}")
        checks[prefix] = observed

    primary_sensitivity = list(
        rows(args.out_dir / "primary_formal_tissue_enriched_with_sensitivity.tsv")
    )
    if len(primary_sensitivity) != int(
        summary_rows[15]["n_formal_tissue_enriched_recurrent_pairs"]
    ):
        errors.append("primary sensitivity table row count mismatch")

    report = {
        "status": "PASS" if not errors else "FAIL",
        "version": VERSION,
        "checks": checks,
        "errors": errors,
    }
    (args.out_dir / "validation.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
