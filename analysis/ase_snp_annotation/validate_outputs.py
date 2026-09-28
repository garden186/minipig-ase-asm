#!/usr/bin/env python3
"""Validate v0.2.5-compatible ASE SNP feature annotation outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Dict, Iterator, Sequence


class ValidationError(RuntimeError):
    pass


def iter_tsv(path: Path) -> Iterator[Dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle = opener("r", encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {key: value if value is not None else "" for key, value in row.items()}


def validate(out_dir: Path) -> Dict[str, object]:
    required = [
        "all_tested_gene_tissue_snp_feature_annotations.tsv.gz",
        "recurrent_gene_tissue_snp_feature_annotations.tsv",
        "recurrent_snp_transcript_annotations.tsv.gz",
        "annotation_category_counts.tsv",
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

    observed = {
        "n_all_tested_annotation_rows": 0,
        "n_recurrent_annotation_rows": 0,
        "n_recurrent_transcript_annotation_rows": 0,
        "n_recurrent_unique_exonic": 0,
        "n_recurrent_ambiguous_shared_exon": 0,
    }
    recurrent_keys = set()
    for row in iter_tsv(out_dir / "all_tested_gene_tissue_snp_feature_annotations.tsv.gz"):
        observed["n_all_tested_annotation_rows"] += 1
        if row["annotation_status"] != "PASS":
            raise ValidationError(f"non-PASS tested annotation: {row}")
        if row["original_ase_assignment_class"] not in {"UNIQUE_EXONIC", "AMBIGUOUS_SHARED_EXON"}:
            raise ValidationError(f"unexpected tested assignment class: {row}")
        if int(row["target_gene_in_protein_coding_exon_union"]) != 1:
            raise ValidationError(f"target gene is not in original exon union: {row}")

    for row in iter_tsv(out_dir / "recurrent_gene_tissue_snp_feature_annotations.tsv"):
        observed["n_recurrent_annotation_rows"] += 1
        key = (row["gene_id"], row["tissue"], row["variant_id"])
        if key in recurrent_keys:
            raise ValidationError(f"duplicate recurrent annotation key: {key}")
        recurrent_keys.add(key)
        assignment = row["original_ase_assignment_class"]
        observed["n_recurrent_unique_exonic"] += int(assignment == "UNIQUE_EXONIC")
        observed["n_recurrent_ambiguous_shared_exon"] += int(assignment == "AMBIGUOUS_SHARED_EXON")

    transcript_hypotheses = set()
    for row in iter_tsv(out_dir / "recurrent_snp_transcript_annotations.tsv.gz"):
        observed["n_recurrent_transcript_annotation_rows"] += 1
        transcript_hypotheses.add((row["gene_id"], row["tissue"], row["variant_id"]))
        if row["transcript_feature"] not in {
            "CDS",
            "FIVE_PRIME_UTR",
            "THREE_PRIME_UTR",
            "START_CODON",
            "STOP_CODON",
            "EXON_NONCODING_OR_UNSPECIFIED",
            "INTRON",
            "OUTSIDE_TRANSCRIPT",
        }:
            raise ValidationError(f"unexpected transcript feature: {row}")
    if transcript_hypotheses != recurrent_keys:
        missing = sorted(recurrent_keys - transcript_hypotheses)[:10]
        raise ValidationError(f"recurrent hypotheses missing transcript annotations: {missing}")

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
            "all_tested_targets_reproduced_in_original_exon_union": True,
            "only_unique_or_shared_exonic_assignment_present": True,
            "all_recurrent_hypotheses_have_transcript_annotations": True,
            "transcript_features_use_original_v025_vocabulary": True,
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
    except (ValidationError, KeyError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(f"[PASS] validation: {result['observed_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
