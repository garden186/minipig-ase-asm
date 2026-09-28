#!/usr/bin/env python3
"""Audit exact PigGTEx gene-eQTL and fine-mapping support for all candidates.

The candidate table is expected to contain every mapped gene--tissue--SNP
context from the independent GATK recurrent SNP analysis. The two PigGTEx
archives are streamed without first restricting candidates to the canonical
haplotype-level recurrence set.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import csv
from collections import defaultdict
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
import time
from typing import Iterable


VARIANT_RE = re.compile(r"^(?:chr)?([^_]+)_(\d+)_([ACGTN]+)_([ACGTN]+)$", re.I)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-input", type=Path, required=True)
    parser.add_argument("--significant-eqtl", type=Path, required=True)
    parser.add_argument("--finemapped-eqtl", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--archive-workers",
        type=int,
        choices=(1, 2),
        default=2,
        help="Use 2 to scan the significant and fine-mapping archives concurrently.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def open_gz_member(tar: tarfile.TarFile, member: tarfile.TarInfo) -> io.TextIOWrapper:
    raw = tar.extractfile(member)
    if raw is None:
        raise FileNotFoundError(member.name)
    return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8-sig")


def clean_gene(value: str) -> str:
    return str(value).split(".", 1)[0]


def unordered_variant_key(value: str) -> str | None:
    match = VARIANT_RE.match(str(value))
    if match is None:
        return None
    chrom, position, allele_1, allele_2 = match.groups()
    first, second = sorted((allele_1.upper(), allele_2.upper()))
    return f"{chrom}_{int(position)}_{first}_{second}"


def infer_tissue(member_name: str, marker: str) -> str:
    return Path(member_name).name.split(marker, 1)[0]


def true_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_candidates(path: Path):
    rows: list[dict[str, str]] = []
    exact: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    pairs: set[tuple[str, str]] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "internal_gene_id",
            "internal_gene_name",
            "internal_tissue",
            "piggtex_gene_id",
            "piggtex_tissue",
            "internal_variant_unordered_key",
            "is_primary_tissue_enriched_pair",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Candidate input is missing fields: {sorted(missing)}")
        for index, source in enumerate(reader, start=1):
            row = dict(source)
            row["candidate_row_id"] = f"C{index:07d}"
            tissue = row["piggtex_tissue"]
            gene = clean_gene(row["piggtex_gene_id"])
            variant = row["internal_variant_unordered_key"]
            key = (tissue, gene, variant)
            rows.append(row)
            exact[key].append(index - 1)
            pairs.add((tissue, gene))
    return rows, exact, pairs


def default_evidence() -> dict[str, object]:
    return {
        "exact_significant_eqtl": 0,
        "minimum_exact_significant_eqtl_p": None,
        "significant_eqtl_slopes": set(),
        "significant_eqtl_slope_ses": set(),
        "exact_finemapped_eqtl": 0,
        "maximum_exact_finemapping_pip": None,
        "finemapping_credible_sets": set(),
    }


def stream_significant(
    path: Path,
    exact: dict[tuple[str, str, str], list[int]],
    pairs: set[tuple[str, str]],
    evidence: list[dict[str, object]],
):
    inventories: list[dict[str, object]] = []
    tissues = {tissue for tissue, _ in pairs}
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".gz"):
                continue
            tissue = infer_tissue(member.name, ".cis_qtl")
            used = tissue in tissues
            scanned = 0
            pair_rows = 0
            exact_rows = 0
            if used:
                with open_gz_member(archive, member) as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        scanned += 1
                        gene = clean_gene(row["phenotype_id"])
                        if (tissue, gene) not in pairs:
                            continue
                        pair_rows += 1
                        variant = unordered_variant_key(row["variant_id"])
                        if variant is None:
                            continue
                        indexes = exact.get((tissue, gene, variant), [])
                        if not indexes:
                            continue
                        exact_rows += 1
                        p_value = finite_float(row.get("pval_nominal"))
                        for index in indexes:
                            target = evidence[index]
                            target["exact_significant_eqtl"] = 1
                            current = target["minimum_exact_significant_eqtl_p"]
                            if p_value is not None and (current is None or p_value < current):
                                target["minimum_exact_significant_eqtl_p"] = p_value
                            if row.get("slope") not in {None, ""}:
                                target["significant_eqtl_slopes"].add(row["slope"])
                            if row.get("slope_se") not in {None, ""}:
                                target["significant_eqtl_slope_ses"].add(row["slope_se"])
            inventories.append(
                {
                    "archive": "significant_gene_eQTL",
                    "member": member.name,
                    "tissue": tissue,
                    "used": used,
                    "rows_scanned": scanned if used else None,
                    "candidate_pair_rows": pair_rows if used else None,
                    "exact_candidate_rows": exact_rows if used else None,
                }
            )
            if used:
                print(
                    f"[INFO] significant_eQTL tissue={tissue} "
                    f"rows={scanned} candidate_pair_rows={pair_rows} "
                    f"exact_rows={exact_rows}",
                    flush=True,
                )
    return inventories


def stream_finemapped(
    path: Path,
    exact: dict[tuple[str, str, str], list[int]],
    pairs: set[tuple[str, str]],
    evidence: list[dict[str, object]],
):
    inventories: list[dict[str, object]] = []
    tissues = {tissue for tissue, _ in pairs}
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".gz"):
                continue
            tissue = infer_tissue(member.name, ".susieinf")
            used = tissue in tissues
            scanned = 0
            pair_rows = 0
            exact_rows = 0
            if used:
                with open_gz_member(archive, member) as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        scanned += 1
                        gene = clean_gene(row["gene_id"])
                        if (tissue, gene) not in pairs:
                            continue
                        pair_rows += 1
                        variant = unordered_variant_key(row["variant_id"])
                        if variant is None:
                            continue
                        indexes = exact.get((tissue, gene, variant), [])
                        if not indexes:
                            continue
                        exact_rows += 1
                        pip = finite_float(row.get("prob"))
                        for index in indexes:
                            target = evidence[index]
                            target["exact_finemapped_eqtl"] = 1
                            current = target["maximum_exact_finemapping_pip"]
                            if pip is not None and (current is None or pip > current):
                                target["maximum_exact_finemapping_pip"] = pip
                            if row.get("cs") not in {None, ""}:
                                target["finemapping_credible_sets"].add(str(row["cs"]))
            inventories.append(
                {
                    "archive": "finemapped_gene_eQTL",
                    "member": member.name,
                    "tissue": tissue,
                    "used": used,
                    "rows_scanned": scanned if used else None,
                    "candidate_pair_rows": pair_rows if used else None,
                    "exact_candidate_rows": exact_rows if used else None,
                }
            )
            if used:
                print(
                    f"[INFO] finemapped_eQTL tissue={tissue} "
                    f"rows={scanned} candidate_pair_rows={pair_rows} "
                    f"exact_rows={exact_rows}",
                    flush=True,
                )
    return inventories


def run_significant_job(
    path: Path,
    exact: dict[tuple[str, str, str], list[int]],
    pairs: set[tuple[str, str]],
    row_count: int,
):
    evidence = [default_evidence() for _ in range(row_count)]
    started = time.monotonic()
    inventory = stream_significant(path, exact, pairs, evidence)
    print(
        f"[INFO] significant_eQTL archive complete elapsed_seconds={time.monotonic() - started:.1f}",
        flush=True,
    )
    return evidence, inventory


def run_finemapped_job(
    path: Path,
    exact: dict[tuple[str, str, str], list[int]],
    pairs: set[tuple[str, str]],
    row_count: int,
):
    evidence = [default_evidence() for _ in range(row_count)]
    started = time.monotonic()
    inventory = stream_finemapped(path, exact, pairs, evidence)
    print(
        f"[INFO] finemapped_eQTL archive complete elapsed_seconds={time.monotonic() - started:.1f}",
        flush=True,
    )
    return evidence, inventory


def merge_evidence(
    significant: list[dict[str, object]],
    finemapped: list[dict[str, object]],
) -> list[dict[str, object]]:
    if len(significant) != len(finemapped):
        raise ValueError("Evidence arrays have different lengths.")
    merged: list[dict[str, object]] = []
    for significant_row, finemapped_row in zip(significant, finemapped):
        row = default_evidence()
        for field in (
            "exact_significant_eqtl",
            "minimum_exact_significant_eqtl_p",
            "significant_eqtl_slopes",
            "significant_eqtl_slope_ses",
        ):
            row[field] = significant_row[field]
        for field in (
            "exact_finemapped_eqtl",
            "maximum_exact_finemapping_pip",
            "finemapping_credible_sets",
        ):
            row[field] = finemapped_row[field]
        merged.append(row)
    return merged


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    if path.suffix == ".gz":
        handle = gzip.open(path, "wt", encoding="utf-8", newline="")
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, object]], label: str) -> dict[str, object]:
    external_keys = {
        (row["piggtex_tissue"], row["piggtex_gene_id"], row["internal_variant_unordered_key"])
        for row in rows
    }
    pair_keys = {(row["internal_gene_id"], row["internal_tissue"]) for row in rows}
    return {
        "stratum": label,
        "candidate_rows": len(rows),
        "unique_external_gene_tissue_variant_keys": len(external_keys),
        "internal_gene_tissue_pairs": len(pair_keys),
        "rows_with_exact_significant_eqtl": sum(int(row["exact_significant_eqtl"]) for row in rows),
        "rows_with_exact_finemapped_eqtl": sum(int(row["exact_finemapped_eqtl"]) for row in rows),
        "rows_with_pip_ge_0_1": sum((row["maximum_exact_finemapping_pip"] or 0) >= 0.1 for row in rows),
        "rows_with_pip_ge_0_5": sum((row["maximum_exact_finemapping_pip"] or 0) >= 0.5 for row in rows),
        "rows_with_pip_ge_0_9": sum((row["maximum_exact_finemapping_pip"] or 0) >= 0.9 for row in rows),
    }


def main() -> None:
    args = parse_args()
    run_started = time.monotonic()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, exact, pairs = load_candidates(args.candidate_input)
    print(
        f"[INFO] candidates={len(rows)} exact_keys={len(exact)} "
        f"gene_tissue_pairs={len(pairs)} archive_workers={args.archive_workers}",
        flush=True,
    )
    if args.archive_workers == 2:
        with ProcessPoolExecutor(max_workers=2) as executor:
            significant_future = executor.submit(
                run_significant_job,
                args.significant_eqtl,
                exact,
                pairs,
                len(rows),
            )
            finemapped_future = executor.submit(
                run_finemapped_job,
                args.finemapped_eqtl,
                exact,
                pairs,
                len(rows),
            )
            significant_evidence, significant_inventory = significant_future.result()
            finemapped_evidence, finemapped_inventory = finemapped_future.result()
    else:
        significant_evidence, significant_inventory = run_significant_job(
            args.significant_eqtl, exact, pairs, len(rows)
        )
        finemapped_evidence, finemapped_inventory = run_finemapped_job(
            args.finemapped_eqtl, exact, pairs, len(rows)
        )
    evidence = merge_evidence(significant_evidence, finemapped_evidence)
    inventories = significant_inventory + finemapped_inventory

    output_rows: list[dict[str, object]] = []
    for source, observed in zip(rows, evidence):
        row: dict[str, object] = dict(source)
        row.update(observed)
        row["significant_eqtl_slopes"] = ",".join(sorted(observed["significant_eqtl_slopes"]))
        row["significant_eqtl_slope_ses"] = ",".join(sorted(observed["significant_eqtl_slope_ses"]))
        row["finemapping_credible_sets"] = ",".join(sorted(observed["finemapping_credible_sets"]))
        pip = observed["maximum_exact_finemapping_pip"]
        row["exact_finemapped_pip_ge_0_1"] = int((pip or 0) >= 0.1)
        row["exact_finemapped_pip_ge_0_5"] = int((pip or 0) >= 0.5)
        row["exact_finemapped_pip_ge_0_9"] = int((pip or 0) >= 0.9)
        output_rows.append(row)

    evidence_fields = [
        "exact_significant_eqtl",
        "minimum_exact_significant_eqtl_p",
        "significant_eqtl_slopes",
        "significant_eqtl_slope_ses",
        "exact_finemapped_eqtl",
        "maximum_exact_finemapping_pip",
        "finemapping_credible_sets",
        "exact_finemapped_pip_ge_0_1",
        "exact_finemapped_pip_ge_0_5",
        "exact_finemapped_pip_ge_0_9",
    ]
    output_fields = list(rows[0]) + evidence_fields if rows else evidence_fields
    write_tsv(args.out_dir / "full_candidate_gene_eqtl_finemap_audit.tsv.gz", output_rows, output_fields)
    exact_rows = [
        row for row in output_rows
        if row["exact_significant_eqtl"] or row["exact_finemapped_eqtl"]
    ]
    write_tsv(args.out_dir / "exact_gene_eqtl_or_finemap_matches.tsv.gz", exact_rows, output_fields)

    strata = [summarize(output_rows, "all_mapped_candidates")]
    primary = [row for row in output_rows if true_value(row["is_primary_tissue_enriched_pair"])]
    strata.append(summarize(primary, "primary_39_pair_candidates"))
    summary_fields = list(strata[0])
    write_tsv(args.out_dir / "gene_eqtl_finemap_stage_summary.tsv", strata, summary_fields)
    inventory_fields = list(inventories[0]) if inventories else []
    write_tsv(args.out_dir / "archive_member_audit.tsv", inventories, inventory_fields)

    errors: list[str] = []
    if len(output_rows) != len(rows):
        errors.append("Output row count differs from candidate input row count.")
    if len({row["candidate_row_id"] for row in output_rows}) != len(output_rows):
        errors.append("Candidate row IDs are not unique.")
    for row in output_rows:
        pip = row["maximum_exact_finemapping_pip"]
        if pip is not None and not 0 <= float(pip) <= 1:
            errors.append(f"PIP outside [0,1]: {row['candidate_row_id']}")
        if row["exact_finemapped_eqtl"] and pip is None:
            errors.append(f"Fine-mapped match lacks PIP: {row['candidate_row_id']}")
    validation = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "candidate_rows": len(rows),
        "unique_exact_external_keys": len(exact),
        "candidate_external_gene_tissue_pairs": len(pairs),
        "stage_summary": strata,
        "matching_rule": "same mapped tissue, release-100 Ensembl gene ID, position, and unordered allele pair",
        "input_sha256": sha256(args.candidate_input),
    }
    (args.out_dir / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    with ThreadPoolExecutor(max_workers=2) as executor:
        significant_hash_future = executor.submit(sha256, args.significant_eqtl)
        finemapped_hash_future = executor.submit(sha256, args.finemapped_eqtl)
        significant_archive_sha256 = significant_hash_future.result()
        finemapped_archive_sha256 = finemapped_hash_future.result()
    metadata = {
        "candidate_input": str(args.candidate_input.resolve()),
        "significant_gene_eqtl_archive": str(args.significant_eqtl.resolve()),
        "finemapped_gene_eqtl_archive": str(args.finemapped_eqtl.resolve()),
        "significant_gene_eqtl_archive_sha256": significant_archive_sha256,
        "finemapped_gene_eqtl_archive_sha256": finemapped_archive_sha256,
        "archive_workers": args.archive_workers,
        "elapsed_seconds": time.monotonic() - run_started,
        "archive_members": inventories,
    }
    (args.out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(validation, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
