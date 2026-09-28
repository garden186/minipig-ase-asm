#!/usr/bin/env python3
"""Call independent SNP-level ASE from GATK ASEReadCounter tables.

The caller supports two explicit exclusion scopes. Excluding a WGS animal
removes every RNA-seq unit from that animal. Excluding an RNA unit removes
only the named animal-tissue unit. The effective cohort is frozen in every
output directory so downstream recurrence can reconstruct its tested family.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Iterator, Mapping, Sequence, Tuple


VERSION = "ase_snp_level_gatk_v0.1.0"
REQUIRED_COUNT_COLUMNS = {
    "contig",
    "position",
    "variantID",
    "refAllele",
    "altAllele",
    "refCount",
    "altCount",
    "totalCount",
    "lowMAPQDepth",
    "lowBaseQDepth",
    "rawDepth",
    "otherBases",
    "improperPairs",
}
AUTOSOMES = {str(value) for value in range(1, 19)}
TISSUE_NORMALIZATION = {
    "Tenderlo": "Tenderloin",
    "blood": "Blood",
    "L-blood": "Blood",
}

SITE_HEADER = [
    "sample",
    "tissue",
    "job_id",
    "contig",
    "position",
    "variant_id",
    "ref_allele",
    "alt_allele",
    "ref_count",
    "alt_count",
    "total_count",
    "minor_allele_count",
    "alt_fraction",
    "major_allele_fraction",
    "major_allele",
    "imbalance_direction",
    "p_exact_two_sided",
    "q_bh_within_animal_tissue",
    "bh_scope",
    "bh_family_size",
    "minimum_total_count",
    "minimum_major_allele_fraction",
    "fdr_threshold",
    "is_fdr_significant",
    "is_strong_effect",
    "is_snp_level_ase",
    "low_mapq_depth",
    "low_baseq_depth",
    "raw_depth",
    "other_bases",
    "improper_pairs",
    "count_table",
]

SUMMARY_HEADER = [
    "sample",
    "tissue",
    "job_id",
    "n_raw_count_rows",
    "n_testable_sites",
    "n_fdr_significant_sites",
    "n_strong_effect_sites",
    "n_snp_level_ase_calls",
    "count_table",
]

COHORT_HEADER = [
    "sample",
    "tissue",
    "job_id",
    "bam",
    "vcf",
    "count_table",
    "inclusion_status",
    "exclusion_reason",
]


class AnalysisError(RuntimeError):
    """Raised for invalid inputs or violated analysis invariants."""


@dataclass(frozen=True)
class Unit:
    sample: str
    tissue: str
    job_id: str
    bam: str
    vcf: str
    count_table: Path
    inclusion_status: str
    exclusion_reason: str


def normalize_sample(value: str) -> str:
    text = str(value).strip()
    return text.zfill(4) if text.isdigit() else text


def normalize_tissue(value: str) -> str:
    text = str(value).strip()
    return TISSUE_NORMALIZATION.get(text, text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def fmt_float(value: float) -> str:
    return f"{value:.12g}"


@lru_cache(maxsize=None)
def exact_binom_two_sided_half_counts(total: int, minor: int) -> float:
    """Return the exact two-sided p value for a symmetric binomial test.

    For p=0.5, the probability-ordering and doubled-smaller-tail definitions
    coincide. Summation starts at the observed minor count and moves toward
    zero, avoiding a 2**(-n) starting value that underflows at high depth.
    """

    if total <= 0 or minor < 0 or minor > total // 2:
        raise ValueError("invalid total/minor counts")
    log_probability = (
        math.lgamma(total + 1)
        - math.lgamma(minor + 1)
        - math.lgamma(total - minor + 1)
        - total * math.log(2.0)
    )
    if log_probability < math.log(sys.float_info.min):
        return sys.float_info.min
    probability = math.exp(log_probability)
    cumulative = probability
    current = probability
    index = minor
    while index > 0:
        current *= index / (total - index + 1)
        cumulative += current
        if cumulative >= 0.5:
            return 1.0
        index -= 1
    return max(sys.float_info.min, min(1.0, 2.0 * cumulative))


def exact_binom_two_sided_half(ref_count: int, alt_count: int) -> float:
    total = ref_count + alt_count
    return exact_binom_two_sided_half_counts(total, min(ref_count, alt_count))


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    ranked = sorted(enumerate(pvalues), key=lambda item: (item[1], item[0]))
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(ranked)
    for reverse_index in range(total - 1, -1, -1):
        original_index, pvalue = ranked[reverse_index]
        rank = reverse_index + 1
        running = min(running, pvalue * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def read_excluded_animals(path: Path | None) -> Dict[str, str]:
    if path is None:
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) != {"sample", "reason"}:
            raise AnalysisError(f"invalid WGS exclusion header: {path}")
        result: Dict[str, str] = {}
        for row in reader:
            sample = normalize_sample(row["sample"])
            if not sample or sample in result:
                raise AnalysisError(f"blank or duplicate WGS exclusion: {sample!r}")
            result[sample] = row["reason"].strip()
    return result


def read_excluded_rna_units(path: Path | None) -> Dict[Tuple[str, str], str]:
    if path is None:
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) != {"sample", "tissue", "reason"}:
            raise AnalysisError(f"invalid RNA exclusion header: {path}")
        result: Dict[Tuple[str, str], str] = {}
        for row in reader:
            key = (normalize_sample(row["sample"]), normalize_tissue(row["tissue"]))
            if not key[0] or not key[1] or key in result:
                raise AnalysisError(f"blank or duplicate RNA exclusion: {key}")
            result[key] = row["reason"].strip()
    return result


def build_cohort(
    manifest_path: Path,
    count_dir: Path,
    excluded_animals: Mapping[str, str],
    excluded_rna_units: Mapping[Tuple[str, str], str],
) -> list[Unit]:
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", "tissue", "job_id", "bam", "vcf"}
    if not rows or not required.issubset(rows[0]):
        raise AnalysisError(f"manifest missing required columns: {manifest_path}")

    manifest_animals = {normalize_sample(row["sample"]) for row in rows}
    manifest_units = {
        (normalize_sample(row["sample"]), normalize_tissue(row["tissue"])) for row in rows
    }
    unknown_animals = set(excluded_animals) - manifest_animals
    unknown_units = set(excluded_rna_units) - manifest_units
    if unknown_animals:
        raise AnalysisError(f"unknown WGS exclusions: {sorted(unknown_animals)}")
    if unknown_units:
        raise AnalysisError(f"unknown RNA exclusions: {sorted(unknown_units)}")

    units: list[Unit] = []
    seen_jobs: set[str] = set()
    for row in rows:
        sample = normalize_sample(row["sample"])
        tissue = normalize_tissue(row["tissue"])
        job_id = row["job_id"].strip()
        if job_id in seen_jobs:
            raise AnalysisError(f"duplicate manifest job: {job_id}")
        seen_jobs.add(job_id)
        count_table = count_dir / f"{job_id}.ase_readcounter.tsv"
        if not count_table.is_file() or count_table.stat().st_size == 0:
            raise AnalysisError(f"missing count table: {count_table}")
        if sample in excluded_animals:
            status = "EXCLUDED_WGS_ANIMAL"
            reason = excluded_animals[sample]
        elif (sample, tissue) in excluded_rna_units:
            status = "EXCLUDED_RNA_UNIT"
            reason = excluded_rna_units[(sample, tissue)]
        else:
            status = "INCLUDED"
            reason = ""
        units.append(
            Unit(
                sample=sample,
                tissue=tissue,
                job_id=job_id,
                bam=row["bam"],
                vcf=row["vcf"],
                count_table=count_table,
                inclusion_status=status,
                exclusion_reason=reason,
            )
        )
    if not any(unit.inclusion_status == "INCLUDED" for unit in units):
        raise AnalysisError("the effective cohort is empty")
    return sorted(units, key=lambda unit: (unit.sample, unit.tissue, unit.job_id))


def write_cohort(path: Path, units: Iterable[Unit], included_only: bool) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COHORT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for unit in units:
            if included_only and unit.inclusion_status != "INCLUDED":
                continue
            writer.writerow(
                {
                    "sample": unit.sample,
                    "tissue": unit.tissue,
                    "job_id": unit.job_id,
                    "bam": unit.bam,
                    "vcf": unit.vcf,
                    "count_table": str(unit.count_table.resolve()),
                    "inclusion_status": unit.inclusion_status,
                    "exclusion_reason": unit.exclusion_reason,
                }
            )


def iter_testable_rows(
    unit: Unit,
    minimum_total_count: int,
) -> tuple[int, list[dict[str, object]]]:
    raw_rows = 0
    rows: list[dict[str, object]] = []
    with unit.count_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COUNT_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise AnalysisError(f"{unit.count_table}: missing columns {sorted(missing)}")
        last_order = 0
        last_position = 0
        for row_number, row in enumerate(reader, start=2):
            raw_rows += 1
            contig = row["contig"]
            if contig not in AUTOSOMES:
                raise AnalysisError(f"non-autosomal count at {unit.count_table}:{row_number}")
            try:
                position = int(row["position"])
                ref_count = int(row["refCount"])
                alt_count = int(row["altCount"])
                total_count = int(row["totalCount"])
            except ValueError as error:
                raise AnalysisError(
                    f"invalid core count at {unit.count_table}:{row_number}"
                ) from error
            if ref_count + alt_count != total_count or total_count < 1:
                raise AnalysisError(f"count invariant failed at {unit.count_table}:{row_number}")
            order = int(contig)
            if order < last_order or (order == last_order and position <= last_position):
                raise AnalysisError(
                    f"duplicate or unsorted position at {unit.count_table}:{row_number}"
                )
            last_order = order
            last_position = position
            if total_count < minimum_total_count:
                continue
            ref = row["refAllele"].upper()
            alt = row["altAllele"].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                raise AnalysisError(f"non-biallelic SNV at {unit.count_table}:{row_number}")
            minor = min(ref_count, alt_count)
            alt_fraction = alt_count / total_count
            major_fraction = max(ref_count, alt_count) / total_count
            if ref_count > alt_count:
                major_allele = ref
                direction = "REF"
            elif alt_count > ref_count:
                major_allele = alt
                direction = "ALT"
            else:
                major_allele = ""
                direction = "TIE"
            rows.append(
                {
                    "sample": unit.sample,
                    "tissue": unit.tissue,
                    "job_id": unit.job_id,
                    "contig": contig,
                    "position": position,
                    "variant_id": f"{contig}_{position}_{ref}_{alt}",
                    "ref_allele": ref,
                    "alt_allele": alt,
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                    "total_count": total_count,
                    "minor_allele_count": minor,
                    "alt_fraction": alt_fraction,
                    "major_allele_fraction": major_fraction,
                    "major_allele": major_allele,
                    "imbalance_direction": direction,
                    "p_exact_two_sided": exact_binom_two_sided_half_counts(total_count, minor),
                    "low_mapq_depth": int(row["lowMAPQDepth"]),
                    "low_baseq_depth": int(row["lowBaseQDepth"]),
                    "raw_depth": int(row["rawDepth"]),
                    "other_bases": int(row["otherBases"]),
                    "improper_pairs": int(row["improperPairs"]),
                    "count_table": str(unit.count_table.resolve()),
                }
            )
    return raw_rows, rows


def serializable_site_row(
    row: Mapping[str, object],
    qvalue: float,
    family_size: int,
    minimum_total_count: int,
    minimum_major_fraction: float,
    fdr: float,
) -> dict[str, object]:
    pvalue = float(row["p_exact_two_sided"])
    major_fraction = float(row["major_allele_fraction"])
    fdr_significant = qvalue <= fdr
    strong_effect = major_fraction >= minimum_major_fraction
    output = {key: row.get(key, "") for key in SITE_HEADER}
    output.update(
        {
            "alt_fraction": fmt_float(float(row["alt_fraction"])),
            "major_allele_fraction": fmt_float(major_fraction),
            "p_exact_two_sided": fmt_float(pvalue),
            "q_bh_within_animal_tissue": fmt_float(qvalue),
            "bh_scope": f"{row['sample']}|{row['tissue']}",
            "bh_family_size": family_size,
            "minimum_total_count": minimum_total_count,
            "minimum_major_allele_fraction": fmt_float(minimum_major_fraction),
            "fdr_threshold": fmt_float(fdr),
            "is_fdr_significant": int(fdr_significant),
            "is_strong_effect": int(strong_effect),
            "is_snp_level_ase": int(fdr_significant and strong_effect),
        }
    )
    return output


def validate_args(args: argparse.Namespace) -> None:
    if args.minimum_total_count < 1:
        raise AnalysisError("--minimum-total-count must be >=1")
    if not 0.5 <= args.minimum_major_fraction <= 1.0:
        raise AnalysisError("--minimum-major-fraction must be in [0.5, 1]")
    if not 0.0 < args.fdr < 1.0:
        raise AnalysisError("--fdr must be in (0, 1)")
    if not args.profile_name.strip():
        raise AnalysisError("--profile-name must not be blank")
    for path in (args.manifest, args.count_dir):
        if not path.exists():
            raise AnalysisError(f"missing input: {path}")
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise AnalysisError(f"refusing to overwrite non-empty output: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)


def call(args: argparse.Namespace) -> dict[str, object]:
    validate_args(args)
    excluded_animals = read_excluded_animals(args.exclude_wgs_animals)
    excluded_rna_units = read_excluded_rna_units(args.exclude_rna_units)
    units = build_cohort(args.manifest, args.count_dir, excluded_animals, excluded_rna_units)
    included = [unit for unit in units if unit.inclusion_status == "INCLUDED"]

    cohort_audit_path = args.out_dir / "cohort_audit.tsv"
    effective_cohort_path = args.out_dir / "effective_cohort.tsv"
    test_path = args.out_dir / "site_tests.tsv.gz"
    call_path = args.out_dir / "snp_level_ase_calls.tsv.gz"
    summary_path = args.out_dir / "sample_tissue_summary.tsv"
    metadata_path = args.out_dir / "run_metadata.json"
    write_cohort(cohort_audit_path, units, included_only=False)
    write_cohort(effective_cohort_path, units, included_only=True)

    summaries: list[dict[str, object]] = []
    totals = defaultdict(int)
    with gzip.open(test_path, "wt", encoding="utf-8", newline="") as test_handle, gzip.open(
        call_path, "wt", encoding="utf-8", newline=""
    ) as call_handle:
        test_writer = csv.DictWriter(
            test_handle, fieldnames=SITE_HEADER, delimiter="\t", lineterminator="\n"
        )
        call_writer = csv.DictWriter(
            call_handle, fieldnames=SITE_HEADER, delimiter="\t", lineterminator="\n"
        )
        test_writer.writeheader()
        call_writer.writeheader()
        for index, unit in enumerate(included, start=1):
            print(
                f"[{index}/{len(included)}] {unit.job_id}: reading {unit.count_table}",
                flush=True,
            )
            raw_rows, rows = iter_testable_rows(unit, args.minimum_total_count)
            pvalues = [float(row["p_exact_two_sided"]) for row in rows]
            qvalues = bh_adjust(pvalues)
            n_fdr = 0
            n_effect = 0
            n_calls = 0
            for row, qvalue in zip(rows, qvalues):
                serialized = serializable_site_row(
                    row,
                    qvalue,
                    len(rows),
                    args.minimum_total_count,
                    args.minimum_major_fraction,
                    args.fdr,
                )
                test_writer.writerow(serialized)
                n_fdr += int(serialized["is_fdr_significant"])
                n_effect += int(serialized["is_strong_effect"])
                if serialized["is_snp_level_ase"]:
                    call_writer.writerow(serialized)
                    n_calls += 1
            summaries.append(
                {
                    "sample": unit.sample,
                    "tissue": unit.tissue,
                    "job_id": unit.job_id,
                    "n_raw_count_rows": raw_rows,
                    "n_testable_sites": len(rows),
                    "n_fdr_significant_sites": n_fdr,
                    "n_strong_effect_sites": n_effect,
                    "n_snp_level_ase_calls": n_calls,
                    "count_table": str(unit.count_table.resolve()),
                }
            )
            totals["raw_count_rows"] += raw_rows
            totals["testable_sites"] += len(rows)
            totals["fdr_significant_sites"] += n_fdr
            totals["strong_effect_sites"] += n_effect
            totals["snp_level_ase_calls"] += n_calls

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SUMMARY_HEADER, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summaries)

    metadata: dict[str, object] = {
        "status": "PASS",
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profile_name": args.profile_name,
        "analysis_scope": (
            "Independent autosomal SNP-level ASE calls from GATK ASEReadCounter "
            "REF/ALT fragment counts"
        ),
        "statistical_unit": "animal-tissue-physical_SNP",
        "null_hypothesis": "REF:ALT fragment probability = 0.5:0.5",
        "multiple_testing": (
            "Benjamini-Hochberg within each retained animal-tissue across all "
            "SNPs with totalCount at least the configured minimum"
        ),
        "call_definition": {
            "minimum_total_count": args.minimum_total_count,
            "exact_two_sided_binomial": True,
            "maximum_q_bh_within_animal_tissue": args.fdr,
            "minimum_major_allele_fraction": args.minimum_major_fraction,
        },
        "cohort_policy": {
            "excluded_wgs_animal": "remove every RNA unit for that animal",
            "excluded_rna_unit": "remove only the named animal-tissue unit",
            "downstream_requirement": (
                "recurrence and its multiple-testing family must be rebuilt from effective_cohort.tsv"
            ),
        },
        "inputs": {
            "manifest": str(args.manifest.resolve()),
            "manifest_sha256": sha256(args.manifest),
            "count_dir": str(args.count_dir.resolve()),
            "exclude_wgs_animals": (
                str(args.exclude_wgs_animals.resolve()) if args.exclude_wgs_animals else None
            ),
            "exclude_wgs_animals_sha256": (
                sha256(args.exclude_wgs_animals) if args.exclude_wgs_animals else None
            ),
            "exclude_rna_units": (
                str(args.exclude_rna_units.resolve()) if args.exclude_rna_units else None
            ),
            "exclude_rna_units_sha256": (
                sha256(args.exclude_rna_units) if args.exclude_rna_units else None
            ),
        },
        "cohort": {
            "manifest_units": len(units),
            "retained_units": len(included),
            "retained_animals": sorted({unit.sample for unit in included}),
            "excluded_wgs_animals": sorted(excluded_animals),
            "excluded_rna_units": [
                {"sample": sample, "tissue": tissue}
                for sample, tissue in sorted(excluded_rna_units)
            ],
        },
        "counts": dict(totals),
        "outputs": {
            "cohort_audit": str(cohort_audit_path.resolve()),
            "effective_cohort": str(effective_cohort_path.resolve()),
            "site_tests": str(test_path.resolve()),
            "snp_level_ase_calls": str(call_path.resolve()),
            "sample_tissue_summary": str(summary_path.resolve()),
        },
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--count-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--profile-name", required=True)
    parser.add_argument("--exclude-wgs-animals", type=Path)
    parser.add_argument("--exclude-rna-units", type=Path)
    parser.add_argument("--minimum-total-count", type=int, default=15)
    parser.add_argument("--minimum-major-fraction", type=float, default=0.65)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        metadata = call(args)
    except (OSError, KeyError, ValueError, AnalysisError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": metadata["status"],
                "version": metadata["version"],
                "profile_name": metadata["profile_name"],
                "cohort": metadata["cohort"],
                "counts": metadata["counts"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
