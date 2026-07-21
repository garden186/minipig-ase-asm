#!/usr/bin/env python3

"""Combine concise Stage 1 QC metrics across all samples."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "sample",
    "mean_autosome_depth",
    "autosome_fraction_ge_10x",
    "primary_mapped_fraction",
    "properly_paired_fraction",
    "duplicate_fraction",
    "heterozygous_snv_before_dp10",
    "final_heterozygous_snv_count",
    "dp10_retention_fraction",
    "phased_snv_count",
    "phased_fraction",
    "multi_variant_blocks",
    "block_span_n50_bp",
    "variants_per_block_n50",
]


def read_samples(path: Path) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            sample = line.rstrip("\n").split("\t", 1)[0]
            if sample == "sample_id":
                continue
            if not sample:
                continue
            if sample in seen:
                raise ValueError(f"Duplicate sample in {path}: {sample}")
            seen.add(sample)
            samples.append(sample)
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def read_metric_value(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"metric", "value"}.issubset(rows[0]):
        raise ValueError(f"Expected metric/value TSV: {path}")
    return {row["metric"]: row["value"] for row in rows}


def read_all_phasing_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("chrom") == "ALL":
                return row
    raise ValueError(f"Missing ALL row in {path}")


def require_keys(data: dict[str, str], keys: list[str], path: Path) -> None:
    missing = [key for key in keys if key not in data or data[key] == ""]
    if missing:
        raise ValueError(f"Missing fields in {path}: {', '.join(missing)}")


def build_dp_summary(
    script_dir: Path,
    project_dir: Path,
    sample: str,
    autosomes: list[str],
    minimum_dp: int,
    output: Path,
) -> None:
    subprocess.run(
        [
            sys.executable,
            str(script_dir / "summarize_dp_retention.py"),
            "--sample",
            sample,
            "--deepvariant-vcf",
            str(project_dir / "results/wgs/variants" / f"{sample}.deepvariant.vcf.gz"),
            "--filtered-vcf",
            str(project_dir / "results/wgs/variants" / f"{sample}.het.snp.vcf.gz"),
            "--autosomes",
            *autosomes,
            "--minimum-dp",
            str(minimum_dp),
            "--output",
            str(output),
        ],
        check=True,
    )


def collect_sample(
    project_dir: Path,
    script_dir: Path,
    sample: str,
    autosomes: list[str],
    minimum_dp: int,
    build_missing_dp_qc: bool,
) -> dict[str, str]:
    qc_dir = project_dir / "results/qc/wgs"
    phasing_dir = project_dir / "results/phasing"
    alignment_path = qc_dir / f"{sample}.stage1_alignment_qc.tsv"
    het_path = qc_dir / f"{sample}.het_snp_qc.tsv"
    dp_path = qc_dir / f"{sample}.het_dp_retention_qc.tsv"
    phasing_path = phasing_dir / f"{sample}.phased.wgs.honest.summary.tsv"

    required = [alignment_path, het_path, phasing_path]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing required summary file(s): " + ", ".join(missing))
    if not dp_path.is_file() or dp_path.stat().st_size == 0:
        if not build_missing_dp_qc:
            raise FileNotFoundError(
                f"Missing {dp_path}; rerun with --build-missing-dp-qc"
            )
        build_dp_summary(
            script_dir, project_dir, sample, autosomes, minimum_dp, dp_path
        )

    alignment = read_metric_value(alignment_path)
    het = read_metric_value(het_path)
    dp = read_metric_value(dp_path)
    phasing = read_all_phasing_row(phasing_path)
    require_keys(
        alignment,
        [
            "mean_autosome_depth",
            "autosome_fraction_ge_10x",
            "primary_mapped_fraction",
            "properly_paired_fraction",
            "duplicate_fraction",
        ],
        alignment_path,
    )
    require_keys(het, ["heterozygous_snv_count"], het_path)
    require_keys(
        dp,
        [
            "minimum_dp",
            "heterozygous_snv_before_dp",
            "heterozygous_snv_after_dp",
            "dp_retention_fraction",
        ],
        dp_path,
    )
    require_keys(
        phasing,
        [
            "phased",
            "phased_fraction",
            "multi_variant_blocks",
            "block_span_n50_bp",
            "variants_per_block_n50",
            "variants",
        ],
        phasing_path,
    )

    if int(dp["minimum_dp"]) != minimum_dp:
        raise ValueError(f"Unexpected DP threshold for {sample}: {dp['minimum_dp']}")
    final_count = int(het["heterozygous_snv_count"])
    if final_count != int(dp["heterozygous_snv_after_dp"]):
        raise ValueError(f"Het-SNV and DP-retention counts differ for {sample}")
    if final_count != int(phasing["variants"]):
        raise ValueError(f"Het-SNV and phasing input counts differ for {sample}")

    return {
        "sample": sample,
        "mean_autosome_depth": alignment["mean_autosome_depth"],
        "autosome_fraction_ge_10x": alignment["autosome_fraction_ge_10x"],
        "primary_mapped_fraction": alignment["primary_mapped_fraction"],
        "properly_paired_fraction": alignment["properly_paired_fraction"],
        "duplicate_fraction": alignment["duplicate_fraction"],
        "heterozygous_snv_before_dp10": dp["heterozygous_snv_before_dp"],
        "final_heterozygous_snv_count": str(final_count),
        "dp10_retention_fraction": dp["dp_retention_fraction"],
        "phased_snv_count": phasing["phased"],
        "phased_fraction": phasing["phased_fraction"],
        "multi_variant_blocks": phasing["multi_variant_blocks"],
        "block_span_n50_bp": phasing["block_span_n50_bp"],
        "variants_per_block_n50": phasing["variants_per_block_n50"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--autosomes", nargs="+", default=[str(i) for i in range(1, 19)])
    parser.add_argument("--minimum-dp", type=int, default=10)
    parser.add_argument("--build-missing-dp-qc", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    samples = read_samples(args.sample_sheet)
    rows = [
        collect_sample(
            args.project_dir,
            script_dir,
            sample,
            args.autosomes,
            args.minimum_dp,
            args.build_missing_dp_qc,
        )
        for sample in samples
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.output)
    print(f"Wrote {len(rows)} samples to {args.output}")


if __name__ == "__main__":
    main()
