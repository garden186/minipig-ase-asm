#!/usr/bin/env python3

"""Report retention of primary heterozygous SNVs after the Stage 1 DP cutoff."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def vcf_sample_count(vcf: Path) -> int:
    result = subprocess.run(
        ["bcftools", "query", "-l", str(vcf)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return len([line for line in result.stdout.splitlines() if line])


def count_before_and_after_dp(
    vcf: Path, autosomes: list[str], minimum_dp: int
) -> tuple[int, int]:
    view = subprocess.Popen(
        [
            "bcftools",
            "view",
            "-r",
            ",".join(autosomes),
            "-m2",
            "-M2",
            "-v",
            "snps",
            "-f",
            "PASS",
            "-g",
            "het",
            str(vcf),
            "-Ou",
        ],
        stdout=subprocess.PIPE,
    )
    assert view.stdout is not None
    query = subprocess.Popen(
        ["bcftools", "query", "-f", "[%DP\\n]", "-"],
        stdin=view.stdout,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    view.stdout.close()
    assert query.stdout is not None

    before = 0
    after = 0
    for line in query.stdout:
        depth_text = line.strip()
        before += 1
        if depth_text != "." and int(depth_text) >= minimum_dp:
            after += 1

    query_return_code = query.wait()
    view_return_code = view.wait()
    if view_return_code != 0 or query_return_code != 0:
        raise SystemExit(
            "bcftools failed while counting pre-DP heterozygous SNVs "
            f"(view={view_return_code}, query={query_return_code})"
        )
    if before == 0:
        raise SystemExit(f"No PASS biallelic heterozygous SNVs found in {vcf}")
    return before, after


def indexed_record_count(vcf: Path) -> int:
    result = subprocess.run(
        ["bcftools", "index", "-n", str(vcf)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return int(result.stdout.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--deepvariant-vcf", required=True, type=Path)
    parser.add_argument("--filtered-vcf", required=True, type=Path)
    parser.add_argument("--autosomes", required=True, nargs="+")
    parser.add_argument("--minimum-dp", type=int, default=10)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if args.minimum_dp < 0:
        raise SystemExit("--minimum-dp must be non-negative")
    for path in (args.deepvariant_vcf, args.filtered_vcf):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Required non-empty VCF not found: {path}")
    if vcf_sample_count(args.deepvariant_vcf) != 1:
        raise SystemExit("Stage 1 DP-retention summary requires a single-sample VCF")

    before, after_from_raw = count_before_and_after_dp(
        args.deepvariant_vcf, args.autosomes, args.minimum_dp
    )
    after_from_filtered = indexed_record_count(args.filtered_vcf)
    if after_from_raw != after_from_filtered:
        raise SystemExit(
            "DP-retained count does not match the filtered VCF: "
            f"{after_from_raw} versus {after_from_filtered}"
        )

    rows = [
        ("sample", args.sample),
        ("minimum_dp", str(args.minimum_dp)),
        ("heterozygous_snv_before_dp", str(before)),
        ("heterozygous_snv_after_dp", str(after_from_filtered)),
        ("dp_retention_fraction", f"{after_from_filtered / before:.6f}"),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(rows)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
