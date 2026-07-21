#!/usr/bin/env python3

"""Stream genotype fields from a het-SNV VCF and report Stage 1 QC."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    command = [
        "bcftools",
        "query",
        "-f",
        "[%GQ\\t%DP\\t%AD\\n]",
        str(args.vcf),
    ]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None

    count = 0
    depth_count = 0
    depth_sum = 0
    gq_count = 0
    gq_lt20 = 0
    gq_lt30 = 0
    ab_count = 0
    ab_outside = 0

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3:
            continue
        gq_text, dp_text, ad_text = fields
        count += 1

        if dp_text != ".":
            depth_sum += int(dp_text)
            depth_count += 1

        if gq_text != ".":
            gq = float(gq_text)
            gq_count += 1
            gq_lt20 += gq < 20
            gq_lt30 += gq < 30

        if ad_text not in {".", ""}:
            allele_depths = ad_text.split(",")
            if len(allele_depths) == 2 and "." not in allele_depths:
                ref_depth, alt_depth = map(int, allele_depths)
                total = ref_depth + alt_depth
                if total > 0:
                    allele_balance = alt_depth / total
                    ab_count += 1
                    ab_outside += allele_balance < 0.25 or allele_balance > 0.75

    return_code = process.wait()
    if return_code != 0:
        raise SystemExit(f"bcftools query failed with exit code {return_code}")
    if count == 0:
        raise SystemExit("No variants were found in the heterozygous-SNV VCF")

    def fraction(numerator: int, denominator: int) -> str:
        return f"{numerator / denominator:.6f}" if denominator else "NA"

    metrics = [
        ("sample", args.sample),
        ("heterozygous_snv_count", str(count)),
        ("mean_dp", f"{depth_sum / depth_count:.6f}" if depth_count else "NA"),
        ("gq_lt_20_fraction", fraction(gq_lt20, gq_count)),
        ("gq_lt_30_fraction", fraction(gq_lt30, gq_count)),
        ("ab_outside_0.25_0.75_fraction", fraction(ab_outside, ab_count)),
        ("gq_nonmissing_count", str(gq_count)),
        ("allele_balance_evaluable_count", str(ab_count)),
    ]

    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(metrics)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
