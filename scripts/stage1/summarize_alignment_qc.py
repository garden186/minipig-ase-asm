#!/usr/bin/env python3

"""Create an autosome-specific Stage 1 alignment QC table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, length = line.rstrip("\n").split("\t")[:2]
            lengths[chrom] = int(length)
    if not lengths:
        raise ValueError("No autosome lengths were found")
    return lengths


def weighted_depth(path: Path, lengths: dict[str, int]) -> float:
    means: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["chrom"] in lengths:
                means[row["chrom"]] = float(row["mean"])
    missing = sorted(set(lengths) - set(means))
    if missing:
        raise ValueError(f"Missing mosdepth summary rows: {', '.join(missing)}")
    total_length = sum(lengths.values())
    return sum(means[c] * lengths[c] for c in lengths) / total_length


def coverage_breadth(
    path: Path, lengths: dict[str, int], thresholds: tuple[int, ...]
) -> dict[int, float]:
    values: dict[tuple[str, int], float] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom, depth_text, fraction_text = fields[:3]
            if chrom in lengths:
                depth = int(depth_text)
                if depth in thresholds:
                    values[(chrom, depth)] = float(fraction_text)

    total_length = sum(lengths.values())
    result: dict[int, float] = {}
    for threshold in thresholds:
        missing = [c for c in lengths if (c, threshold) not in values]
        if missing:
            raise ValueError(
                f"Missing depth {threshold} rows for: {', '.join(sorted(missing))}"
            )
        result[threshold] = (
            sum(values[(c, threshold)] * lengths[c] for c in lengths) / total_length
        )
    return result


def flagstat_metrics(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    patterns = {
        "primary_mapped_fraction": re.compile(r"primary mapped \(([0-9.]+)%"),
        "properly_paired_fraction": re.compile(r"properly paired \(([0-9.]+)%"),
    }
    text = path.read_text(encoding="utf-8")
    for key, pattern in patterns.items():
        match = pattern.search(text)
        if not match:
            raise ValueError(f"Could not parse {key} from {path}")
        result[key] = f"{float(match.group(1)) / 100:.6f}"
    return result


def duplicate_metrics(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if line.startswith("LIBRARY\t") and index + 1 < len(lines):
            header = line.split("\t")
            values = lines[index + 1].split("\t")
            row = dict(zip(header, values))
            return {
                "duplicate_fraction": row["PERCENT_DUPLICATION"],
                "estimated_library_size": row.get("ESTIMATED_LIBRARY_SIZE", "NA"),
            }
    raise ValueError(f"Could not parse duplication metrics from {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mosdepth-summary", required=True, type=Path)
    parser.add_argument("--mosdepth-global-dist", required=True, type=Path)
    parser.add_argument("--autosome-lengths", required=True, type=Path)
    parser.add_argument("--flagstat", required=True, type=Path)
    parser.add_argument("--duplication-metrics", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    lengths = read_lengths(args.autosome_lengths)
    metrics: list[tuple[str, str]] = [
        ("sample", args.sample),
        (
            "mean_autosome_depth",
            f"{weighted_depth(args.mosdepth_summary, lengths):.6f}",
        ),
    ]
    breadth = coverage_breadth(
        args.mosdepth_global_dist, lengths, thresholds=(1, 10, 20, 30, 40)
    )
    metrics.extend(
        (f"autosome_fraction_ge_{depth}x", f"{fraction:.6f}")
        for depth, fraction in breadth.items()
    )
    metrics.extend(flagstat_metrics(args.flagstat).items())
    metrics.extend(duplicate_metrics(args.duplication_metrics).items())

    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(metrics)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
