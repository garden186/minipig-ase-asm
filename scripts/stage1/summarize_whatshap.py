#!/usr/bin/env python3

"""Summarize WhatsHap stats and block-list TSVs by named columns."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def n50(values: list[int]) -> int:
    if not values:
        return 0
    total = sum(values)
    cumulative = 0
    for value in sorted(values, reverse=True):
        cumulative += value
        if cumulative * 2 >= total:
            return value
    return 0


def integer(row: dict[str, str], key: str) -> int:
    value = row.get(key, "0")
    return int(float(value)) if value else 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--blocks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--autosomes", required=True, nargs="+")
    args = parser.parse_args()

    stats_rows = read_tsv(args.stats)
    if not stats_rows:
        raise SystemExit("WhatsHap stats TSV contains no data rows")

    stats_by_chrom = {row["chromosome"]: row for row in stats_rows}
    blocks_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    all_blocks: list[tuple[int, int]] = []
    with args.blocks.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row["chromosome"]
            span = integer(row, "to") - integer(row, "from") + 1
            variants = integer(row, "variants")
            if span <= 0 or variants <= 0:
                raise SystemExit(f"Invalid block row in {args.blocks}: {row}")
            block = (variants, span)
            blocks_by_chrom[chrom].append(block)
            all_blocks.append(block)

    fieldnames = [
        "chrom",
        "variants",
        "phased",
        "unphased",
        "phased_fraction",
        "singletons",
        "multi_variant_blocks",
        "block_rows_including_singletons",
        "variants_per_block_n50",
        "block_span_n50_bp",
        "largest_block_variants",
        "largest_block_span_bp",
        "sum_block_spans_bp",
    ]

    output_rows: list[dict[str, str | int]] = []
    for chrom in [*args.autosomes, "ALL"]:
        if chrom not in stats_by_chrom:
            raise SystemExit(f"Missing chromosome {chrom} in {args.stats}")
        stats = stats_by_chrom[chrom]
        blocks = all_blocks if chrom == "ALL" else blocks_by_chrom.get(chrom, [])
        variants = integer(stats, "variants")
        phased = integer(stats, "phased")
        output_rows.append(
            {
                "chrom": chrom,
                "variants": variants,
                "phased": phased,
                "unphased": integer(stats, "unphased"),
                "phased_fraction": f"{phased / variants:.6f}" if variants else "0.000000",
                "singletons": integer(stats, "singletons"),
                "multi_variant_blocks": integer(stats, "blocks"),
                "block_rows_including_singletons": len(blocks),
                "variants_per_block_n50": n50([item[0] for item in blocks]),
                "block_span_n50_bp": n50([item[1] for item in blocks]),
                "largest_block_variants": max((item[0] for item in blocks), default=0),
                "largest_block_span_bp": max((item[1] for item in blocks), default=0),
                "sum_block_spans_bp": sum(item[1] for item in blocks),
            }
        )

    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
