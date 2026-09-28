#!/usr/bin/env python3
"""Shared contracts and statistics for reference-fixed exact CpG-K ASM."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO


VERSION = "0.3.3"
MAX_K = 5
COUNT_FIELDS = (
    "sample", "tissue", "chrom", "k", "window_id", "first_cpg_index",
    "last_cpg_index", "cpg_positions_1based", "start0", "end0", "span_bp",
    "ps",
    *(f"H1_nM{value}" for value in range(MAX_K + 1)),
    *(f"H2_nM{value}" for value in range(MAX_K + 1)),
    "H1_all", "H2_all",
)


class CpGKError(RuntimeError):
    pass


@dataclass(frozen=True)
class TestabilityRule:
    name: str
    minimum_total_um: int
    minimum_minor_um: int

    def passes(self, h1_um: int, h2_um: int) -> bool:
        return (
            h1_um + h2_um >= self.minimum_total_um
            and min(h1_um, h2_um) >= self.minimum_minor_um
        )


TESTABILITY_RULES = (
    TestabilityRule("M3", 6, 3),
    TestabilityRule("M4", 8, 4),
    TestabilityRule("M5", 10, 5),
    TestabilityRule("T8m3", 8, 3),
    TestabilityRule("T10m3", 10, 3),
    TestabilityRule("T10m4", 10, 4),
    TestabilityRule("T12m4", 12, 4),
)


def parse_positive_ints(text: str, *, label: str) -> list[int]:
    result = sorted({int(token.strip()) for token in text.split(",") if token.strip()})
    if not result or result[0] < 1:
        raise CpGKError(f"{label} must contain positive integers")
    return result


def parse_chromosomes(text: str) -> list[str]:
    result: list[str] = []
    for part in text.split(","):
        token = part.strip().removeprefix("chr")
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            result.extend(str(value) for value in range(int(left), int(right) + 1))
        else:
            result.append(token)
    if not result:
        raise CpGKError("no chromosomes were selected")
    return result


def read_cohort(path: Path) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("enabled", "1").strip().lower() in {"0", "false", "no"}:
                continue
            result.append((row["sample"].strip(), row["tissue"].strip()))
    if not result:
        raise CpGKError(f"no enabled cohort units in {path}")
    return result


def metadata_path(root: Path, sample: str, tissue: str, chrom: str) -> Path:
    directory = root / sample / tissue / "mxf020" / "chromosomes" / f"chr{chrom}"
    for name in (
        f"{sample}.{tissue}.chr{chrom}.discover.run_metadata.json",
        f"{sample}.{tissue}.chr{chrom}.run_metadata.json",
    ):
        candidate = directory / name
        if candidate.exists():
            return candidate
    raise CpGKError(f"source metadata not found below {directory}")


def metadata_input(metadata: Mapping, key: str) -> str:
    value = metadata.get("inputs", {}).get(key)
    if isinstance(value, Mapping):
        value = value.get("path")
    if not value:
        raise CpGKError(f"source metadata is missing inputs.{key}")
    return str(value)


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


@contextmanager
def open_text(path: Path, mode: str = "rt") -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(partial, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fieldnames})
        os.replace(partial, path)
    finally:
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


def write_json(path: Path, value: Mapping) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    try:
        with partial.open("wt", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(partial, path)
    finally:
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


def load_reference_cpgs(path: Path) -> tuple[str, list[int]]:
    chrom: str | None = None
    positions: list[int] = []
    with open_text(path) as handle:
        for expected, row in enumerate(csv.DictReader(handle, delimiter="\t")):
            row_chrom = row["chrom"].removeprefix("chr")
            if chrom is None:
                chrom = row_chrom
            elif row_chrom != chrom:
                raise CpGKError(f"multiple chromosomes in {path}")
            if int(row["cpg_index"]) != expected:
                raise CpGKError(f"non-consecutive CpG index in {path}")
            positions.append(int(row["cpg_pos_1based"]) - 1)
    if chrom is None or not positions:
        raise CpGKError(f"empty reference CpG index: {path}")
    return chrom, positions


def exact_window_id(chrom: str, positions0: Sequence[int], first: int, k: int) -> str:
    coords = ",".join(str(positions0[index] + 1) for index in range(first, first + k))
    return f"chr{chrom}:K{k}:{coords}"


def state_for_methylated_count(n_methylated: int, k: int, scheme: str) -> str:
    if not 0 <= n_methylated <= k:
        raise CpGKError("methylated count is outside the CpG-K window")
    if scheme == "strict":
        return "U" if n_methylated == 0 else "M" if n_methylated == k else "X"
    if scheme == "q25":
        fraction = n_methylated / k
        return "U" if fraction <= 0.25 else "M" if fraction >= 0.75 else "X"
    if scheme == "q35":
        fraction = n_methylated / k
        return "U" if fraction <= 0.35 else "M" if fraction >= 0.65 else "X"
    raise CpGKError(f"unknown state scheme {scheme!r}")


def collapse_state_counts(row: Mapping[str, object], scheme: str) -> dict[str, int]:
    k = int(row["k"])
    result = {f"{hap}_{state}": 0 for hap in ("H1", "H2") for state in "UXM"}
    for hap in ("H1", "H2"):
        for n_methylated in range(k + 1):
            state = state_for_methylated_count(n_methylated, k, scheme)
            result[f"{hap}_{state}"] += int(row[f"{hap}_nM{n_methylated}"])
    return result


def get_rule(name: str) -> TestabilityRule:
    for rule in TESTABILITY_RULES:
        if rule.name == name:
            return rule
    raise CpGKError(f"unknown testability rule {name!r}")


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    if not pvalues:
        return []
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in pvalues):
        raise CpGKError("BH received a non-finite or invalid p-value")
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    result = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank_index in range(total - 1, -1, -1):
        index = order[rank_index]
        running = min(running, pvalues[index] * total / (rank_index + 1))
        result[index] = min(1.0, running)
    return result


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    if min(a, b, c, d) < 0:
        raise CpGKError("Fisher counts must be non-negative")
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2
    if row1 == 0 or row2 == 0 or total == 0:
        return math.nan

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return -math.inf
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def probability(x: int) -> float:
        return math.exp(
            log_choose(col1, x)
            + log_choose(total - col1, row1 - x)
            - log_choose(total, row1)
        )

    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    observed = probability(a)
    tolerance = max(1e-15, observed * 1e-12)
    return min(
        1.0,
        sum(
            probability(value)
            for value in range(lower, upper + 1)
            if probability(value) <= observed + tolerance
        ),
    )


def methylated_fraction_difference(h1_u: int, h1_m: int, h2_u: int, h2_m: int) -> float:
    h1_total, h2_total = h1_u + h1_m, h2_u + h2_m
    if h1_total == 0 or h2_total == 0:
        return math.nan
    return h1_m / h1_total - h2_m / h2_total


def fmt(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.12g}"
    return value
