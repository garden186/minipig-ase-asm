#!/usr/bin/env python3
"""Shared dependency-free I/O and frozen constants for ASM annotation v0.1.1."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

VERSION = "0.1.1"
NO_DIRECT_FEATURE = "NO_DIRECT_PROTEIN_CODING_ANNOTATION"
PRIMARY_FEATURE_ORDER = (
    "PROMOTER",
    "FIVE_PRIME_UTR",
    "CDS",
    "THREE_PRIME_UTR",
    "INTRON",
    "OTHER_GENE_BODY",
)
SUMMARY_FEATURE_ORDER = (*PRIMARY_FEATURE_ORDER, NO_DIRECT_FEATURE)
FEATURE_RANK = {feature: index for index, feature in enumerate(SUMMARY_FEATURE_ORDER)}


class AnnotationError(RuntimeError):
    """Raised when an input or output violates the frozen annotation contract."""


def normalize_chrom(value: str) -> str:
    return value.strip().removeprefix("chr")


def chrom_key(value: str) -> tuple[int, str]:
    normalized = normalize_chrom(value)
    return (int(normalized), "") if normalized.isdigit() else (10**9, normalized)


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_signature(path: Path, *, include_sha256: bool = True) -> dict[str, object]:
    if not path.is_file():
        raise AnnotationError(f"required input does not exist: {path}")
    stat = path.stat()
    result: dict[str, object] = {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_sha256:
        result["sha256"] = sha256_file(path)
    return result


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    try:
        partial.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(partial, path)
    finally:
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


@contextmanager
def atomic_tsv_writer(path: Path, fieldnames: Sequence[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(partial, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            yield writer
        os.replace(partial, path)
    finally:
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with atomic_tsv_writer(path, fieldnames) as writer:
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_tsv(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def require_fields(fieldnames: Sequence[str] | None, required: Iterable[str], path: Path) -> None:
    missing = sorted(set(required).difference(fieldnames or ()))
    if missing:
        raise AnnotationError(f"{path}: missing fields: {','.join(missing)}")
