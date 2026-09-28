#!/usr/bin/env python3
"""Opportunity-aware recurrence test for block-concordant SNP-level ASE.

The unit of inference is gene x tissue x SNP.  Eligibility is defined from the
complete MAIN/technical-PASS SNP-level audit, independently of whether the
parent block is core ASE.  A success is a site-level ASE call whose major
allele agrees with the major haplotype of a core parent block for that gene.

For every evaluable hypothesis (at least two eligible animals), the null uses
leave-one-gene-out sample x tissue burdens and an exact Poisson-binomial upper
tail.  One global BH correction is applied across all evaluable hypotheses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sqlite3
import sys
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple


VERSION = "ase_snp_recurrence_v0.1.0"

SITE_REQUIRED = {
    "sample",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "ref_fragment_count",
    "alt_fragment_count",
    "total_informative_fragments",
    "major_allele_fraction",
    "major_allele",
    "site_p_exact",
    "site_q_bh",
    "site_testable",
    "site_level_ase_signal",
    "core_concordant_site_level_ase_signal",
    "all_parent_gene_ids",
    "all_parent_gene_names",
}

CORE_SUMMARY_REQUIRED = {
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_samples_core_direction_concordant_site_ase",
    "samples_core_direction_concordant_site_ase",
}

MEASUREMENT_HEADER = [
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "ref_fragment_count",
    "alt_fragment_count",
    "total_informative_fragments",
    "major_allele_fraction",
    "major_allele",
    "site_p_exact",
    "site_q_bh",
    "site_level_ase_signal",
    "block_concordant_snp_level_ase",
]

RECURRENCE_HEADER = [
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_tissue_animals_total",
    "n_animals_testable",
    "animals_testable",
    "n_animals_with_block_concordant_snp_ase",
    "animals_with_block_concordant_snp_ase",
    "recurrence_fraction",
    "expected_animals_leave_one_gene_out",
    "null_probabilities_leave_one_gene_out",
    "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
]

SAMPLE_TISSUE_BURDEN_HEADER = [
    "sample",
    "tissue",
    "n_testable_gene_snp_measurements",
    "n_block_concordant_snp_ase_measurements",
    "block_concordant_snp_ase_probability_inclusive",
]

SAMPLE_TISSUE_GENE_BURDEN_HEADER = [
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "n_testable_snp_measurements_for_gene",
    "n_block_concordant_snp_ase_measurements_for_gene",
]


class AnalysisError(RuntimeError):
    """Raised when an input or an analysis invariant is violated."""


@contextmanager
def open_text(path: Path, mode: str = "rt"):
    if "b" in mode:
        raise ValueError("open_text supports text mode only")
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def iter_tsv(path: Path, required: Set[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise AnalysisError(f"{path}: missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value if value is not None else "" for key, value in row.items()}


def split_values(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_gene_id(value: str) -> str:
    return value.strip().split(".")[0]


def parse_int(value: str, field: str, context: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{context}: invalid integer in {field}: {value!r}") from error


def parse_bool(value: str, field: str, context: str) -> int:
    normalized = str(value).strip().upper()
    if normalized in {"1", "TRUE", "YES"}:
        return 1
    if normalized in {"0", "FALSE", "NO", ""}:
        return 0
    raise AnalysisError(f"{context}: invalid boolean in {field}: {value!r}")


def parse_optional_float(value: str, field: str, context: str) -> float | None:
    if value == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{context}: invalid float in {field}: {value!r}") from error
    if not math.isfinite(result):
        raise AnalysisError(f"{context}: non-finite float in {field}: {value!r}")
    return result


def fmt_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def poisson_binomial_tail(probabilities: Sequence[float], observed: int) -> float:
    """Return P(X >= observed) for independent Bernoulli probabilities."""
    if observed <= 0:
        return 1.0
    if observed > len(probabilities):
        return 0.0
    pmf = [1.0]
    for probability in probabilities:
        if not 0.0 <= probability <= 1.0:
            raise AnalysisError(f"invalid null probability: {probability}")
        updated = [0.0] * (len(pmf) + 1)
        for count, mass in enumerate(pmf):
            updated[count] += mass * (1.0 - probability)
            updated[count + 1] += mass * probability
        pmf = updated
    return min(1.0, max(0.0, math.fsum(pmf[observed:])))


def bh_values(values: Sequence[float]) -> List[float]:
    """Reference in-memory BH implementation used by unit tests."""
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 1.0
    for reverse_index in range(len(order) - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, min(1.0, values[index] * len(values) / rank))
        adjusted[index] = running
    return adjusted


def register_annotation(annotations: MutableMapping[str, str], gene_id: str, gene_name: str) -> None:
    gene_id = normalize_gene_id(gene_id)
    gene_name = gene_name.strip()
    if not gene_id or not gene_name:
        return
    current = annotations.get(gene_id, "")
    if current and current != gene_name:
        raise AnalysisError(
            f"conflicting gene names for {gene_id}: {current!r} versus {gene_name!r}"
        )
    annotations[gene_id] = gene_name


def load_success_keys(
    path: Path,
) -> Tuple[Set[Tuple[str, str, str, str]], Dict[str, str], Dict[str, int]]:
    success_keys: Set[Tuple[str, str, str, str]] = set()
    annotations: Dict[str, str] = {}
    counts = defaultdict(int)
    seen_hypotheses: Set[Tuple[str, str, str]] = set()

    for row in iter_tsv(path, CORE_SUMMARY_REQUIRED):
        gene_id = normalize_gene_id(row["gene_id"])
        tissue = row["tissue"].strip()
        variant_id = row["variant_id"].strip()
        context = f"{gene_id}|{tissue}|{variant_id}"
        if not gene_id or not tissue or not variant_id:
            raise AnalysisError(f"{path}: blank hypothesis key: {context}")
        hypothesis = (gene_id, tissue, variant_id)
        if hypothesis in seen_hypotheses:
            raise AnalysisError(f"{path}: duplicate hypothesis row: {context}")
        seen_hypotheses.add(hypothesis)
        register_annotation(annotations, gene_id, row["gene_name"])

        samples = split_values(row["samples_core_direction_concordant_site_ase"])
        n_reported = parse_int(
            row["n_samples_core_direction_concordant_site_ase"],
            "n_samples_core_direction_concordant_site_ase",
            context,
        )
        if len(samples) != n_reported or len(samples) != len(set(samples)):
            raise AnalysisError(
                f"{path}: concordant sample count/list mismatch for {context}: "
                f"reported={n_reported}, listed={samples}"
            )
        counts["n_core_summary_rows"] += 1
        counts["n_core_summary_rows_with_at_least_two_observed"] += int(n_reported >= 2)
        counts["n_core_summary_rows_with_at_least_three_observed"] += int(n_reported >= 3)
        for sample in samples:
            key = (sample, tissue, gene_id, variant_id)
            if key in success_keys:
                raise AnalysisError(f"{path}: duplicate success key: {key}")
            success_keys.add(key)

    counts["n_gene_linked_concordant_success_keys"] = len(success_keys)
    return success_keys, annotations, dict(counts)


def configure_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-200000")
    connection.executescript(
        """
        CREATE TABLE annotation (
            gene_id TEXT PRIMARY KEY,
            gene_name TEXT NOT NULL
        ) WITHOUT ROWID;

        CREATE TABLE measurement (
            sample TEXT NOT NULL,
            tissue TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            variant_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            ref_count INTEGER NOT NULL,
            alt_count INTEGER NOT NULL,
            total_count INTEGER NOT NULL,
            major_fraction REAL,
            major_allele TEXT NOT NULL,
            site_p REAL,
            site_q REAL,
            site_signal INTEGER NOT NULL,
            success INTEGER NOT NULL,
            PRIMARY KEY (sample, tissue, gene_id, variant_id)
        ) WITHOUT ROWID;
        """
    )
    return connection


def insert_measurement_batch(
    connection: sqlite3.Connection,
    batch: List[Tuple[object, ...]],
) -> None:
    if not batch:
        return
    try:
        connection.executemany(
            """
            INSERT INTO measurement VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            )
            """,
            batch,
        )
    except sqlite3.IntegrityError as error:
        raise AnalysisError(
            "duplicate sample-tissue-gene-SNP measurement encountered while loading site input"
        ) from error
    batch.clear()


def load_measurements(
    connection: sqlite3.Connection,
    path: Path,
    success_keys: Set[Tuple[str, str, str, str]],
    annotations: MutableMapping[str, str],
    batch_size: int,
) -> Tuple[Dict[str, int], Dict[str, Set[str]]]:
    counts = defaultdict(int)
    tissue_animals: Dict[str, Set[str]] = defaultdict(set)
    seen_success: Set[Tuple[str, str, str, str]] = set()
    batch: List[Tuple[object, ...]] = []

    connection.execute("BEGIN")
    try:
        for row in iter_tsv(path, SITE_REQUIRED):
            counts["n_physical_site_rows"] += 1
            sample = row["sample"].strip()
            tissue = row["tissue"].strip()
            variant_id = row["variant_id"].strip()
            context = f"{sample}|{tissue}|{variant_id}"
            if not sample or not tissue or not variant_id:
                raise AnalysisError(f"{path}: blank site key: {context}")
            tissue_animals[tissue].add(sample)

            physical_concordant = parse_bool(
                row["core_concordant_site_level_ase_signal"],
                "core_concordant_site_level_ase_signal",
                context,
            )
            counts["n_physical_block_concordant_site_observations"] += physical_concordant

            site_testable = parse_bool(row["site_testable"], "site_testable", context)
            if not site_testable:
                continue
            counts["n_testable_physical_site_rows"] += 1

            gene_ids = sorted({normalize_gene_id(value) for value in split_values(row["all_parent_gene_ids"])} - {""})
            if not gene_ids:
                raise AnalysisError(f"{path}: testable site has no parent gene: {context}")
            gene_names = sorted(set(split_values(row["all_parent_gene_names"])))
            if len(gene_ids) == 1 and len(gene_names) == 1:
                register_annotation(annotations, gene_ids[0], gene_names[0])

            chrom = row["chrom"].strip()
            position = parse_int(row["position"], "position", context)
            ref = row["ref"].strip().upper()
            alt = row["alt"].strip().upper()
            ref_count = parse_int(row["ref_fragment_count"], "ref_fragment_count", context)
            alt_count = parse_int(row["alt_fragment_count"], "alt_fragment_count", context)
            total_count = parse_int(
                row["total_informative_fragments"], "total_informative_fragments", context
            )
            if ref_count + alt_count != total_count:
                raise AnalysisError(f"{path}: REF+ALT != total for {context}")
            major_fraction = parse_optional_float(
                row["major_allele_fraction"], "major_allele_fraction", context
            )
            site_p = parse_optional_float(row["site_p_exact"], "site_p_exact", context)
            site_q = parse_optional_float(row["site_q_bh"], "site_q_bh", context)
            site_signal = parse_bool(row["site_level_ase_signal"], "site_level_ase_signal", context)

            for gene_id in gene_ids:
                key = (sample, tissue, gene_id, variant_id)
                success = int(key in success_keys)
                if success and not site_signal:
                    raise AnalysisError(
                        f"{path}: concordant success is not a site-level ASE signal: {key}"
                    )
                if success:
                    seen_success.add(key)
                batch.append(
                    (
                        sample,
                        tissue,
                        gene_id,
                        variant_id,
                        chrom,
                        position,
                        ref,
                        alt,
                        ref_count,
                        alt_count,
                        total_count,
                        major_fraction,
                        row["major_allele"].strip().upper(),
                        site_p,
                        site_q,
                        site_signal,
                        success,
                    )
                )
                counts["n_testable_gene_snp_measurements"] += 1
                counts["n_gene_linked_block_concordant_successes"] += success
                if len(batch) >= batch_size:
                    insert_measurement_batch(connection, batch)
        insert_measurement_batch(connection, batch)
        connection.commit()
    except Exception:
        connection.rollback()
        raise

    missing_success = success_keys - seen_success
    if missing_success:
        examples = sorted(missing_success)[:10]
        raise AnalysisError(
            f"{len(missing_success)} gene-linked concordant successes were not found as testable "
            f"measurements; examples={examples}"
        )
    if counts["n_gene_linked_block_concordant_successes"] != len(success_keys):
        raise AnalysisError(
            "gene-linked concordant success count does not equal the core summary success-key count"
        )
    return dict(counts), tissue_animals


def build_burdens_and_candidates(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE INDEX measurement_hypothesis_idx
        ON measurement (gene_id, tissue, variant_id, sample);

        CREATE INDEX measurement_sample_tissue_gene_idx
        ON measurement (sample, tissue, gene_id);

        CREATE TABLE sample_tissue_burden AS
        SELECT
            sample,
            tissue,
            COUNT(*) AS n_eligible,
            SUM(success) AS n_success
        FROM measurement
        GROUP BY sample, tissue;

        CREATE UNIQUE INDEX sample_tissue_burden_key
        ON sample_tissue_burden (sample, tissue);

        CREATE TABLE sample_tissue_gene_burden AS
        SELECT
            sample,
            tissue,
            gene_id,
            COUNT(*) AS n_eligible,
            SUM(success) AS n_success
        FROM measurement
        GROUP BY sample, tissue, gene_id;

        CREATE UNIQUE INDEX sample_tissue_gene_burden_key
        ON sample_tissue_gene_burden (sample, tissue, gene_id);

        CREATE TABLE candidate (
            hypothesis_id INTEGER PRIMARY KEY,
            gene_id TEXT NOT NULL,
            tissue TEXT NOT NULL,
            variant_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            n_testable INTEGER NOT NULL,
            n_success INTEGER NOT NULL
        );

        INSERT INTO candidate (
            gene_id, tissue, variant_id, chrom, position, ref, alt, n_testable, n_success
        )
        SELECT
            gene_id,
            tissue,
            variant_id,
            MIN(chrom),
            MIN(position),
            MIN(ref),
            MIN(alt),
            COUNT(*),
            SUM(success)
        FROM measurement
        GROUP BY gene_id, tissue, variant_id
        HAVING COUNT(*) >= 2
        ORDER BY gene_id, tissue, variant_id;

        CREATE UNIQUE INDEX candidate_key
        ON candidate (gene_id, tissue, variant_id);

        CREATE TABLE hypothesis (
            hypothesis_id INTEGER PRIMARY KEY,
            gene_id TEXT NOT NULL,
            gene_name TEXT NOT NULL,
            tissue TEXT NOT NULL,
            variant_id TEXT NOT NULL,
            chrom TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            n_tissue_animals_total INTEGER NOT NULL,
            n_testable INTEGER NOT NULL,
            samples_testable TEXT NOT NULL,
            n_success INTEGER NOT NULL,
            samples_success TEXT NOT NULL,
            recurrence_fraction REAL NOT NULL,
            expected_loo REAL NOT NULL,
            probabilities_loo TEXT NOT NULL,
            recurrence_p REAL NOT NULL,
            recurrence_q REAL,
            is_recurrent INTEGER NOT NULL DEFAULT 0
        );
        """
    )
    connection.commit()


def build_hypotheses(
    connection: sqlite3.Connection,
    tissue_animals: Mapping[str, Set[str]],
) -> Dict[str, int]:
    query = """
        SELECT
            c.hypothesis_id,
            c.gene_id,
            COALESCE(a.gene_name, ''),
            c.tissue,
            c.variant_id,
            c.chrom,
            c.position,
            c.ref,
            c.alt,
            c.n_testable,
            c.n_success,
            m.sample,
            m.success,
            st.n_eligible,
            st.n_success,
            stg.n_eligible,
            stg.n_success
        FROM candidate AS c
        JOIN measurement AS m
          ON m.gene_id = c.gene_id
         AND m.tissue = c.tissue
         AND m.variant_id = c.variant_id
        JOIN sample_tissue_burden AS st
          ON st.sample = m.sample
         AND st.tissue = m.tissue
        JOIN sample_tissue_gene_burden AS stg
          ON stg.sample = m.sample
         AND stg.tissue = m.tissue
         AND stg.gene_id = m.gene_id
        LEFT JOIN annotation AS a
          ON a.gene_id = c.gene_id
        ORDER BY c.hypothesis_id, m.sample
    """
    insert_sql = """
        INSERT INTO hypothesis VALUES (
            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL, 0
        )
    """
    counts = defaultdict(int)
    current_id: int | None = None
    metadata: Tuple[object, ...] | None = None
    samples: List[str] = []
    success_samples: List[str] = []
    probabilities: List[float] = []

    def flush() -> None:
        nonlocal current_id, metadata, samples, success_samples, probabilities
        if current_id is None or metadata is None:
            return
        (
            gene_id,
            gene_name,
            tissue,
            variant_id,
            chrom,
            position,
            ref,
            alt,
            n_testable,
            n_success,
        ) = metadata
        if len(samples) != int(n_testable) or len(success_samples) != int(n_success):
            raise AnalysisError(f"candidate aggregation mismatch for {gene_id}|{tissue}|{variant_id}")
        recurrence_p = poisson_binomial_tail(probabilities, int(n_success))
        expected = math.fsum(probabilities)
        tissue_total = len(tissue_animals.get(str(tissue), set()))
        if tissue_total < int(n_testable):
            raise AnalysisError(f"tissue animal count is smaller than testable count for {metadata}")
        connection.execute(
            insert_sql,
            (
                current_id,
                gene_id,
                gene_name,
                tissue,
                variant_id,
                chrom,
                position,
                ref,
                alt,
                tissue_total,
                n_testable,
                ",".join(samples),
                n_success,
                ",".join(success_samples),
                int(n_success) / int(n_testable),
                expected,
                ",".join(fmt_float(value) for value in probabilities),
                recurrence_p,
            ),
        )
        counts["n_evaluable_hypotheses"] += 1
        counts["n_hypotheses_with_at_least_two_observed"] += int(int(n_success) >= 2)
        counts["n_hypotheses_with_at_least_three_observed"] += int(int(n_success) >= 3)

    connection.execute("BEGIN")
    try:
        for row in connection.execute(query):
            hypothesis_id = int(row[0])
            if current_id != hypothesis_id:
                flush()
                current_id = hypothesis_id
                metadata = tuple(row[1:11])
                samples = []
                success_samples = []
                probabilities = []
            sample = str(row[11])
            success = int(row[12])
            st_eligible = int(row[13])
            st_success = int(row[14])
            stg_eligible = int(row[15])
            stg_success = int(row[16])
            denominator = st_eligible - stg_eligible
            numerator = st_success - stg_success
            if denominator <= 0:
                raise AnalysisError(
                    f"cannot compute leave-one-gene-out null for {sample}|{row[3]}|{row[1]}: "
                    f"denominator={denominator}"
                )
            if not 0 <= numerator <= denominator:
                raise AnalysisError(
                    f"invalid leave-one-gene-out burden for {sample}|{row[3]}|{row[1]}: "
                    f"{numerator}/{denominator}"
                )
            samples.append(sample)
            if success:
                success_samples.append(sample)
            probabilities.append(numerator / denominator)
        flush()
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return dict(counts)


def apply_global_bh(connection: sqlite3.Connection, fdr: float, batch_size: int = 50000) -> int:
    n_hypotheses = int(connection.execute("SELECT COUNT(*) FROM hypothesis").fetchone()[0])
    if n_hypotheses == 0:
        return 0
    rank = n_hypotheses
    running = 1.0
    updates: List[Tuple[float, int]] = []
    connection.execute("BEGIN")
    try:
        for hypothesis_id, p_value in connection.execute(
            "SELECT hypothesis_id, recurrence_p FROM hypothesis "
            "ORDER BY recurrence_p DESC, hypothesis_id DESC"
        ):
            running = min(running, min(1.0, float(p_value) * n_hypotheses / rank))
            updates.append((running, int(hypothesis_id)))
            rank -= 1
            if len(updates) >= batch_size:
                connection.executemany(
                    "UPDATE hypothesis SET recurrence_q = ? WHERE hypothesis_id = ?", updates
                )
                updates.clear()
        if updates:
            connection.executemany(
                "UPDATE hypothesis SET recurrence_q = ? WHERE hypothesis_id = ?", updates
            )
        connection.execute(
            """
            UPDATE hypothesis
            SET is_recurrent = CASE
                WHEN n_success >= 2 AND recurrence_q <= ? THEN 1 ELSE 0 END
            """,
            (fdr,),
        )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return int(connection.execute("SELECT SUM(is_recurrent) FROM hypothesis").fetchone()[0] or 0)


def write_annotations(connection: sqlite3.Connection, annotations: Mapping[str, str]) -> None:
    connection.executemany(
        "INSERT INTO annotation (gene_id, gene_name) VALUES (?, ?)",
        sorted(annotations.items()),
    )
    connection.commit()


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle = opener(path, "wt", encoding="utf-8", newline="")
    else:
        handle = opener("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def measurement_rows(connection: sqlite3.Connection) -> Iterator[Sequence[object]]:
    query = """
        SELECT
            m.sample, m.tissue, m.gene_id, COALESCE(a.gene_name, ''),
            m.variant_id, m.chrom, m.position, m.ref, m.alt,
            m.ref_count, m.alt_count, m.total_count, m.major_fraction,
            m.major_allele, m.site_p, m.site_q, m.site_signal, m.success
        FROM measurement AS m
        LEFT JOIN annotation AS a ON a.gene_id = m.gene_id
        ORDER BY m.sample, m.tissue, m.gene_id, m.chrom, m.position, m.variant_id
    """
    for row in connection.execute(query):
        mutable = list(row)
        for index in (12, 14, 15):
            mutable[index] = fmt_float(mutable[index])
        yield mutable


def recurrence_rows(connection: sqlite3.Connection, significant_only: bool) -> Iterator[Sequence[object]]:
    where = "WHERE is_recurrent = 1" if significant_only else ""
    order = (
        "recurrence_q, recurrence_p, n_success DESC, tissue, gene_id, chrom, position"
        if significant_only
        else "tissue, gene_id, chrom, position, variant_id"
    )
    query = f"""
        SELECT
            gene_id, gene_name, tissue, variant_id, chrom, position, ref, alt,
            n_tissue_animals_total, n_testable, samples_testable,
            n_success, samples_success, recurrence_fraction, expected_loo,
            probabilities_loo, recurrence_p, recurrence_q, is_recurrent
        FROM hypothesis
        {where}
        ORDER BY {order}
    """
    for row in connection.execute(query):
        mutable = list(row)
        for index in (13, 14, 16, 17):
            mutable[index] = fmt_float(mutable[index])
        yield mutable


def sample_tissue_burden_rows(connection: sqlite3.Connection) -> Iterator[Sequence[object]]:
    for sample, tissue, n_eligible, n_success in connection.execute(
        "SELECT sample, tissue, n_eligible, n_success FROM sample_tissue_burden "
        "ORDER BY sample, tissue"
    ):
        yield (
            sample,
            tissue,
            n_eligible,
            n_success,
            fmt_float(n_success / n_eligible),
        )


def sample_tissue_gene_burden_rows(connection: sqlite3.Connection) -> Iterator[Sequence[object]]:
    query = """
        SELECT b.sample, b.tissue, b.gene_id, COALESCE(a.gene_name, ''),
               b.n_eligible, b.n_success
        FROM sample_tissue_gene_burden AS b
        LEFT JOIN annotation AS a ON a.gene_id = b.gene_id
        ORDER BY b.sample, b.tissue, b.gene_id
    """
    yield from connection.execute(query)


def file_metadata(path: Path) -> Dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
    }


def prepare_output_directory(path: Path) -> None:
    if path.exists():
        entries = list(path.iterdir())
        if entries:
            raise AnalysisError(
                f"output directory is not empty: {path}; use a new directory to avoid overwriting"
            )
    else:
        path.mkdir(parents=True)


def run_analysis(args: argparse.Namespace) -> Dict[str, object]:
    if not args.site_input.is_file():
        raise AnalysisError(f"site input not found: {args.site_input}")
    if not args.core_summary.is_file():
        raise AnalysisError(f"core summary not found: {args.core_summary}")
    if not 0.0 < args.fdr < 1.0:
        raise AnalysisError(f"FDR must be between zero and one: {args.fdr}")
    if args.min_testable_animals != 2:
        raise AnalysisError(
            "v0.1.0 evaluates hypotheses that are testable in at least 2 animals"
        )
    if args.min_recurrent_animals != 2:
        raise AnalysisError(
            "v0.1.0 formally defines recurrence as at least 2 observed animals; "
            "do not replace the exact observed count with a >=3 category"
        )
    prepare_output_directory(args.out_dir)

    database_path = args.out_dir / "_snp_recurrence_work.sqlite3"
    success_keys, annotations, core_counts = load_success_keys(args.core_summary)
    connection = configure_database(database_path)
    try:
        measurement_counts, tissue_animals = load_measurements(
            connection,
            args.site_input,
            success_keys,
            annotations,
            args.batch_size,
        )
        write_annotations(connection, annotations)
        build_burdens_and_candidates(connection)
        hypothesis_counts = build_hypotheses(connection, tissue_animals)
        n_recurrent = apply_global_bh(connection, args.fdr)

        if (
            hypothesis_counts.get("n_hypotheses_with_at_least_two_observed", 0)
            != core_counts.get("n_core_summary_rows_with_at_least_two_observed", 0)
        ):
            raise AnalysisError(
                "pre-statistical >=2-animal hypothesis count does not reproduce the core summary"
            )
        if (
            hypothesis_counts.get("n_hypotheses_with_at_least_three_observed", 0)
            != core_counts.get("n_core_summary_rows_with_at_least_three_observed", 0)
        ):
            raise AnalysisError(
                "pre-statistical >=3-animal hypothesis count does not reproduce the core summary"
            )

        output_paths = {
            "eligible_measurements": args.out_dir / "eligible_gene_variant_measurements.tsv.gz",
            "sample_tissue_burden": args.out_dir / "sample_tissue_snp_burden.tsv",
            "sample_tissue_gene_burden": args.out_dir / "sample_tissue_gene_snp_burden.tsv.gz",
            "recurrence_all": args.out_dir / "snp_recurrence_all.tsv.gz",
            "statistically_recurrent": args.out_dir
            / "statistically_recurrent_block_concordant_snp_ase.tsv",
        }
        write_tsv(
            output_paths["eligible_measurements"],
            MEASUREMENT_HEADER,
            measurement_rows(connection),
        )
        write_tsv(
            output_paths["sample_tissue_burden"],
            SAMPLE_TISSUE_BURDEN_HEADER,
            sample_tissue_burden_rows(connection),
        )
        write_tsv(
            output_paths["sample_tissue_gene_burden"],
            SAMPLE_TISSUE_GENE_BURDEN_HEADER,
            sample_tissue_gene_burden_rows(connection),
        )
        write_tsv(
            output_paths["recurrence_all"],
            RECURRENCE_HEADER,
            recurrence_rows(connection, significant_only=False),
        )
        write_tsv(
            output_paths["statistically_recurrent"],
            RECURRENCE_HEADER,
            recurrence_rows(connection, significant_only=True),
        )

        counts = {
            **core_counts,
            **measurement_counts,
            **hypothesis_counts,
            "n_statistically_recurrent_gene_tissue_snp_hypotheses": n_recurrent,
            "n_statistically_recurrent_unique_physical_snps": int(
                connection.execute(
                    "SELECT COUNT(DISTINCT variant_id) FROM hypothesis WHERE is_recurrent = 1"
                ).fetchone()[0]
            ),
            "n_statistically_recurrent_genes": int(
                connection.execute(
                    "SELECT COUNT(DISTINCT gene_id) FROM hypothesis WHERE is_recurrent = 1"
                ).fetchone()[0]
            ),
        }
        qc = {
            "version": VERSION,
            "status": "PASS",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "counts": counts,
            "invariants": {
                "success_is_subset_of_testable_measurements": True,
                "gene_linked_success_count_matches_core_summary": True,
                "prestatistical_ge2_count_matches_core_summary": True,
                "prestatistical_ge3_count_matches_core_summary": True,
                "all_evaluable_hypotheses_have_at_least_two_testable_animals": True,
                "recurrence_uses_exact_counts_not_replication_categories": True,
            },
        }
        with (args.out_dir / "analysis_qc.json").open("w", encoding="utf-8") as handle:
            json.dump(qc, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")

        metadata = {
            "version": VERSION,
            "status": "PASS",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "unit": "gene-tissue-SNP",
            "definitions": {
                "testable_measurement": (
                    "sample-tissue-gene-SNP represented in a MAIN/technical-PASS parent block "
                    "and site-testable with at least 15 informative fragments in the upstream audit"
                ),
                "success": (
                    "site-level ASE (sample-tissue BH q<=0.05 and major-allele fraction>=0.65) "
                    "whose major allele agrees with a core ASE parent block for the same gene"
                ),
                "evaluable_hypothesis": (
                    f"the same gene-tissue-SNP is testable in at least {args.min_testable_animals} animals"
                ),
                "null_probability": (
                    "within each eligible sample-tissue, successful gene-SNP measurements divided "
                    "by testable gene-SNP measurements after removing all measurements of the target gene"
                ),
                "statistically_recurrent": (
                    f"observed in at least {args.min_recurrent_animals} animals and global BH q<={args.fdr}"
                ),
                "reporting": (
                    "exact n_animals_with_block_concordant_snp_ase / n_animals_testable; "
                    "the observed animal count is not converted to a replication category"
                ),
            },
            "statistics": {
                "tail_test": "exact Poisson-binomial upper tail",
                "background": "leave-one-gene-out sample-tissue burden",
                "multiple_testing": "one global Benjamini-Hochberg correction",
                "fdr": args.fdr,
            },
            "inputs": {
                "site_level_ase_measurements": file_metadata(args.site_input),
                "core_gene_variant_tissue_summary": file_metadata(args.core_summary),
            },
            "outputs": {key: str(value.resolve()) for key, value in output_paths.items()},
            "counts": counts,
            "software": {
                "python": sys.version,
                "sqlite": sqlite3.sqlite_version,
            },
        }
        with (args.out_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
    finally:
        connection.close()

    if not args.keep_database:
        for suffix in ("", "-wal", "-shm"):
            candidate = Path(str(database_path) + suffix)
            if candidate.exists():
                candidate.unlink()
    return metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-input", required=True, type=Path)
    parser.add_argument("--core-summary", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-testable-animals", type=int, default=2)
    parser.add_argument("--min-recurrent-animals", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument("--keep-database", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        metadata = run_analysis(args)
    except AnalysisError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 2
    print(
        "[PASS] statistically recurrent gene-tissue-SNP hypotheses: "
        f"{metadata['counts']['n_statistically_recurrent_gene_tissue_snp_hypotheses']}"
    )
    print(f"[PASS] output: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
