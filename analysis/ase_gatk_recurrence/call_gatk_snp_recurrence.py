#!/usr/bin/env python3
"""Call opportunity-aware recurrence of independent GATK SNP-level ASE.

The inference unit is tissue x physical biallelic SNP. Only animals in which
the SNP passed the upstream depth threshold are eligible. The null probability
for an eligible animal is its animal-tissue SNP ASE burden after excluding the
target SNP itself. Exact Poisson-binomial upper-tail P values are corrected by
one global Benjamini-Hochberg procedure across all evaluable hypotheses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Set, Tuple


VERSION = "ase_gatk_snp_recurrence_v0.1.0"

SITE_REQUIRED = {
    "sample",
    "tissue",
    "contig",
    "position",
    "variant_id",
    "ref_allele",
    "alt_allele",
    "ref_count",
    "alt_count",
    "total_count",
    "major_allele_fraction",
    "major_allele",
    "imbalance_direction",
    "p_exact_two_sided",
    "q_bh_within_animal_tissue",
    "is_snp_level_ase",
}

COHORT_REQUIRED = {"sample", "tissue", "inclusion_status"}

RECURRENCE_HEADER = [
    "tissue",
    "physical_snp_id",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "n_tissue_animals_total",
    "n_animals_testable",
    "animals_testable",
    "n_animals_with_snp_level_ase",
    "animals_with_snp_level_ase",
    "recurrence_fraction",
    "expected_animals_leave_one_site_out",
    "null_probabilities_leave_one_site_out",
    "recurrence_p_poisson_binomial",
    "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase",
]

BURDEN_HEADER = [
    "sample",
    "tissue",
    "n_testable_physical_snps",
    "n_snp_level_ase_calls",
    "snp_level_ase_probability_inclusive",
]

PHYSICAL_SUMMARY_HEADER = [
    "physical_snp_id",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "n_testable_animal_tissue_measurements",
    "n_snp_level_ase_measurements",
    "n_testable_animals",
    "n_ase_positive_animals",
    "n_testable_tissues",
    "n_ase_positive_tissues",
]


class AnalysisError(RuntimeError):
    """Raised when an input or analysis invariant is violated."""


@contextmanager
def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def iter_tsv(path: Path, required: Set[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AnalysisError(f"{path}: missing columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise AnalysisError(f"{path}:{line_number}: malformed TSV row")
            yield {key: value or "" for key, value in row.items()}


def parse_int(value: str, field: str, context: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise AnalysisError(f"{context}: invalid integer {field}={value!r}") from error


def parse_float(value: str, field: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise AnalysisError(f"{context}: invalid float {field}={value!r}") from error
    if not math.isfinite(result):
        raise AnalysisError(f"{context}: non-finite {field}={value!r}")
    return result


def parse_bool(value: str, field: str, context: str) -> int:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return 1
    if normalized in {"0", "false", "no"}:
        return 0
    raise AnalysisError(f"{context}: invalid Boolean {field}={value!r}")


def fmt_float(value: object) -> str:
    if value is None:
        return ""
    return format(float(value), ".12g")


def physical_snp_id(contig: str, position: int, ref: str, alt: str) -> str:
    return f"{contig}_{position}_{ref}_{alt}"


def poisson_binomial_tail(probabilities: Sequence[float], observed: int) -> float:
    if observed <= 0:
        return 1.0
    if observed > len(probabilities):
        return 0.0
    distribution = [1.0] + [0.0] * len(probabilities)
    populated = 0
    for probability in probabilities:
        if probability < 0.0 or probability > 1.0 or not math.isfinite(probability):
            raise AnalysisError(f"invalid Poisson-binomial probability: {probability}")
        for count in range(populated + 1, 0, -1):
            distribution[count] = (
                distribution[count] * (1.0 - probability)
                + distribution[count - 1] * probability
            )
        distribution[0] *= 1.0 - probability
        populated += 1
    return min(1.0, max(0.0, math.fsum(distribution[observed:])))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_metadata(path: Path) -> Dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def load_effective_cohort(path: Path) -> Tuple[Set[Tuple[str, str]], Dict[str, Set[str]]]:
    units: Set[Tuple[str, str]] = set()
    tissue_animals: Dict[str, Set[str]] = defaultdict(set)
    for row in iter_tsv(path, COHORT_REQUIRED):
        if row["inclusion_status"].strip().upper() != "INCLUDED":
            continue
        sample = row["sample"].strip()
        tissue = row["tissue"].strip()
        key = (sample, tissue)
        if key in units:
            raise AnalysisError(f"duplicate effective cohort unit: {sample}|{tissue}")
        units.add(key)
        tissue_animals[tissue].add(sample)
    if not units:
        raise AnalysisError("effective cohort is empty")
    return units, dict(tissue_animals)


def create_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(str(path))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.executescript(
        """
        CREATE TABLE measurement (
            sample TEXT NOT NULL,
            tissue TEXT NOT NULL,
            contig TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            variant_id TEXT NOT NULL,
            ref_count INTEGER NOT NULL,
            alt_count INTEGER NOT NULL,
            total_count INTEGER NOT NULL,
            major_fraction REAL NOT NULL,
            major_allele TEXT NOT NULL,
            direction TEXT NOT NULL,
            site_p REAL NOT NULL,
            site_q REAL NOT NULL,
            success INTEGER NOT NULL,
            PRIMARY KEY (sample, tissue, contig, position, ref, alt)
        ) WITHOUT ROWID;
        """
    )
    return connection


def load_measurements(
    connection: sqlite3.Connection,
    path: Path,
    cohort_units: Set[Tuple[str, str]],
    batch_size: int,
) -> Dict[str, int]:
    counts = Counter()
    observed_units: Set[Tuple[str, str]] = set()
    batch: List[Tuple[object, ...]] = []
    connection.execute("BEGIN")
    try:
        for row in iter_tsv(path, SITE_REQUIRED):
            sample = row["sample"].strip()
            tissue = row["tissue"].strip()
            unit = (sample, tissue)
            context = f"{sample}|{tissue}|{row['variant_id']}"
            if unit not in cohort_units:
                raise AnalysisError(f"site row outside effective cohort: {context}")
            observed_units.add(unit)
            contig = row["contig"].strip()
            position = parse_int(row["position"], "position", context)
            ref = row["ref_allele"].strip().upper()
            alt = row["alt_allele"].strip().upper()
            if not contig or not ref or not alt or ref == alt:
                raise AnalysisError(f"invalid SNP identity: {context}")
            ref_count = parse_int(row["ref_count"], "ref_count", context)
            alt_count = parse_int(row["alt_count"], "alt_count", context)
            total_count = parse_int(row["total_count"], "total_count", context)
            if ref_count + alt_count != total_count:
                raise AnalysisError(f"REF+ALT != total: {context}")
            major_fraction = parse_float(
                row["major_allele_fraction"], "major_allele_fraction", context
            )
            site_p = parse_float(row["p_exact_two_sided"], "p_exact_two_sided", context)
            site_q = parse_float(
                row["q_bh_within_animal_tissue"],
                "q_bh_within_animal_tissue",
                context,
            )
            success = parse_bool(row["is_snp_level_ase"], "is_snp_level_ase", context)
            batch.append(
                (
                    sample,
                    tissue,
                    contig,
                    position,
                    ref,
                    alt,
                    row["variant_id"].strip(),
                    ref_count,
                    alt_count,
                    total_count,
                    major_fraction,
                    row["major_allele"].strip().upper(),
                    row["imbalance_direction"].strip().upper(),
                    site_p,
                    site_q,
                    success,
                )
            )
            counts["n_testable_measurements"] += 1
            counts["n_snp_level_ase_calls"] += success
            if len(batch) >= batch_size:
                connection.executemany(
                    "INSERT INTO measurement VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO measurement VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    missing_units = cohort_units - observed_units
    if missing_units:
        raise AnalysisError(f"effective units without site rows: {sorted(missing_units)[:10]}")
    counts["n_effective_units"] = len(cohort_units)
    return dict(counts)


def build_aggregates(connection: sqlite3.Connection, minimum_testable_animals: int) -> None:
    connection.executescript(
        f"""
        CREATE INDEX measurement_hypothesis_idx
        ON measurement (tissue, contig, position, ref, alt, sample);

        CREATE INDEX measurement_unit_idx
        ON measurement (sample, tissue);

        CREATE TABLE burden AS
        SELECT sample, tissue, COUNT(*) AS n_testable, SUM(success) AS n_success
        FROM measurement
        GROUP BY sample, tissue;

        CREATE UNIQUE INDEX burden_key ON burden (sample, tissue);

        CREATE TABLE candidate AS
        SELECT
            tissue, contig, position, ref, alt,
            COUNT(*) AS n_testable,
            SUM(success) AS n_success
        FROM measurement
        GROUP BY tissue, contig, position, ref, alt
        HAVING COUNT(*) >= {int(minimum_testable_animals)};

        CREATE UNIQUE INDEX candidate_key
        ON candidate (tissue, contig, position, ref, alt);

        CREATE TABLE hypothesis (
            tissue TEXT NOT NULL,
            contig TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            n_testable INTEGER NOT NULL,
            samples_testable TEXT NOT NULL,
            n_success INTEGER NOT NULL,
            samples_success TEXT NOT NULL,
            recurrence_fraction REAL NOT NULL,
            expected_loo REAL NOT NULL,
            probabilities_loo TEXT NOT NULL,
            recurrence_p REAL NOT NULL,
            recurrence_q REAL,
            is_recurrent INTEGER NOT NULL DEFAULT 0,
            PRIMARY KEY (tissue, contig, position, ref, alt)
        ) WITHOUT ROWID;
        """
    )
    connection.commit()


def build_hypotheses(connection: sqlite3.Connection) -> Dict[str, int]:
    query = """
        SELECT
            c.tissue, c.contig, c.position, c.ref, c.alt,
            c.n_testable, c.n_success,
            m.sample, m.success, b.n_testable, b.n_success
        FROM candidate AS c
        JOIN measurement AS m
          ON m.tissue=c.tissue AND m.contig=c.contig AND m.position=c.position
         AND m.ref=c.ref AND m.alt=c.alt
        JOIN burden AS b ON b.sample=m.sample AND b.tissue=m.tissue
        ORDER BY c.tissue, c.contig, c.position, c.ref, c.alt, m.sample
    """
    insert = "INSERT INTO hypothesis VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,NULL,0)"
    counts = Counter()
    current_key: Tuple[str, str, int, str, str] | None = None
    expected_counts: Tuple[int, int] | None = None
    samples: List[str] = []
    successes: List[str] = []
    probabilities: List[float] = []

    def flush() -> None:
        nonlocal current_key, expected_counts, samples, successes, probabilities
        if current_key is None or expected_counts is None:
            return
        n_testable, n_success = expected_counts
        if len(samples) != n_testable or len(successes) != n_success:
            raise AnalysisError(f"hypothesis aggregation mismatch: {current_key}")
        p_value = poisson_binomial_tail(probabilities, n_success)
        connection.execute(
            insert,
            (
                *current_key,
                n_testable,
                ",".join(samples),
                n_success,
                ",".join(successes),
                n_success / n_testable,
                math.fsum(probabilities),
                ",".join(fmt_float(value) for value in probabilities),
                p_value,
            ),
        )
        counts["n_evaluable_hypotheses"] += 1
        counts["n_hypotheses_observed_ge2"] += int(n_success >= 2)
        counts["n_hypotheses_observed_ge3"] += int(n_success >= 3)

    connection.execute("BEGIN")
    try:
        for row in connection.execute(query):
            key = (str(row[0]), str(row[1]), int(row[2]), str(row[3]), str(row[4]))
            if key != current_key:
                flush()
                current_key = key
                expected_counts = (int(row[5]), int(row[6]))
                samples = []
                successes = []
                probabilities = []
            sample = str(row[7])
            success = int(row[8])
            denominator = int(row[9]) - 1
            numerator = int(row[10]) - success
            if denominator <= 0 or not 0 <= numerator <= denominator:
                raise AnalysisError(
                    f"invalid leave-one-site-out burden for {sample}|{row[0]}: "
                    f"{numerator}/{denominator}"
                )
            samples.append(sample)
            if success:
                successes.append(sample)
            probabilities.append(numerator / denominator)
        flush()
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return dict(counts)


def apply_global_bh(
    connection: sqlite3.Connection,
    fdr: float,
    minimum_observed_animals: int,
) -> int:
    n = int(connection.execute("SELECT COUNT(*) FROM hypothesis").fetchone()[0])
    if n == 0:
        return 0
    rows = connection.execute(
        "SELECT tissue,contig,position,ref,alt,recurrence_p FROM hypothesis "
        "ORDER BY recurrence_p DESC,tissue DESC,contig DESC,position DESC,ref DESC,alt DESC"
    )
    rank = n
    running = 1.0
    updates: List[Tuple[float, str, str, int, str, str]] = []
    connection.execute("BEGIN")
    try:
        for tissue, contig, position, ref, alt, p_value in rows:
            running = min(running, min(1.0, float(p_value) * n / rank))
            updates.append((running, tissue, contig, position, ref, alt))
            rank -= 1
            if len(updates) >= 50000:
                connection.executemany(
                    "UPDATE hypothesis SET recurrence_q=? WHERE "
                    "tissue=? AND contig=? AND position=? AND ref=? AND alt=?",
                    updates,
                )
                updates.clear()
        if updates:
            connection.executemany(
                "UPDATE hypothesis SET recurrence_q=? WHERE "
                "tissue=? AND contig=? AND position=? AND ref=? AND alt=?",
                updates,
            )
        connection.execute(
            "UPDATE hypothesis SET is_recurrent=CASE WHEN n_success>=? "
            "AND recurrence_q<=? THEN 1 ELSE 0 END",
            (minimum_observed_animals, fdr),
        )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return int(connection.execute("SELECT SUM(is_recurrent) FROM hypothesis").fetchone()[0] or 0)


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with open_text(path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def recurrence_rows(
    connection: sqlite3.Connection,
    tissue_animals: Mapping[str, Set[str]],
    significant_only: bool,
) -> Iterator[Sequence[object]]:
    where = "WHERE is_recurrent=1" if significant_only else ""
    order = (
        "recurrence_q,recurrence_p,n_success DESC,tissue,contig,position,ref,alt"
        if significant_only
        else "tissue,contig,position,ref,alt"
    )
    query = f"""
        SELECT tissue,contig,position,ref,alt,n_testable,samples_testable,
               n_success,samples_success,recurrence_fraction,expected_loo,
               probabilities_loo,recurrence_p,recurrence_q,is_recurrent
        FROM hypothesis {where} ORDER BY {order}
    """
    for row in connection.execute(query):
        tissue, contig, position, ref, alt = row[:5]
        yield (
            tissue,
            physical_snp_id(contig, int(position), ref, alt),
            contig,
            position,
            ref,
            alt,
            len(tissue_animals[str(tissue)]),
            row[5],
            row[6],
            row[7],
            row[8],
            fmt_float(row[9]),
            fmt_float(row[10]),
            row[11],
            fmt_float(row[12]),
            fmt_float(row[13]),
            row[14],
        )


def burden_rows(connection: sqlite3.Connection) -> Iterator[Sequence[object]]:
    for sample, tissue, n_testable, n_success in connection.execute(
        "SELECT sample,tissue,n_testable,n_success FROM burden ORDER BY sample,tissue"
    ):
        yield sample, tissue, n_testable, n_success, fmt_float(n_success / n_testable)


def physical_summary_rows(connection: sqlite3.Connection) -> Iterator[Sequence[object]]:
    query = """
        SELECT contig,position,ref,alt,COUNT(*),SUM(success),
               COUNT(DISTINCT sample),COUNT(DISTINCT CASE WHEN success=1 THEN sample END),
               COUNT(DISTINCT tissue),COUNT(DISTINCT CASE WHEN success=1 THEN tissue END)
        FROM measurement
        GROUP BY contig,position,ref,alt
        ORDER BY contig,position,ref,alt
    """
    for row in connection.execute(query):
        yield physical_snp_id(row[0], int(row[1]), row[2], row[3]), *row


def prepare_output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise AnalysisError(f"output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def run(args: argparse.Namespace) -> None:
    site_tests = args.site_tests.resolve()
    effective_cohort = args.effective_cohort.resolve()
    output = args.out_dir.resolve()
    prepare_output_directory(output)

    cohort_units, tissue_animals = load_effective_cohort(effective_cohort)
    temporary_parent = args.tmp_dir.resolve() if args.tmp_dir else output
    temporary_parent.mkdir(parents=True, exist_ok=True)
    database_fd, database_name = tempfile.mkstemp(
        prefix="gatk_snp_recurrence_", suffix=".sqlite", dir=temporary_parent
    )
    os.close(database_fd)
    database = Path(database_name)
    connection: sqlite3.Connection | None = None
    started = datetime.now(timezone.utc)
    try:
        connection = create_database(database)
        counts = load_measurements(connection, site_tests, cohort_units, args.batch_size)
        build_aggregates(connection, args.minimum_testable_animals)
        counts.update(build_hypotheses(connection))
        counts["n_statistically_recurrent_hypotheses"] = apply_global_bh(
            connection, args.fdr, args.minimum_observed_animals
        )
        counts["n_unique_testable_physical_snps"] = int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT 1 FROM measurement GROUP BY contig,position,ref,alt)"
            ).fetchone()[0]
        )
        counts["n_unique_called_physical_snps"] = int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT 1 FROM measurement WHERE success=1 "
                "GROUP BY contig,position,ref,alt)"
            ).fetchone()[0]
        )

        all_path = output / "gatk_snp_recurrence_all.tsv.gz"
        recurrent_path = output / "statistically_recurrent_gatk_snp_ase.tsv"
        burden_path = output / "sample_tissue_snp_burden.tsv"
        physical_path = output / "physical_snp_summary.tsv.gz"
        write_tsv(all_path, RECURRENCE_HEADER, recurrence_rows(connection, tissue_animals, False))
        write_tsv(recurrent_path, RECURRENCE_HEADER, recurrence_rows(connection, tissue_animals, True))
        write_tsv(burden_path, BURDEN_HEADER, burden_rows(connection))
        write_tsv(physical_path, PHYSICAL_SUMMARY_HEADER, physical_summary_rows(connection))

        metadata = {
            "version": VERSION,
            "status": "PASS",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": (datetime.now(timezone.utc) - started).total_seconds(),
            "statistical_unit": "tissue x physical biallelic SNP",
            "eligibility": "SNP present in upstream GATK site_tests for an effective animal-tissue unit",
            "success": "upstream independent GATK SNP-level ASE call",
            "null_probability": "animal-tissue SNP ASE burden after removing the target physical SNP",
            "test": "exact one-sided Poisson-binomial upper tail",
            "multiple_testing": "one global Benjamini-Hochberg correction across all evaluable tissue-SNP hypotheses",
            "parameters": {
                "minimum_testable_animals": args.minimum_testable_animals,
                "minimum_observed_animals": args.minimum_observed_animals,
                "fdr": args.fdr,
            },
            "counts": counts,
            "inputs": {
                "site_tests": file_metadata(site_tests),
                "effective_cohort": file_metadata(effective_cohort),
            },
            "outputs": {},
        }
        for path in (all_path, recurrent_path, burden_path, physical_path):
            metadata["outputs"][path.name] = file_metadata(path)
        with (output / "run_metadata.json").open("w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")
    finally:
        if connection is not None:
            connection.close()
        if database.exists():
            if args.keep_database:
                shutil.move(str(database), str(output / "analysis.sqlite"))
            else:
                database.unlink()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-tests", type=Path, required=True)
    parser.add_argument("--effective-cohort", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tmp-dir", type=Path)
    parser.add_argument("--minimum-testable-animals", type=int, default=2)
    parser.add_argument("--minimum-observed-animals", type=int, default=2)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument("--keep-database", action="store_true")
    args = parser.parse_args(argv)
    if args.minimum_testable_animals < 2:
        parser.error("--minimum-testable-animals must be at least 2")
    if args.minimum_observed_animals < 2:
        parser.error("--minimum-observed-animals must be at least 2")
    if not 0 < args.fdr <= 1:
        parser.error("--fdr must be in (0,1]")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    try:
        run(parse_args(argv))
    except (AnalysisError, OSError, sqlite3.Error) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
