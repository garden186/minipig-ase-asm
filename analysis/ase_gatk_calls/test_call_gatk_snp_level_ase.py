#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
CALLER = HERE / "call_gatk_snp_level_ase.py"
VALIDATOR = HERE / "validate_outputs.py"
SPEC = importlib.util.spec_from_file_location("gatk_snp_caller", CALLER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


COUNT_HEADER = [
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
]


def write_tsv(path: Path, header: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle_context = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle_context = opener("r", encoding="utf-8", newline="")
    with handle_context as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_row(position: int, ref_count: int, alt_count: int) -> dict[str, object]:
    total = ref_count + alt_count
    return {
        "contig": "1",
        "position": position,
        "variantID": ".",
        "refAllele": "A",
        "altAllele": "G",
        "refCount": ref_count,
        "altCount": alt_count,
        "totalCount": total,
        "lowMAPQDepth": 0,
        "lowBaseQDepth": 0,
        "rawDepth": total,
        "otherBases": 0,
        "improperPairs": 0,
    }


class GatkSnpCallerTest(unittest.TestCase):
    def test_exact_binomial_matches_brute_force(self) -> None:
        for total in range(1, 81):
            for minor in range(total // 2 + 1):
                expected = min(
                    1.0,
                    2.0
                    * sum(math.comb(total, value) for value in range(minor + 1))
                    / (2**total),
                )
                observed = MODULE.exact_binom_two_sided_half_counts(total, minor)
                self.assertAlmostEqual(expected, observed, places=13)

    def test_bh_adjustment(self) -> None:
        observed = MODULE.bh_adjust([0.01, 0.04, 0.03])
        expected = [0.03, 0.04, 0.04]
        for left, right in zip(observed, expected):
            self.assertAlmostEqual(left, right)

    def test_wgs_and_rna_exclusions_are_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            count_dir = root / "counts"
            manifest = root / "manifest.tsv"
            excluded_animals = root / "exclude_animals.tsv"
            excluded_rna = root / "exclude_rna.tsv"
            out_dir = root / "out"
            manifest_rows = []
            units = [
                ("0001", "Brain"),
                ("0001", "Liver"),
                ("0002", "Liver"),
                ("0003", "Heart"),
            ]
            for sample, tissue in units:
                job = f"{sample}_{tissue}"
                manifest_rows.append(
                    {
                        "sample": sample,
                        "tissue": tissue,
                        "job_id": job,
                        "bam": f"/{job}.bam",
                        "vcf": f"/{sample}.vcf.gz",
                    }
                )
                write_tsv(
                    count_dir / f"{job}.ase_readcounter.tsv",
                    COUNT_HEADER,
                    [
                        count_row(100, 18, 2),
                        count_row(200, 10, 10),
                        count_row(300, 9, 1),
                    ],
                )
            write_tsv(
                manifest,
                ["sample", "tissue", "job_id", "bam", "vcf"],
                manifest_rows,
            )
            write_tsv(
                excluded_animals,
                ["sample", "reason"],
                [{"sample": "0002", "reason": "WGS QC"}],
            )
            write_tsv(
                excluded_rna,
                ["sample", "tissue", "reason"],
                [{"sample": "0001", "tissue": "Brain", "reason": "RNA QC"}],
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(CALLER),
                    "--manifest",
                    str(manifest),
                    "--count-dir",
                    str(count_dir),
                    "--out-dir",
                    str(out_dir),
                    "--profile-name",
                    "test_exclusions",
                    "--exclude-wgs-animals",
                    str(excluded_animals),
                    "--exclude-rna-units",
                    str(excluded_rna),
                    "--minimum-total-count",
                    "15",
                    "--minimum-major-fraction",
                    "0.65",
                    "--fdr",
                    "0.05",
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn('"status": "PASS"', result.stdout)
            effective = read_tsv(out_dir / "effective_cohort.tsv")
            self.assertEqual(
                {"0001_Liver", "0003_Heart"},
                {row["job_id"] for row in effective},
            )
            audit = {row["job_id"]: row for row in read_tsv(out_dir / "cohort_audit.tsv")}
            self.assertEqual("EXCLUDED_RNA_UNIT", audit["0001_Brain"]["inclusion_status"])
            self.assertEqual("EXCLUDED_WGS_ANIMAL", audit["0002_Liver"]["inclusion_status"])
            tests = read_tsv(out_dir / "site_tests.tsv.gz")
            calls = read_tsv(out_dir / "snp_level_ase_calls.tsv.gz")
            self.assertEqual(4, len(tests))
            self.assertEqual(2, len(calls))
            self.assertEqual({"0001_Liver", "0003_Heart"}, {row["job_id"] for row in calls})
            metadata = json.loads((out_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(4, metadata["cohort"]["manifest_units"])
            self.assertEqual(2, metadata["cohort"]["retained_units"])
            self.assertEqual(["0002"], metadata["cohort"]["excluded_wgs_animals"])
            self.assertEqual(4, metadata["counts"]["testable_sites"])
            validation = subprocess.run(
                [
                    sys.executable,
                    str(VALIDATOR),
                    "--out-dir",
                    str(out_dir),
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("[PASS]", validation.stdout)
            validation_json = json.loads(
                (out_dir / "validation.json").read_text(encoding="utf-8")
            )
            self.assertEqual("PASS", validation_json["status"])


if __name__ == "__main__":
    unittest.main()
