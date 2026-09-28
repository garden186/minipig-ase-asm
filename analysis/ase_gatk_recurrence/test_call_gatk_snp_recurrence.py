#!/usr/bin/env python3
"""Focused tests for independent GATK SNP-level ASE recurrence."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("call_gatk_snp_recurrence.py")
SPEC = importlib.util.spec_from_file_location("gatk_recurrence", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class GATKSNPRecurrenceTests(unittest.TestCase):
    def test_poisson_binomial_tail(self) -> None:
        self.assertAlmostEqual(MODULE.poisson_binomial_tail([0.5, 0.5], 2), 0.25)
        self.assertAlmostEqual(MODULE.poisson_binomial_tail([0.2, 0.3], 1), 0.44)
        self.assertEqual(MODULE.poisson_binomial_tail([0.2], 0), 1.0)

    def test_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cohort = root / "effective_cohort.tsv"
            sites = root / "site_tests.tsv.gz"
            output = root / "output"

            with cohort.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["sample", "tissue", "inclusion_status"])
                for sample in ("A", "B", "C"):
                    writer.writerow([sample, "Liver", "INCLUDED"])

            header = [
                "sample", "tissue", "contig", "position", "variant_id",
                "ref_allele", "alt_allele", "ref_count", "alt_count",
                "total_count", "major_allele_fraction", "major_allele",
                "imbalance_direction", "p_exact_two_sided",
                "q_bh_within_animal_tissue", "is_snp_level_ase",
            ]
            with gzip.open(sites, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for sample in ("A", "B", "C"):
                    writer.writerow(
                        [sample, "Liver", "1", 100, "1_100_A_G", "A", "G",
                         15, 0, 15, 1.0, "A", "REF", 0.0001, 0.001, 1]
                    )
                    for position in (200, 300, 400):
                        writer.writerow(
                            [sample, "Liver", "1", position, f"1_{position}_A_G",
                             "A", "G", 8, 7, 15, 8 / 15, "A", "REF", 1.0, 1.0, 0]
                        )

            args = MODULE.parse_args(
                [
                    "--site-tests", str(sites),
                    "--effective-cohort", str(cohort),
                    "--out-dir", str(output),
                ]
            )
            MODULE.run(args)
            metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["counts"]["n_testable_measurements"], 12)
            self.assertEqual(metadata["counts"]["n_evaluable_hypotheses"], 4)
            self.assertEqual(metadata["counts"]["n_statistically_recurrent_hypotheses"], 1)

            with (output / "statistically_recurrent_gatk_snp_ase.tsv").open(
                encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["physical_snp_id"], "1_100_A_G")
            self.assertEqual(rows[0]["n_animals_with_snp_level_ase"], "3")
            self.assertEqual(rows[0]["recurrence_p_poisson_binomial"], "0")


if __name__ == "__main__":
    unittest.main()
