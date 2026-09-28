#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from call_snp_recurrence import bh_values, poisson_binomial_tail, run_analysis


SITE_HEADER = [
    "sample", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "ref_fragment_count", "alt_fragment_count", "total_informative_fragments",
    "major_allele_fraction", "major_allele", "site_p_exact", "site_q_bh",
    "site_testable", "site_level_ase_signal", "core_concordant_site_level_ase_signal",
    "all_parent_gene_ids", "all_parent_gene_names",
]

CORE_HEADER = [
    "gene_id", "gene_name", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "n_samples_core_direction_concordant_site_ase",
    "samples_core_direction_concordant_site_ase",
]


class StatisticsTests(unittest.TestCase):
    def test_poisson_binomial(self) -> None:
        self.assertAlmostEqual(poisson_binomial_tail([0.1, 0.2], 2), 0.02)
        self.assertAlmostEqual(poisson_binomial_tail([0.1, 0.2], 1), 0.28)

    def test_bh(self) -> None:
        observed = bh_values([0.01, 0.04, 0.03])
        expected = [0.03, 0.04, 0.04]
        for left, right in zip(observed, expected):
            self.assertAlmostEqual(left, right)


class IntegrationTest(unittest.TestCase):
    def test_exact_counts_are_reported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            site_path = root / "sites.tsv.gz"
            core_path = root / "core.tsv.gz"
            out_dir = root / "out"

            with gzip.open(site_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=SITE_HEADER, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for sample in ("S1", "S2", "S3"):
                    # Target G1-v1: testable in three animals and successful in exactly two.
                    writer.writerow(self.site_row(sample, "v1", "G1", "Gene1", sample != "S3"))
                    # Background opportunities outside G1 ensure a valid leave-one-gene-out null.
                    for index in range(1, 5):
                        writer.writerow(
                            self.site_row(sample, f"b{index}", "G2", "Gene2", False, position=100 + index)
                        )

            with gzip.open(core_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=CORE_HEADER, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(
                    {
                        "gene_id": "G1",
                        "gene_name": "Gene1",
                        "tissue": "Liver",
                        "variant_id": "v1",
                        "chrom": "1",
                        "position": "1",
                        "ref": "A",
                        "alt": "G",
                        "n_samples_core_direction_concordant_site_ase": "2",
                        "samples_core_direction_concordant_site_ase": "S1,S2",
                    }
                )

            args = argparse.Namespace(
                site_input=site_path,
                core_summary=core_path,
                out_dir=out_dir,
                fdr=0.05,
                min_testable_animals=2,
                min_recurrent_animals=2,
                batch_size=10,
                keep_database=False,
            )
            metadata = run_analysis(args)
            self.assertEqual(metadata["counts"]["n_hypotheses_with_at_least_two_observed"], 1)

            with gzip.open(out_dir / "snp_recurrence_all.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            target = next(row for row in rows if row["gene_id"] == "G1")
            self.assertEqual(target["n_animals_testable"], "3")
            self.assertEqual(target["n_animals_with_block_concordant_snp_ase"], "2")
            self.assertEqual(target["animals_with_block_concordant_snp_ase"], "S1,S2")
            self.assertEqual(target["recurrence_fraction"], "0.666666666667")

    @staticmethod
    def site_row(
        sample: str,
        variant: str,
        gene_id: str,
        gene_name: str,
        success: bool,
        position: int = 1,
    ) -> dict[str, str]:
        return {
            "sample": sample,
            "tissue": "Liver",
            "variant_id": variant,
            "chrom": "1",
            "position": str(position),
            "ref": "A",
            "alt": "G",
            "ref_fragment_count": "18" if success else "8",
            "alt_fragment_count": "2" if success else "12",
            "total_informative_fragments": "20",
            "major_allele_fraction": "0.9" if success else "0.6",
            "major_allele": "A" if success else "G",
            "site_p_exact": "0.0004" if success else "0.5",
            "site_q_bh": "0.004" if success else "1",
            "site_testable": "1",
            "site_level_ase_signal": "1" if success else "0",
            "core_concordant_site_level_ase_signal": "1" if success else "0",
            "all_parent_gene_ids": gene_id,
            "all_parent_gene_names": gene_name,
        }


if __name__ == "__main__":
    unittest.main()
