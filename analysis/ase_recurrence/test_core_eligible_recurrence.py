#!/usr/bin/env python3
"""Small deterministic tests for the recurrence implementation."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("core_eligible_recurrence.py")
SPEC = importlib.util.spec_from_file_location("core_eligible_recurrence", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def example_row(**updates: str):
    row = {
        "sample": "A",
        "tissue": "Brain",
        "gene_id": "GENE1",
        "gene_name": "Gene1",
        "gene_biotype": "protein_coding",
        "haplotype_block_id": "BLOCK1",
        "measurement_id": "MEAS1",
        "analysis_set": "MAIN",
        "count_validation": "MATCH",
        "n_wgs_balanced_gene_snps": "2",
        "fragment_hap_conflict_fraction": "0.01",
        "ase_a_count": "12",
        "ase_b_count": "8",
        "ase_total_count": "20",
        "major_haplotype_fraction": "0.6",
        "ase_p_exact": "0.5",
        "ase_q_bh": "0.5",
        "n_gene_variants_with_reads": "2",
        "snp_support_status": "SINGLE_SNP_DOMINANT",
        "local_phase_status": "PASS",
    }
    row.update(updates)
    return row


class RecurrenceTests(unittest.TestCase):
    def test_bh_values(self):
        observed = MODULE.bh_values([0.01, 0.04, 0.03])
        self.assertEqual([round(value, 8) for value in observed], [0.03, 0.04, 0.04])

    def test_poisson_binomial_tail(self):
        self.assertAlmostEqual(MODULE.poisson_binomial_tail([0.1, 0.1], 1), 0.19)
        self.assertAlmostEqual(MODULE.poisson_binomial_tail([0.1, 0.1], 2), 0.01)
        self.assertEqual(MODULE.poisson_binomial_tail([0.1, 0.1], 0), 1.0)

    def test_eligibility_is_outcome_independent(self):
        row = example_row()
        # Eligibility must not depend on significance, imbalance, or the LOO
        # support class.  This row is eligible but is not a core ASE call.
        self.assertTrue(MODULE.eligibility_pass(row, threshold=15, max_conflict=0.10))
        self.assertFalse(
            MODULE.core_pass(
                row,
                q_value=0.5,
                threshold=15,
                max_conflict=0.10,
                fdr=0.05,
                min_major_fraction=0.65,
            )
        )

    def test_eligibility_requires_two_contributing_snps(self):
        row = example_row(n_gene_variants_with_reads="1")
        self.assertFalse(MODULE.eligibility_pass(row, threshold=15, max_conflict=0.10))

    def test_core_is_nested_within_eligibility(self):
        row = example_row(
            major_haplotype_fraction="0.8",
            snp_support_status="MULTI_SNP_SUPPORTED",
        )
        self.assertTrue(
            MODULE.core_pass(
                row,
                q_value=0.01,
                threshold=15,
                max_conflict=0.10,
                fdr=0.05,
                min_major_fraction=0.65,
            )
        )
        low_coverage = dict(row, ase_total_count="14")
        self.assertFalse(
            MODULE.core_pass(
                low_coverage,
                q_value=0.01,
                threshold=15,
                max_conflict=0.10,
                fdr=0.05,
                min_major_fraction=0.65,
            )
        )

    def test_bins_cover_boundaries(self):
        self.assertEqual(MODULE.block_bin(3), "3_PLUS")
        self.assertEqual(MODULE.coverage_bin(29), "15_29")
        self.assertEqual(MODULE.coverage_bin(30), "30_59")
        self.assertEqual(MODULE.coverage_bin(240), "240_PLUS")
        self.assertEqual(MODULE.snp_bin(5), "4_5")
        self.assertEqual(MODULE.snp_bin(6), "6_PLUS")

    def test_tissue_aliases_are_normalized(self):
        self.assertEqual(MODULE.normalize_tissue("Tenderlo"), "Tenderloin")
        self.assertEqual(MODULE.normalize_tissue("L-blood"), "Blood")

    def test_sample_inference_for_flat_and_nested_sources(self):
        flat = Path("/input/0235.ase_haplotype_block_results.tsv.gz")
        nested = Path("/input/0235/ase_haplotype_block_results.tsv")
        self.assertEqual(MODULE.infer_sample(flat), "0235")
        self.assertEqual(MODULE.infer_sample(nested), "0235")


if __name__ == "__main__":
    unittest.main(verbosity=2)
