#!/usr/bin/env python3
"""Focused unit tests for the canonical SNP/PigGTEx evidence audit."""

from __future__ import annotations

import unittest

from build_evidence_audit import is_full_high_pip, normalize_tissue, parse_variant_id, support_signature, unordered_site


class EvidenceAuditTests(unittest.TestCase):
    def test_tissue_normalization(self):
        self.assertEqual(normalize_tissue("blood"), "Blood")
        self.assertEqual(normalize_tissue("L-blood"), "Blood")
        self.assertEqual(normalize_tissue("Tenderlo"), "Tenderloin")

    def test_unordered_site_key(self):
        self.assertEqual(unordered_site("chr9", 100, "T", "G"), ("9", 100, "G", "T"))
        self.assertEqual(parse_variant_id("9_100_T_G"), ("9", 100, "G", "T"))

    def test_full_high_pip_definition(self):
        row = {
            "canonical_recurrent_gene_tissue_snp": 1,
            "n_primary_block_direction_concordant_samples": 1,
            "piggtex_significant_site_rows_same_mapped_tissue": 3,
            "piggtex_exact_significant_eqtl": 1,
            "piggtex_exact_finemapped_eqtl_pip_ge_0_5": 1,
        }
        self.assertTrue(is_full_high_pip(row))
        row["n_primary_block_direction_concordant_samples"] = 0
        self.assertFalse(is_full_high_pip(row))

    def test_signature_keeps_independent_recurrence_distinct(self):
        row = {
            "canonical_recurrent_gene_tissue_snp": 1,
            "n_primary_block_direction_concordant_samples": 2,
            "independent_gatk_statistically_recurrent": 0,
            "piggtex_significant_site_rows_same_mapped_tissue": 1,
            "piggtex_exact_significant_eqtl": 1,
            "piggtex_exact_finemapped_eqtl_pip_ge_0_5": 1,
            "piggtex_exact_finemapped_eqtl_pip_ge_0_9": 1,
        }
        signature = support_signature(row)
        self.assertIn("CANONICAL_SNP_RECURRENCE", signature)
        self.assertIn("INDEPENDENT_GATK_DIRECTION_CONCORDANT", signature)
        self.assertNotIn("INDEPENDENT_GATK_RECURRENCE;", signature)


if __name__ == "__main__":
    unittest.main()
