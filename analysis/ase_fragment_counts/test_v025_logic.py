#!/usr/bin/env python3

import copy
import types
import unittest
from collections import defaultdict

import call_phaser_ase_candidates as caller
import wgs_snp_qc


def gene_count(a_by_snp, b_by_snp):
    record = caller.empty_strict_gene_count()
    for snp, qnames in a_by_snp.items():
        record["variant_qnames"]["A"][snp].update(qnames)
        record["strict"]["A"].update(qnames)
    for snp, qnames in b_by_snp.items():
        record["variant_qnames"]["B"][snp].update(qnames)
        record["strict"]["B"].update(qnames)
    return record


class V025LogicTest(unittest.TestCase):
    def setUp(self):
        self.args = types.SimpleNamespace(
            min_total=15,
            min_major_fraction=0.65,
            single_variant_driver_fraction=0.80,
        )

    def test_wgs_balance_is_snp_specific(self):
        requested = {
            "1_10_A_G": ("1", 10, "A", "G"),
            "1_20_C_T": ("1", 20, "C", "T"),
            "1_30_G_A": ("1", 30, "G", "A"),
        }
        vcf = {
            snp: {
                "filter": "PASS", "gt": "0|1", "gq": 50.0,
                "ps": "10", "ref": values[2], "alt": values[3],
            }
            for snp, values in requested.items()
        }
        counts = {
            "1_10_A_G": {
                "refCount": "20", "altCount": "20", "totalCount": "40"
            },
            "1_20_C_T": {
                "refCount": "34", "altCount": "6", "totalCount": "40"
            },
            "1_30_G_A": {
                "refCount": "5", "altCount": "4", "totalCount": "9"
            },
        }
        rows, balanced = wgs_snp_qc.evaluate_wgs_snps(
            requested, vcf, counts, 20, 15, 5, 0.30, 0.70
        )
        status = {row["variant_id"]: row["wgs_snp_status"] for row in rows}
        self.assertEqual(balanced, {"1_10_A_G"})
        self.assertEqual(status["1_10_A_G"], "WGS_BALANCED")
        self.assertEqual(status["1_20_C_T"], "WGS_UNBALANCED")
        self.assertEqual(status["1_30_G_A"], "WGS_INSUFFICIENT")

    def test_leave_one_snp_out_support_classes(self):
        robust = gene_count(
            {
                "1_10_A_G": {"a{}".format(i) for i in range(12)},
                "1_20_C_T": {"c{}".format(i) for i in range(12)},
                "1_30_G_A": {"e{}".format(i) for i in range(12)},
            },
            {
                "1_10_A_G": {"b{}".format(i) for i in range(3)},
                "1_20_C_T": {"d{}".format(i) for i in range(3)},
                "1_30_G_A": {"f{}".format(i) for i in range(3)},
            },
        )
        metrics = caller.calculate_strict_variant_metrics(
            robust, self.args
        )
        self.assertEqual(metrics["support_class"], "MULTI_SNP_LOO_ROBUST")
        self.assertTrue(metrics["loo_all_robust"])

        dominated = gene_count(
            {
                "1_10_A_G": {"a{}".format(i) for i in range(20)},
                "1_20_C_T": {"shared"},
            },
            {
                "1_10_A_G": {"b0", "b1", "b2"},
                "1_20_C_T": set(),
            },
        )
        metrics = caller.calculate_strict_variant_metrics(
            dominated, self.args
        )
        self.assertEqual(metrics["support_class"], "SINGLE_SNP_DOMINANT")

    def test_primary_class_is_independent_of_wgs_phase_orientation(self):
        base = {
            "analysis_set": "MAIN",
            "count_validation": "MATCH",
            "ase_test_status": "TESTED",
            "ase_q_bh": 0.01,
            "major_haplotype_fraction": 0.75,
            "fragment_hap_conflict_fraction": 0.0,
            "minor_support_class": "ADEQUATE_MINOR",
            "n_wgs_balanced_gene_snps": 3,
            "snp_support_status": "MULTI_SNP_LOO_ROBUST",
            "local_phase_status": "PASS",
            "wgs_phase_status": "DISCORDANT",
        }
        thresholds = {
            "fdr": 0.05,
            "min_major_fraction": 0.65,
            "max_fragment_hap_conflict_fraction": 0.10,
        }
        candidate_class, reasons = caller.classify_candidate(
            base, None, thresholds
        )
        self.assertEqual(candidate_class, "PRIMARY_ROBUST")
        self.assertEqual(reasons, "")

        single = copy.deepcopy(base)
        single["snp_support_status"] = "SINGLE_SNP_ONLY"
        candidate_class, reasons = caller.classify_candidate(
            single, None, thresholds
        )
        self.assertEqual(candidate_class, "SECONDARY_SINGLE_SNP")
        self.assertIn("SINGLE_SNP_ONLY", reasons)

        extreme = copy.deepcopy(base)
        extreme["minor_support_class"] = "EXTREME_MINOR"
        candidate_class, reasons = caller.classify_candidate(
            extreme, None, thresholds
        )
        self.assertEqual(candidate_class, "NOT_PRIMARY")
        self.assertIn("EXTREME_MINOR", reasons)

    def test_local_phase_requires_connected_retained_edges(self):
        haplotypes = {
            "B1": {"edges_supporting": "2", "edges_total": "2"}
        }
        edges = defaultdict(list)
        edges["B1"] = [
            {
                "variant_a": "v1", "variant_b": "v2",
                "pair": caller.pair_key("v1", "v2"),
                "supporting_connections": 5, "total_connections": 5,
                "conflicting_reads": 0, "retained": True,
            },
            {
                "variant_a": "v2", "variant_b": "v3",
                "pair": caller.pair_key("v2", "v3"),
                "supporting_connections": 4, "total_connections": 4,
                "conflicting_reads": 0, "retained": True,
            },
        ]
        result = caller.assess_local_phase_for_snps(
            "B1", ["v1", "v2", "v3"], edges, haplotypes, 3, 0.10
        )
        self.assertEqual(result["local_phase_status"], "PASS")

        result = caller.assess_local_phase_for_snps(
            "B1", ["v1", "v3"], edges, haplotypes, 3, 0.10
        )
        self.assertEqual(result["local_phase_status"], "FAIL")
        self.assertEqual(
            result["local_phase_reason"],
            "WGS_FILTERED_SNP_GRAPH_DISCONNECTED",
        )


if __name__ == "__main__":
    unittest.main()
