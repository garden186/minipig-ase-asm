#!/usr/bin/env python3
import csv
import json
import os
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


PKG_DIR = Path(__file__).resolve().parent
if str(PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PKG_DIR))

import reclassify_ase_v026 as reclassifier
import validate_ase_v026_outputs as validator


SOURCE_FIELDS = sorted(reclassifier.REQUIRED_BLOCK_FIELDS)


def make_row(
    block,
    gene,
    a_count,
    b_count,
    qvalue,
    pvalue,
    support,
    phase,
    old_class="NOT_PRIMARY",
    old_reasons="",
    balanced=2,
    conflict=0.0,
    test_status="TESTED",
):
    total = a_count + b_count
    row = {field: "" for field in SOURCE_FIELDS}
    row.update(
        {
            "sample": "TEST1",
            "tissue": "Liver",
            "gene_id": gene,
            "gene_name": gene,
            "gene_biotype": "protein_coding",
            "haplotype_block_id": block,
            "measurement_id": "TEST1|Liver|{}".format(block),
            "analysis_set": "MAIN",
            "count_validation": "MATCH",
            "n_wgs_balanced_gene_snps": str(balanced),
            "fragment_hap_conflict_fraction": str(conflict),
            "ase_a_count": str(a_count),
            "ase_b_count": str(b_count),
            "ase_total_count": str(total),
            "major_haplotype_fraction": str(max(a_count, b_count) / float(total)),
            "ase_test_status": test_status,
            "ase_p_exact": "" if pvalue is None else str(pvalue),
            "ase_q_bh": "" if qvalue is None else str(qvalue),
            "snp_support_status": support,
            "local_phase_status": phase,
            "wgs_phase_status": "CONCORDANT",
            "wgs_phase_orientation": "DIRECT",
            "wgs_phase_set": "1000",
            "candidate_class": old_class,
            "is_primary_candidate": "1" if old_class.startswith("PRIMARY") else "0",
            "exclusion_reasons": old_reasons,
        }
    )
    return row


class ClassificationTests(unittest.TestCase):
    def setUp(self):
        self.thresholds = dict(reclassifier.DEFAULT_THRESHOLDS)

    def classify(self, row):
        return reclassifier.classify_row(row, self.thresholds, {})

    def test_monoallelic_multi_snp_is_core(self):
        row = make_row(
            "b1",
            "SGCE",
            20,
            0,
            0.01,
            1e-5,
            "MULTI_SNP_SUPPORTED",
            "PASS",
            old_reasons="EXTREME_MINOR",
        )
        result = self.classify(row)
        self.assertEqual(result["v026_imbalance_shape"], "MONOALLELIC_LIKE")
        self.assertEqual(result["v026_evidence_class"], "ASE_CORE_SUPPORTED")
        self.assertIn(
            "V025_MINOR_HARD_FAILURE_REMOVED",
            result["v026_reclassification_reasons"],
        )

    def test_monoallelic_single_snp_stays_supported(self):
        row = make_row(
            "b2",
            "PEG10",
            19,
            1,
            0.01,
            1e-5,
            "SINGLE_SNP_ONLY",
            "SINGLE_SNP",
        )
        result = self.classify(row)
        self.assertEqual(
            result["v026_evidence_class"], "ASE_SUPPORTED_SINGLE_SNP"
        )
        self.assertEqual(result["v026_is_core_ase"], "0")

    def test_phase_review_is_not_core(self):
        row = make_row(
            "b3",
            "GENE3",
            18,
            2,
            0.01,
            1e-4,
            "MULTI_SNP_LOO_ROBUST",
            "REVIEW",
        )
        result = self.classify(row)
        self.assertEqual(
            result["v026_evidence_class"],
            "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED",
        )

    def test_minor_count_boundaries_are_shape_only(self):
        for minor, expected in [
            (1, "MONOALLELIC_LIKE"),
            (2, "LOW_MINOR"),
            (3, "TWO_HAPLOTYPE_SUPPORTED"),
        ]:
            row = make_row(
                "minor{}".format(minor),
                "GENE",
                20 - minor,
                minor,
                0.01,
                1e-4,
                "MULTI_SNP_SUPPORTED",
                "PASS",
            )
            result = self.classify(row)
            self.assertEqual(result["v026_imbalance_shape"], expected)
            self.assertEqual(result["v026_evidence_class"], "ASE_CORE_SUPPORTED")

    def test_signal_boundaries(self):
        cases = [
            (13, 7, 0.05, "FDR_STRONG_EFFECT"),
            (12, 8, 0.05, "FDR_MODERATE_EFFECT"),
            (11, 9, 0.05, "FDR_SMALL_EFFECT"),
            (13, 7, 0.10, "NEAR_FDR_STRONG_EFFECT"),
        ]
        for index, (a_count, b_count, qvalue, expected) in enumerate(cases):
            row = make_row(
                "signal{}".format(index),
                "GENE",
                a_count,
                b_count,
                qvalue,
                0.01,
                "MULTI_SNP_SUPPORTED",
                "PASS",
            )
            self.assertEqual(
                reclassifier.statistical_signal(row, self.thresholds), expected
            )

    def test_hard_qc_overrides_signal(self):
        row = make_row(
            "b4",
            "GENE4",
            20,
            0,
            0.001,
            1e-8,
            "MULTI_SNP_LOO_ROBUST",
            "PASS",
            balanced=0,
        )
        result = self.classify(row)
        self.assertEqual(result["v026_evidence_class"], "TECHNICAL_FAIL")
        self.assertIn("NO_WGS_BALANCED_SNP", result["v026_technical_qc_reasons"])


class IntegrationTests(unittest.TestCase):
    def test_independent_bh_recalculation(self):
        with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as handle:
            db_path = handle.name
        try:
            connection = validator.create_db(db_path)
            observations = [
                ("m1", "TEST1", "Liver", 0.01, 0.03),
                ("m2", "TEST1", "Liver", 0.02, 0.03),
                ("m3", "TEST1", "Liver", 0.20, 0.20),
            ]
            connection.executemany(
                "INSERT INTO measurement_obs VALUES (?,?,?,?,?)", observations
            )
            connection.commit()
            errors = []
            n_measurements = validator.validate_bh(connection, errors, 1e-12)
            self.assertEqual(n_measurements, 3)
            self.assertEqual(errors, [])
            connection.close()
        finally:
            try:
                os.remove(db_path)
            except FileNotFoundError:
                pass

    def test_reclassifier_and_validator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            audit = source / "audit"
            audit.mkdir(parents=True)

            rows = [
                make_row(
                    "core",
                    "SHARED_GENE",
                    18,
                    2,
                    0.01,
                    0.001,
                    "MULTI_SNP_LOO_ROBUST",
                    "PASS",
                    old_class="PRIMARY_ROBUST",
                ),
                make_row(
                    "single",
                    "SHARED_GENE",
                    20,
                    0,
                    0.001,
                    1e-5,
                    "SINGLE_SNP_ONLY",
                    "SINGLE_SNP",
                    old_reasons="EXTREME_MINOR",
                ),
                make_row(
                    "phase",
                    "PHASE_GENE",
                    18,
                    2,
                    0.01,
                    0.001,
                    "MULTI_SNP_SUPPORTED",
                    "FAIL",
                ),
                make_row(
                    "moderate",
                    "MODERATE_GENE",
                    12,
                    8,
                    0.01,
                    0.01,
                    "MULTI_SNP_SUPPORTED",
                    "PASS",
                ),
                make_row(
                    "near",
                    "NEAR_GENE",
                    13,
                    7,
                    0.08,
                    0.01,
                    "MULTI_SNP_SUPPORTED",
                    "PASS",
                ),
            ]
            with open(
                source / "ase_haplotype_block_results.tsv",
                "w",
                encoding="utf-8",
                newline="",
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=SOURCE_FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            with open(source / "ase_v025_validation.json", "w", encoding="utf-8") as handle:
                json.dump({"status": "PASS", "n_errors": 0, "version": "0.2.5"}, handle)
            with open(source / "ase_v025_qc.json", "w", encoding="utf-8") as handle:
                json.dump({"version": "0.2.5", "n_rows": len(rows)}, handle)
            with open(
                audit / "ase_measurement_unique.tsv",
                "w",
                encoding="utf-8",
                newline="",
            ) as handle:
                fields = [
                    "measurement_id",
                    "direct_total_count",
                    "strand_resolved_total",
                    "rescue_fraction",
                    "rescue_dependency",
                ]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(
                        {
                            "measurement_id": row["measurement_id"],
                            "direct_total_count": row["ase_total_count"],
                            "strand_resolved_total": "0",
                            "rescue_fraction": "0",
                            "rescue_dependency": "NO_STRAND_RESOLVED_FRAGMENTS",
                        }
                    )

            args = Namespace(
                input_dir=str(source),
                block_results=None,
                validation_json=None,
                qc_json=None,
                assignment_audit=None,
                sample="TEST1",
                out_dir=str(output),
                allow_unvalidated_input=False,
                force=False,
                fdr=0.05,
                near_fdr=0.10,
                nominal_p=0.05,
                strong_major_fraction=0.65,
                moderate_major_fraction=0.60,
                max_fragment_hap_conflict_fraction=0.10,
            )
            metadata = reclassifier.run_reclassification(args)
            self.assertEqual(metadata["n_block_gene_rows"], len(rows))
            self.assertEqual(metadata["n_gene_summary_rows"], 4)

            with open(output / "ase_gene_summary.tsv", "r", encoding="utf-8") as handle:
                summaries = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(
                summaries["SHARED_GENE"]["top_haplotype_block_id"], "core"
            )
            self.assertEqual(
                summaries["SHARED_GENE"]["gene_ase_status_v026"],
                "ASE_GENE_CORE",
            )
            validation_code = validator.main(
                ["--out-dir", str(output), "--expect-sample", "TEST1", "--skip-bh-recompute"]
            )
            self.assertEqual(validation_code, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
