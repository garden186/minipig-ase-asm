#!/usr/bin/env python3
import csv
import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


PKG_DIR = Path(__file__).resolve().parent
if str(PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PKG_DIR))

import build_ase_v026_cohort as builder
import validate_ase_v026_cohort as validator


def gene_row(sample, tissue, gene_id, status):
    row = {field: "" for field in builder.GENE_HEADER}
    row.update(
        {
            "sample": sample,
            "tissue": tissue,
            "gene_id": gene_id,
            "gene_name": gene_id,
            "gene_biotype": "protein_coding",
            "n_informative_haplotype_blocks": "1",
            "n_testable_haplotype_blocks": "1",
            "n_fdr_strong_blocks": "0",
            "n_core_robust_blocks": "0",
            "n_core_supported_blocks": "0",
            "n_supported_single_snp_blocks": "0",
            "n_phase_unresolved_blocks": "0",
            "n_support_unresolved_blocks": "0",
            "n_fdr_moderate_effect_blocks": "0",
            "n_near_fdr_blocks": "0",
            "n_nominal_directional_blocks": "0",
            "n_statistical_small_effect_blocks": "0",
            "n_monoallelic_like_blocks": "0",
            "n_low_minor_blocks": "0",
            "n_two_haplotype_supported_blocks": "1",
            "gene_phase_scope": "SINGLE_BLOCK",
            "gene_ase_status_v026": status,
            "top_haplotype_block_id": f"{sample}:{gene_id}:B1",
            "top_measurement_id": f"{sample}|{tissue}|{gene_id}",
            "top_total_count": "20",
            "core_haplotype_block_ids": "",
            "supported_haplotype_block_ids": "",
            "sensitivity_haplotype_block_ids": "",
        }
    )
    if status == builder.CORE_STATUS:
        row["n_fdr_strong_blocks"] = "1"
        row["n_core_robust_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_CORE_ROBUST"
        row["core_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_SUPPORTED_SINGLE_SNP":
        row["n_fdr_strong_blocks"] = "1"
        row["n_supported_single_snp_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_SUPPORTED_SINGLE_SNP"
        row["supported_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED":
        row["n_fdr_strong_blocks"] = "1"
        row["n_phase_unresolved_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_LOCAL_SIGNAL_PHASE_UNRESOLVED"
        row["supported_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_FDR_MODERATE":
        row["n_fdr_moderate_effect_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_FDR_MODERATE_EFFECT"
        row["sensitivity_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_NEAR_FDR":
        row["n_near_fdr_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_NEAR_FDR"
        row["sensitivity_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_NOMINAL":
        row["n_nominal_directional_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_NOMINAL_DIRECTIONAL"
        row["sensitivity_haplotype_block_ids"] = row["top_haplotype_block_id"]
    elif status == "ASE_GENE_STATISTICAL_SMALL_EFFECT":
        row["n_statistical_small_effect_blocks"] = "1"
        row["top_evidence_class_v026"] = "ASE_STATISTICAL_SMALL_EFFECT"
        row["sensitivity_haplotype_block_ids"] = row["top_haplotype_block_id"]
    else:
        row["top_evidence_class_v026"] = "NO_ASE_SIGNAL"
    return row


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=builder.GENE_HEADER,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(root, tissues=None):
    samples = ["S1", "S2", "S3"]
    tissues = tissues or {sample: "Heart" for sample in samples}
    for sample in samples:
        sample_dir = root / sample
        sample_dir.mkdir(parents=True)
        statuses = [builder.CORE_STATUS]
        if sample == "S1":
            statuses.append("ASE_GENE_SUPPORTED_SINGLE_SNP")
        elif sample == "S2":
            statuses.append("ASE_GENE_PHASE_OR_SUPPORT_UNRESOLVED")
        else:
            statuses.append("ASE_GENE_FDR_MODERATE")
        statuses.extend(["NOT_ASE_GENE"] * 8)
        rows = [
            gene_row(sample, tissues[sample], f"G{index + 1}", status)
            for index, status in enumerate(statuses)
        ]
        core = [row for row in rows if row["gene_ase_status_v026"] == builder.CORE_STATUS]
        supported = [
            row
            for row in rows
            if row["gene_ase_status_v026"] in builder.SUPPORTED_STATUSES
        ]
        sensitivity = [
            row
            for row in rows
            if row["gene_ase_status_v026"] in builder.SENSITIVITY_STATUSES
        ]
        write_tsv(sample_dir / "ase_gene_summary.tsv", rows)
        write_tsv(sample_dir / "ase_core_genes.tsv", core)
        write_tsv(sample_dir / "ase_supported_genes.tsv", supported)
        write_tsv(sample_dir / "ase_sensitivity_genes.tsv", sensitivity)
        counts = dict(sorted(Counter(
            row["gene_ase_status_v026"] for row in rows
        ).items()))
        validation = {
            "version": "0.2.6",
            "status": "PASS",
            "n_errors": 0,
            "n_warnings": 0,
            "n_gene_summary_rows": len(rows),
            "gene_status_counts": counts,
        }
        qc = {
            "version": "0.2.6",
            "source_validation_status": "PASS",
            "source_qc_loaded": True,
            "source_block_results_sha256": sample * 16,
            "thresholds": {
                "fdr": 0.05,
                "near_fdr": 0.10,
                "nominal_p": 0.05,
                "strong_major_fraction": 0.65,
                "moderate_major_fraction": 0.60,
                "max_fragment_hap_conflict_fraction": 0.10,
            },
            "gene_status_counts": counts,
        }
        with (sample_dir / "ase_v026_validation.json").open("w", encoding="utf-8") as handle:
            json.dump(validation, handle)
        with (sample_dir / "ase_v026_qc.json").open("w", encoding="utf-8") as handle:
            json.dump(qc, handle)
    return samples


class CohortLogicTests(unittest.TestCase):
    def test_exact_core_recurrence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = make_fixture(root)
            outputs, metadata = builder.build_expected(
                root, samples, set(), recurrence_q_threshold=0.05
            )
            recurrence = {
                (row["gene_id"], row["tissue"]): row
                for row in outputs["gene_tissue_core_recurrence_test_results.tsv"]
            }
            g1 = recurrence[("G1", "Heart")]
            self.assertAlmostEqual(float(g1["recurrence_p"]), 0.001)
            self.assertAlmostEqual(float(g1["recurrence_q_bh"]), 0.01)
            self.assertEqual(g1["n_samples_core"], 3)
            self.assertEqual(g1["is_statistically_recurrent_core"], 1)
            self.assertEqual(metadata["n_statistically_recurrent_core_gene_tissue_pairs"], 1)

    def test_supported_and_sensitivity_are_descriptive_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = make_fixture(root)
            outputs, metadata = builder.build_expected(root, samples, set())
            burdens = outputs["sample_tissue_burden.tsv"]
            self.assertTrue(all(row["n_core_genes"] == 1 for row in burdens))
            self.assertEqual(len(outputs["supported_gene_measurements.tsv"]), 2)
            self.assertEqual(len(outputs["sensitivity_gene_measurements.tsv"]), 1)
            self.assertIn(
                "descriptive only",
                metadata["statistical_recurrence"]["supported_and_sensitivity_policy"],
            )

    def test_blood_aliases_merge_before_recurrence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = make_fixture(
                root, {"S1": "blood", "S2": "L-blood", "S3": "blood"}
            )
            outputs, metadata = builder.build_expected(
                root,
                samples,
                set(),
                tissue_aliases=builder.DEFAULT_TISSUE_ALIASES,
                expect_sample_tissue_combinations=3,
                expect_tissues=1,
            )
            row = next(
                row
                for row in outputs["gene_tissue_core_recurrence_test_results.tsv"]
                if row["gene_id"] == "G1"
            )
            self.assertEqual(row["tissue"], "Blood")
            self.assertEqual(row["n_samples_core"], 3)
            self.assertEqual(metadata["n_tissues"], 1)

    def test_modified_subset_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = make_fixture(root)
            path = root / "S1" / "ase_core_genes.tsv"
            write_tsv(path, [])
            with self.assertRaises(builder.CohortError):
                builder.build_expected(root, samples, set())


class CohortIntegrationTests(unittest.TestCase):
    def test_build_materialize_and_validate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            project_tmp = root / "tmp"
            samples = make_fixture(source)
            expected, metadata = builder.build_expected(source, samples, set())
            builder.materialize_outputs(
                expected, metadata, output, project_tmp, force=False
            )
            code = validator.main(
                [
                    "--input-root", str(source),
                    "--out-dir", str(output),
                    "--samples", ",".join(samples),
                    "--exclude-tissues", "",
                    "--tissue-aliases", "",
                    "--expect-sample-tissue-combinations", "3",
                    "--expect-tissues", "1",
                ]
            )
            self.assertEqual(code, 0)
            with (output / "cohort_validation.json").open(encoding="utf-8") as handle:
                report = json.load(handle)
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["n_errors"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
