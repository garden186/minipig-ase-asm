#!/usr/bin/env python3
"""Focused tests for GATK-to-primary-block linkage."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("validate_gatk_against_primary_blocks.py")
SPEC = importlib.util.spec_from_file_location("block_validation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class PrimaryBlockValidationTests(unittest.TestCase):
    def test_tissue_normalization(self) -> None:
        self.assertEqual(MODULE.normalize_tissue("L-blood"), "Blood")
        self.assertEqual(MODULE.normalize_tissue("blood"), "Blood")
        self.assertEqual(MODULE.normalize_tissue("Tenderlo"), "Tenderloin")

    def test_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            measurements = root / "core_gene_measurements.tsv"
            archive = root / "blocks.tar.gz"
            sites = root / "site_tests.tsv.gz"
            recurrent = root / "recurrent.tsv"
            output = root / "output"

            measurement_header = sorted(MODULE.CORE_MEASUREMENT_REQUIRED)
            measurement_row = {
                "sample": "A", "tissue": "Blood", "gene_id": "GENE1",
                "gene_name": "GeneOne", "top_haplotype_block_id": "BLOCK1",
                "core_haplotype_block_ids": "BLOCK1", "gene_ase_status_v026": "ASE_GENE_CORE",
            }
            with measurements.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=measurement_header, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerow(measurement_row)

            block_header = sorted(MODULE.BLOCK_REQUIRED)
            block_row = {
                "sample": "A", "tissue": "L-blood", "gene_id": "GENE1",
                "gene_name": "GeneOne", "haplotype_block_id": "BLOCK1",
                "measurement_id": "M1", "block_snps": "1_100_A_G,1_200_C_T",
                "gene_variants": "1_100_A_G,1_200_C_T", "haplotypeA": "A,C",
                "haplotypeB": "G,T", "v026_ase_direction": "HAP_A",
                "v026_evidence_class": "ASE_CORE_ROBUST", "v026_is_core_ase": "1",
            }
            with tarfile.open(archive, "w:gz") as handle:
                for index in range(10):
                    buffer = io.StringIO()
                    writer = csv.DictWriter(buffer, fieldnames=block_header, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    if index == 0:
                        writer.writerow(block_row)
                    data = buffer.getvalue().encode("utf-8")
                    info = tarfile.TarInfo(f"S{index}/ase_core_haplotype_blocks.tsv")
                    info.size = len(data)
                    handle.addfile(info, io.BytesIO(data))

            site_header = sorted(MODULE.SITE_REQUIRED)
            common = {
                "sample": "A", "tissue": "Blood", "contig": "1",
                "ref_count": "15", "alt_count": "0", "total_count": "15",
                "major_allele_fraction": "1", "imbalance_direction": "REF",
                "p_exact_two_sided": "0.0001", "q_bh_within_animal_tissue": "0.001",
            }
            rows = []
            row = dict(common); row.update({"position":"100","variant_id":"1_100_A_G","ref_allele":"A","alt_allele":"G","major_allele":"A","is_snp_level_ase":"1"}); rows.append(row)
            row = dict(common); row.update({"position":"200","variant_id":"1_200_C_T","ref_allele":"C","alt_allele":"T","major_allele":"C","is_snp_level_ase":"0"}); rows.append(row)
            with gzip.open(sites, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=site_header, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(rows)

            recurrent_header = sorted(MODULE.RECURRENCE_REQUIRED)
            recurrent_row = {
                "tissue":"Blood", "physical_snp_id":"1_100_A_G",
                "n_animals_testable":"3", "n_animals_with_snp_level_ase":"3",
                "recurrence_q_bh_global":"0.01", "is_statistically_recurrent_snp_ase":"1",
            }
            with recurrent.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=recurrent_header, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerow(recurrent_row)

            args = MODULE.parse_args([
                "--site-tests", str(sites), "--core-measurements", str(measurements),
                "--block-archive", str(archive), "--recurrent-snps", str(recurrent),
                "--out-dir", str(output),
            ])
            MODULE.run(args)
            with (output / "primary_gene_measurement_validation_summary.tsv").open(encoding="utf-8") as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["n_gene_variant_contexts"], "2")
            self.assertEqual(summary["n_gatk_snp_level_ase_variants"], "1")
            self.assertEqual(summary["n_direction_concordant_variants"], "1")
            self.assertEqual(summary["gatk_measurement_validation_status"], "DIRECTION_CONCORDANT_ONLY")


if __name__ == "__main__":
    unittest.main()
