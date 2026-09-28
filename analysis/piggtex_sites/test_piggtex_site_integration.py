#!/usr/bin/env python3
"""Focused tests for GATK–PigGTEx exact-site integration."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("integrate_gatk_with_piggtex_sites.py")
SPEC = importlib.util.spec_from_file_location("piggtex_site", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_table(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        handle = opener(path, "wt", encoding="utf-8", newline="")
    else:
        handle = opener("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


class PigGTExSiteIntegrationTests(unittest.TestCase):
    def test_unordered_key(self) -> None:
        self.assertEqual(
            MODULE.unordered_key("chr1", 100, "G", "A"),
            MODULE.unordered_key("1", 100, "A", "G"),
        )

    def test_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calls = root / "calls.tsv.gz"
            recurrent = root / "recurrent.tsv"
            blocks = root / "blocks.tsv.gz"
            pairs = root / "pairs.tsv"
            tissue_map = root / "tissue_map.tsv"
            pig = root / "pig.tsv.gz"
            output = root / "output"

            call_header = sorted(MODULE.GATK_REQUIRED)
            call_rows = []
            for sample in ("A", "B"):
                call_rows.append({
                    "sample":sample, "tissue":"Blood", "contig":"1", "position":"100",
                    "ref_allele":"A", "alt_allele":"G", "major_allele":"A",
                    "q_bh_within_animal_tissue":"0.01", "is_snp_level_ase":"1",
                })
            write_table(calls, call_header, call_rows)

            recurrence_row = {
                "tissue":"Blood", "contig":"1", "position":"100", "ref_allele":"A",
                "alt_allele":"G", "n_animals_testable":"2",
                "n_animals_with_snp_level_ase":"2", "recurrence_q_bh_global":"0.01",
                "is_statistically_recurrent_snp_ase":"1",
            }
            write_table(recurrent, sorted(MODULE.RECURRENCE_REQUIRED), [recurrence_row])

            block_rows = []
            for sample in ("A", "B"):
                block_rows.append({
                    "sample":sample, "tissue":"Blood", "gene_id":"GENE1",
                    "gene_name":"GeneOne", "haplotype_block_id":f"BLOCK_{sample}",
                    "variant_id":"1_100_A_G", "contig":"1", "position":"100",
                    "ref_allele":"A", "alt_allele":"G", "gatk_snp_level_ase":"1",
                    "gatk_block_direction_status":"DIRECTION_CONCORDANT",
                })
            write_table(blocks, sorted(MODULE.BLOCK_REQUIRED), block_rows)

            pair_rows = []
            for index in range(39):
                pair_rows.append({
                    "gene_id":"GENE1" if index == 0 else f"GENE{index + 1}",
                    "gene_name":"GeneOne" if index == 0 else f"Gene{index + 1}",
                    "target_tissue":"Blood", "pair_level_primary":"1",
                })
            write_table(pairs, sorted(MODULE.TISSUE_PAIR_REQUIRED), pair_rows)
            write_table(tissue_map, sorted(MODULE.TISSUE_MAP_REQUIRED), [{"internal_tissue":"Blood","piggtex_tissue":"Blood"}])
            pig_row = {
                "chr":"1", "position":"100", "refAllele":"G", "altAllele":"A",
                "refCount":"2", "altCount":"18", "FDR":"0.01", "SampleID":"P1",
                "Breed":"Cross", "Tissue":"Blood", "pCADD":"12.5",
            }
            write_table(pig, sorted(MODULE.PIGGTEX_REQUIRED), [pig_row])

            args = MODULE.parse_args([
                "--gatk-calls",str(calls), "--gatk-recurrent",str(recurrent),
                "--block-validation",str(blocks), "--tissue-enriched-pairs",str(pairs),
                "--tissue-map",str(tissue_map), "--piggtex-sites",str(pig),
                "--out-dir",str(output),
            ])
            MODULE.run(args)
            with (output / "recurrent_gatk_snp_piggtex_site_support.tsv").open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["n_internal_gatk_ase_measurements"], "2")
            self.assertEqual(row["piggtex_significant_rows_same_mapped_tissue"], "1")
            self.assertEqual(row["piggtex_same_tissue_major_internal_ref_rows"], "1")
            self.assertEqual(row["piggtex_site_support_status"], "SIGNIFICANT_ASE_SITE_SAME_MAPPED_TISSUE")


if __name__ == "__main__":
    unittest.main()
