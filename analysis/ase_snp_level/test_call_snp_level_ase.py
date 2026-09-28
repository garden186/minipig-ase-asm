#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
CALLER = HERE / "call_snp_level_ase.py"


def write_tsv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        handle_context = gzip.open(path, "wt", newline="")
    else:
        handle_context = path.open("w", encoding="utf-8", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        handle_context = gzip.open(path, "rt", newline="")
    else:
        handle_context = path.open("r", encoding="utf-8", newline="")
    with handle_context as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


BLOCK_HEADER = [
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id", "measurement_id",
    "analysis_set", "block_snps", "haplotypeA", "haplotypeB", "gene_variants_with_reads",
    "v026_technical_qc_status", "v026_is_core_ase", "v026_ase_direction",
]
FRAGMENT_HEADER = [
    "sample", "tissue", "haplotype_block_id", "qname", "gene_id", "haplotype",
    "final_status", "gene_supported_variants",
]


class CallerTest(unittest.TestCase):
    def test_site_counts_bh_and_core_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            block_root = root / "v026"
            fragment_root = root / "fragment"
            v025_root = root / "v025"
            out = root / "out"
            sample = "0001"
            variants = ["1_100_A_G", "1_110_C_T"]
            write_tsv(
                block_root / sample / "ase_haplotype_block_reclassified.tsv.gz",
                BLOCK_HEADER,
                [
                    {
                        "sample": sample,
                        "tissue": "Liver",
                        "gene_id": "G1",
                        "gene_name": "GENE1",
                        "haplotype_block_id": "B1",
                        "measurement_id": "M1",
                        "analysis_set": "MAIN",
                        "block_snps": ",".join(variants),
                        "haplotypeA": "A,C",
                        "haplotypeB": "G,T",
                        "gene_variants_with_reads": ",".join(variants),
                        "v026_technical_qc_status": "PASS",
                        "v026_is_core_ase": "1",
                        "v026_ase_direction": "HAP_A",
                    }
                ],
            )
            rows = []
            for index in range(18):
                rows.append(
                    {
                        "sample": sample,
                        "tissue": "Liver",
                        "haplotype_block_id": "B1",
                        "qname": f"a{index}",
                        "gene_id": "G1",
                        "haplotype": "A",
                        "final_status": "ASSIGNED",
                        "gene_supported_variants": "1_100_A_G" if index >= 10 else ",".join(variants),
                    }
                )
            for index in range(2):
                rows.append(
                    {
                        "sample": sample,
                        "tissue": "Liver",
                        "haplotype_block_id": "B1",
                        "qname": f"b{index}",
                        "gene_id": "G1",
                        "haplotype": "B",
                        "final_status": "ASSIGNED",
                        "gene_supported_variants": ",".join(variants),
                    }
                )
            for index in range(8):
                rows.append(
                    {
                        "sample": sample,
                        "tissue": "Liver",
                        "haplotype_block_id": "B1",
                        "qname": f"b_extra{index}",
                        "gene_id": "G1",
                        "haplotype": "B",
                        "final_status": "ASSIGNED",
                        "gene_supported_variants": "1_110_C_T",
                    }
                )
            write_tsv(fragment_root / sample / "audit/ase_fragment_assignment.tsv.gz", FRAGMENT_HEADER, rows)
            (v025_root / sample).mkdir(parents=True)

            result = subprocess.run(
                [
                    sys.executable,
                    str(CALLER),
                    "--block-root", str(block_root),
                    "--fragment-root", str(fragment_root),
                    "--v025-root", str(v025_root),
                    "--out-dir", str(out),
                    "--samples", sample,
                    "--min-total", "10",
                    "--min-major-fraction", "0.65",
                    "--fdr", "0.05",
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn('"status": "PASS"', result.stdout)
            site_rows = read_tsv(out / "core_block_site_level_ase_measurements.tsv.gz")
            self.assertEqual(2, len(site_rows))
            by_variant = {row["variant_id"]: row for row in site_rows}
            self.assertEqual("18", by_variant["1_100_A_G"]["ref_fragment_count"])
            self.assertEqual("2", by_variant["1_100_A_G"]["alt_fragment_count"])
            self.assertEqual("1", by_variant["1_100_A_G"]["site_level_ase_signal"])
            self.assertEqual("CONSISTENT", by_variant["1_100_A_G"]["core_parent_direction_status"])
            self.assertEqual("10", by_variant["1_110_C_T"]["ref_fragment_count"])
            self.assertEqual("10", by_variant["1_110_C_T"]["alt_fragment_count"])
            self.assertEqual("0", by_variant["1_110_C_T"]["site_level_ase_signal"])
            gene_rows = read_tsv(out / "core_gene_variant_tissue_summary.tsv.gz")
            self.assertEqual(2, len(gene_rows))
            gene_by_variant = {row["variant_id"]: row for row in gene_rows}
            self.assertEqual("1", gene_by_variant["1_100_A_G"]["n_samples_site_level_ase_signal"])
            self.assertEqual("0", gene_by_variant["1_110_C_T"]["n_samples_site_level_ase_signal"])
            metadata = json.loads((out / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual("PASS", metadata["status"])


if __name__ == "__main__":
    unittest.main()
