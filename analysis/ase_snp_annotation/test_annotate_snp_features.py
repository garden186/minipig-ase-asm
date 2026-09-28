#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from annotate_snp_features import run_analysis


INPUT_HEADER = [
    "gene_id", "gene_name", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "n_animals_testable", "n_animals_with_block_concordant_snp_ase",
    "recurrence_q_bh_global", "is_statistically_recurrent_snp_ase",
]


class IntegrationTest(unittest.TestCase):
    def test_original_exon_union_and_transcript_features(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "test.gtf"
            all_path = root / "all.tsv.gz"
            recurrent_path = root / "recurrent.tsv"
            out_dir = root / "out"
            gtf.write_text(self.gtf_text(), encoding="utf-8")

            rows = [
                self.row("G1", "Gene1", "v1", 120, "1"),
                self.row("G1", "Gene1", "v2", 160, "1"),
            ]
            with gzip.open(all_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=INPUT_HEADER, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with recurrent_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=INPUT_HEADER, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

            args = argparse.Namespace(
                all_hypotheses=all_path,
                recurrent=recurrent_path,
                gtf=gtf,
                out_dir=out_dir,
            )
            metadata = run_analysis(args)
            self.assertEqual(metadata["counts"]["n_recurrent_annotation_rows"], 2)

            with (out_dir / "recurrent_gene_tissue_snp_feature_annotations.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                annotated = {row["variant_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(annotated["v1"]["original_ase_assignment_class"], "UNIQUE_EXONIC")
            self.assertEqual(annotated["v1"]["canonical_transcript_features"], "CDS")
            self.assertEqual(annotated["v2"]["original_ase_assignment_class"], "AMBIGUOUS_SHARED_EXON")
            self.assertEqual(
                annotated["v2"]["site_gene_assignment_requirement"],
                "STRAND_RESOLUTION_REQUIRED_FOR_SITE",
            )

    @staticmethod
    def row(gene_id: str, gene_name: str, variant_id: str, position: int, recurrent: str) -> dict[str, str]:
        return {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "tissue": "Liver",
            "variant_id": variant_id,
            "chrom": "1",
            "position": str(position),
            "ref": "A",
            "alt": "G",
            "n_animals_testable": "3",
            "n_animals_with_block_concordant_snp_ase": "2",
            "recurrence_q_bh_global": "0.01",
            "is_statistically_recurrent_snp_ase": recurrent,
        }

    @staticmethod
    def gtf_text() -> str:
        return """1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; gene_biotype \"protein_coding\";
1\ttest\ttranscript\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_name \"T1\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\";
1\ttest\texon\t100\t170\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\";
1\ttest\tCDS\t110\t140\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\";
1\ttest\tthree_prime_utr\t141\t170\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\";
1\ttest\tgene\t150\t250\t.\t-\t.\tgene_id \"G2\"; gene_name \"Gene2\"; gene_biotype \"protein_coding\";
1\ttest\ttranscript\t150\t250\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\"; transcript_name \"T2\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\";
1\ttest\texon\t150\t180\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\";
"""


if __name__ == "__main__":
    unittest.main()
