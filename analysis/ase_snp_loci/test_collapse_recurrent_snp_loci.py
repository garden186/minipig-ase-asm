#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import io
import tarfile
import tempfile
import unittest
from pathlib import Path

from collapse_recurrent_snp_loci import run_analysis


RECURRENT_FIELDS = [
    "gene_id", "gene_name", "tissue", "variant_id", "chrom", "position", "ref", "alt",
    "n_animals_testable", "animals_testable", "n_animals_with_block_concordant_snp_ase",
    "animals_with_block_concordant_snp_ase", "recurrence_fraction",
    "recurrence_p_poisson_binomial", "recurrence_q_bh_global",
    "is_statistically_recurrent_snp_ase", "original_ase_assignment_class",
    "canonical_transcript_features", "annotation_status",
]
SITE_FIELDS = [
    "sample", "tissue", "variant_id", "ref_fragment_count", "alt_fragment_count",
    "total_informative_fragments", "major_allele", "site_level_ase_signal",
    "core_parent_measurement_ids",
]
BLOCK_FIELDS = [
    "sample", "tissue", "gene_id", "gene_name", "haplotype_block_id", "block_start",
    "block_end", "block_snps", "haplotypeA", "haplotypeB", "gene_variants_with_reads",
    "measurement_id", "analysis_set", "v026_technical_qc_status", "v026_is_core_ase",
    "v026_ase_direction",
]


class IntegrationTest(unittest.TestCase):
    def test_shared_core_block_connects_only_supported_snps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recurrent = root / "recurrent.tsv"
            sites = root / "sites.tsv.gz"
            archive = root / "blocks.tar.gz"
            out_dir = root / "out"

            recurrent_rows = [
                self.recurrent("1_100_A_G", 100, "S1,S2", "0.01"),
                self.recurrent("1_120_C_T", 120, "S1,S2", "0.02"),
                self.recurrent("1_140_G_A", 140, "S2,S3", "0.03"),
            ]
            self.write_tsv(recurrent, RECURRENT_FIELDS, recurrent_rows)

            site_rows = [
                self.site("S1", "1_100_A_G", "M1", "G", 20),
                self.site("S1", "1_120_C_T", "M1", "T", 22),
                self.site("S2", "1_100_A_G", "M2", "G", 24),
                self.site("S2", "1_120_C_T", "M2", "T", 26),
                self.site("S2", "1_140_G_A", "M3", "A", 28),
                self.site("S3", "1_140_G_A", "M4", "A", 30),
            ]
            with gzip.open(sites, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=SITE_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(site_rows)

            block_rows = [
                self.block("S1", "M1", "B1", "1_100_A_G,1_120_C_T", "G,T"),
                self.block("S2", "M2", "B2", "1_100_A_G,1_120_C_T", "G,T"),
                self.block("S2", "M3", "B3", "1_140_G_A", "A"),
                self.block("S3", "M4", "B4", "1_140_G_A", "A"),
            ]
            self.write_archive(archive, block_rows)

            metadata = run_analysis(
                argparse.Namespace(
                    recurrent_annotation=recurrent,
                    core_site_measurements=sites,
                    block_archive=archive,
                    out_dir=out_dir,
                )
            )
            self.assertEqual(metadata["counts"]["n_recurrent_gene_tissue_snp_hypotheses"], 3)
            self.assertEqual(metadata["counts"]["n_recurrent_ase_loci"], 2)
            self.assertEqual(metadata["counts"]["n_direct_core_block_edges"], 1)

            with (out_dir / "recurrent_ase_loci.tsv").open(encoding="utf-8", newline="") as handle:
                loci = list(csv.DictReader(handle, delimiter="\t"))
            multi = next(row for row in loci if row["n_recurrent_snps"] == "2")
            singleton = next(row for row in loci if row["n_recurrent_snps"] == "1")
            self.assertEqual(multi["recurrent_snp_ids"], "1_100_A_G,1_120_C_T")
            self.assertEqual(multi["lead_snp_id"], "1_100_A_G")
            self.assertEqual(singleton["recurrent_snp_ids"], "1_140_G_A")

    @staticmethod
    def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def recurrent(variant: str, position: int, samples: str, qvalue: str) -> dict[str, str]:
        _, _, ref, alt = variant.split("_")
        return {
            "gene_id": "G1", "gene_name": "Gene1", "tissue": "Liver", "variant_id": variant,
            "chrom": "1", "position": str(position), "ref": ref, "alt": alt,
            "n_animals_testable": "3", "animals_testable": "S1,S2,S3",
            "n_animals_with_block_concordant_snp_ase": "2",
            "animals_with_block_concordant_snp_ase": samples, "recurrence_fraction": "0.6666666667",
            "recurrence_p_poisson_binomial": qvalue, "recurrence_q_bh_global": qvalue,
            "is_statistically_recurrent_snp_ase": "1", "original_ase_assignment_class": "UNIQUE_EXONIC",
            "canonical_transcript_features": "CDS", "annotation_status": "PASS",
        }

    @staticmethod
    def site(sample: str, variant: str, measurement: str, major: str, total: int) -> dict[str, str]:
        return {
            "sample": sample, "tissue": "Liver", "variant_id": variant,
            "ref_fragment_count": "2", "alt_fragment_count": str(total - 2),
            "total_informative_fragments": str(total), "major_allele": major,
            "site_level_ase_signal": "1", "core_parent_measurement_ids": measurement,
        }

    @staticmethod
    def block(sample: str, measurement: str, block_id: str, variants: str, major: str) -> dict[str, str]:
        variant_values = variants.split(",")
        major_values = major.split(",")
        minor_values = []
        for variant, allele in zip(variant_values, major_values):
            ref, alt = variant.split("_")[2:]
            minor_values.append(ref if allele == alt else alt)
        return {
            "sample": sample, "tissue": "Liver", "gene_id": "G1", "gene_name": "Gene1",
            "haplotype_block_id": block_id, "block_start": "90", "block_end": "150",
            "block_snps": variants, "haplotypeA": major, "haplotypeB": ",".join(minor_values),
            "gene_variants_with_reads": variants, "measurement_id": measurement,
            "analysis_set": "MAIN", "v026_technical_qc_status": "PASS", "v026_is_core_ase": "1",
            "v026_ase_direction": "HAP_A",
        }

    @staticmethod
    def write_archive(path: Path, rows: list[dict[str, str]]) -> None:
        by_sample: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            by_sample.setdefault(row["sample"], []).append(row)
        with tarfile.open(path, "w:gz") as archive:
            for sample, sample_rows in by_sample.items():
                buffer = io.StringIO()
                writer = csv.DictWriter(buffer, fieldnames=BLOCK_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(sample_rows)
                payload = buffer.getvalue().encode("utf-8")
                info = tarfile.TarInfo(f"{sample}/ase_core_haplotype_blocks.tsv")
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))


if __name__ == "__main__":
    unittest.main()
