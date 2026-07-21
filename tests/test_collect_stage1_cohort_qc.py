#!/usr/bin/env python3

from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "stage1" / "collect_stage1_cohort_qc.py"


def write_metric_value(path: Path, values: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(values.items())


class CollectStage1CohortQcTest(unittest.TestCase):
    def test_collects_concise_consistent_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            project = Path(temporary_directory)
            sample_sheet = project / "config" / "samples.tsv"
            sample_sheet.parent.mkdir(parents=True)
            sample_sheet.write_text(
                "sample_id\twgs_r1\twgs_r2\nS1\ta\tb\nS2\tc\td\n",
                encoding="utf-8",
            )
            for index, sample in enumerate(("S1", "S2"), start=1):
                final = 1000 + index
                write_metric_value(
                    project / "results/qc/wgs" / f"{sample}.stage1_alignment_qc.tsv",
                    {
                        "sample": sample,
                        "mean_autosome_depth": str(30 + index),
                        "autosome_fraction_ge_10x": "0.95",
                        "primary_mapped_fraction": "0.99",
                        "properly_paired_fraction": "0.97",
                        "duplicate_fraction": "0.08",
                    },
                )
                write_metric_value(
                    project / "results/qc/wgs" / f"{sample}.het_snp_qc.tsv",
                    {"sample": sample, "heterozygous_snv_count": str(final)},
                )
                write_metric_value(
                    project / "results/qc/wgs" / f"{sample}.het_dp_retention_qc.tsv",
                    {
                        "sample": sample,
                        "minimum_dp": "10",
                        "heterozygous_snv_before_dp": str(final + 10),
                        "heterozygous_snv_after_dp": str(final),
                        "dp_retention_fraction": f"{final / (final + 10):.6f}",
                    },
                )
                phasing = project / "results/phasing" / f"{sample}.phased.wgs.honest.summary.tsv"
                phasing.parent.mkdir(parents=True, exist_ok=True)
                with phasing.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=[
                            "chrom", "variants", "phased", "unphased",
                            "phased_fraction", "singletons", "multi_variant_blocks",
                            "block_rows_including_singletons", "variants_per_block_n50",
                            "block_span_n50_bp", "largest_block_variants",
                            "largest_block_span_bp", "sum_block_spans_bp",
                        ],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow(
                        {
                            "chrom": "ALL", "variants": final, "phased": final - 20,
                            "unphased": 20, "phased_fraction": f"{(final - 20) / final:.6f}",
                            "singletons": 10, "multi_variant_blocks": 100,
                            "block_rows_including_singletons": 110,
                            "variants_per_block_n50": 60, "block_span_n50_bp": 5000,
                            "largest_block_variants": 200, "largest_block_span_bp": 20000,
                            "sum_block_spans_bp": 500000,
                        }
                    )

            output = project / "results/qc/stage1_cohort_qc.tsv"
            subprocess.run(
                [
                    "python3", str(SCRIPT), "--project-dir", str(project),
                    "--sample-sheet", str(sample_sheet), "--output", str(output),
                ],
                check=True,
            )
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual([row["sample"] for row in rows], ["S1", "S2"])
        self.assertEqual(rows[0]["final_heterozygous_snv_count"], "1001")
        self.assertEqual(rows[1]["phased_snv_count"], "982")
        self.assertEqual(len(rows[0]), 14)


if __name__ == "__main__":
    unittest.main()
