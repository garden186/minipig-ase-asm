#!/usr/bin/env python3

from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "stage1" / "summarize_whatshap.py"
FIXTURES = ROOT / "tests" / "fixtures"


class SummarizeWhatsHapTest(unittest.TestCase):
    def test_named_column_summary_and_n50(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "summary.tsv"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--stats",
                    str(FIXTURES / "whatshap.stats.tsv"),
                    "--blocks",
                    str(FIXTURES / "whatshap.blocks.tsv"),
                    "--output",
                    str(output),
                    "--autosomes",
                    "1",
                    "2",
                ],
                check=True,
            )
            with output.open(encoding="utf-8", newline="") as handle:
                rows = {row["chrom"]: row for row in csv.DictReader(handle, delimiter="\t")}

        self.assertEqual(rows["1"]["block_rows_including_singletons"], "2")
        self.assertEqual(rows["1"]["multi_variant_blocks"], "1")
        self.assertEqual(rows["1"]["variants_per_block_n50"], "2")
        self.assertEqual(rows["1"]["block_span_n50_bp"], "10")
        self.assertEqual(rows["ALL"]["variants_per_block_n50"], "3")
        self.assertEqual(rows["ALL"]["block_span_n50_bp"], "101")
        self.assertEqual(rows["ALL"]["sum_block_spans_bp"], "112")


if __name__ == "__main__":
    unittest.main()
