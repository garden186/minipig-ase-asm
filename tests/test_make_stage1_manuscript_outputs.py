#!/usr/bin/env python3

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "stage1" / "make_stage1_manuscript_outputs.py"
FIXTURE = ROOT / "tests" / "fixtures" / "stage1_cohort_qc.tsv"


class MakeStage1ManuscriptOutputsTest(unittest.TestCase):
    def test_writes_vector_raster_and_table_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory)
            subprocess.run(
                [
                    "python3", str(SCRIPT), "--cohort-qc", str(FIXTURE),
                    "--output-dir", str(output),
                ],
                check=True,
            )
            expected = [
                "Supplementary_Figure_S1_Stage1_workflow.pdf",
                "Supplementary_Figure_S1_Stage1_workflow.svg",
                "Supplementary_Figure_S1_Stage1_workflow.png",
                "Supplementary_Figure_S2_Stage1_cohort_QC.pdf",
                "Supplementary_Figure_S2_Stage1_cohort_QC.svg",
                "Supplementary_Figure_S2_Stage1_cohort_QC.png",
                "Supplementary_Table_S1_Stage1_cohort_QC.tsv",
                "Supplementary_Table_S1_metric_dictionary.tsv",
                "stage1_cohort_descriptive_statistics.tsv",
                "stage1_supplementary_legends.md",
            ]
            for name in expected:
                path = output / name
                self.assertTrue(path.is_file(), name)
                self.assertGreater(path.stat().st_size, 0, name)


if __name__ == "__main__":
    unittest.main()
