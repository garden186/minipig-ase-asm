#!/usr/bin/env python3
"""Focused tests for the exclusive pair-support hierarchy."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("summarize_pair_support.py")
SPEC = importlib.util.spec_from_file_location("pair_support", MODULE_PATH)
assert SPEC and SPEC.loader
pair_support = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pair_support)


class SupportLogicTest(unittest.TestCase):
    def test_strongest_evidence_wins(self):
        self.assertEqual(
            pair_support.best_support_class(True, True, True, True),
            "DIRECTION_CONCORDANT_PIGGTEX_EXACT_SITE",
        )

    def test_exact_support_without_direction(self):
        self.assertEqual(
            pair_support.best_support_class(True, True, True, False),
            "PIGGTEX_EXACT_SITE_WITHOUT_CONCORDANT_DIRECTION",
        )

    def test_no_assigned_site(self):
        self.assertEqual(
            pair_support.best_support_class(False, False, False, False),
            "NO_COORDINATE_ASSIGNED_INDEPENDENT_RECURRENT_SNP",
        )


if __name__ == "__main__":
    unittest.main()
