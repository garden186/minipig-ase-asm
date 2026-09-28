#!/usr/bin/env python3
"""Focused tests for full-universe eQTL input preparation."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("prepare_full_universe.py")
SPEC = importlib.util.spec_from_file_location("preflight", MODULE_PATH)
assert SPEC and SPEC.loader
preflight = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(preflight)


class HelperTest(unittest.TestCase):
    def test_unordered_variant_key(self):
        self.assertEqual(preflight.variant_key("chr9", 100, "T", "A"), "9_100_A_T")
        self.assertEqual(preflight.variant_key_from_id("9_100_T_A"), "9_100_A_T")

    def test_versioned_gene_id(self):
        self.assertEqual(preflight.clean_gene("ENSSSCG0001.5"), "ENSSSCG0001")


if __name__ == "__main__":
    unittest.main()
