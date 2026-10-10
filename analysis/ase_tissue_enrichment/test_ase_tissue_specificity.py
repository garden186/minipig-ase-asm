#!/usr/bin/env python3
"""Deterministic unit tests for the ASE tissue-specificity engine."""

from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("ase_tissue_specificity.py")
SPEC = importlib.util.spec_from_file_location("ase_tissue_specificity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class SpecificityTests(unittest.TestCase):
    def test_equal_weight_conditional_marginals(self):
        tissues = ["A", "B", "C"]
        one = MODULE.conditional_marginals(tissues, [1.0, 1.0, 1.0], 1)
        two = MODULE.conditional_marginals(tissues, [1.0, 1.0, 1.0], 2)
        for tissue in tissues:
            self.assertAlmostEqual(one[tissue], 1.0 / 3.0)
            self.assertAlmostEqual(two[tissue], 2.0 / 3.0)

    def test_weighted_two_tissue_marginal(self):
        observed = MODULE.conditional_marginals(["A", "B"], [3.0, 1.0], 1)
        self.assertAlmostEqual(observed["A"], 0.75)
        self.assertAlmostEqual(observed["B"], 0.25)

    def test_elementary_symmetric(self):
        self.assertAlmostEqual(MODULE.elementary_symmetric([2.0, 3.0, 5.0], 1), 10.0)
        self.assertAlmostEqual(MODULE.elementary_symmetric([2.0, 3.0, 5.0], 2), 31.0)
        self.assertAlmostEqual(MODULE.elementary_symmetric([2.0, 3.0, 5.0], 3), 30.0)

    def test_poisson_binomial_tail(self):
        self.assertAlmostEqual(
            MODULE.poisson_binomial_tail([1.0 / 3.0] * 3, 3),
            1.0 / 27.0,
        )

    def test_profile_odds_ratio_is_one_for_balanced_outcome(self):
        estimate, lower, upper = MODULE.profile_log_odds_ratio(
            [1, 0, 1, 0], [0.5, 0.5, 0.5, 0.5]
        )
        self.assertAlmostEqual(estimate, 1.0, places=8)
        self.assertLess(lower, 1.0)
        self.assertGreater(upper, 1.0)

    def test_profile_boundary(self):
        estimate, lower, upper = MODULE.profile_log_odds_ratio(
            [1, 1, 1], [1.0 / 3.0] * 3
        )
        self.assertTrue(math.isinf(estimate))
        self.assertGreater(lower, 0.0)
        self.assertTrue(math.isinf(upper))

    def test_bh(self):
        observed = MODULE.bh_values([0.01, 0.04, 0.03])
        self.assertEqual([round(value, 8) for value in observed], [0.03, 0.04, 0.04])

    def test_tissue_alias(self):
        self.assertEqual(MODULE.normalize_tissue("Tenderlo"), "Tenderloin")
        self.assertEqual(MODULE.normalize_tissue("L-blood"), "Blood")


if __name__ == "__main__":
    unittest.main(verbosity=2)
