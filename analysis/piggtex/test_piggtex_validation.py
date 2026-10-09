#!/usr/bin/env python3

import math
import unittest

from piggtex_validation import (
    bh_values,
    fisher_upper_tail,
    hypergeom_pmf,
    stratified_exact_tail,
)


class StatisticsTests(unittest.TestCase):
    def test_hypergeom_normalizes(self):
        distribution = hypergeom_pmf(10, 4, 3)
        self.assertAlmostEqual(sum(distribution.values()), 1.0, places=12)

    def test_fisher_known_perfect_split(self):
        # R: fisher.test(matrix(c(5,0,0,5),2), alternative="greater")
        self.assertAlmostEqual(fisher_upper_tail(5, 0, 0, 5), 1 / 252, places=12)

    def test_stratified_single_equals_fisher(self):
        table = (4, 1, 2, 7)
        self.assertAlmostEqual(
            stratified_exact_tail([table]), fisher_upper_tail(*table), places=12
        )

    def test_bh(self):
        observed = bh_values([0.01, 0.04, 0.03, 0.002])
        expected = [0.02, 0.04, 0.04, 0.008]
        for left, right in zip(observed, expected):
            self.assertTrue(math.isclose(left, right, abs_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
