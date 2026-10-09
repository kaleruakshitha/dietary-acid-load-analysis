import unittest

import numpy as np

from src.reliability import (
    cohen_kappa_from_matrix,
    icc_two_way_random_absolute,
    weighted_mean,
    weighted_pearson,
    weighted_quantile,
    weighted_standard_deviation,
)


class ReliabilityTests(unittest.TestCase):
    def test_weighted_summary_functions(self):
        values = np.array([1.0, 2.0, 3.0, 4.0])
        weights = np.array([1.0, 1.0, 1.0, 5.0])
        self.assertAlmostEqual(weighted_mean(values, weights), 3.25)
        self.assertAlmostEqual(
            weighted_standard_deviation(values, weights),
            np.sqrt(np.average((values - 3.25) ** 2, weights=weights)),
        )
        np.testing.assert_array_equal(
            weighted_quantile(values, weights, [0.25, 0.5, 0.75]),
            np.array([2.0, 4.0, 4.0]),
        )

    def test_weighted_pearson_is_one_for_linear_match(self):
        first = np.array([1.0, 2.0, 4.0, 8.0])
        second = 3 * first + 7
        weights = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(weighted_pearson(first, second, weights), 1.0)

    def test_icc_distinguishes_consistency_from_absolute_agreement(self):
        first = np.arange(1.0, 11.0)
        measurements = np.column_stack([first, first + 5])
        absolute, consistency = icc_two_way_random_absolute(measurements)
        self.assertLess(absolute, consistency)
        self.assertAlmostEqual(consistency, 1.0)

    def test_kappa_is_one_for_perfect_agreement(self):
        matrix = np.diag([12.0, 15.0, 11.0, 14.0])
        for scheme in ("unweighted", "linear", "quadratic"):
            self.assertAlmostEqual(cohen_kappa_from_matrix(matrix, scheme), 1.0)

    def test_rejects_nonpositive_weights(self):
        with self.assertRaises(ValueError):
            weighted_mean([1.0, 2.0], [1.0, 0.0])


if __name__ == "__main__":
    unittest.main()