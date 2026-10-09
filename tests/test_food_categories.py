import unittest

import numpy as np
import pandas as pd

from scripts.analyze_food_categories import (
    add_pral_components,
    assign_weighted_quartiles,
)


class FoodCategoryTests(unittest.TestCase):
    def test_pral_components_reproduce_formula(self):
        foods = pd.DataFrame(
            {
                "protein_g": [80.0],
                "phosphorus_mg": [1200.0],
                "potassium_mg": [3000.0],
                "magnesium_mg": [350.0],
                "calcium_mg": [900.0],
            }
        )
        observed = add_pral_components(foods).iloc[0]
        expected_positive = 0.49 * 80 + 0.037 * 1200
        expected_negative = -(0.021 * 3000 + 0.026 * 350 + 0.013 * 900)
        self.assertAlmostEqual(observed["positive_component"], expected_positive)
        self.assertAlmostEqual(observed["negative_component"], expected_negative)
        self.assertAlmostEqual(
            observed["net_pral"], expected_positive + expected_negative
        )

    def test_structural_missing_nutrient_is_treated_as_zero(self):
        foods = pd.DataFrame(
            {
                "protein_g": [np.nan],
                "phosphorus_mg": [10.0],
                "potassium_mg": [20.0],
                "magnesium_mg": [np.nan],
                "calcium_mg": [30.0],
            }
        )
        observed = add_pral_components(foods).iloc[0]
        self.assertAlmostEqual(
            observed["net_pral"], 0.037 * 10 - 0.021 * 20 - 0.013 * 30
        )

    def test_pral_components_reject_negative_nutrients(self):
        foods = pd.DataFrame(
            {
                "protein_g": [-1.0],
                "phosphorus_mg": [10.0],
                "potassium_mg": [20.0],
                "magnesium_mg": [5.0],
                "calcium_mg": [30.0],
            }
        )
        with self.assertRaises(ValueError):
            add_pral_components(foods)

    def test_weighted_quartiles_are_valid(self):
        values = np.arange(1.0, 9.0)
        weights = np.ones(8)
        quartiles, cutpoints = assign_weighted_quartiles(values, weights)
        np.testing.assert_array_equal(cutpoints, [2.0, 4.0, 6.0])
        np.testing.assert_array_equal(quartiles, [1, 1, 2, 2, 3, 3, 4, 4])


if __name__ == "__main__":
    unittest.main()
