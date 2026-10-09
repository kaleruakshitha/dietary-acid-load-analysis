"""Unit checks for the standalone reference calculator only.

These are new illustrative tests and not the original project's 14 tests.
"""
import unittest
from src.pral import pral


class TestPRAL(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(pral(0, 0, 0, 0, 0), 0)

    def test_protein_contribution(self):
        self.assertAlmostEqual(pral(10, 0, 0, 0, 0), 4.9)

    def test_component_sum(self):
        expected = 0.49 * 70 + 0.037 * 1000 - 0.021 * 2500 - 0.026 * 300 - 0.013 * 800
        self.assertAlmostEqual(pral(70, 1000, 2500, 300, 800), expected)

    def test_negative_input_rejected(self):
        with self.assertRaises(ValueError):
            pral(-1, 0, 0, 0, 0)


if __name__ == "__main__":
    unittest.main()
