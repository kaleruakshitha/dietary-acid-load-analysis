import math
import unittest
from src.features import calculate_egfr_2021, calculate_pral

class FeatureTests(unittest.TestCase):
    def test_pral_matches_hand_calculation(self):
        observed=calculate_pral(80,1200,3000,350,900)
        expected=0.49*80+0.037*1200-0.021*3000-0.026*350-0.013*900
        self.assertAlmostEqual(observed,expected)
    def test_pral_rejects_negative_nutrient(self):
        with self.assertRaises(ValueError): calculate_pral(80,1200,-1,350,900)
    def test_egfr_matches_published_equation_structure(self):
        observed=calculate_egfr_2021(1.0,50,"female")
        ratio=1.0/0.7
        expected=142*min(ratio,1)**-0.241*max(ratio,1)**-1.2*0.9938**50*1.012
        self.assertAlmostEqual(observed,expected)
        self.assertTrue(math.isfinite(observed))
    def test_egfr_rejects_invalid_sex(self):
        with self.assertRaises(ValueError): calculate_egfr_2021(1.0,50,"unknown")
    def test_egfr_rejects_nonpositive_creatinine(self):
        with self.assertRaises(ValueError): calculate_egfr_2021(0,50,"male")
if __name__=="__main__": unittest.main()
