"""Reference PRAL calculator reconstructed from the published project report.

This standalone function does not reproduce the full NHANES analysis.
All quantities are daily totals in the specified units.
"""
from math import isfinite


def pral(protein_g, phosphorus_mg, potassium_mg, magnesium_mg, calcium_mg):
    """Estimate potential renal acid load (mEq/day).

    Remer and Manz (1995) nutrient-based equation.
    """
    values = (protein_g, phosphorus_mg, potassium_mg, magnesium_mg, calcium_mg)
    if any(not isinstance(x, (int, float)) or not isfinite(x) or x < 0 for x in values):
        raise ValueError("All nutrient inputs must be finite, nonnegative numbers")
    return (0.49 * protein_g + 0.037 * phosphorus_mg
            - 0.021 * potassium_mg - 0.026 * magnesium_mg
            - 0.013 * calcium_mg)
