# Reproducing the analysis

The original full NHANES/USDA processing and R model source files are not yet in this GitHub repository. The repository currently provides a transparent nutrient-equation reference and reported aggregate results, **not a full rerunnable pipeline**.

## Reference calculator

Using Python 3 from the repository root:

```bash
python -m unittest discover -s tests -v
```

Example:

```python
from src.pral import pral
print(pral(70, 1000, 2500, 300, 800))
```

Inputs must be daily nutrient totals in grams (protein) and milligrams (minerals).

## Full analysis requirements

The original project used NHANES August 2021–August 2023 data, two-day dietary sampling weights (`WTDR2D`), strata and PSUs, and an R `survey` v4.5 analysis with survey-weighted regression. Reproduction requires reobtaining the original public-use source files, documenting complete variable definitions and exclusion criteria, retrieving cycle-matched USDA/FNDDS/WWEIA mapping files, and restoring the original source code and environment. Do not assume unweighted calculations reproduce survey-weighted results.

## Data safety

Do not commit local participant-level extracts or intermediary tables. Public inputs should be downloaded from source agencies, and only reviewed aggregate outputs should be published.

## Source

*Dietary Acid Load and Serum Bicarbonate*, independent project report (September 2026).
