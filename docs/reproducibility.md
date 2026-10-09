# Reproducibility and verification

## Files now available
The original scripts for NHANES downloads, participant-cohort construction, recall reliability, survey-weighted regression (R), food-category processing and food-source validation are published in `scripts/`. Supporting Python functions are in `src/` and the original tests in `tests/`.

## What was actually verified
- **Python unit tests:** All 14 tests in the supplied original project ZIP passed locally with `python -m unittest discover -s tests -v` (October 2026).
- **Reported results:** Saved aggregate output tables in the supplied project package are consistent, after rounding, with the reported two-day-mean PRAL fully adjusted coefficient (−0.0316), confidence interval (−0.0935 to 0.0303) and p-value (0.293), n=3,321.
- **Not yet verified:** A fresh end-to-end run from original raw NHANES and USDA files, including execution of the R `survey` regression. Rscript was not available in the review environment. Passing unit tests and reading existing output tables do not independently reproduce the regression estimates.

## Run the tests
Create a Python environment and install the dependencies from `requirements.txt`. Then run from the repository root:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

## Reproduce the analysis using your original saved data
Follow the sequence below with the precise NHANES August 2021–August 2023 and USDA/WWEIA input files used in the project. This is a **run plan**, not a newly validated fresh pipeline run.

```bash
python scripts/download_nhanes.py
python scripts/build_cohort.py
python scripts/analyze_reliability.py
Rscript scripts/install_r_dependencies.R
Rscript scripts/analyze_regression.R
python scripts/validate_food_sources.py
python scripts/analyze_food_categories.py
```

The USDA category crosswalk is a separate required input and must be placed under `data/raw/usda/` with the filename expected by `validate_food_sources.py`. The raw NHANES and participant-level output folders are intentionally not published. Source files downloaded today can differ from the project's September 2026 snapshot and source websites may change.

## Data handling
Do not publish individual participant records, processed participant-level cohorts, or data containing direct identifiers. The portfolio only needs source code, reproducible commands, non-disclosive aggregate results, and documented methods.

## Interpretation
NHANES is observational and uses a complex sampling design. The reported small adjusted coefficient is statistically uncertain; no causal or clinical effect is established by this analysis.
