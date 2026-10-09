# Dietary Acid Load and Serum Bicarbonate

**Measurement reliability, survey-weighted regression, and food-category contributions in NHANES August 2021–August 2023**

Independent public-health analytics project | Akshitha Kaleru | September 2026

## Research question

How reliably does potential renal acid load (PRAL) estimated from two dietary recalls characterize intake, and how is PRAL associated with serum bicarbonate after adjustment for other factors?

## Findings from the completed project report

| Measure | Result |
| --- | --- |
| NHANES participants released | 11,933 |
| Adults age 20+ | 7,809 |
| Core analytic cohort | 3,848 |
| Complete-case regression cohort | 3,321 |
| Day-to-day absolute-agreement ICC | 0.336 (95% bootstrap CI: 0.297–0.374) |
| Participants who changed PRAL quartiles | 66.7% |
| Fully adjusted bicarbonate estimate per 10 mEq/day increase in two-day mean PRAL | −0.0316 mmol/L |
| 95% confidence interval | −0.0935 to 0.0303 |
| p-value | 0.293 |
| Food records mapped | 120,956 across 163 food categories |

**Interpretation:** Day-to-day PRAL estimates varied substantially. The fully adjusted association with serum bicarbonate was small and statistically uncertain; the confidence interval included zero. This observational analysis **does not establish a causal effect** of dietary acid load on bicarbonate.

## Methods

- **Data:** NHANES August 2021–August 2023 public-use dietary, demographic, laboratory, examination, kidney, smoking, and diabetes information; cycle-matched USDA FNDDS/WWEIA food categories.
- **Exposure:** PRAL calculated separately for each 24-hour dietary recall and averaged across the two days.
- **Calculation:** PRAL (mEq/day) = 0.49 × protein (g) + 0.037 × phosphorus (mg) − 0.021 × potassium (mg) − 0.026 × magnesium (mg) − 0.013 × calcium (mg).
- **Reliability:** Pearson/Spearman correlations, absolute-agreement ICC, Bland–Altman agreement, and changes between day-specific quartiles.
- **Model:** Survey-weighted linear regression of serum bicarbonate, adjusting for demographic, socioeconomic, behavioral, dietary energy, and kidney-function measures. Used NHANES two-day dietary weights, strata and PSUs.
- **Food contributions:** Summarized food-level PRAL components by USDA/WWEIA categories.

## Programming and analytical skills

| Technology | Role |
| --- | --- |
| Python | Data cleaning, nutrient calculations, checks, and food-category summaries |
| R (`survey` package v4.5) | Complex-survey regression and adjusted estimates |
| SQL | Structured data querying in the broader project workflow |
| Statistics | Reliability assessment, interval estimation, covariate adjustment, and interpretation |

The original project reported 14 passing validation tests. I reran the original ZIP's 14 Python unit tests in a local Python environment and all passed; a fresh full-data end-to-end run of the complete pipeline, including R regression, has not been performed.

## My motivation

I studied Nutrition Science and wanted to use analytical tools to investigate a question with real public-health relevance. The project helped me connect subject knowledge with data quality checks, meaningful statistical comparisons, and communicating results without overstating what observational data can show. It supports my transition toward graduate study in analytics.

## Start here: reviewer navigation

| What to inspect | File |
| --- | --- |
| Key statistical findings | [Results summary](docs/results.md) |
| Food-category interpretation | [Food-category summary](docs/food-category-summary.md) |
| Python analytical feature engineering | [PRAL/eGFR functions](src/features.py) |
| Full cohort build | [build_cohort.py](scripts/build_cohort.py) |
| Recall-day reliability analysis | [analyze_reliability.py](scripts/analyze_reliability.py) |
| Survey-weighted R regression | [analyze_regression.R](scripts/analyze_regression.R) |
| Food-category analysis | [analyze_food_categories.py](scripts/analyze_food_categories.py) |
| Python dependencies | [requirements.txt](requirements.txt) |
| Python unit testing | [Feature tests](tests/test_features.py) |
| SQL data integration | [Core participant SQL](sql/01_core_participant.sql) |
| Reproducibility and limitations | [Reproducibility guide](docs/reproducibility.md) |

The original cohort, reliability, regression (R), and food-category analysis scripts are now in `scripts/`, with their supporting calculations in `src/` and original unit tests in `tests/`. The redundant example calculator and its separate example tests have been removed; the core original feature code and tests are under `src/` and `tests/`. The original data and full end-to-end regression results have not been rerun here.

## Repository contents

- [Reported results and limitations](docs/results.md)
- [Original PRAL and eGFR feature implementation](src/features.py)
- [Original unit tests](tests/test_features.py)
- [Survey-weighted regression script](scripts/analyze_regression.R)
- [Data handling and workflow notes](docs/reproducibility.md)

**Reproducibility status:** The original Python cohort, reliability, food-category, validation and R regression scripts have been restored. The complete analysis was not rerun in this environment: the saved source-data snapshot is not published and R was unavailable. The values above are checked against saved aggregate project output tables, not independently re-estimated from the raw source files in this review.

## References

- [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/)
- Remer T, Manz F. Potential renal acid load of foods and its influence on urine pH. *Journal of the American Dietetic Association* (1995).
- The project's September 2026 report, *Dietary Acid Load and Serum Bicarbonate*.

