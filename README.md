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

The original project report records 14 passing validation tests. This refers to the reported working analysis, **not a newly verified test run from the files currently in this GitHub repository**.

## My motivation

I studied Nutrition Science and wanted to use analytical tools to investigate a question with real public-health relevance. The project helped me connect subject knowledge with data quality checks, meaningful statistical comparisons, and communicating results without overstating what observational data can show. It supports my transition toward graduate study in analytics.

## Start here: reviewer navigation

| What to inspect | File |
| --- | --- |
| Key statistical findings | [Results summary](docs/results.md) |
| Food-category interpretation | [Food-category summary](docs/food-category-summary.md) |
| Python analytical feature engineering | [Original PRAL/eGFR functions](src/features.py) |
| Python unit testing | [Feature tests](tests/test_features.py) |
| SQL data integration | [Core participant SQL](sql/01_core_participant.sql) |
| Reproducibility and limitations | [Reproducibility guide](docs/reproducibility.md) |

The original full NHANES pipeline (including the complete R regression and food-category analysis scripts) is **still being transferred** from the provided ZIP. The short `src/pral.py` is an earlier example and is not the main original implementation; use `src/features.py` for the recovered project feature code. Reviewer-facing summaries distinguish reported results from independently reproduced runs.

## Repository contents

- [Reported results and limitations](docs/results.md)
- [PRAL calculation reference implementation](src/pral.py)
- [Calculator tests](tests/test_pral.py)
- [Data handling and workflow notes](docs/reproducibility.md)

**Reproducibility status:** The accompanying calculation implementation and tests were reconstructed from the final report for transparent demonstration. The **original full NHANES extraction, merging, food mapping, and R modeling scripts are not yet in this repository**. The published study estimates should be treated as *reported results*, not as reproduced by the small reference implementation alone.

## References

- [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/)
- Remer T, Manz F. Potential renal acid load of foods and its influence on urine pH. *Journal of the American Dietetic Association* (1995).
- The project's September 2026 report, *Dietary Acid Load and Serum Bicarbonate*.

