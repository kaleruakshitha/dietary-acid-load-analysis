# Reported project results

These figures are transcribed from *Dietary Acid Load and Serum Bicarbonate: Measurement reliability and food-category contributions in NHANES 2021–2023* (September 2026). They have **not yet been reproduced from original analysis scripts in this repository**.

## Sample construction

| Stage | N |
| --- | ---: |
| Released participants | 11,933 |
| Adults 20+ | 7,809 |
| Core analytic cohort | 3,848 |
| Regression complete cases | 3,321 |

## Dietary recall reliability

- Pearson correlation: 0.336; Spearman correlation: 0.326.
- Absolute-agreement ICC[A,1]: 0.336 (95% bootstrap CI 0.297–0.374).
- Day 1 minus Day 2 mean difference: −0.53 mEq/day (95% CI −1.38 to 0.32).
- Bland–Altman 95% limits of agreement: −53.24 to 52.19 mEq/day.
- Weighted quartile movement: same quartile 33.3%, one quartile 40.4%, two quartiles 19.7%, lowest-to-highest 6.6%. A total of 66.7% moved at least one quartile.

## Survey-weighted linear regression

Fully adjusted regression sample: n=3,321; estimates per 10 mEq/day higher PRAL, mmol/L serum bicarbonate.

| PRAL definition | Coefficient | 95% CI | p |
| --- | ---: | --- | ---: |
| Day 1 | −0.0159 | −0.0628 to 0.0310 | 0.480 |
| Day 2 | −0.0223 | −0.0809 to 0.0363 | 0.430 |
| Mean of two days | −0.0316 | −0.0935 to 0.0303 | 0.293 |

Reported adjusted model sequence for mean PRAL:

| Model | Coefficient (95% CI) | p |
| --- | --- | ---: |
| PRAL alone | −0.0595 (−0.1177 to −0.0013) | 0.046 |
| + Age and sex | −0.0552 (−0.1061 to −0.0043) | 0.035 |
| + Race/ethnicity, education, income | −0.0448 (−0.0886 to −0.0011) | 0.045 |
| + BMI and mean energy | −0.0376 (−0.0898 to 0.0147) | 0.147 |
| + Smoking and diabetes | −0.0351 (−0.0913 to 0.0210) | 0.203 |
| + Kidney function | −0.0316 (−0.0935 to 0.0303) | 0.293 |

## Food sources

The report describes 120,956 matched food records across 163 USDA/WWEIA categories with no unmatched records. Chicken, eggs and cheese were among the largest positive net contributors and coffee, bananas and bottled water among the largest negative contributors. Category values are components of calculated PRAL, not causal effects on bicarbonate.

**Limitations:** short-term recall variability, possible confounding, model specification, and cross-sectional design. None of the fully adjusted PRAL comparisons excludes a zero association.
