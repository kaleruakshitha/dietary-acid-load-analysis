-- Original project: one row per participant after NHANES XPT-to-Parquet conversion.
CREATE OR REPLACE VIEW core_participant AS
SELECT d.SEQN AS seqn,d.RIDAGEYR AS age_years,d.RIAGENDR AS sex_code,
d.RIDRETH3 AS race_ethnicity_code,d.DMDEDUC2 AS education_code,
d.INDFMPIR AS income_poverty_ratio,d.SDMVSTRA AS survey_stratum,
d.SDMVPSU AS survey_psu,r1.WTDR2D AS two_day_diet_weight,
r1.DR1DRSTZ AS day1_recall_status,r2.DR2DRSTZ AS day2_recall_status,
r1.DR1TKCAL AS energy_day1,r2.DR2TKCAL AS energy_day2,
r1.DR1TPROT AS protein_day1,r1.DR1TPHOS AS phosphorus_day1,
r1.DR1TPOTA AS potassium_day1,r1.DR1TMAGN AS magnesium_day1,
r1.DR1TCALC AS calcium_day1,r2.DR2TPROT AS protein_day2,
r2.DR2TPHOS AS phosphorus_day2,r2.DR2TPOTA AS potassium_day2,
r2.DR2TMAGN AS magnesium_day2,r2.DR2TCALC AS calcium_day2,
b.LBXSC3SI AS bicarbonate,b.LBXSCR AS serum_creatinine,
x.BMXBMI AS bmi,x.BMXWT AS weight_kg,a.URDACT AS uacr
FROM read_parquet('data/interim/DEMO_L.parquet') d
INNER JOIN read_parquet('data/interim/DR1TOT_L.parquet') r1 USING (SEQN)
INNER JOIN read_parquet('data/interim/DR2TOT_L.parquet') r2 USING (SEQN)
INNER JOIN read_parquet('data/interim/BIOPRO_L.parquet') b USING (SEQN)
LEFT JOIN read_parquet('data/interim/BMX_L.parquet') x USING (SEQN)
LEFT JOIN read_parquet('data/interim/ALB_CR_L.parquet') a USING (SEQN);
