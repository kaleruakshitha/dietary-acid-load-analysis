"""Build and validate the Milestone 2 NHANES participant cohort."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.features import calculate_egfr_2021, calculate_pral


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "nhanes"
USDA_RAW = ROOT / "data" / "raw" / "usda"
PROCESSED = ROOT / "data" / "processed"
DOCS = ROOT / "docs"

REQUIRED_FILES = [
    "DEMO_L",
    "DR1TOT_L",
    "DR2TOT_L",
    "BIOPRO_L",
    "BMX_L",
    "ALB_CR_L",
    "SMQ_L",
    "DIQ_L",
    "KIQ_U_L",
]

PRAL_DAY1 = ["DR1TPROT", "DR1TPHOS", "DR1TPOTA", "DR1TMAGN", "DR1TCALC"]
PRAL_DAY2 = ["DR2TPROT", "DR2TPHOS", "DR2TPOTA", "DR2TMAGN", "DR2TCALC"]


def load_xpt(stem: str) -> pd.DataFrame:
    path = RAW / f"{stem}.xpt"
    if not path.exists():
        raise FileNotFoundError(f"Missing required source file: {path}")
    frame = pd.read_sas(path, format="xport")
    numeric_columns = frame.select_dtypes(include="number").columns
    for column in numeric_columns:
        values = frame[column]
        sas_missing = values.ne(0) & values.abs().lt(1e-50)
        if sas_missing.any():
            frame.loc[sas_missing, column] = np.nan
    frame["SEQN"] = frame["SEQN"].astype("int64")
    if frame["SEQN"].duplicated().any():
        raise ValueError(f"{stem} contains duplicate participant identifiers")
    return frame


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_smoking(frame: pd.DataFrame) -> pd.Series:
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    result.loc[frame["SMQ020"].eq(2)] = "never"
    result.loc[frame["SMQ020"].eq(1) & frame["SMQ040"].eq(3)] = "former"
    result.loc[
        frame["SMQ020"].eq(1) & frame["SMQ040"].isin([1, 2])
    ] = "current"
    return result


def classify_diabetes(frame: pd.DataFrame) -> pd.Series:
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    result.loc[frame["DIQ010"].eq(1)] = "yes"
    result.loc[frame["DIQ010"].eq(2)] = "no"
    result.loc[frame["DIQ010"].eq(3)] = "borderline"
    return result


def build_joined_source(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frame = frames["DEMO_L"].copy()
    selections = {
        "DR1TOT_L": None,
        "DR2TOT_L": None,
        "BIOPRO_L": ["SEQN", "LBXSC3SI", "LBXSCR"],
        "BMX_L": ["SEQN", "BMXBMI", "BMXWT", "BMXHT"],
        "ALB_CR_L": ["SEQN", "URDACT"],
        "SMQ_L": ["SEQN", "SMQ020", "SMQ040"],
        "DIQ_L": ["SEQN", "DIQ010"],
        "KIQ_U_L": ["SEQN", "KIQ022", "KIQ025"],
    }
    for name, columns in selections.items():
        right = frames[name] if columns is None else frames[name][columns]
        frame = frame.merge(
            right,
            on="SEQN",
            how="left",
            validate="one_to_one",
            suffixes=("", f"_{name.lower()}"),
        )
    return frame


def make_attrition(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    records = [
        {
            "stage": 0,
            "rule": "Total released NHANES participants",
            "removed_at_stage": 0,
            "remaining": len(frame),
        }
    ]
    mask = pd.Series(True, index=frame.index)

    def apply_rule(label: str, condition: pd.Series) -> None:
        nonlocal mask
        before = int(mask.sum())
        mask = mask & condition.fillna(False)
        after = int(mask.sum())
        records.append(
            {
                "stage": len(records),
                "rule": label,
                "removed_at_stage": before - after,
                "remaining": after,
            }
        )

    apply_rule("Adults age 20 years or older", frame["RIDAGEYR"].ge(20))
    apply_rule("Reliable Day 1 dietary recall", frame["DR1DRSTZ"].eq(1))
    apply_rule("Reliable Day 2 dietary recall", frame["DR2DRSTZ"].eq(1))
    apply_rule("Positive two-day dietary weight", frame["WTDR2D"].gt(0))
    apply_rule(
        "All Day 1 and Day 2 PRAL nutrients present",
        frame[PRAL_DAY1 + PRAL_DAY2].notna().all(axis=1),
    )
    apply_rule("Serum bicarbonate present", frame["LBXSC3SI"].notna())
    apply_rule("Serum creatinine present", frame["LBXSCR"].notna())
    apply_rule(
        "Survey strata and PSU present",
        frame[["SDMVSTRA", "SDMVPSU"]].notna().all(axis=1),
    )
    apply_rule("Not identified as pregnant", ~frame["RIDEXPRG"].eq(1))
    apply_rule("No dialysis in the past 12 months", ~frame["KIQ025"].eq(1))
    core_mask = mask.copy()

    apply_rule("BMI present", frame["BMXBMI"].notna())
    apply_rule("Income-to-poverty ratio present", frame["INDFMPIR"].notna())
    apply_rule("Adult education present", frame["DMDEDUC2"].notna())
    apply_rule("Smoking status classifiable", frame["smoking_status"].notna())
    apply_rule("Diabetes status classifiable", frame["diabetes_status"].notna())
    fully_adjusted_without_uacr = mask.copy()
    apply_rule("Urine albumin-creatinine ratio present", frame["URDACT"].notna())

    attrition = pd.DataFrame.from_records(records)
    return attrition, core_mask, fully_adjusted_without_uacr


def add_features(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["pral_day1"] = result.apply(
        lambda row: calculate_pral(
            row["DR1TPROT"],
            row["DR1TPHOS"],
            row["DR1TPOTA"],
            row["DR1TMAGN"],
            row["DR1TCALC"],
        ),
        axis=1,
    )
    result["pral_day2"] = result.apply(
        lambda row: calculate_pral(
            row["DR2TPROT"],
            row["DR2TPHOS"],
            row["DR2TPOTA"],
            row["DR2TMAGN"],
            row["DR2TCALC"],
        ),
        axis=1,
    )
    result["pral_mean"] = result[["pral_day1", "pral_day2"]].mean(axis=1)
    result["pral_difference"] = result["pral_day1"] - result["pral_day2"]
    result["pral_absolute_difference"] = result["pral_difference"].abs()
    result["mean_energy"] = result[["DR1TKCAL", "DR2TKCAL"]].mean(axis=1)
    result["sex_label"] = result["RIAGENDR"].map({1.0: "male", 2.0: "female"})
    result["egfr"] = result.apply(
        lambda row: calculate_egfr_2021(
            row["LBXSCR"], row["RIDAGEYR"], row["sex_label"]
        ),
        axis=1,
    )
    result["protein_g_per_kg"] = (
        result[["DR1TPROT", "DR2TPROT"]].mean(axis=1) / result["BMXWT"]
    )
    return result


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    frames = {name: load_xpt(name) for name in REQUIRED_FILES}
    joined = build_joined_source(frames)
    joined["smoking_status"] = classify_smoking(joined)
    joined["diabetes_status"] = classify_diabetes(joined)
    attrition, core_mask, full_without_uacr = make_attrition(joined)

    cohort = add_features(joined.loc[core_mask].copy())
    cohort["fully_adjusted_without_uacr"] = full_without_uacr.loc[cohort.index]
    cohort["fully_adjusted_with_uacr"] = (
        cohort["fully_adjusted_without_uacr"] & cohort["URDACT"].notna()
    )

    output_columns = [
        "SEQN",
        "RIDAGEYR",
        "RIAGENDR",
        "sex_label",
        "RIDRETH3",
        "DMDEDUC2",
        "INDFMPIR",
        "SDMVSTRA",
        "SDMVPSU",
        "WTDR2D",
        "DR1TKCAL",
        "DR2TKCAL",
        "mean_energy",
        *PRAL_DAY1,
        *PRAL_DAY2,
        "pral_day1",
        "pral_day2",
        "pral_mean",
        "pral_difference",
        "pral_absolute_difference",
        "LBXSC3SI",
        "LBXSCR",
        "egfr",
        "BMXBMI",
        "BMXWT",
        "protein_g_per_kg",
        "URDACT",
        "smoking_status",
        "diabetes_status",
        "KIQ022",
        "KIQ025",
        "fully_adjusted_without_uacr",
        "fully_adjusted_with_uacr",
    ]
    cohort[output_columns].to_csv(PROCESSED / "analysis_cohort.csv", index=False)
    attrition.to_csv(DOCS / "cohort_attrition.csv", index=False)

    manifest_rows = []
    for path in sorted(RAW.glob("*.xpt")):
        manifest_rows.append(
            {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        )
    with (RAW / "sha256_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file_name", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    usda_manifest_rows = []
    for path in sorted(USDA_RAW.glob("*.xlsx")):
        usda_manifest_rows.append(
            {"file_name": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        )
    with (USDA_RAW / "sha256_manifest.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["file_name", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(usda_manifest_rows)

    summary = {
        "released_participants": int(len(joined)),
        "core_analytic_cohort": int(len(cohort)),
        "fully_adjusted_without_uacr": int(cohort["fully_adjusted_without_uacr"].sum()),
        "fully_adjusted_with_uacr": int(cohort["fully_adjusted_with_uacr"].sum()),
        "participant_id_unique": bool(cohort["SEQN"].is_unique),
        "missing_pral_day1": int(cohort["pral_day1"].isna().sum()),
        "missing_pral_day2": int(cohort["pral_day2"].isna().sum()),
        "missing_bicarbonate": int(cohort["LBXSC3SI"].isna().sum()),
        "missing_egfr": int(cohort["egfr"].isna().sum()),
    }
    (DOCS / "cohort_validation.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()