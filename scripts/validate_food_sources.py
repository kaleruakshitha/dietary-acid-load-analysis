"""Validate NHANES individual-food files against USDA food categories."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NHANES = ROOT / "data" / "raw" / "nhanes"
USDA = ROOT / "data" / "raw" / "usda"
PROCESSED = ROOT / "data" / "processed"
DOCS = ROOT / "docs"


def read_xpt(name: str) -> pd.DataFrame:
    frame = pd.read_sas(NHANES / f"{name}.xpt", format="xport")
    numeric_columns = frame.select_dtypes(include="number").columns
    for column in numeric_columns:
        values = frame[column]
        sas_missing = values.ne(0) & values.abs().lt(1e-50)
        if sas_missing.any():
            frame.loc[sas_missing, column] = np.nan
    if "SEQN" in frame:
        frame["SEQN"] = frame["SEQN"].astype("int64")
    return frame


def prepare_day(day: int, core_ids: set[int]) -> pd.DataFrame:
    prefix = f"DR{day}"
    frame = read_xpt(f"DR{day}IFF_L")
    columns = [
        "SEQN",
        f"{prefix}ILINE",
        f"{prefix}IFDCD",
        f"{prefix}IPROT",
        f"{prefix}IPHOS",
        f"{prefix}IPOTA",
        f"{prefix}IMAGN",
        f"{prefix}ICALC",
    ]
    frame = frame.loc[frame["SEQN"].isin(core_ids), columns].copy()
    frame.columns = [
        "seqn",
        "food_line",
        "food_code",
        "protein_g",
        "phosphorus_mg",
        "potassium_mg",
        "magnesium_mg",
        "calcium_mg",
    ]
    frame["food_code"] = frame["food_code"].astype("int64")
    frame["day"] = day
    if frame.duplicated(["seqn", "food_line"]).any():
        raise ValueError(f"Day {day} has duplicate participant-food line keys")
    return frame


def main() -> None:
    cohort = pd.read_csv(PROCESSED / "analysis_cohort.csv")
    core_ids = set(cohort["SEQN"].astype("int64"))
    day1 = prepare_day(1, core_ids)
    day2 = prepare_day(2, core_ids)
    foods = pd.concat([day1, day2], ignore_index=True)

    categories = pd.read_excel(
        USDA / "WWEIA_August2021_August2023_foodcat_FNDDS.xlsx",
        sheet_name="Aug2021-Aug2023_FNDDS_foodcat",
    )
    categories["food_code"] = categories["food_code"].astype("int64")
    if categories["food_code"].duplicated().any():
        raise ValueError("USDA category crosswalk has duplicate food codes")

    joined = foods.merge(
        categories[["food_code", "category_number", "category_description"]],
        on="food_code",
        how="left",
        validate="many_to_one",
    )

    nutrient_columns = [
        "protein_g",
        "phosphorus_mg",
        "potassium_mg",
        "magnesium_mg",
        "calcium_mg",
    ]
    totals = joined.groupby(["seqn", "day"], as_index=False)[nutrient_columns].sum()

    comparison = cohort[[
        "SEQN",
        "DR1TPROT",
        "DR1TPHOS",
        "DR1TPOTA",
        "DR1TMAGN",
        "DR1TCALC",
        "DR2TPROT",
        "DR2TPHOS",
        "DR2TPOTA",
        "DR2TMAGN",
        "DR2TCALC",
    ]].copy()
    long_totals = []
    for day in (1, 2):
        piece = comparison[[
            "SEQN",
            f"DR{day}TPROT",
            f"DR{day}TPHOS",
            f"DR{day}TPOTA",
            f"DR{day}TMAGN",
            f"DR{day}TCALC",
        ]].copy()
        piece.columns = ["seqn", *nutrient_columns]
        piece["day"] = day
        long_totals.append(piece)
    expected = pd.concat(long_totals, ignore_index=True)
    check = expected.merge(
        totals,
        on=["seqn", "day"],
        how="left",
        validate="one_to_one",
        suffixes=("_total", "_food_sum"),
    )
    max_differences = {
        nutrient: float(
            (check[f"{nutrient}_total"] - check[f"{nutrient}_food_sum"])
            .abs()
            .max()
        )
        for nutrient in nutrient_columns
    }

    result = {
        "core_participants": len(core_ids),
        "participants_with_day1_food_records": int(day1["seqn"].nunique()),
        "participants_with_day2_food_records": int(day2["seqn"].nunique()),
        "day1_food_records": int(len(day1)),
        "day2_food_records": int(len(day2)),
        "usda_food_codes": int(categories["food_code"].nunique()),
        "unmatched_food_records": int(joined["category_number"].isna().sum()),
        "unmatched_unique_food_codes": int(
            joined.loc[joined["category_number"].isna(), "food_code"].nunique()
        ),
        "maximum_absolute_nutrient_sum_differences": max_differences,
    }
    (DOCS / "food_source_validation.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()