"""Complete Milestone 5: USDA/WWEIA food-category PRAL contributions."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.validate_food_sources import prepare_day
from src.reliability import weighted_mean, weighted_quantile


PROCESSED = ROOT / "data" / "processed"
USDA = ROOT / "data" / "raw" / "usda"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
DOCS = ROOT / "docs"
EXPECTED_COHORT_SIZE = 3848
EXPECTED_FOOD_RECORDS = 120956
NUTRIENTS = [
    "protein_g",
    "phosphorus_mg",
    "potassium_mg",
    "magnesium_mg",
    "calcium_mg",
]
COMPONENTS = ["positive_component", "negative_component", "net_pral"]
COLOR_POSITIVE = "#8F3A4C"
COLOR_NEGATIVE = "#2F7D66"
COLOR_VARIABILITY = "#2F6B8A"
COLOR_DARK = "#263238"
COLOR_GRID = "#D9E1E5"


def add_pral_components(frame: pd.DataFrame) -> pd.DataFrame:
    """Add the positive, negative, and net terms of the PRAL equation."""

    missing = set(NUTRIENTS).difference(frame.columns)
    if missing:
        raise ValueError(f"Food records are missing nutrient columns: {sorted(missing)}")
    result = frame.copy()
    # NHANES SAS special missing values in individual-food nutrient fields are
    # structural zeroes; participant-level nutrient sums verify this treatment.
    result[NUTRIENTS] = result[NUTRIENTS].fillna(0.0)
    values = result[NUTRIENTS].to_numpy(dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Food-level PRAL nutrients must be finite and nonnegative")
    result["positive_component"] = (
        0.49 * result["protein_g"] + 0.037 * result["phosphorus_mg"]
    )
    result["negative_component"] = -(
        0.021 * result["potassium_mg"]
        + 0.026 * result["magnesium_mg"]
        + 0.013 * result["calcium_mg"]
    )
    result["net_pral"] = (
        result["positive_component"] + result["negative_component"]
    )
    return result


def assign_weighted_quartiles(
    values: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Assign quartiles using weighted empirical inverse-CDF cut points."""

    cutpoints = weighted_quantile(values, weights, [0.25, 0.5, 0.75])
    if not np.all(np.diff(cutpoints) > 0):
        raise ValueError("Weighted PRAL quartile cut points must be increasing")
    quartiles = np.digitize(values, cutpoints, right=True) + 1
    return quartiles, cutpoints


def weighted_component_columns(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for component in COMPONENTS:
        result[f"weighted_{component}"] = result[component] * result["WTDR2D"]
    return result


def summarize_overall(
    category_day: pd.DataFrame, cohort: pd.DataFrame
) -> pd.DataFrame:
    denominator = 2 * cohort["WTDR2D"].sum()
    rows = weighted_component_columns(category_day)
    summary = rows.groupby(
        ["category_number", "category_description"], as_index=False
    ).agg(
        participant_day_records=("seqn", "size"),
        weighted_presence=("WTDR2D", "sum"),
        weighted_positive=("weighted_positive_component", "sum"),
        weighted_negative=("weighted_negative_component", "sum"),
        weighted_net=("weighted_net_pral", "sum"),
    )
    summary["unweighted_participant_day_prevalence_percent"] = (
        summary["participant_day_records"] / (2 * len(cohort)) * 100
    )
    summary["survey_weighted_participant_day_prevalence_percent"] = (
        summary["weighted_presence"] / denominator * 100
    )
    summary["survey_weighted_mean_positive_component_mEq_day"] = (
        summary["weighted_positive"] / denominator
    )
    summary["survey_weighted_mean_negative_component_mEq_day"] = (
        summary["weighted_negative"] / denominator
    )
    summary["survey_weighted_mean_net_pral_mEq_day"] = (
        summary["weighted_net"] / denominator
    )
    positive_total = summary[
        "survey_weighted_mean_positive_component_mEq_day"
    ].sum()
    negative_total = -summary[
        "survey_weighted_mean_negative_component_mEq_day"
    ].sum()
    summary["positive_component_share_percent"] = (
        summary["survey_weighted_mean_positive_component_mEq_day"]
        / positive_total
        * 100
    )
    summary["negative_component_share_percent"] = (
        -summary["survey_weighted_mean_negative_component_mEq_day"]
        / negative_total
        * 100
    )
    summary["positive_net_rank"] = pd.Series(pd.NA, index=summary.index, dtype="Int64")
    positive = summary["survey_weighted_mean_net_pral_mEq_day"].gt(0)
    summary.loc[positive, "positive_net_rank"] = (
        summary.loc[positive, "survey_weighted_mean_net_pral_mEq_day"]
        .rank(method="first", ascending=False)
        .astype("int64")
    )
    summary["negative_net_rank"] = pd.Series(pd.NA, index=summary.index, dtype="Int64")
    negative = summary["survey_weighted_mean_net_pral_mEq_day"].lt(0)
    summary.loc[negative, "negative_net_rank"] = (
        summary.loc[negative, "survey_weighted_mean_net_pral_mEq_day"]
        .rank(method="first", ascending=True)
        .astype("int64")
    )
    return summary.drop(
        columns=[
            "weighted_presence",
            "weighted_positive",
            "weighted_negative",
            "weighted_net",
        ]
    ).sort_values("survey_weighted_mean_net_pral_mEq_day", ascending=False)


def summarize_by_quartile(
    category_day: pd.DataFrame, cohort: pd.DataFrame
) -> pd.DataFrame:
    categories = category_day[
        ["category_number", "category_description"]
    ].drop_duplicates()
    outputs = []
    for quartile in range(1, 5):
        cohort_piece = cohort.loc[cohort["pral_mean_weighted_quartile"].eq(quartile)]
        denominator = 2 * cohort_piece["WTDR2D"].sum()
        piece = weighted_component_columns(
            category_day.loc[
                category_day["pral_mean_weighted_quartile"].eq(quartile)
            ]
        )
        grouped = piece.groupby(
            ["category_number", "category_description"], as_index=False
        ).agg(
            participant_day_records=("seqn", "size"),
            weighted_presence=("WTDR2D", "sum"),
            weighted_positive=("weighted_positive_component", "sum"),
            weighted_negative=("weighted_negative_component", "sum"),
            weighted_net=("weighted_net_pral", "sum"),
        )
        grouped = categories.merge(
            grouped,
            on=["category_number", "category_description"],
            how="left",
            validate="one_to_one",
        ).fillna(0)
        grouped["pral_mean_weighted_quartile"] = quartile
        grouped["participants_in_quartile"] = len(cohort_piece)
        grouped["survey_weighted_participant_day_prevalence_percent"] = (
            grouped["weighted_presence"] / denominator * 100
        )
        grouped["survey_weighted_mean_positive_component_mEq_day"] = (
            grouped["weighted_positive"] / denominator
        )
        grouped["survey_weighted_mean_negative_component_mEq_day"] = (
            grouped["weighted_negative"] / denominator
        )
        grouped["survey_weighted_mean_net_pral_mEq_day"] = (
            grouped["weighted_net"] / denominator
        )
        grouped["net_contribution_rank_within_quartile"] = (
            grouped["survey_weighted_mean_net_pral_mEq_day"]
            .rank(method="first", ascending=False)
            .astype("int64")
        )
        outputs.append(grouped)
    result = pd.concat(outputs, ignore_index=True)
    return result.drop(
        columns=[
            "weighted_presence",
            "weighted_positive",
            "weighted_negative",
            "weighted_net",
        ]
    ).sort_values(
        ["pral_mean_weighted_quartile", "net_contribution_rank_within_quartile"]
    )


def summarize_day_differences(
    category_day: pd.DataFrame, cohort: pd.DataFrame
) -> pd.DataFrame:
    pivot = category_day.pivot(
        index=["seqn", "category_number", "category_description"],
        columns="day",
        values="net_pral",
    ).fillna(0.0)
    pivot = pivot.rename(columns={1: "day1_net_pral", 2: "day2_net_pral"}).reset_index()
    for column in ("day1_net_pral", "day2_net_pral"):
        if column not in pivot:
            pivot[column] = 0.0
    pivot = pivot.merge(
        cohort[["seqn", "WTDR2D"]], on="seqn", how="left", validate="many_to_one"
    )
    pivot["signed_difference"] = pivot["day1_net_pral"] - pivot["day2_net_pral"]
    pivot["absolute_difference"] = pivot["signed_difference"].abs()
    for column in (
        "day1_net_pral",
        "day2_net_pral",
        "signed_difference",
        "absolute_difference",
    ):
        pivot[f"weighted_{column}"] = pivot[column] * pivot["WTDR2D"]
    denominator = cohort["WTDR2D"].sum()
    summary = pivot.groupby(
        ["category_number", "category_description"], as_index=False
    ).agg(
        participants_consuming_either_day=("seqn", "size"),
        weighted_presence=("WTDR2D", "sum"),
        weighted_day1=("weighted_day1_net_pral", "sum"),
        weighted_day2=("weighted_day2_net_pral", "sum"),
        weighted_signed_difference=("weighted_signed_difference", "sum"),
        weighted_absolute_difference=("weighted_absolute_difference", "sum"),
    )
    summary["survey_weighted_prevalence_any_day_percent"] = (
        summary["weighted_presence"] / denominator * 100
    )
    summary["survey_weighted_mean_day1_net_pral_mEq_day"] = (
        summary["weighted_day1"] / denominator
    )
    summary["survey_weighted_mean_day2_net_pral_mEq_day"] = (
        summary["weighted_day2"] / denominator
    )
    summary["survey_weighted_mean_signed_difference_mEq_day"] = (
        summary["weighted_signed_difference"] / denominator
    )
    summary["survey_weighted_mean_absolute_difference_mEq_day"] = (
        summary["weighted_absolute_difference"] / denominator
    )
    total_variability = summary[
        "survey_weighted_mean_absolute_difference_mEq_day"
    ].sum()
    summary["category_variability_share_percent"] = (
        summary["survey_weighted_mean_absolute_difference_mEq_day"]
        / total_variability
        * 100
    )
    summary["absolute_difference_rank"] = (
        summary["survey_weighted_mean_absolute_difference_mEq_day"]
        .rank(method="first", ascending=False)
        .astype("int64")
    )
    return summary.drop(
        columns=[
            "weighted_presence",
            "weighted_day1",
            "weighted_day2",
            "weighted_signed_difference",
            "weighted_absolute_difference",
        ]
    ).sort_values("absolute_difference_rank")


def write_figures(overall: pd.DataFrame, differences: pd.DataFrame) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/dietary-acid-load-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titleweight": "bold",
            "axes.edgecolor": COLOR_DARK,
            "axes.labelcolor": COLOR_DARK,
            "xtick.color": COLOR_DARK,
            "ytick.color": COLOR_DARK,
        }
    )

    value = "survey_weighted_mean_net_pral_mEq_day"
    positive = overall.loc[overall[value].gt(0)].nlargest(10, value)
    negative = overall.loc[overall[value].lt(0)].nsmallest(10, value)
    plot_data = pd.concat([negative, positive]).sort_values(value)
    colors = np.where(plot_data[value].gt(0), COLOR_POSITIVE, COLOR_NEGATIVE)
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    ax.barh(plot_data["category_description"], plot_data[value], color=colors)
    ax.axvline(0, color=COLOR_DARK, linewidth=1)
    ax.set_title("Leading Food-Category Contributions to Two-Day Mean PRAL")
    ax.set_xlabel("Survey-weighted mean net PRAL contribution (mEq/day)")
    ax.set_ylabel("")
    ax.grid(axis="x", color=COLOR_GRID, linewidth=0.7, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.35, right=0.97, top=0.92, bottom=0.1)
    fig.savefig(
        FIGURES / "food_category_net_pral_contributions.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    variability = "survey_weighted_mean_absolute_difference_mEq_day"
    plot_data = differences.nsmallest(15, "absolute_difference_rank").sort_values(
        variability
    )
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    ax.barh(
        plot_data["category_description"],
        plot_data[variability],
        color=COLOR_VARIABILITY,
    )
    ax.set_title("Food Categories Driving Day-to-Day PRAL Variation")
    ax.set_xlabel(
        "Survey-weighted mean absolute Day 1-Day 2 contribution difference (mEq/day)"
    )
    ax.set_ylabel("")
    ax.grid(axis="x", color=COLOR_GRID, linewidth=0.7, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.35, right=0.97, top=0.92, bottom=0.1)
    fig.savefig(
        FIGURES / "food_category_day_to_day_variation.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(PROCESSED / "analysis_cohort.csv").rename(
        columns={"SEQN": "seqn"}
    )
    required_cohort = {
        "seqn",
        "WTDR2D",
        "pral_day1",
        "pral_day2",
        "pral_mean",
        "pral_difference",
    }
    missing_cohort = required_cohort.difference(cohort.columns)
    if missing_cohort:
        raise ValueError(f"Cohort is missing columns: {sorted(missing_cohort)}")
    if len(cohort) != EXPECTED_COHORT_SIZE or not cohort["seqn"].is_unique:
        raise ValueError("The core cohort must contain 3,848 unique participants")
    if not cohort["WTDR2D"].gt(0).all():
        raise ValueError("All two-day dietary weights must be positive")

    quartiles, cutpoints = assign_weighted_quartiles(
        cohort["pral_mean"].to_numpy(dtype=float),
        cohort["WTDR2D"].to_numpy(dtype=float),
    )
    cohort["pral_mean_weighted_quartile"] = quartiles
    core_ids = set(cohort["seqn"].astype("int64"))
    foods = pd.concat(
        [prepare_day(1, core_ids), prepare_day(2, core_ids)], ignore_index=True
    )
    if len(foods) != EXPECTED_FOOD_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_FOOD_RECORDS:,} food records, found {len(foods):,}"
        )
    if foods.duplicated(["seqn", "day", "food_line"]).any():
        raise ValueError("Food records contain duplicate participant-day-line keys")

    categories = pd.read_excel(
        USDA / "WWEIA_August2021_August2023_foodcat_FNDDS.xlsx",
        sheet_name="Aug2021-Aug2023_FNDDS_foodcat",
    )
    categories["food_code"] = categories["food_code"].astype("int64")
    if categories["food_code"].duplicated().any():
        raise ValueError("USDA category crosswalk contains duplicate food codes")
    joined = foods.merge(
        categories[
            [
                "food_code",
                "food_code_description",
                "category_number",
                "category_description",
            ]
        ],
        on="food_code",
        how="left",
        validate="many_to_one",
    )
    unmatched = joined["category_number"].isna()
    if unmatched.any():
        raise ValueError(f"{int(unmatched.sum())} food records lack a USDA category")
    joined["category_number"] = joined["category_number"].astype("int64")
    joined = add_pral_components(joined)

    category_day = joined.groupby(
        ["seqn", "day", "category_number", "category_description"],
        as_index=False,
    )[COMPONENTS].sum()
    if category_day.duplicated(["seqn", "day", "category_number"]).any():
        raise ValueError("Participant-day-category output contains duplicate keys")
    category_day = category_day.merge(
        cohort[["seqn", "WTDR2D", "pral_mean_weighted_quartile"]],
        on="seqn",
        how="left",
        validate="many_to_one",
    )
    category_day.to_csv(
        PROCESSED / "food_category_pral_contributions.csv", index=False
    )

    category_totals = category_day.groupby(["seqn", "day"], as_index=False)[
        COMPONENTS
    ].sum()
    expected = pd.concat(
        [
            cohort[["seqn", "pral_day1"]]
            .rename(columns={"pral_day1": "expected_pral"})
            .assign(day=1),
            cohort[["seqn", "pral_day2"]]
            .rename(columns={"pral_day2": "expected_pral"})
            .assign(day=2),
        ],
        ignore_index=True,
    )
    total_check = expected.merge(
        category_totals, on=["seqn", "day"], how="left", validate="one_to_one"
    )
    participant_total_max_error = float(
        (total_check["expected_pral"] - total_check["net_pral"]).abs().max()
    )
    component_identity_max_error = float(
        (
            joined["net_pral"]
            - joined["positive_component"]
            - joined["negative_component"]
        )
        .abs()
        .max()
    )

    overall = summarize_overall(category_day, cohort)
    by_quartile = summarize_by_quartile(category_day, cohort)
    day_differences = summarize_day_differences(category_day, cohort)
    overall.to_csv(TABLES / "food_category_pral_overall.csv", index=False)
    by_quartile.to_csv(TABLES / "food_category_pral_by_quartile.csv", index=False)
    day_differences.to_csv(
        TABLES / "food_category_pral_day_differences.csv", index=False
    )
    pd.DataFrame(
        {
            "percentile": [25, 50, 75],
            "cutpoint_mEq_per_day": cutpoints,
            "method": "WTDR2D-weighted empirical inverse CDF of two-day mean PRAL",
        }
    ).to_csv(TABLES / "pral_mean_weighted_quartile_cutpoints.csv", index=False)
    write_figures(overall, day_differences)

    overall_net_sum = float(
        overall["survey_weighted_mean_net_pral_mEq_day"].sum()
    )
    expected_weighted_mean = weighted_mean(
        cohort["pral_mean"].to_numpy(dtype=float),
        cohort["WTDR2D"].to_numpy(dtype=float),
    )
    category_signed_difference_sum = float(
        day_differences[
            "survey_weighted_mean_signed_difference_mEq_day"
        ].sum()
    )
    expected_weighted_difference = weighted_mean(
        cohort["pral_difference"].to_numpy(dtype=float),
        cohort["WTDR2D"].to_numpy(dtype=float),
    )
    quartile_category_means = by_quartile.groupby(
        "pral_mean_weighted_quartile"
    )["survey_weighted_mean_net_pral_mEq_day"].sum()
    quartile_expected_means = cohort.groupby("pral_mean_weighted_quartile").apply(
        lambda group: weighted_mean(
            group["pral_mean"].to_numpy(dtype=float),
            group["WTDR2D"].to_numpy(dtype=float),
        ),
        include_groups=False,
    )
    quartile_mean_max_error = float(
        (quartile_category_means - quartile_expected_means).abs().max()
    )
    positive_share_total = float(overall["positive_component_share_percent"].sum())
    negative_share_total = float(overall["negative_component_share_percent"].sum())
    variability_share_total = float(
        day_differences["category_variability_share_percent"].sum()
    )
    identity_tolerance = 1e-8
    identity_errors = {
        "participant PRAL totals": participant_total_max_error,
        "food component formula": component_identity_max_error,
        "overall weighted mean": abs(overall_net_sum - expected_weighted_mean),
        "weighted signed day difference": abs(
            category_signed_difference_sum - expected_weighted_difference
        ),
        "weighted quartile means": quartile_mean_max_error,
        "positive contribution shares": abs(positive_share_total - 100),
        "negative contribution shares": abs(negative_share_total - 100),
        "variability shares": abs(variability_share_total - 100),
    }
    failed_identities = {
        name: error
        for name, error in identity_errors.items()
        if error > identity_tolerance
    }
    if failed_identities:
        raise ValueError(f"Food-category validation failed: {failed_identities}")
    validation = {
        "source_food_records": int(len(joined)),
        "participant_day_category_rows": int(len(category_day)),
        "participants": int(len(cohort)),
        "categories_observed": int(overall["category_number"].nunique()),
        "unmatched_food_records": int(unmatched.sum()),
        "duplicate_food_keys": int(
            foods.duplicated(["seqn", "day", "food_line"]).sum()
        ),
        "duplicate_participant_day_category_keys": int(
            category_day.duplicated(["seqn", "day", "category_number"]).sum()
        ),
        "participant_pral_total_max_absolute_error": participant_total_max_error,
        "food_component_identity_max_absolute_error": component_identity_max_error,
        "overall_category_net_sum": overall_net_sum,
        "expected_survey_weighted_mean_pral": expected_weighted_mean,
        "overall_weighted_mean_identity_absolute_error": abs(
            overall_net_sum - expected_weighted_mean
        ),
        "category_signed_difference_sum": category_signed_difference_sum,
        "expected_survey_weighted_day1_minus_day2_difference": expected_weighted_difference,
        "signed_difference_identity_absolute_error": abs(
            category_signed_difference_sum - expected_weighted_difference
        ),
        "quartile_weighted_mean_identity_max_absolute_error": quartile_mean_max_error,
        "positive_component_share_total_percent": positive_share_total,
        "negative_component_share_total_percent": negative_share_total,
        "variability_share_total_percent": variability_share_total,
        "identity_tolerance": identity_tolerance,
        "pral_mean_quartile_cutpoints": cutpoints.tolist(),
        "quartile_values_within_1_to_4": bool(
            np.isin(quartiles, [1, 2, 3, 4]).all()
        ),
        "weighting_note": (
            "WTDR2D-weighted descriptive point estimates; zero contribution is "
            "included when a category was not consumed; no design-based CIs"
        ),
    }
    (DOCS / "food_category_analysis_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )

    net_column = "survey_weighted_mean_net_pral_mEq_day"
    positive_top = overall.loc[overall[net_column].gt(0)].nlargest(5, net_column)
    negative_top = overall.loc[overall[net_column].lt(0)].nsmallest(5, net_column)
    variability_column = "survey_weighted_mean_absolute_difference_mEq_day"
    variability_top = day_differences.nsmallest(5, "absolute_difference_rank")
    console_summary = {
        "participants": int(len(cohort)),
        "food_records": int(len(joined)),
        "categories_observed": int(overall["category_number"].nunique()),
        "top_positive_net_categories": [
            {
                "category": row.category_description,
                "mean_net_mEq_day": float(getattr(row, net_column)),
            }
            for row in positive_top.itertuples()
        ],
        "top_negative_net_categories": [
            {
                "category": row.category_description,
                "mean_net_mEq_day": float(getattr(row, net_column)),
            }
            for row in negative_top.itertuples()
        ],
        "top_day_to_day_variability_categories": [
            {
                "category": row.category_description,
                "mean_absolute_difference_mEq_day": float(
                    getattr(row, variability_column)
                ),
            }
            for row in variability_top.itertuples()
        ],
    }
    print(json.dumps(console_summary, indent=2))


if __name__ == "__main__":
    main()