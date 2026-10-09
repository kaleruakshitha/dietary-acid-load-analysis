"""Complete Milestone 3: Day 1 versus Day 2 PRAL reliability and agreement."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/dietary-acid-load-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.reliability import (
    cohen_kappa_from_matrix,
    icc_two_way_random_absolute,
    weighted_mean,
    weighted_pearson,
    weighted_quantile,
    weighted_standard_deviation,
)


ROOT = Path(__file__).resolve().parents[1]
COHORT_PATH = ROOT / "data" / "processed" / "analysis_cohort.csv"
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
DOCS = ROOT / "docs"
EXPECTED_COHORT_SIZE = 3848
BOOTSTRAP_REPLICATES = 2000
RANDOM_SEED = 20260819
DAY1 = "pral_day1"
DAY2 = "pral_day2"
WEIGHT = "WTDR2D"
COLOR_DAY1 = "#2F6B8A"
COLOR_DAY2 = "#D9822B"
COLOR_ACCENT = "#8F3A4C"
COLOR_DARK = "#263238"
COLOR_GRID = "#D9E1E5"


def validate_cohort(frame: pd.DataFrame) -> None:
    required = {
        "SEQN",
        DAY1,
        DAY2,
        "pral_mean",
        "pral_difference",
        "pral_absolute_difference",
        WEIGHT,
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Analysis cohort is missing columns: {sorted(missing)}")
    if len(frame) != EXPECTED_COHORT_SIZE:
        raise ValueError(
            f"Expected {EXPECTED_COHORT_SIZE} core participants, found {len(frame)}"
        )
    if not frame["SEQN"].is_unique:
        raise ValueError("Analysis cohort contains duplicate participant identifiers")
    if frame[list(required)].isna().any().any():
        raise ValueError("Required reliability variables contain missing values")
    if not np.isfinite(frame[[DAY1, DAY2, WEIGHT]].to_numpy(dtype=float)).all():
        raise ValueError("Required reliability variables contain nonfinite values")
    if not frame[WEIGHT].gt(0).all():
        raise ValueError("All two-day dietary weights must be positive")


def descriptive_row(
    name: str, values: np.ndarray, weights: np.ndarray | None
) -> dict[str, float | int | str]:
    if weights is None:
        quantiles = np.quantile(values, [0.25, 0.5, 0.75])
        mean = float(np.mean(values))
        standard_deviation = float(np.std(values, ddof=1))
        weighting = "unweighted"
    else:
        quantiles = weighted_quantile(values, weights, [0.25, 0.5, 0.75])
        mean = weighted_mean(values, weights)
        standard_deviation = weighted_standard_deviation(values, weights)
        weighting = "survey_weighted"
    return {
        "variable": name,
        "weighting": weighting,
        "n": int(values.size),
        "mean": mean,
        "standard_deviation": standard_deviation,
        "q1": float(quantiles[0]),
        "median": float(quantiles[1]),
        "q3": float(quantiles[2]),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
    }


def fisher_confidence_interval(correlation: float, n: int) -> tuple[float, float]:
    z = np.arctanh(correlation)
    margin = stats.norm.ppf(0.975) / np.sqrt(n - 3)
    return float(np.tanh(z - margin)), float(np.tanh(z + margin))


def bootstrap_icc(measurements: np.ndarray) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(RANDOM_SEED)
    n = measurements.shape[0]
    absolute = np.empty(BOOTSTRAP_REPLICATES)
    consistency = np.empty(BOOTSTRAP_REPLICATES)
    for replicate in range(BOOTSTRAP_REPLICATES):
        sample = measurements[rng.integers(0, n, size=n)]
        absolute[replicate], consistency[replicate] = (
            icc_two_way_random_absolute(sample)
        )
    return {
        "absolute": tuple(np.percentile(absolute, [2.5, 97.5])),
        "consistency": tuple(np.percentile(consistency, [2.5, 97.5])),
    }


def assign_quartiles(values: np.ndarray, cutpoints: np.ndarray) -> np.ndarray:
    if cutpoints.shape != (3,) or not np.all(np.diff(cutpoints) > 0):
        raise ValueError("Quartile cut points must contain three increasing values")
    return np.digitize(values, cutpoints, right=True) + 1


def transition_matrices(
    day1_quartile: np.ndarray, day2_quartile: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    counts = np.zeros((4, 4), dtype=int)
    weighted = np.zeros((4, 4), dtype=float)
    np.add.at(counts, (day1_quartile - 1, day2_quartile - 1), 1)
    np.add.at(weighted, (day1_quartile - 1, day2_quartile - 1), weights)
    return counts, weighted


def reliability_row(
    statistic: str,
    weighting: str,
    estimate: float,
    ci_lower: float | None = None,
    ci_upper: float | None = None,
    unit: str = "coefficient",
    notes: str = "",
) -> dict[str, float | str | None]:
    return {
        "statistic": statistic,
        "weighting": weighting,
        "estimate": float(estimate),
        "ci_lower": None if ci_lower is None else float(ci_lower),
        "ci_upper": None if ci_upper is None else float(ci_upper),
        "unit": unit,
        "notes": notes,
    }


def write_figures(
    day1: np.ndarray,
    day2: np.ndarray,
    differences: np.ndarray,
    quartile_summary: pd.DataFrame,
    bias: float,
    lower_limit: float,
    upper_limit: float,
    pearson: float,
    icc_absolute: float,
) -> None:
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

    participant_mean = (day1 + day2) / 2
    fig, ax = plt.subplots(figsize=(9, 6.2))
    ax.scatter(
        participant_mean,
        differences,
        s=13,
        alpha=0.28,
        color=COLOR_DAY1,
        edgecolors="none",
        rasterized=True,
    )
    ax.axhline(bias, color=COLOR_DARK, linewidth=2, label=f"Mean difference: {bias:.2f}")
    ax.axhline(
        upper_limit,
        color=COLOR_ACCENT,
        linestyle="--",
        linewidth=1.8,
        label=f"Upper limit: {upper_limit:.2f}",
    )
    ax.axhline(
        lower_limit,
        color=COLOR_ACCENT,
        linestyle="--",
        linewidth=1.8,
        label=f"Lower limit: {lower_limit:.2f}",
    )
    ax.set_title("Agreement of Day 1 and Day 2 PRAL")
    ax.set_xlabel("Participant mean PRAL across two recalls (mEq/day)")
    ax.set_ylabel("Day 1 − Day 2 PRAL (mEq/day)")
    ax.grid(color=COLOR_GRID, linewidth=0.7, alpha=0.75)
    ax.legend(frameon=False, loc="upper right")
    ax.text(
        0.01,
        0.01,
        f"n = {day1.size:,}; all prespecified eligible observations retained",
        transform=ax.transAxes,
        fontsize=9,
        color=COLOR_DARK,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "pral_bland_altman.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "pral_bland_altman.pdf", bbox_inches="tight")
    plt.close(fig)

    combined_min = float(min(day1.min(), day2.min()))
    combined_max = float(max(day1.max(), day2.max()))
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    density = ax.hexbin(
        day1,
        day2,
        gridsize=48,
        mincnt=1,
        cmap="Blues",
        linewidths=0,
    )
    ax.plot(
        [combined_min, combined_max],
        [combined_min, combined_max],
        color=COLOR_ACCENT,
        linestyle="--",
        linewidth=1.6,
        label="Perfect day-to-day agreement",
    )
    ax.set_xlim(combined_min, combined_max)
    ax.set_ylim(combined_min, combined_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Day 1 versus Day 2 PRAL")
    ax.set_xlabel("Day 1 PRAL (mEq/day)")
    ax.set_ylabel("Day 2 PRAL (mEq/day)")
    ax.grid(color=COLOR_GRID, linewidth=0.6, alpha=0.55)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.98,
        0.03,
        f"Pearson r = {pearson:.3f}\nICC(A,1) = {icc_absolute:.3f}",
        ha="right",
        va="bottom",
        transform=ax.transAxes,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    colorbar = fig.colorbar(density, ax=ax, pad=0.02)
    colorbar.set_label("Participant count per hexagon")
    fig.tight_layout()
    fig.savefig(FIGURES / "pral_day1_day2_hexbin.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    plot_data = quartile_summary.copy()
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    bars = ax.bar(
        plot_data["category"],
        plot_data["survey_weighted_percent"],
        color=[COLOR_DAY1, "#6E9F76", COLOR_DAY2, COLOR_ACCENT],
        width=0.68,
    )
    ax.bar_label(bars, fmt="%.1f%%", padding=3, color=COLOR_DARK, fontweight="bold")
    ax.set_title("Movement Between PRAL Quartiles Across Recall Days")
    ax.set_ylabel("Survey-weighted percentage")
    ax.set_xlabel("")
    ax.set_ylim(0, max(plot_data["survey_weighted_percent"]) * 1.18)
    ax.grid(axis="y", color=COLOR_GRID, linewidth=0.7, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", rotation=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "pral_quartile_reclassification.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(COHORT_PATH)
    validate_cohort(frame)

    day1 = frame[DAY1].to_numpy(dtype=float)
    day2 = frame[DAY2].to_numpy(dtype=float)
    weights = frame[WEIGHT].to_numpy(dtype=float)
    measurements = np.column_stack([day1, day2])
    differences = day1 - day2
    absolute_differences = np.abs(differences)
    participant_mean = (day1 + day2) / 2
    n = day1.size

    descriptives = []
    for name, values in (
        ("PRAL Day 1", day1),
        ("PRAL Day 2", day2),
        ("PRAL two-day mean", participant_mean),
        ("Day 1 minus Day 2", differences),
        ("Absolute day difference", absolute_differences),
    ):
        descriptives.append(descriptive_row(name, values, None))
        descriptives.append(descriptive_row(name, values, weights))
    pd.DataFrame(descriptives).to_csv(
        TABLES / "pral_descriptive_statistics.csv", index=False
    )

    pearson_result = stats.pearsonr(day1, day2)
    pearson = float(pearson_result.statistic)
    pearson_ci = fisher_confidence_interval(pearson, n)
    spearman_result = stats.spearmanr(day1, day2)
    survey_weighted_pearson = weighted_pearson(day1, day2, weights)

    icc_absolute, icc_consistency = icc_two_way_random_absolute(measurements)
    icc_intervals = bootstrap_icc(measurements)

    bias = float(np.mean(differences))
    difference_sd = float(np.std(differences, ddof=1))
    lower_limit = bias - 1.96 * difference_sd
    upper_limit = bias + 1.96 * difference_sd
    t_critical = float(stats.t.ppf(0.975, n - 1))
    bias_margin = t_critical * difference_sd / np.sqrt(n)
    limit_standard_error = difference_sd * np.sqrt(1 / n + 1.96**2 / (2 * (n - 1)))
    limit_margin = t_critical * limit_standard_error

    weighted_bias = weighted_mean(differences, weights)
    weighted_difference_sd = weighted_standard_deviation(differences, weights)
    weighted_lower_limit = weighted_bias - 1.96 * weighted_difference_sd
    weighted_upper_limit = weighted_bias + 1.96 * weighted_difference_sd

    proportional_bias = stats.linregress(participant_mean, differences)

    day1_cutpoints = weighted_quantile(day1, weights, [0.25, 0.5, 0.75])
    day2_cutpoints = weighted_quantile(day2, weights, [0.25, 0.5, 0.75])
    day1_quartile = assign_quartiles(day1, day1_cutpoints)
    day2_quartile = assign_quartiles(day2, day2_cutpoints)
    quartile_distance = np.abs(day1_quartile - day2_quartile)
    count_matrix, weighted_matrix = transition_matrices(
        day1_quartile, day2_quartile, weights
    )

    reliability = [
        reliability_row(
            "Pearson correlation",
            "unweighted",
            pearson,
            pearson_ci[0],
            pearson_ci[1],
            notes="Ranking statistic; 95% CI uses Fisher transformation",
        ),
        reliability_row(
            "Spearman rank correlation",
            "unweighted",
            float(spearman_result.statistic),
            notes=f"Ranking statistic; two-sided p={spearman_result.pvalue:.6g}",
        ),
        reliability_row(
            "Pearson correlation",
            "survey_weighted_point_estimate",
            survey_weighted_pearson,
            notes="WTDR2D-weighted descriptive point estimate; no design-based CI",
        ),
        reliability_row(
            "ICC(A,1): two-way random absolute agreement, single measure",
            "unweighted",
            icc_absolute,
            icc_intervals["absolute"][0],
            icc_intervals["absolute"][1],
            notes=f"Percentile bootstrap CI; {BOOTSTRAP_REPLICATES} participant resamples",
        ),
        reliability_row(
            "ICC(C,1): two-way consistency, single measure",
            "unweighted_sensitivity",
            icc_consistency,
            icc_intervals["consistency"][0],
            icc_intervals["consistency"][1],
            notes=f"Percentile bootstrap CI; {BOOTSTRAP_REPLICATES} participant resamples",
        ),
        reliability_row(
            "Bland-Altman mean difference (Day 1 - Day 2)",
            "unweighted",
            bias,
            bias - bias_margin,
            bias + bias_margin,
            "mEq/day",
        ),
        reliability_row(
            "Bland-Altman lower 95% limit of agreement",
            "unweighted",
            lower_limit,
            lower_limit - limit_margin,
            lower_limit + limit_margin,
            "mEq/day",
        ),
        reliability_row(
            "Bland-Altman upper 95% limit of agreement",
            "unweighted",
            upper_limit,
            upper_limit - limit_margin,
            upper_limit + limit_margin,
            "mEq/day",
        ),
        reliability_row(
            "Bland-Altman mean difference (Day 1 - Day 2)",
            "survey_weighted_point_estimate",
            weighted_bias,
            unit="mEq/day",
            notes="WTDR2D-weighted descriptive estimate; no design-based CI",
        ),
        reliability_row(
            "Bland-Altman lower 95% limit of agreement",
            "survey_weighted_point_estimate",
            weighted_lower_limit,
            unit="mEq/day",
            notes="WTDR2D-weighted descriptive estimate; no design-based CI",
        ),
        reliability_row(
            "Bland-Altman upper 95% limit of agreement",
            "survey_weighted_point_estimate",
            weighted_upper_limit,
            unit="mEq/day",
            notes="WTDR2D-weighted descriptive estimate; no design-based CI",
        ),
        reliability_row(
            "Linearly weighted Cohen kappa",
            "unweighted",
            cohen_kappa_from_matrix(count_matrix, "linear"),
            notes="Day-specific survey-weighted quartile cut points",
        ),
        reliability_row(
            "Linearly weighted Cohen kappa",
            "survey_weighted_point_estimate",
            cohen_kappa_from_matrix(weighted_matrix, "linear"),
            notes="WTDR2D-weighted table and day-specific weighted quartile cut points",
        ),
        reliability_row(
            "Quadratically weighted Cohen kappa",
            "survey_weighted_sensitivity",
            cohen_kappa_from_matrix(weighted_matrix, "quadratic"),
            notes="WTDR2D-weighted table and day-specific weighted quartile cut points",
        ),
        reliability_row(
            "Bland-Altman proportional-bias slope",
            "unweighted_diagnostic",
            float(proportional_bias.slope),
            float(proportional_bias.slope - t_critical * proportional_bias.stderr),
            float(proportional_bias.slope + t_critical * proportional_bias.stderr),
            "mEq difference per mEq mean",
            notes=f"Exploratory diagnostic; two-sided p={proportional_bias.pvalue:.6g}",
        ),
    ]
    pd.DataFrame(reliability).to_csv(
        TABLES / "reliability_statistics.csv", index=False
    )

    cutpoint_rows = []
    for day, cutpoints in (("Day 1", day1_cutpoints), ("Day 2", day2_cutpoints)):
        for percentile, value in zip((25, 50, 75), cutpoints):
            cutpoint_rows.append(
                {
                    "recall_day": day,
                    "percentile": percentile,
                    "cutpoint_mEq_per_day": float(value),
                    "method": "WTDR2D-weighted empirical inverse CDF",
                }
            )
    pd.DataFrame(cutpoint_rows).to_csv(
        TABLES / "pral_weighted_quartile_cutpoints.csv", index=False
    )

    category_definitions = [
        ("Same quartile", 0),
        ("Adjacent quartile", 1),
        ("Two-quartile change", 2),
        ("Extreme Q1-Q4 change", 3),
    ]
    quartile_rows = []
    for category, distance in category_definitions:
        mask = quartile_distance == distance
        quartile_rows.append(
            {
                "category": category,
                "quartile_distance": distance,
                "participants": int(mask.sum()),
                "unweighted_percent": float(mask.mean() * 100),
                "survey_weighted_percent": float(weights[mask].sum() / weights.sum() * 100),
            }
        )
    quartile_summary = pd.DataFrame(quartile_rows)
    quartile_summary.to_csv(
        TABLES / "quartile_reclassification_summary.csv", index=False
    )

    transition_rows = []
    weighted_row_totals = weighted_matrix.sum(axis=1)
    for row in range(4):
        for column in range(4):
            transition_rows.append(
                {
                    "day1_quartile": row + 1,
                    "day2_quartile": column + 1,
                    "unweighted_count": int(count_matrix[row, column]),
                    "unweighted_overall_percent": float(
                        count_matrix[row, column] / n * 100
                    ),
                    "survey_weighted_overall_percent": float(
                        weighted_matrix[row, column] / weighted_matrix.sum() * 100
                    ),
                    "survey_weighted_row_percent": float(
                        weighted_matrix[row, column] / weighted_row_totals[row] * 100
                    ),
                }
            )
    pd.DataFrame(transition_rows).to_csv(
        TABLES / "quartile_transition_matrix.csv", index=False
    )

    threshold_rows = []
    for threshold in (5, 10, 20):
        mask = absolute_differences > threshold
        threshold_rows.append(
            {
                "absolute_difference_threshold_mEq_per_day": threshold,
                "comparison": "greater_than",
                "participants": int(mask.sum()),
                "unweighted_percent": float(mask.mean() * 100),
                "survey_weighted_percent": float(weights[mask].sum() / weights.sum() * 100),
            }
        )
    pd.DataFrame(threshold_rows).to_csv(
        TABLES / "absolute_difference_thresholds.csv", index=False
    )

    participant_output = pd.DataFrame(
        {
            "SEQN": frame["SEQN"].astype("int64"),
            "WTDR2D": weights,
            "pral_day1": day1,
            "pral_day2": day2,
            "pral_mean": participant_mean,
            "pral_difference": differences,
            "pral_absolute_difference": absolute_differences,
            "pral_day1_weighted_quartile": day1_quartile,
            "pral_day2_weighted_quartile": day2_quartile,
            "quartile_distance": quartile_distance,
        }
    )
    participant_output.to_csv(
        PROCESSED / "pral_reliability_participants.csv", index=False
    )

    write_figures(
        day1,
        day2,
        differences,
        quartile_summary,
        bias,
        lower_limit,
        upper_limit,
        pearson,
        icc_absolute,
    )

    validation = {
        "source_file": str(COHORT_PATH.relative_to(ROOT)),
        "participants": int(n),
        "participant_id_unique": bool(frame["SEQN"].is_unique),
        "required_values_complete": bool(
            frame[[DAY1, DAY2, WEIGHT]].notna().all().all()
        ),
        "all_weights_positive": bool(frame[WEIGHT].gt(0).all()),
        "no_outliers_removed": True,
        "day1_quartile_cutpoints_increasing": bool(np.all(np.diff(day1_cutpoints) > 0)),
        "day2_quartile_cutpoints_increasing": bool(np.all(np.diff(day2_cutpoints) > 0)),
        "quartile_values_within_1_to_4": bool(
            np.isin(day1_quartile, [1, 2, 3, 4]).all()
            and np.isin(day2_quartile, [1, 2, 3, 4]).all()
        ),
        "transition_count_total": int(count_matrix.sum()),
        "weighted_transition_percent_total": float(
            weighted_matrix.sum() / weights.sum() * 100
        ),
        "reclassification_percent_total": float(
            quartile_summary["survey_weighted_percent"].sum()
        ),
        "difference_identity_max_absolute_error": float(
            np.max(np.abs(differences - frame["pral_difference"].to_numpy(dtype=float)))
        ),
        "mean_identity_max_absolute_error": float(
            np.max(np.abs(participant_mean - frame["pral_mean"].to_numpy(dtype=float)))
        ),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "random_seed": RANDOM_SEED,
        "primary_icc_specification": "ICC(A,1), two-way random-effects, absolute agreement, single measure",
        "bland_altman_difference_order": "Day 1 minus Day 2",
        "quartile_method": "Day-specific WTDR2D-weighted empirical inverse CDF cut points",
    }
    (DOCS / "reliability_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )

    console_summary = {
        "participants": int(n),
        "pearson_r": pearson,
        "spearman_rho": float(spearman_result.statistic),
        "survey_weighted_pearson_r": survey_weighted_pearson,
        "icc_absolute": icc_absolute,
        "icc_absolute_95_ci": list(icc_intervals["absolute"]),
        "bland_altman_bias_day1_minus_day2": bias,
        "bland_altman_limits": [lower_limit, upper_limit],
        "survey_weighted_same_quartile_percent": float(
            quartile_summary.loc[
                quartile_summary["category"].eq("Same quartile"),
                "survey_weighted_percent",
            ].iloc[0]
        ),
        "survey_weighted_extreme_reclassification_percent": float(
            quartile_summary.loc[
                quartile_summary["category"].eq("Extreme Q1-Q4 change"),
                "survey_weighted_percent",
            ].iloc[0]
        ),
        "survey_weighted_linear_kappa": cohen_kappa_from_matrix(
            weighted_matrix, "linear"
        ),
    }
    print(json.dumps(console_summary, indent=2))


if __name__ == "__main__":
    main()