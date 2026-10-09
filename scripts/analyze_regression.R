#!/usr/bin/env Rscript

# Milestone 4: survey-weighted PRAL and serum bicarbonate regression.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1) {
  stop("Run this script with Rscript from the repository checkout.")
}

script_path <- normalizePath(sub("^--file=", "", script_argument))
root <- dirname(dirname(script_path))
cohort_path <- file.path(root, "data", "processed", "analysis_cohort.csv")
table_directory <- file.path(root, "outputs", "tables")
figure_directory <- file.path(root, "outputs", "figures")
document_directory <- file.path(root, "docs")

if (!requireNamespace("survey", quietly = TRUE)) {
  stop("Missing R package 'survey'. Run: Rscript scripts/install_r_dependencies.R")
}

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

frame <- read.csv(cohort_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "SEQN", "LBXSC3SI", "pral_day1", "pral_day2", "pral_mean",
  "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
  "BMXBMI", "mean_energy", "smoking_status", "diabetes_status",
  "egfr", "URDACT", "WTDR2D", "SDMVSTRA", "SDMVPSU",
  "fully_adjusted_with_uacr"
)
missing_columns <- setdiff(required_columns, names(frame))
if (length(missing_columns) > 0) {
  stop(paste("Analysis cohort is missing columns:", paste(missing_columns, collapse = ", ")))
}
if (anyDuplicated(frame$SEQN)) {
  stop("Analysis cohort contains duplicate participant identifiers.")
}

frame$eligible <- frame$fully_adjusted_with_uacr %in% c(TRUE, "TRUE", "True", 1, "1")
frame$sex <- factor(frame$RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
frame$race_ethnicity <- factor(
  frame$RIDRETH3,
  levels = c(3, 4, 1, 2, 6, 7),
  labels = c(
    "Non-Hispanic White", "Non-Hispanic Black", "Mexican American",
    "Other Hispanic", "Non-Hispanic Asian", "Other or multiracial"
  )
)
frame$education <- factor(
  frame$DMDEDUC2,
  levels = c(5, 4, 3, 2, 1),
  labels = c(
    "College graduate or above", "Some college or associate degree",
    "High school graduate", "Grades 9-11", "Less than grade 9"
  )
)
frame$smoking <- factor(
  frame$smoking_status,
  levels = c("never", "former", "current"),
  labels = c("Never", "Former", "Current")
)
frame$diabetes <- factor(
  frame$diabetes_status,
  levels = c("no", "borderline", "yes"),
  labels = c("No", "Borderline", "Yes")
)

# Scaling makes the coefficients clinically readable without changing model fit.
frame$pral_day1_10 <- frame$pral_day1 / 10
frame$pral_day2_10 <- frame$pral_day2 / 10
frame$pral_mean_10 <- frame$pral_mean / 10
frame$mean_energy_1000 <- frame$mean_energy / 1000
frame$egfr_10 <- frame$egfr / 10
frame$uacr_100 <- frame$URDACT / 100

model_variables <- c(
  "LBXSC3SI", "pral_day1_10", "pral_day2_10", "pral_mean_10",
  "RIDAGEYR", "sex", "race_ethnicity", "education", "INDFMPIR",
  "BMXBMI", "mean_energy_1000", "smoking", "diabetes", "egfr_10",
  "uacr_100", "WTDR2D", "SDMVSTRA", "SDMVPSU"
)
analysis_frame <- frame[frame$eligible, model_variables]
if (nrow(analysis_frame) != 3321) {
  stop(sprintf("Expected 3,321 complete cases, found %s.", nrow(analysis_frame)))
}
if (anyNA(analysis_frame)) {
  stop("Complete-case regression variables contain missing values.")
}
numeric_analysis_values <- as.matrix(
  analysis_frame[vapply(analysis_frame, is.numeric, logical(1))]
)
if (any(!is.finite(numeric_analysis_values))) {
  stop("Regression variables contain nonfinite numeric values.")
}
if (any(analysis_frame$WTDR2D <= 0)) {
  stop("All survey weights must be positive.")
}

options(survey.lonely.psu = "adjust")
core_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  nest = TRUE,
  data = frame
)
analysis_design <- subset(core_design, eligible)
design_degrees_freedom <- survey::degf(analysis_design)

variables_added <- list(
  M1 = "PRAL only",
  M2 = "Add age and sex",
  M3 = "Add race/ethnicity, education, and income",
  M4 = "Add BMI and mean energy intake",
  M5 = "Add smoking and diabetes",
  M6 = "Add eGFR and UACR"
)
model_formulas <- list(
  M1 = LBXSC3SI ~ pral_mean_10,
  M2 = LBXSC3SI ~ pral_mean_10 + RIDAGEYR + sex,
  M3 = LBXSC3SI ~ pral_mean_10 + RIDAGEYR + sex + race_ethnicity + education + INDFMPIR,
  M4 = LBXSC3SI ~ pral_mean_10 + RIDAGEYR + sex + race_ethnicity + education + INDFMPIR + BMXBMI + mean_energy_1000,
  M5 = LBXSC3SI ~ pral_mean_10 + RIDAGEYR + sex + race_ethnicity + education + INDFMPIR + BMXBMI + mean_energy_1000 + smoking + diabetes,
  M6 = LBXSC3SI ~ pral_mean_10 + RIDAGEYR + sex + race_ethnicity + education + INDFMPIR + BMXBMI + mean_energy_1000 + smoking + diabetes + egfr_10 + uacr_100
)

extract_exposure <- function(fit, term, model_name, added_variables, exposure_label) {
  estimate <- unname(stats::coef(fit)[term])
  standard_error <- sqrt(stats::vcov(fit)[term, term])
  critical_value <- stats::qt(0.975, df = design_degrees_freedom)
  test_statistic <- estimate / standard_error
  data.frame(
    model = model_name,
    cumulative_model = paste(deparse(stats::formula(fit)), collapse = " "),
    variables_added = added_variables,
    exposure = exposure_label,
    participants = nrow(analysis_frame),
    design_degrees_freedom = design_degrees_freedom,
    estimate_per_10_mEq_day = estimate,
    standard_error = standard_error,
    ci_lower = estimate - critical_value * standard_error,
    ci_upper = estimate + critical_value * standard_error,
    t_statistic = test_statistic,
    p_value = 2 * stats::pt(-abs(test_statistic), df = design_degrees_freedom),
    outcome_unit = "mmol/L serum bicarbonate",
    stringsAsFactors = FALSE
  )
}

model_fits <- lapply(
  model_formulas,
  survey::svyglm,
  design = analysis_design,
  family = stats::gaussian()
)
model_rows <- Map(
  function(fit, model_name, added_variables) {
    extract_exposure(
      fit, "pral_mean_10", model_name, added_variables,
      "Two-day mean PRAL"
    )
  },
  model_fits,
  names(model_fits),
  unname(variables_added)
)
model_sequence <- do.call(rbind, model_rows)
write.csv(
  model_sequence,
  file.path(table_directory, "pral_bicarbonate_model_sequence.csv"),
  row.names = FALSE
)

full_covariates <- paste(
  "RIDAGEYR + sex + race_ethnicity + education + INDFMPIR +",
  "BMXBMI + mean_energy_1000 + smoking + diabetes + egfr_10 + uacr_100"
)
exposure_terms <- c(
  "pral_day1_10" = "Day 1 PRAL",
  "pral_day2_10" = "Day 2 PRAL",
  "pral_mean_10" = "Two-day mean PRAL"
)
comparison_fits <- lapply(names(exposure_terms), function(term) {
  formula <- stats::as.formula(paste("LBXSC3SI ~", term, "+", full_covariates))
  survey::svyglm(formula, design = analysis_design, family = stats::gaussian())
})
names(comparison_fits) <- names(exposure_terms)
comparison_rows <- Map(function(fit, term) {
  extract_exposure(
    fit, term, "M6", "All prespecified covariates", exposure_terms[[term]]
  )
}, comparison_fits, names(comparison_fits))
exposure_comparison <- do.call(rbind, comparison_rows)
write.csv(
  exposure_comparison,
  file.path(table_directory, "pral_bicarbonate_exposure_comparison.csv"),
  row.names = FALSE
)

weighted_means <- survey::svymean(
  ~LBXSC3SI + pral_day1 + pral_day2 + pral_mean,
  analysis_design,
  na.rm = TRUE
)
descriptive_statistics <- data.frame(
  variable = c("Serum bicarbonate", "PRAL Day 1", "PRAL Day 2", "Two-day mean PRAL"),
  survey_weighted_mean = unname(stats::coef(weighted_means)),
  standard_error = unname(survey::SE(weighted_means)),
  participants = nrow(analysis_frame),
  stringsAsFactors = FALSE
)
write.csv(
  descriptive_statistics,
  file.path(table_directory, "regression_sample_descriptive_statistics.csv"),
  row.names = FALSE
)

plot_rows <- exposure_comparison[c(1, 2, 3), ]
plot_y <- c(3, 2, 1)
grDevices::png(
  file.path(figure_directory, "pral_bicarbonate_exposure_comparison.png"),
  width = 2100,
  height = 1350,
  res = 300
)
graphics::par(mar = c(5, 9, 3, 2))
plot_range <- range(c(plot_rows$ci_lower, plot_rows$ci_upper, 0))
graphics::plot(
  plot_rows$estimate_per_10_mEq_day,
  plot_y,
  xlim = plot_range,
  ylim = c(0.5, 3.5),
  yaxt = "n",
  ylab = "",
  xlab = "Difference in serum bicarbonate per 10 mEq/day higher PRAL (mmol/L)",
  pch = 19,
  col = "#2F6B8A",
  cex = 1.2
)
graphics::segments(
  plot_rows$ci_lower,
  plot_y,
  plot_rows$ci_upper,
  plot_y,
  col = "#2F6B8A",
  lwd = 2
)
graphics::abline(v = 0, lty = 2, col = "#8F3A4C")
graphics::axis(2, at = plot_y, labels = plot_rows$exposure, las = 1)
graphics::title("Fully Adjusted PRAL-Bicarbonate Association")
graphics::box()
grDevices::dev.off()

validation <- data.frame(
  check = c(
    "source_file", "core_participants", "analysis_participants",
    "participant_id_unique", "regression_values_complete", "all_weights_positive",
    "survey_strata", "survey_psus", "design_degrees_freedom",
    "all_sequence_models_same_sample", "all_exposure_models_same_sample",
    "survey_package_version"
  ),
  value = c(
    "data/processed/analysis_cohort.csv", nrow(frame), nrow(analysis_frame),
    !anyDuplicated(frame$SEQN), !anyNA(analysis_frame), all(analysis_frame$WTDR2D > 0),
    length(unique(analysis_frame$SDMVSTRA)),
    nrow(unique(analysis_frame[c("SDMVSTRA", "SDMVPSU")])),
    design_degrees_freedom,
    all(vapply(model_fits, function(fit) length(stats::residuals(fit)), integer(1)) == nrow(analysis_frame)),
    all(vapply(comparison_fits, function(fit) length(stats::residuals(fit)), integer(1)) == nrow(analysis_frame)),
    as.character(utils::packageVersion("survey"))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  validation,
  file.path(document_directory, "regression_validation.csv"),
  row.names = FALSE
)

primary <- model_sequence[model_sequence$model == "M6", ]
cat(sprintf("Participants: %d\n", primary$participants))
cat(sprintf("Design degrees of freedom: %d\n", primary$design_degrees_freedom))
cat(sprintf(
  "Fully adjusted two-day mean PRAL estimate: %.4f mmol/L per 10 mEq/day\n",
  primary$estimate_per_10_mEq_day
))
cat(sprintf("95%% CI: %.4f to %.4f\n", primary$ci_lower, primary$ci_upper))
cat(sprintf("p-value: %.6g\n", primary$p_value))