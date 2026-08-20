# 03_descriptives.R
# Produce auditable descriptive summaries before estimating any model.

if (!exists("panel_main", inherits = FALSE) || !exists("panel_wb", inherits = FALSE)) {
  stop("Run 01_import.R and 02_cleaning.R first.", call. = FALSE)
}

describe_frame <- function(data, variables) {
  selected <- intersect(variables, names(data))
  description <- psych::describe(data[selected], fast = FALSE)
  data.frame(
    Variable = rownames(description),
    description,
    row.names = NULL,
    check.names = FALSE
  )
}

main_vars <- c(
  "Applicants", "log_app", "EU", "Post", "event_time",
  "regression_weight"
)

wb_vars <- c(
  "Applicants", "log_app", "EU", "Post", "did", "event_time",
  "regression_weight", "GDPpc", "Pop15_65", "YouthUnempl",
  "EduSpendGDP", "c_log_GDPpc", "c_log_Pop15_65", "c_YouthUnempl"
)

desc_main <- describe_frame(panel_main, main_vars)
desc_wb <- describe_frame(panel_wb, wb_vars)

sample_overview <- dplyr::bind_rows(
  panel_main |>
    dplyr::summarise(
      sample = "Baseline panel",
      observations = dplyr::n(),
      countries = dplyr::n_distinct(Domicile_named_country),
      treated_countries = dplyr::n_distinct(Domicile_named_country[EU == 1L]),
      control_countries = dplyr::n_distinct(Domicile_named_country[EU == 0L]),
      first_year = min(Year),
      last_year = max(Year)
    ),
  panel_wb |>
    dplyr::summarise(
      sample = "World Bank complete-case panel",
      observations = dplyr::n(),
      countries = dplyr::n_distinct(Domicile_named_country),
      treated_countries = dplyr::n_distinct(Domicile_named_country[EU == 1L]),
      control_countries = dplyr::n_distinct(Domicile_named_country[EU == 0L]),
      first_year = min(Year),
      last_year = max(Year)
    )
)

group_year_summary <- panel_main |>
  dplyr::group_by(EU_group, Year) |>
  dplyr::summarise(
    countries = dplyr::n_distinct(Domicile_named_country),
    applicants_mean = mean(Applicants),
    applicants_sd = stats::sd(Applicants),
    log_app_mean = stats::weighted.mean(log_app, regression_weight),
    .groups = "drop"
  )

donor_weight_summary <- weights |>
  dplyr::summarise(
    donors = dplyr::n(),
    sum_normalized_weights = sum(normalized_weight),
    sum_regression_weights = sum(regression_weight),
    minimum_weight = min(regression_weight),
    maximum_weight = max(regression_weight)
  )

openxlsx::write.xlsx(
  list(
    sample_overview = sample_overview,
    baseline_descriptives = desc_main,
    controls_descriptives = desc_wb,
    group_by_year = group_year_summary,
    donor_weights = weights,
    donor_weight_checks = donor_weight_summary
  ),
  file.path(OUTPUT_DIR, "03_descriptives.xlsx"),
  overwrite = TRUE
)

descriptive_checkpoint <- list(
  desc_main = desc_main,
  desc_wb = desc_wb,
  sample_overview = sample_overview,
  group_year_summary = group_year_summary,
  donor_weight_summary = donor_weight_summary
)

saveRDS(descriptive_checkpoint, file.path(CHECKPOINT_DIR, "03_descriptives.rds"))

print(sample_overview)
message("03_descriptives.R complete")

