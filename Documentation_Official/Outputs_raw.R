# Outputs_raw.R
# Print numerical results, figure inputs, model summaries, and tests to console.

output_script_directory <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]

  if (length(frame_files) > 0L) {
    return(dirname(normalizePath(tail(frame_files, 1L), mustWork = TRUE)))
  }

  command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(command_file) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", command_file[1L]), mustWork = TRUE)))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

PROJECT_DIR <- output_script_directory()
CHECKPOINT_DIR <- file.path(PROJECT_DIR, "checkpoints")

required_packages <- c("broom", "dplyr", "stringr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the missing packages before running Outputs_raw.R\n",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

checkpoint_files <- c(
  import = "01_import.rds",
  cleaning = "02_cleaning.rds",
  descriptives = "03_descriptives.rds",
  baseline = "04_baseline.rds",
  robustness = "05_robustness.rds"
)
checkpoint_paths <- stats::setNames(
  file.path(CHECKPOINT_DIR, unname(checkpoint_files)),
  names(checkpoint_files)
)

missing_checkpoints <- checkpoint_paths[!file.exists(checkpoint_paths)]
if (length(missing_checkpoints) > 0L) {
  stop(
    "Run scripts 01 through 05 before Outputs_raw.R. Missing checkpoints\n",
    paste(basename(missing_checkpoints), collapse = ", "),
    call. = FALSE
  )
}

import_checkpoint <- readRDS(checkpoint_paths[["import"]])
cleaning <- readRDS(checkpoint_paths[["cleaning"]])
descriptives <- readRDS(checkpoint_paths[["descriptives"]])
baseline <- readRDS(checkpoint_paths[["baseline"]])
robustness <- readRDS(checkpoint_paths[["robustness"]])

as_plain_data_frame <- function(data) {
  output <- as.data.frame(data, stringsAsFactors = FALSE)
  factor_columns <- vapply(output, is.factor, logical(1))
  output[factor_columns] <- lapply(output[factor_columns], as.character)
  rownames(output) <- NULL
  output
}

tidy_one_model <- function(model, model_name) {
  output <- tryCatch(
    broom::tidy(model, conf.int = TRUE),
    error = function(error) {
      data.frame(
        term = NA_character_,
        error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    }
  )
  output$model <- model_name
  output <- output[c("model", setdiff(names(output), "model"))]
  as_plain_data_frame(output)
}

tidy_named_models <- function(models, prefix = NULL) {
  dplyr::bind_rows(lapply(names(models), function(model_name) {
    label <- if (is.null(prefix)) {
      model_name
    } else {
      paste(prefix, model_name, sep = "_")
    }
    tidy_one_model(models[[model_name]], label)
  }))
}

extract_saturated_raw <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE) |>
    dplyr::filter(grepl("^treated_factor::.*:Post$", term)) |>
    dplyr::transmute(
      model = model_name,
      Country_UCAS = term |>
        stringr::str_remove("^treated_factor::") |>
        stringr::str_remove(":Post$"),
      term,
      raw_log_effect = estimate,
      std_error = std.error,
      statistic,
      p_value = p.value,
      conf_low = conf.low,
      conf_high = conf.high
    )
}

panel_main <- cleaning$panel_main
eu_real <- import_checkpoint$eu_real
loan_parameters <- robustness$loan_parameters

fee_pre <- 9250
fee_post <- 22200
log_fee_change <- log(fee_post) - log(fee_pre)
arithmetic_fee_change <- (fee_post - fee_pre) / fee_pre

pooled_coefficients <- stats::coef(baseline$models$pooled_DK)
if (!"EU:Post" %in% names(pooled_coefficients)) {
  stop("The pooled baseline model has no EU:Post coefficient.", call. = FALSE)
}
beta_pooled <- unname(pooled_coefficients[["EU:Post"]])

# Figure 1. Pooled event-study coefficients on the log-applicant scale.
figure_01_event <- baseline$event_coefficients |>
  dplyr::mutate(scale = "raw log-applicant coefficient")

# Figure 2. Leads-only coefficients on the log-applicant scale.
figure_02_leads <- broom::tidy(
  robustness$placebo_models$leads_only,
  conf.int = TRUE
) |>
  dplyr::filter(grepl("^event_time::", term)) |>
  dplyr::mutate(
    event_time = as.integer(sub(
      "^event_time::(-?[0-9]+):EU$",
      "\\1",
      term
    )),
    scale = "raw log-applicant coefficient"
  ) |>
  dplyr::arrange(event_time)

# Figure 3. Country-saturated event-study coefficients.
figure_03_saturated_event <- baseline$saturated_event_coefficients |>
  dplyr::mutate(scale = "raw log-applicant coefficient")

# Figure 4. Raw country effects plus each plotted elasticity transformation.
fee_mechanical <- panel_main |>
  dplyr::filter(Domicile_named_country %in% eu_real) |>
  dplyr::group_by(Domicile_named_country) |>
  dplyr::summarise(
    applicants_2020 = mean(Applicants[Year == 2020L], na.rm = TRUE),
    applicants_2021 = mean(Applicants[Year == 2021L], na.rm = TRUE),
    raw_application_change = applicants_2021 - applicants_2020,
    proportional_application_change = raw_application_change / applicants_2020,
    mechanical_fee_elasticity = proportional_application_change /
      arithmetic_fee_change,
    .groups = "drop"
  ) |>
  dplyr::rename(Country_UCAS = Domicile_named_country)

figure_04_decomposition <- baseline$baseline_country_effects |>
  dplyr::mutate(
    raw_log_effect = Estimate_TWFE,
    log_fee_change = log_fee_change,
    arithmetic_fee_change = arithmetic_fee_change,
    structural_log_elasticity = Estimate_TWFE / log_fee_change,
    structural_log_elasticity_low = conf_low / log_fee_change,
    structural_log_elasticity_high = conf_high / log_fee_change,
    did_arc_elasticity = (exp(Estimate_TWFE) - 1) /
      arithmetic_fee_change
  ) |>
  dplyr::left_join(fee_mechanical, by = "Country_UCAS")

# Figure 5. Robust income-gradient coefficients and their country inputs.
figure_05_gradient_coefficients <- robustness$gradient_coefficients
figure_05_country_inputs <- robustness$elas_df |>
  dplyr::select(
    dplyr::any_of(c(
      "Country", "iso2c", "GDP_PPP", "log_gdp", "Estimate_TWFE",
      "Estimate_controls", "Estimate_Pooled", "k_i",
      "elasticity_TWFE", "elasticity_new_TWFE",
      "elasticity_controls", "elasticity_controls_new",
      "elasticity_Pooled", "elasticity_new_Pooled"
    ))
  )

# Figure 6. Raw COVID-split log effects and plotted elasticities.
figure_06_covid <- baseline$covid_coefficients |>
  dplyr::mutate(
    raw_log_effect = estimate,
    log_fee_change = log_fee_change,
    plotted_elasticity = estimate / log_fee_change,
    plotted_elasticity_low = conf.low / log_fee_change,
    plotted_elasticity_high = conf.high / log_fee_change
  )

# Figure 7. Raw country effects from every saturated control specification.
figure_07_controlled <- dplyr::bind_rows(lapply(
  names(robustness$saturated_control_models),
  function(model_name) {
    extract_saturated_raw(
      robustness$saturated_control_models[[model_name]],
      model_name
    )
  }
)) |>
  dplyr::mutate(
    log_fee_change = log_fee_change,
    plotted_elasticity = raw_log_effect / log_fee_change,
    plotted_elasticity_low = conf_low / log_fee_change,
    plotted_elasticity_high = conf_high / log_fee_change
  )

# Figure 8. Pooled robustness effects on the original log-applicant scale.
figure_08_pooled <- dplyr::bind_rows(lapply(
  names(robustness$pooled_robustness_models),
  function(sample_name) {
    dplyr::bind_rows(lapply(
      names(robustness$pooled_robustness_models[[sample_name]]),
      function(weighting) {
        broom::tidy(
          robustness$pooled_robustness_models[[sample_name]][[weighting]],
          conf.int = TRUE
        ) |>
          dplyr::filter(term == "EU:Post") |>
          dplyr::mutate(
            sample = sample_name,
            weighting = weighting,
            model = paste(sample_name, weighting, sep = "_")
          )
      }
    ))
  }
)) |>
  dplyr::bind_rows(
    broom::tidy(robustness$m_ireland_only, conf.int = TRUE) |>
      dplyr::filter(term == "EU:Post") |>
      dplyr::mutate(
        sample = "ireland",
        weighting = "unweighted",
        model = "ireland_only"
      )
  ) |>
  dplyr::mutate(scale = "raw log-applicant coefficient")

# Figure 9. Country-saturated robustness effects on their original scale.
figure_09_saturated <- dplyr::bind_rows(lapply(
  names(robustness$saturated_robustness_models),
  function(sample_name) {
    dplyr::bind_rows(lapply(
      names(robustness$saturated_robustness_models[[sample_name]]),
      function(weighting) {
        extract_saturated_raw(
          robustness$saturated_robustness_models[[sample_name]][[weighting]],
          paste(sample_name, weighting, sep = "_")
        ) |>
          dplyr::mutate(sample = sample_name, weighting = weighting)
      }
    ))
  }
)) |>
  dplyr::bind_rows(
    extract_saturated_raw(robustness$sat_ireland_only, "ireland_only") |>
      dplyr::mutate(sample = "ireland", weighting = "unweighted")
  )

# Figure 10. Calibration values, slopes, p-values, and sample sizes.
figure_10_calibration <- robustness$sweep_all

# Figure 11. Raw pooled effect and every loan-value denominator in the curve.
pv_grid <- seq(500, 9250, by = 50)
figure_11_pv_curve <- data.frame(
  PV_loan = pv_grid,
  raw_pooled_log_effect = beta_pooled,
  post_brexit_price = as.numeric(loan_parameters$Ppost),
  effective_log_price_change = log(as.numeric(loan_parameters$Ppost)) -
    log(pv_grid),
  plotted_elasticity_magnitude = abs(beta_pooled) /
    (log(as.numeric(loan_parameters$Ppost)) - log(pv_grid)),
  stringsAsFactors = FALSE
)
figure_11_salary_points <- robustness$aggregate_pv_table |>
  dplyr::mutate(raw_pooled_log_effect = beta_pooled)

# Figure 12. Raw country effects and both common-k and corrected elasticities.
figure_12_gradient_shift <- robustness$elas_df |>
  dplyr::select(
    dplyr::any_of(c(
      "Country", "iso2c", "Wage_EURO", "wage_gbp", "PV_loan", "k_i",
      "GDP_PPP", "log_gdp", "Estimate_TWFE", "Estimate_controls",
      "Estimate_Pooled", "elasticity_TWFE", "elasticity_new_TWFE",
      "elasticity_controls", "elasticity_controls_new",
      "elasticity_Pooled", "elasticity_new_Pooled"
    ))
  ) |>
  dplyr::mutate(common_log_price_change = as.numeric(loan_parameters$k_common)) |>
  dplyr::arrange(log_gdp)

# Figure 13. Robust standardized slopes and the country-level inputs.
figure_13_correlations <- robustness$standardized_correlations
figure_13_country_inputs <- robustness$elas_df |>
  dplyr::select(
    dplyr::any_of(c(
      "Country", "elasticity_new_TWFE", "GDP_PPP", "gini",
      "youth_unemp", "educ_spend_gdp", "tertiary_enroll",
      "net_migration"
    ))
  )

# Additional diagnostic figures produced by 06_figures.R.
extra_sdid_comparison <- robustness$sdid_comparison
extra_placebo_denominator <- robustness$placebo_summary
extra_price_compression <- robustness$elas_df |>
  dplyr::select(dplyr::any_of(c("Country", "k_i"))) |>
  dplyr::mutate(
    k_common = as.numeric(loan_parameters$k_common),
    difference_from_common_k = k_i - k_common
  )

# Raw model tables and numbered manuscript table inputs.
baseline_model_coefficients <- tidy_named_models(baseline$models)
control_model_coefficients <- tidy_named_models(c(
  robustness$control_models,
  robustness$saturated_control_models
))
gradient_model_coefficients <- tidy_named_models(robustness$gradient_models)
placebo_model_coefficients <- tidy_named_models(
  robustness$placebo_models,
  prefix = "placebo"
)

table_01_descriptives <- descriptives$desc_wb |>
  dplyr::select(Variable, n, mean, sd, median, min, max)
table_02_donor_weights <- cleaning$weights
table_03_control_models <- control_model_coefficients
table_04_control_groups <- robustness$control_group_composition
table_05_estimator_comparison <- robustness$sdid_comparison

parameter_table <- data.frame(
  parameter = c(
    "fee_pre", "fee_post", "log_fee_change", "arithmetic_fee_change",
    paste0("loan_", names(loan_parameters))
  ),
  value = c(
    fee_pre,
    fee_post,
    log_fee_change,
    arithmetic_fee_change,
    as.numeric(unlist(loan_parameters))
  ),
  stringsAsFactors = FALSE
)

manifest <- data.frame(
  dataset = c(
    "F01_event_raw", "F02_leads_raw", "F03_sat_event_raw",
    "F04_decomposition", "F05_gradient_coefs", "F05_country_inputs",
    "F06_covid_raw", "F07_controls_raw", "F08_pooled_raw",
    "F09_saturated_raw", "F10_calibration", "F11_pv_curve",
    "F11_salary_points", "F12_gradient_shift", "F13_correlations",
    "F13_country_inputs", "X_sdid_comparison", "X_placebo_denominator",
    "X_price_compression", "M_baseline_all", "M_controls_all",
    "M_gradient_all", "M_placebo_all", "T01_descriptives",
    "T02_donor_weights", "T03_control_models", "T04_control_groups",
    "T05_estimator_compare", "Parameters", "Data_source_roles"
  ),
  description = c(
    "Figure 1 pooled event-study coefficients",
    "Figure 2 leads-only coefficients",
    "Figure 3 country-saturated event-study coefficients",
    "Figure 4 raw effects and elasticity transformations",
    "Figure 5 robust income-gradient slopes",
    "Figure 5 country-level regression inputs",
    "Figure 6 raw COVID-split effects and plotted elasticities",
    "Figure 7 raw country effects from control specifications",
    "Figure 8 pooled robustness coefficients",
    "Figure 9 country-saturated robustness coefficients",
    "Figure 10 calibration slopes, p-values, and sample sizes",
    "Figure 11 loan present-value curve inputs",
    "Figure 11 selected starting-salary points",
    "Figure 12 raw effects and common-k or corrected elasticities",
    "Figure 13 standardized robust correlations",
    "Figure 13 country-level correlation inputs",
    "Additional synthetic DiD comparison figure inputs",
    "Additional placebo-denominator figure inputs",
    "Additional price-compression figure inputs",
    "All baseline-model coefficient output",
    "All macroeconomic control-model coefficient output",
    "All income-gradient model coefficient output",
    "All placebo and leads model coefficient output",
    "Manuscript Table 1 raw input",
    "Manuscript Table 2 raw input",
    "Manuscript Table 3 raw coefficient input",
    "Manuscript Table 4 raw input",
    "Manuscript Table 5 raw input",
    "Fee and loan assumptions used in transformations",
    "Applicant file roles and years used"
  ),
  stringsAsFactors = FALSE
)

raw_tables <- list(
  README = manifest,
  F01_event_raw = figure_01_event,
  F02_leads_raw = figure_02_leads,
  F03_sat_event_raw = figure_03_saturated_event,
  F04_decomposition = figure_04_decomposition,
  F05_gradient_coefs = figure_05_gradient_coefficients,
  F05_country_inputs = figure_05_country_inputs,
  F06_covid_raw = figure_06_covid,
  F07_controls_raw = figure_07_controlled,
  F08_pooled_raw = figure_08_pooled,
  F09_saturated_raw = figure_09_saturated,
  F10_calibration = figure_10_calibration,
  F11_pv_curve = figure_11_pv_curve,
  F11_salary_points = figure_11_salary_points,
  F12_gradient_shift = figure_12_gradient_shift,
  F13_correlations = figure_13_correlations,
  F13_country_inputs = figure_13_country_inputs,
  X_sdid_comparison = extra_sdid_comparison,
  X_placebo_denominator = extra_placebo_denominator,
  X_price_compression = extra_price_compression,
  M_baseline_all = baseline_model_coefficients,
  M_controls_all = control_model_coefficients,
  M_gradient_all = gradient_model_coefficients,
  M_placebo_all = placebo_model_coefficients,
  T01_descriptives = table_01_descriptives,
  T02_donor_weights = table_02_donor_weights,
  T03_control_models = table_03_control_models,
  T04_control_groups = table_04_control_groups,
  T05_estimator_compare = table_05_estimator_comparison,
  Parameters = parameter_table,
  Data_source_roles = cleaning$data_source_roles
)

raw_tables <- lapply(raw_tables, as_plain_data_frame)

print_console_section <- function(title, object) {
  separator <- paste(rep("=", 78L), collapse = "")
  cat("\n\n", separator, "\n", title, "\n", separator, "\n", sep = "")
  print(object)
  invisible(object)
}

print_model_collection <- function(section, models) {
  for (model_name in names(models)) {
    print_console_section(
      paste(section, model_name, sep = " | "),
      summary(models[[model_name]])
    )
  }
  invisible(models)
}

model_objects <- list(
  baseline = baseline$models,
  pooled_robustness = robustness$pooled_robustness_models,
  saturated_robustness = robustness$saturated_robustness_models,
  ireland_pooled = robustness$m_ireland_only,
  ireland_saturated = robustness$sat_ireland_only,
  controls = robustness$control_models,
  saturated_controls = robustness$saturated_control_models,
  gradients = robustness$gradient_models,
  placebo = robustness$placebo_models,
  honest_did = robustness$honest_did
)

test_objects <- list(
  pretrend = baseline$pretrend_test,
  leads = robustness$leads_test,
  interactions = robustness$interaction_test
)

aggregate_sdid <- data.frame(
  estimate = as.numeric(robustness$tau_hat),
  std_error = as.numeric(robustness$tau_se),
  conf_low = as.numeric(robustness$tau_hat) - 1.96 * as.numeric(robustness$tau_se),
  conf_high = as.numeric(robustness$tau_hat) + 1.96 * as.numeric(robustness$tau_se)
)

Outputs_raw <- list(
  manifest = manifest,
  tables = raw_tables,
  model_objects = model_objects,
  tests = test_objects,
  aggregate_sdid = aggregate_sdid
)

print_console_section("OUTPUT MANIFEST", manifest)
for (table_name in setdiff(names(raw_tables), "README")) {
  print_console_section(table_name, raw_tables[[table_name]])
}

print_console_section("AGGREGATE SYNTHETIC DID", aggregate_sdid)
print_model_collection("BASELINE MODEL", baseline$models)
print_model_collection("CONTROL MODEL", robustness$control_models)
print_model_collection(
  "SATURATED CONTROL MODEL",
  robustness$saturated_control_models
)
print_model_collection("GRADIENT MODEL", robustness$gradient_models)
print_model_collection("PLACEBO MODEL", robustness$placebo_models)

for (sample_name in names(robustness$pooled_robustness_models)) {
  print_model_collection(
    paste("POOLED ROBUSTNESS", sample_name),
    robustness$pooled_robustness_models[[sample_name]]
  )
}
print_console_section(
  "POOLED ROBUSTNESS | ireland_only",
  summary(robustness$m_ireland_only)
)

for (sample_name in names(robustness$saturated_robustness_models)) {
  print_model_collection(
    paste("SATURATED ROBUSTNESS", sample_name),
    robustness$saturated_robustness_models[[sample_name]]
  )
}
print_console_section(
  "SATURATED ROBUSTNESS | ireland_only",
  summary(robustness$sat_ireland_only)
)

print_console_section("PRETREND TEST", baseline$pretrend_test)
print_console_section("LEADS TEST", robustness$leads_test)
print_console_section("CONTROL INTERACTION TEST", robustness$interaction_test)
print_console_section("HONEST DID", robustness$honest_did)

cat("\n\nOutputs_raw.R complete. All objects remain available in `Outputs_raw`.\n")
invisible(Outputs_raw)
