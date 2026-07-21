# 04_baseline.R
# Estimate the pooled and country-saturated baseline designs.

if (!exists("panel_main", inherits = FALSE)) {
  stop("Run 01_import.R, 02_cleaning.R, and 03_descriptives.R first.", call. = FALSE)
}

m_pooled_DK <- fixest::feols(
  log_app ~ EU:Post | Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m_sat_DK <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") |
    Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m_sat_DK1 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") |
    Domicile_named_country,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m_sat_DK2 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL"),
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m_sat_DK3 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") |
    Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK"
)

m_event <- fixest::feols(
  log_app ~ i(event_time, EU, ref = 0) |
    Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = ~ Domicile_named_country,
  weights = ~ regression_weight
)

m_sat_event <- fixest::feols(
  log_app ~ i(event_time, treated_factor, ref = 0, ref2 = "CONTROL") |
    Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = ~ Domicile_named_country,
  weights = ~ regression_weight
)

panel_covid <- panel_main |>
  dplyr::mutate(
    Post = as.integer(Year %in% c(2021L, 2022L, 2023L)),
    Post_COVID = as.integer(Year >= 2024L)
  )

m_sat_DK_2 <- fixest::feols(
  log_app ~
    i(treated_factor, Post, ref = "CONTROL") +
    i(treated_factor, Post_COVID, ref = "CONTROL") |
    Domicile_named_country + Year,
  data = panel_covid,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

extract_saturated <- function(model, estimate_name = "estimate") {
  output <- broom::tidy(model, conf.int = TRUE) |>
    dplyr::filter(
      grepl("^treated_factor::", term),
      grepl(":Post$", term)
    ) |>
    dplyr::transmute(
      Country_UCAS = sub(
        ":Post$",
        "",
        sub("^treated_factor::", "", term)
      ),
      estimate = estimate,
      std_error = std.error,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value
    )

  names(output)[names(output) == "estimate"] <- estimate_name
  output
}

baseline_country_effects <- extract_saturated(m_sat_DK, "Estimate_TWFE") |>
  dplyr::left_join(
    extract_saturated(m_sat_DK2, "Estimate_Pooled") |>
      dplyr::select(Country_UCAS, Estimate_Pooled),
    by = "Country_UCAS"
  )
pooled_effect <- unname(stats::coef(m_pooled_DK)[["EU:Post"]])

event_coefficients <- broom::tidy(m_event, conf.int = TRUE) |>
  dplyr::filter(grepl("^event_time::", term)) |>
  dplyr::mutate(
    event_time = as.integer(sub(
      "^event_time::(-?[0-9]+):EU$",
      "\\1",
      term
    ))
  ) |>
  dplyr::arrange(event_time)

saturated_event_coefficients <- broom::tidy(m_sat_event, conf.int = TRUE) |>
  dplyr::filter(grepl("^event_time::", term)) |>
  dplyr::mutate(
    event_time = as.integer(stringr::str_extract(term, "-?[0-9]+(?=:)") ),
    Country_UCAS = stringr::str_extract(term, "(?<=treated_factor::).*$")
  ) |>
  dplyr::filter(!is.na(Country_UCAS), Country_UCAS != "CONTROL")

covid_coefficients <- broom::tidy(m_sat_DK_2, conf.int = TRUE) |>
  dplyr::filter(grepl("^treated_factor::", term)) |>
  dplyr::mutate(
    period = dplyr::case_when(
      grepl(":Post_COVID$", term) ~ "Post-COVID",
      grepl(":Post$", term) ~ "Post-Brexit",
      TRUE ~ NA_character_
    ),
    Country_UCAS = term |>
      stringr::str_remove("^treated_factor::") |>
      stringr::str_remove(":(Post|Post_COVID)$")
  ) |>
  dplyr::filter(!is.na(period))

pretrend_terms <- paste0("event_time::", -4:-1, ":EU")
pretrend_terms <- intersect(pretrend_terms, names(stats::coef(m_event)))

pretrend_test <- if (length(pretrend_terms) > 0L) {
  tryCatch(
    fixest::wald(m_event, pretrend_terms),
    error = function(error) error
  )
} else {
  NULL
}

baseline_models <- list(
  pooled_DK = m_pooled_DK,
  saturated_DK = m_sat_DK,
  saturated_unit_FE = m_sat_DK1,
  saturated_no_FE = m_sat_DK2,
  saturated_unweighted = m_sat_DK3,
  event_study = m_event,
  saturated_event_study = m_sat_event,
  covid_split = m_sat_DK_2
)

baseline_model_tables <- lapply(baseline_models, broom::tidy, conf.int = TRUE)
names(baseline_model_tables) <- substr(names(baseline_model_tables), 1L, 31L)

openxlsx::write.xlsx(
  c(
    baseline_model_tables,
    list(
      country_effects = baseline_country_effects,
      event_coefficients = event_coefficients,
      covid_coefficients = covid_coefficients
    )
  ),
  file.path(OUTPUT_DIR, "04_baseline_models.xlsx"),
  overwrite = TRUE
)

baseline_checkpoint <- list(
  models = baseline_models,
  baseline_country_effects = baseline_country_effects,
  event_coefficients = event_coefficients,
  saturated_event_coefficients = saturated_event_coefficients,
  covid_coefficients = covid_coefficients,
  pretrend_test = pretrend_test
)

saveRDS(baseline_checkpoint, file.path(CHECKPOINT_DIR, "04_baseline.rds"))

print(summary(m_pooled_DK))
print(summary(m_sat_DK))
message("04_baseline.R complete")
