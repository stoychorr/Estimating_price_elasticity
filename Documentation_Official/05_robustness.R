# 05_robustness.R
# Run alternative controls, sensitivity checks, synthetic DiD, and loan corrections.

if (!exists("baseline_models", inherits = FALSE) || !exists("panel_wb", inherits = FALSE)) {
  stop("Run scripts 01 through 04 first.", call. = FALSE)
}

m_pooled_DK <- baseline_models$pooled_DK
m_sat_DK <- baseline_models$saturated_DK
m_event <- baseline_models$event_study

# Placebo and leads-only checks
m_placebo_2018 <- fixest::feols(
  log_app ~ i(I(Year - 2018), EU, ref = -1) |
    Domicile_named_country + Year,
  data = panel_main,
  panel.id = ~ Domicile_named_country + Year,
  vcov = ~ Domicile_named_country,
  weights = ~ regression_weight
)

panel_pre <- panel_main |>
  dplyr::filter(Year < TREATMENT_YEAR)

m_leads_only <- fixest::feols(
  log_app ~ i(event_time, EU, ref = 0) |
    Domicile_named_country + Year,
  data = panel_pre,
  panel.id = ~ Domicile_named_country + Year,
  vcov = ~ Domicile_named_country,
  weights = ~ regression_weight
)

leads_terms <- paste0("event_time::", -4:-1, ":EU")
leads_terms <- intersect(leads_terms, names(stats::coef(m_leads_only)))
leads_test <- if (length(leads_terms) > 0L) {
  tryCatch(
    fixest::wald(m_leads_only, leads_terms),
    error = function(error) error
  )
} else {
  NULL
}

# Alternative donor pools
synth_poor <- c(
  "Pakistan", "Nepal", "Bangladesh", "Ethiopia", "Uganda",
  "Tanzania, United Republic of", "Rwanda", "Mozambique", "Malawi",
  "Sierra Leone", "Ghana", "Nigeria", "Kenya", "Zimbabwe"
)

synth_rich <- c(
  "United States of America", "Canada", "Australia", "Norway", "Switzerland",
  "Singapore", "Qatar", "United Arab Emirates", "Kuwait", "Saudi Arabia",
  "Israel", "Japan", "Hong Kong"
)

oecd <- c(
  "Australia", "Canada", "Chile", "Colombia", "Iceland", "Israel", "Japan",
  "Korea, Republic of", "Mexico", "New Zealand", "Norway", "Switzerland",
  "Turkey", "United States of America"
)

make_sample <- function(control_countries, label) {
  apps_clean |>
    dplyr::filter(Domicile_named_country %in% c(eu_real, control_countries)) |>
    dplyr::mutate(control_group = label)
}

alternative_samples <- list(
  baseline = make_sample(eu_synth, "Baseline synthetic"),
  poor = make_sample(synth_poor, "Very poor controls"),
  rich = make_sample(synth_rich, "Very rich controls"),
  oecd = make_sample(oecd, "OECD controls")
)

check_data <- function(data) {
  data.frame(
    observations = nrow(data),
    countries = dplyr::n_distinct(data$Domicile_named_country),
    eu_countries = dplyr::n_distinct(data$Domicile_named_country[data$EU == 1L]),
    control_countries = dplyr::n_distinct(data$Domicile_named_country[data$EU == 0L]),
    first_year = min(data$Year, na.rm = TRUE),
    last_year = max(data$Year, na.rm = TRUE)
  )
}

data_checks <- dplyr::bind_rows(
  lapply(alternative_samples, check_data),
  .id = "sample"
)

control_group_composition <- data.frame(
  group = c("Poor", "Rich", "OECD"),
  control_countries = c(
    paste(synth_poor, collapse = ", "),
    paste(synth_rich, collapse = ", "),
    paste(oecd, collapse = ", ")
  ),
  eu_countries = data_checks$eu_countries[match(c("poor", "rich", "oecd"), data_checks$sample)],
  observations = data_checks$observations[match(c("poor", "rich", "oecd"), data_checks$sample)],
  stringsAsFactors = FALSE
)

if (any(data_checks$control_countries == 0L)) {
  stop("At least one alternative donor pool has no observed controls.", call. = FALSE)
}

make_distance_weights <- function(data) {
  eu_pre <- data |>
    dplyr::filter(EU == 1L, Year < TREATMENT_YEAR) |>
    dplyr::group_by(Year) |>
    dplyr::summarise(
      eu_mean_log_app = mean(log_app, na.rm = TRUE),
      .groups = "drop"
    )

  control_distances <- data |>
    dplyr::filter(EU == 0L, Year < TREATMENT_YEAR) |>
    dplyr::select(Domicile_named_country, Year, log_app) |>
    dplyr::left_join(eu_pre, by = "Year") |>
    dplyr::group_by(Domicile_named_country) |>
    dplyr::summarise(
      distance = sqrt(mean((log_app - eu_mean_log_app)^2, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      donor_weight = 1 / (distance + 0.001),
      donor_weight = donor_weight / mean(donor_weight, na.rm = TRUE)
    )

  data |>
    dplyr::left_join(control_distances, by = "Domicile_named_country") |>
    dplyr::mutate(
      robustness_weight = dplyr::if_else(EU == 1L, 1, donor_weight),
      robustness_weight = dplyr::coalesce(robustness_weight, 1)
    )
}

weighted_samples <- lapply(alternative_samples, make_distance_weights)

run_weighted_unweighted <- function(data) {
  list(
    weighted = fixest::feols(
      log_app ~ EU:Post | Domicile_named_country + Year,
      data = data,
      weights = ~ robustness_weight,
      panel.id = ~ Domicile_named_country + Year,
      vcov = "DK"
    ),
    unweighted = fixest::feols(
      log_app ~ EU:Post | Domicile_named_country + Year,
      data = data,
      panel.id = ~ Domicile_named_country + Year,
      vcov = "DK"
    )
  )
}

pooled_robustness_models <- lapply(weighted_samples, run_weighted_unweighted)

panel_ireland_pooled <- apps_clean |>
  dplyr::filter(Domicile_named_country %in% c(eu_real, "Ireland")) |>
  dplyr::mutate(
    EU = as.integer(Domicile_named_country %in% eu_real),
    Post = as.integer(Year >= TREATMENT_YEAR),
    control_group = "Ireland"
  )

m_ireland_only <- fixest::feols(
  log_app ~ EU:Post | Domicile_named_country + Year,
  data = panel_ireland_pooled,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK"
)

prepare_saturated_data <- function(data) {
  data |>
    dplyr::mutate(
      treated_factor = factor(dplyr::if_else(
        Domicile_named_country %in% eu_real,
        Domicile_named_country,
        "CONTROL"
      ))
    )
}

saturated_samples <- lapply(weighted_samples, prepare_saturated_data)

run_saturated_models <- function(data) {
  list(
    weighted = fixest::feols(
      log_app ~ i(treated_factor, Post, ref = "CONTROL") |
        Domicile_named_country + Year,
      data = data,
      weights = ~ robustness_weight,
      panel.id = ~ Domicile_named_country + Year,
      vcov = "DK"
    ),
    unweighted = fixest::feols(
      log_app ~ i(treated_factor, Post, ref = "CONTROL") |
        Domicile_named_country + Year,
      data = data,
      panel.id = ~ Domicile_named_country + Year,
      vcov = "DK"
    )
  )
}

saturated_robustness_models <- lapply(saturated_samples, run_saturated_models)

data_ireland_sat <- panel_ireland_pooled |>
  dplyr::mutate(
    treated_factor = factor(dplyr::if_else(
      Domicile_named_country %in% eu_real,
      Domicile_named_country,
      "IRELAND"
    ))
  )

sat_ireland_only <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "IRELAND") |
    Domicile_named_country + Year,
  data = data_ireland_sat,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK"
)

# World Bank control specifications
m1 <- m_pooled_DK
m2 <- fixest::feols(
  log_app ~ did + c_log_GDPpc | Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m4 <- fixest::feols(
  log_app ~ did + c_log_GDPpc + c_log_Pop15_65 |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m5 <- fixest::feols(
  log_app ~ did + c_log_GDPpc + c_log_Pop15_65 + c_YouthUnempl |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

m6 <- fixest::feols(
  log_app ~ did * (c_log_GDPpc + c_log_Pop15_65 + c_YouthUnempl) |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

s1 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

s2 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") + c_log_GDPpc |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

s3 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") +
    c_log_GDPpc + c_log_Pop15_65 |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

s4 <- fixest::feols(
  log_app ~ i(treated_factor, Post, ref = "CONTROL") +
    c_log_GDPpc + c_log_Pop15_65 + c_YouthUnempl |
    Domicile_named_country + Year,
  data = panel_wb,
  panel.id = ~ Domicile_named_country + Year,
  vcov = "DK",
  weights = ~ regression_weight
)

control_models <- list(m1 = m1, m2 = m2, m4 = m4, m5 = m5, m6 = m6)
saturated_control_models <- list(s1 = s1, s2 = s2, s3 = s3, s4 = s4)

interaction_test <- tryCatch(
  fixest::wald(
    m6,
    c("did:c_log_GDPpc", "did:c_log_Pop15_65", "did:c_YouthUnempl")
  ),
  error = function(error) error
)

controlled_country_effects <- extract_saturated(s4, "Estimate_controls") |>
  dplyr::select(Country_UCAS, Estimate_controls)

# Synthetic difference-in-differences on the balanced full donor pool
all_years <- sort(unique(apps_clean$Year))
n_years <- length(all_years)

sdid_long <- apps_clean |>
  dplyr::select(Domicile_named_country, Year, log_app, EU) |>
  dplyr::filter(
    !is.na(Domicile_named_country),
    !is.na(Year),
    !is.na(EU)
  ) |>
  dplyr::group_by(Domicile_named_country, Year) |>
  dplyr::summarise(
    log_app = mean(log_app, na.rm = TRUE),
    EU = dplyr::first(EU),
    n_EU_values = dplyr::n_distinct(EU),
    .groups = "drop"
  )

if (!all(sdid_long$n_EU_values == 1L)) {
  stop("A country-year has conflicting EU classifications.", call. = FALSE)
}

sdid_long <- sdid_long |>
  dplyr::select(-n_EU_values) |>
  dplyr::group_by(Domicile_named_country) |>
  dplyr::filter(
    all(is.finite(log_app)),
    dplyr::n() == n_years,
    dplyr::n_distinct(Year) == n_years,
    all(all_years %in% Year)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(treated = as.integer(EU == 1L & Year >= TREATMENT_YEAR)) |>
  dplyr::arrange(Domicile_named_country, Year)

if (anyDuplicated(sdid_long[c("Domicile_named_country", "Year")])) {
  stop("The synthetic DiD panel is not unique by country and year.", call. = FALSE)
}

sdid_setup <- synthdid::panel.matrices(
  as.data.frame(sdid_long),
  unit = "Domicile_named_country",
  time = "Year",
  outcome = "log_app",
  treatment = "treated"
)

tau_hat <- synthdid::synthdid_estimate(
  sdid_setup$Y,
  sdid_setup$N0,
  sdid_setup$T0
)

tau_se <- sqrt(stats::vcov(tau_hat, method = "jackknife"))

eu_countries_sdid <- sort(unique(
  sdid_long$Domicile_named_country[sdid_long$EU == 1L]
))

sdid_by_country <- vapply(eu_countries_sdid, function(country) {
  panel_j <- sdid_long |>
    dplyr::filter(Domicile_named_country == country | EU == 0L) |>
    dplyr::mutate(
      treated_j = as.integer(
        Domicile_named_country == country & Year >= TREATMENT_YEAR
      )
    )

  setup_j <- synthdid::panel.matrices(
    as.data.frame(panel_j),
    unit = "Domicile_named_country",
    time = "Year",
    outcome = "log_app",
    treatment = "treated_j"
  )

  as.numeric(synthdid::synthdid_estimate(
    setup_j$Y,
    setup_j$N0,
    setup_j$T0
  ))
}, numeric(1))

sdid_comparison <- baseline_country_effects |>
  dplyr::select(Country_UCAS, Estimate_TWFE) |>
  dplyr::left_join(
    data.frame(
      Country_UCAS = names(sdid_by_country),
      Estimate_synthetic_DiD = as.numeric(sdid_by_country),
      stringsAsFactors = FALSE
    ),
    by = "Country_UCAS"
  ) |>
  dplyr::mutate(difference = Estimate_TWFE - Estimate_synthetic_DiD) |>
  dplyr::arrange(Estimate_TWFE)

# Optional HonestDiD sensitivity. The rest of the pipeline does not depend on it.
honest_did <- NULL
if (requireNamespace("HonestDiD", quietly = TRUE)) {
  event_beta <- stats::coef(m_event)
  event_vcov <- stats::vcov(m_event)
  event_index <- grep("event_time::-?[0-9]+:EU", names(event_beta))
  event_values <- as.integer(sub(
    "event_time::(-?[0-9]+):EU",
    "\\1",
    names(event_beta)[event_index]
  ))
  event_order <- order(event_values)
  event_index <- event_index[event_order]
  event_values <- event_values[event_order]

  honest_did <- tryCatch(
    list(
      sensitivity = HonestDiD::createSensitivityResults_relativeMagnitudes(
        betahat = event_beta[event_index],
        sigma = event_vcov[event_index, event_index],
        numPrePeriods = sum(event_values < 0),
        numPostPeriods = sum(event_values > 0),
        Mbarvec = seq(0, 2, by = 0.5)
      ),
      original = HonestDiD::constructOriginalCS(
        betahat = event_beta[event_index],
        sigma = event_vcov[event_index, event_index],
        numPrePeriods = sum(event_values < 0),
        numPostPeriods = sum(event_values > 0)
      )
    ),
    error = function(error) {
      message("HonestDiD check failed without stopping the main pipeline\n", conditionMessage(error))
      NULL
    }
  )
} else {
  message("HonestDiD is not installed, so its optional sensitivity check was skipped.")
}

# Present values and loan-corrected elasticities
L <- 9250
Ppost <- 22200
eur_gbp <- 0.85
k_common <- log(Ppost) - log(L)
beta_pooled <- pooled_effect
loan_alpha <- 0.09
loan_tau <- 27295
loan_g <- 0.02
loan_r <- 0.03
loan_T <- 30L

pv_loan <- function(
  Y0,
  L = 9250,
  alpha = 0.09,
  tau = 27295,
  g = 0.02,
  r = 0.03,
  horizon = 30L,
  floor_pv = TRUE
) {
  balance <- L
  present_value <- 0

  for (year in seq_len(horizon)) {
    income <- Y0 * (1 + g)^(year - 1L)
    scheduled_repayment <- alpha * max(0, income - tau)
    payment <- min(scheduled_repayment, balance)
    present_value <- present_value + payment / (1 + r)^year
    balance <- balance - payment
    if (balance <= 0) break
  }

  if (floor_pv) present_value <- max(present_value, 1)
  present_value
}

pv_of <- function(Y0) {
  pv_loan(
    Y0,
    L = L,
    alpha = loan_alpha,
    tau = loan_tau,
    g = loan_g,
    r = loan_r,
    horizon = loan_T
  )
}

robust_se <- function(model) {
  sqrt(diag(sandwich::vcovHC(model, type = "HC1")))
}

robust_tidy <- function(model, conf.int = TRUE) {
  robust_vcov <- sandwich::vcovHC(model, type = "HC1")
  coefficient_test <- lmtest::coeftest(model, vcov. = robust_vcov)
  output <- data.frame(
    term = rownames(coefficient_test),
    estimate = coefficient_test[, "Estimate"],
    std.error = coefficient_test[, "Std. Error"],
    statistic = coefficient_test[, "t value"],
    p.value = coefficient_test[, "Pr(>|t|)"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    confidence_interval <- lmtest::coefci(model, vcov. = robust_vcov)
    output$conf.low <- confidence_interval[, 1L]
    output$conf.high <- confidence_interval[, 2L]
  }

  output
}

country_name_crosswalk <- c(
  "Cyprus (European Union)" = "Cyprus",
  "Czech Republic" = "Czechia"
)

betas_saturated_model <- baseline_country_effects |>
  dplyr::select(Country_UCAS, Estimate_TWFE, Estimate_Pooled) |>
  dplyr::left_join(controlled_country_effects, by = "Country_UCAS") |>
  dplyr::mutate(
    Country = dplyr::recode(Country_UCAS, !!!country_name_crosswalk)
  ) |>
  dplyr::select(Country, Estimate_TWFE, Estimate_controls, Estimate_Pooled)

openxlsx::write.xlsx(
  list(main = betas_saturated_model),
  file.path(OUTPUT_DIR, "betas_saturated_model.xlsx"),
  overwrite = TRUE
)

wages_analysis <- wages |>
  dplyr::mutate(
    wage_gbp = Wage_EURO * eur_gbp,
    PV_loan = vapply(wage_gbp, pv_of, numeric(1)),
    k_i = log(Ppost) - log(PV_loan),
    iso2c = countrycode::countrycode(Country, "country.name", "iso2c")
  ) |>
  dplyr::left_join(betas_saturated_model, by = "Country") |>
  dplyr::mutate(
    dplyr::across(
      c(Estimate_TWFE, Estimate_controls, Estimate_Pooled),
      as.numeric
    ),
    elasticity_new_TWFE = Estimate_TWFE / k_i,
    elasticity_controls_new = Estimate_controls / k_i,
    elasticity_new_Pooled = Estimate_Pooled / k_i,
    elasticity_TWFE = Estimate_TWFE / k_common,
    elasticity_controls = Estimate_controls / k_common,
    elasticity_Pooled = Estimate_Pooled / k_common
  )

unmatched_wage_effects <- wages_analysis |>
  dplyr::filter(is.na(Estimate_TWFE)) |>
  dplyr::pull(Country)

if (length(unmatched_wage_effects) > 0L) {
  message(
    "No treated-country estimate for ",
    paste(unmatched_wage_effects, collapse = ", ")
  )
}

results_country <- wages_analysis |>
  dplyr::select(
    Country, Wage_EURO, PV_loan, k_i,
    elasticity_new_TWFE, elasticity_controls_new, elasticity_new_Pooled,
    elasticity_TWFE, elasticity_controls, elasticity_Pooled
  )

elas_df <- wages_analysis |>
  dplyr::left_join(wb_eu_2022, by = "iso2c") |>
  dplyr::rename(GDP_PPP = gdp_ppp_pc) |>
  dplyr::filter(
    !is.na(GDP_PPP),
    !is.na(Estimate_TWFE),
    !is.na(Estimate_controls),
    !is.na(Estimate_Pooled)
  ) |>
  dplyr::mutate(log_gdp = log(GDP_PPP))

if (nrow(elas_df) < 8L) {
  stop("Fewer than eight countries remain for the elasticity regressions.", call. = FALSE)
}

gradient_models <- list(
  old_TWFE = stats::lm(elasticity_TWFE ~ log_gdp, data = elas_df),
  new_TWFE = stats::lm(elasticity_new_TWFE ~ log_gdp, data = elas_df),
  old_controls = stats::lm(elasticity_controls ~ log_gdp, data = elas_df),
  new_controls = stats::lm(elasticity_controls_new ~ log_gdp, data = elas_df),
  old_pooled = stats::lm(elasticity_Pooled ~ log_gdp, data = elas_df),
  new_pooled = stats::lm(elasticity_new_Pooled ~ log_gdp, data = elas_df)
)

gradient_coefficients <- dplyr::bind_rows(
  lapply(names(gradient_models), function(model_name) {
    robust_tidy(gradient_models[[model_name]]) |>
      dplyr::mutate(model = model_name)
  })
) |>
  dplyr::filter(term == "log_gdp")

gradient_at <- function(r_value, g_value) {
  data <- wages_analysis |>
    dplyr::mutate(
      PV_k = vapply(
        wage_gbp,
        function(wage) {
          pv_loan(
            wage,
            L = L,
            alpha = loan_alpha,
            tau = loan_tau,
            g = g_value,
            r = r_value,
            horizon = loan_T
          )
        },
        numeric(1)
      ),
      k_i_s = log(Ppost) - log(PV_k),
      eta_s = Estimate_TWFE / k_i_s
    ) |>
    dplyr::left_join(
      elas_df |>
        dplyr::select(Country, log_gdp),
      by = "Country"
    ) |>
    dplyr::filter(
      !is.na(eta_s),
      !is.na(log_gdp),
      is.finite(eta_s)
    )

  model <- stats::lm(eta_s ~ log_gdp, data = data)
  coefficient <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC1")
  )["log_gdp", ]

  c(
    slope = unname(coefficient[["Estimate"]]),
    p = unname(coefficient[["Pr(>|t|)"]]),
    n = nrow(data)
  )
}

r_grid <- seq(0.01, 0.06, by = 0.005)
g_grid <- seq(0.00, 0.04, by = 0.005)

sweep_r <- as.data.frame(t(vapply(
  r_grid,
  function(value) gradient_at(value, loan_g),
  numeric(3)
))) |>
  dplyr::mutate(param = "Discount rate r", value = r_grid)

sweep_g <- as.data.frame(t(vapply(
  g_grid,
  function(value) gradient_at(loan_r, value),
  numeric(3)
))) |>
  dplyr::mutate(param = "Wage growth g", value = g_grid)

sweep_all <- dplyr::bind_rows(sweep_r, sweep_g) |>
  dplyr::mutate(significant_5pct = p < 0.05)

placebo_df <- elas_df |>
  dplyr::filter(
    !is.na(k_i),
    !is.na(Estimate_TWFE),
    !is.na(log_gdp)
  ) |>
  dplyr::mutate(
    k_income_raw = log(wage_gbp),
    k_placebo = (k_income_raw - min(k_income_raw)) /
      (max(k_income_raw) - min(k_income_raw)) *
      (max(k_i) - min(k_i)) + min(k_i),
    eta_placebo = Estimate_TWFE / k_placebo
  )

placebo_models <- list(
  equal_k = stats::lm(elasticity_TWFE ~ log_gdp, data = placebo_df),
  real_k = stats::lm(elasticity_new_TWFE ~ log_gdp, data = placebo_df),
  income_only_k = stats::lm(eta_placebo ~ log_gdp, data = placebo_df)
)

placebo_labels <- c(
  equal_k = "Equal k (common price change)",
  real_k = "Real k_i (loan-adjusted)",
  income_only_k = "Placebo k (income only)"
)

placebo_summary <- dplyr::bind_rows(lapply(names(placebo_models), function(name) {
  robust_tidy(placebo_models[[name]]) |>
    dplyr::filter(term == "log_gdp") |>
    dplyr::mutate(denominator = unname(placebo_labels[[name]]))
})) |>
  dplyr::transmute(
    denominator,
    slope = estimate,
    conf.low,
    conf.high,
    p = p.value,
    significant_5pct = p.value < 0.05
  )

correlation_labels <- c(
  gdp_ppp_pc = "GDP per capita (PPP)",
  gini = "Gini",
  youth_unemp = "Youth unemployment",
  educ_spend_gdp = "Education spending (% GDP)",
  tertiary_enroll = "Tertiary enrollment",
  net_migration = "Net migration"
)

correlation_data <- elas_df |>
  dplyr::rename(gdp_ppp_pc = GDP_PPP)

other_correlation_models <- list()
standardized_correlations <- list()

for (variable in names(correlation_labels)) {
  data_subset <- correlation_data |>
    dplyr::select(
      y = elasticity_new_TWFE,
      x = dplyr::all_of(variable)
    ) |>
    dplyr::filter(stats::complete.cases(y, x))

  if (nrow(data_subset) < 8L || stats::sd(data_subset$x) == 0) next

  raw_model <- stats::lm(y ~ x, data = data_subset)
  other_correlation_models[[variable]] <- raw_model

  standardized_data <- data_subset |>
    dplyr::mutate(
      x = as.numeric(scale(x)),
      y = as.numeric(scale(y))
    )
  standardized_model <- stats::lm(y ~ x, data = standardized_data)

  standardized_correlations[[variable]] <- robust_tidy(standardized_model) |>
    dplyr::filter(term == "x") |>
    dplyr::mutate(
      variable = unname(correlation_labels[[variable]]),
      r2 = summary(standardized_model)$r.squared,
      n = nrow(standardized_data)
    )
}

standardized_correlations <- dplyr::bind_rows(standardized_correlations) |>
  dplyr::mutate(significant_5pct = p.value < 0.05)

salary_levels <- c(20000, 26000, 32000, 35000, 42000, 50000)
aggregate_pv_table <- data.frame(Starting_Salary = salary_levels) |>
  dplyr::mutate(
    PV_Loan = vapply(Starting_Salary, pv_of, numeric(1)),
    k = log(Ppost) - log(PV_Loan),
    Elasticity_magnitude = abs(beta_pooled) / k
  )

pooled_model_tables <- list()
for (sample_name in names(pooled_robustness_models)) {
  for (weighting in names(pooled_robustness_models[[sample_name]])) {
    table_name <- paste(sample_name, weighting, sep = "_")
    pooled_model_tables[[table_name]] <- broom::tidy(
      pooled_robustness_models[[sample_name]][[weighting]],
      conf.int = TRUE
    )
  }
}
pooled_model_tables$ireland <- broom::tidy(m_ireland_only, conf.int = TRUE)

control_model_tables <- lapply(
  c(control_models, saturated_control_models),
  broom::tidy,
  conf.int = TRUE
)

output_tables <- c(
  list(
    data_checks = data_checks,
    sdid_comparison = sdid_comparison,
    country_elasticities = results_country,
    gradient_coefficients = gradient_coefficients,
    calibration_sweep = sweep_all,
    placebo_denominator = placebo_summary,
    standardized_correlations = standardized_correlations,
    aggregate_PV = aggregate_pv_table
  ),
  pooled_model_tables,
  control_model_tables
)

names(output_tables) <- make.unique(substr(names(output_tables), 1L, 31L))
openxlsx::write.xlsx(
  output_tables,
  file.path(OUTPUT_DIR, "05_robustness_results.xlsx"),
  overwrite = TRUE
)

robustness_checkpoint <- list(
  placebo_models = list(placebo_2018 = m_placebo_2018, leads_only = m_leads_only),
  leads_test = leads_test,
  data_checks = data_checks,
  control_group_composition = control_group_composition,
  weighted_samples = weighted_samples,
  pooled_robustness_models = pooled_robustness_models,
  m_ireland_only = m_ireland_only,
  saturated_robustness_models = saturated_robustness_models,
  sat_ireland_only = sat_ireland_only,
  control_models = control_models,
  saturated_control_models = saturated_control_models,
  controlled_country_effects = controlled_country_effects,
  interaction_test = interaction_test,
  tau_hat = tau_hat,
  tau_se = tau_se,
  sdid_by_country = sdid_by_country,
  sdid_comparison = sdid_comparison,
  honest_did = honest_did,
  betas_saturated_model = betas_saturated_model,
  wages_analysis = wages_analysis,
  results_country = results_country,
  elas_df = elas_df,
  gradient_models = gradient_models,
  gradient_coefficients = gradient_coefficients,
  sweep_all = sweep_all,
  placebo_summary = placebo_summary,
  standardized_correlations = standardized_correlations,
  aggregate_pv_table = aggregate_pv_table,
  loan_parameters = list(
    L = L,
    Ppost = Ppost,
    eur_gbp = eur_gbp,
    k_common = k_common,
    alpha = loan_alpha,
    tau = loan_tau,
    g = loan_g,
    r = loan_r,
    horizon = loan_T
  )
)

saveRDS(robustness_checkpoint, file.path(CHECKPOINT_DIR, "05_robustness.rds"))

print(data_checks)
print(sdid_comparison)
message("05_robustness.R complete")
