# 02_cleaning.R
# Construct the synthetic donor weights and all cleaned analysis panels.

if (!exists("PROJECT_DIR", inherits = FALSE)) {
  stop("Run 01_import.R first.", call. = FALSE)
}

solve_weights <- function(y, Z) {
  n <- ncol(Z)
  fit <- Rsolnp::solnp(
    pars = rep(1 / n, n),
    fun = function(w) sum((y - Z %*% w)^2),
    eqfun = function(w) sum(w),
    eqB = 1,
    LB = rep(0, n),
    UB = rep(1, n),
    control = list(trace = 0)
  )

  weights <- as.numeric(fit$pars)
  if (any(!is.finite(weights)) || sum(weights) <= 0) {
    stop("The synthetic-control optimizer returned invalid weights.", call. = FALSE)
  }

  weights / sum(weights)
}

normalize_wdi_columns <- function(
  data,
  indicators,
  required_indicators,
  required_keys,
  source_label
) {
  for (alias in names(indicators)) {
    indicator_code <- unname(indicators[[alias]])
    if (!alias %in% names(data) && indicator_code %in% names(data)) {
      names(data)[names(data) == indicator_code] <- alias
    }
  }

  missing_indicators <- setdiff(names(indicators), names(data))
  missing_required_indicators <- intersect(
    required_indicators,
    missing_indicators
  )
  if (length(missing_required_indicators) > 0L) {
    missing_labels <- paste0(
      missing_required_indicators,
      " (",
      unname(indicators[missing_required_indicators]),
      ")"
    )
    stop(
      source_label,
      " is missing World Bank indicators\n",
      paste(missing_labels, collapse = ", "),
      "\nDelete the corresponding world_bank_*_raw.rds cache and retry the download.",
      call. = FALSE
    )
  }

  optional_missing_indicators <- setdiff(
    missing_indicators,
    missing_required_indicators
  )
  for (alias in optional_missing_indicators) {
    data[[alias]] <- NA_real_
    warning(
      source_label,
      " did not return optional indicator ",
      alias,
      "; the pipeline will retain it as missing.",
      call. = FALSE
    )
  }

  if (
    "iso2c" %in% required_keys &&
    !"iso2c" %in% names(data) &&
    "iso3c" %in% names(data)
  ) {
    data$iso2c <- countrycode::countrycode(
      data$iso3c,
      origin = "iso3c",
      destination = "iso2c",
      warn = FALSE
    )
  }

  missing_keys <- setdiff(required_keys, names(data))
  if (length(missing_keys) > 0L) {
    stop(
      source_label,
      " is missing identifier columns\n",
      paste(missing_keys, collapse = ", "),
      call. = FALSE
    )
  }

  data
}

# Donor trajectories come only from the historical UCAS workbook. The newer
# workbook remains the sole source for the 2016-onward econometric panel.
pre_treatment_years <- 2014L:REFERENCE_YEAR
main_econometric_years <- sort(unique(apps_raw$Year[
  apps_raw$Year >= ANALYSIS_START_YEAR
]))

if (length(main_econometric_years) == 0L) {
  stop(
    "data_filtered.xlsx has no observations from 2016 onward.",
    call. = FALSE
  )
}

missing_pre_treatment_years <- setdiff(
  pre_treatment_years,
  unique(apps_historical_raw$Year)
)

if (length(missing_pre_treatment_years) > 0L) {
  stop(
    "data_filtered_2014_2024.xlsx is missing donor-construction years\n",
    paste(missing_pre_treatment_years, collapse = ", "),
    call. = FALSE
  )
}

data_source_roles <- data.frame(
  analysis_role = c("donor construction", "main econometric panel"),
  source_file = c("data_filtered_2014_2024.xlsx", "data_filtered.xlsx"),
  years_used = c(
    paste(pre_treatment_years, collapse = ", "),
    paste(main_econometric_years, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  data_source_roles,
  file.path(OUTPUT_DIR, "02_data_source_roles.csv"),
  row.names = FALSE
)

control_candidates <- apps_historical_raw |>
  dplyr::filter(
    !Domicile_named_country %in% non_countries,
    Year %in% pre_treatment_years
  )

trajectory_eligibility <- control_candidates |>
  dplyr::group_by(Domicile_named_country) |>
  dplyr::summarise(
    observed_years = dplyr::n_distinct(Year),
    minimum_applicants = min(Applicants),
    complete_trajectory = observed_years == length(pre_treatment_years),
    meets_minimum_flow = minimum_applicants >= 100,
    eligible = complete_trajectory & meets_minimum_flow,
    .groups = "drop"
  )

control_source <- control_candidates |>
  dplyr::semi_join(
    trajectory_eligibility |>
      dplyr::filter(eligible),
    by = "Domicile_named_country"
  )

pre_treatment_window_audit <- data.frame(
  source_file = "data_filtered_2014_2024.xlsx",
  requested_start_year = 2014L,
  reference_year = REFERENCE_YEAR,
  years_used = paste(pre_treatment_years, collapse = ", "),
  eligible_countries = sum(trajectory_eligibility$eligible),
  stringsAsFactors = FALSE
)

utils::write.csv(
  pre_treatment_window_audit,
  file.path(OUTPUT_DIR, "02_pre_treatment_window.csv"),
  row.names = FALSE
)

utils::write.csv(
  trajectory_eligibility,
  file.path(OUTPUT_DIR, "02_trajectory_eligibility.csv"),
  row.names = FALSE
)

message(
  "Using pre-treatment years ",
  paste(pre_treatment_years, collapse = ", "),
  " for donor construction."
)

traj <- control_source |>
  dplyr::select(Domicile_named_country, Year, Applicants) |>
  tidyr::pivot_wider(names_from = Year, values_from = Applicants) |>
  as.data.frame()

rownames(traj) <- traj$Domicile_named_country
traj$Domicile_named_country <- NULL

logmat <- log(as.matrix(traj))
valid_log_rows <- apply(logmat, 1L, function(x) all(is.finite(x)))
logmat <- logmat[valid_log_rows, , drop = FALSE]

row_log_sd <- apply(logmat, 1L, stats::sd)
valid_scale_rows <- is.finite(row_log_sd) & row_log_sd > 0
logmat <- logmat[valid_scale_rows, , drop = FALSE]

if (nrow(logmat) == 0L) {
  stop("No complete positive nonconstant pre-treatment trajectories remain.", call. = FALSE)
}

traj_scaled <- t(scale(t(logmat)))
dist_mat <- as.matrix(stats::dist(traj_scaled, method = "euclidean"))

treated_rows <- rownames(traj_scaled) %in% eu_real
donor_rows <- !treated_rows
dist_treated_donors <- dist_mat[treated_rows, donor_rows, drop = FALSE]
eu_present <- rownames(dist_treated_donors)

if (length(eu_present) == 0L) {
  stop("No treated EU countries remain in the pre-treatment panel.", call. = FALSE)
}

if (ncol(dist_treated_donors) < N_NEAREST_DONORS) {
  stop("Too few eligible donor countries remain.", call. = FALSE)
}

donor_selection <- lapply(eu_present, function(country) {
  distances <- sort(dist_treated_donors[country, ])
  selected <- head(distances, N_NEAREST_DONORS)
  data.frame(
    EU_country = country,
    donor_rank = seq_along(selected),
    donor = names(selected),
    distance = as.numeric(selected),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

weight_results <- lapply(eu_present, function(country) {
  donor_names <- donor_selection$donor[donor_selection$EU_country == country]
  y <- as.numeric(logmat[country, ])
  Z <- t(logmat[donor_names, , drop = FALSE])
  fitted_weights <- solve_weights(y, Z)

  data.frame(
    EU_country = country,
    donor = donor_names,
    synth_weight = fitted_weights,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

fit_results <- lapply(eu_present, function(country) {
  donor_names <- donor_selection$donor[donor_selection$EU_country == country]
  y <- as.numeric(logmat[country, ])
  Z <- t(logmat[donor_names, , drop = FALSE])
  w <- weight_results$synth_weight[weight_results$EU_country == country]
  fitted <- as.numeric(Z %*% w)
  errors <- y - fitted

  data.frame(
    EU_country = country,
    year = as.integer(colnames(logmat)),
    observed = y,
    synthetic = fitted,
    error = errors,
    squared_error = errors^2,
    rmspe = sqrt(mean(errors^2)),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

donor_pool <- donor_selection |>
  dplyr::left_join(weight_results, by = c("EU_country", "donor")) |>
  dplyr::left_join(
    fit_results |>
      dplyr::distinct(EU_country, rmspe),
    by = "EU_country"
  ) |>
  dplyr::arrange(EU_country, donor_rank)

J <- length(eu_present)
weights <- donor_pool |>
  dplyr::group_by(donor) |>
  dplyr::summarise(global_weight = sum(synth_weight) / J, .groups = "drop") |>
  dplyr::mutate(
    normalized_weight = global_weight / sum(global_weight),
    regression_weight = J * normalized_weight
  ) |>
  dplyr::rename(donor_country = donor) |>
  dplyr::arrange(dplyr::desc(normalized_weight))

writexl::write_xlsx(
  list(donor_pool = donor_pool, weights = weights),
  file.path(OUTPUT_DIR, "control_group_output.xlsx")
)

eu_synth <- weights$donor_country
regression_weights <- weights |>
  dplyr::select(
    Domicile_named_country = donor_country,
    regression_weight
  ) |>
  dplyr::bind_rows(
    data.frame(
      Domicile_named_country = eu_real,
      regression_weight = 1,
      stringsAsFactors = FALSE
    )
  )

apps_clean <- apps_raw |>
  dplyr::filter(Year >= ANALYSIS_START_YEAR) |>
  dplyr::mutate(
    EU = as.integer(Domicile_named_country %in% eu_real),
    Post = as.integer(Year >= TREATMENT_YEAR),
    log_app = log(Applicants)
  )

panel_main <- apps_clean |>
  dplyr::filter(Domicile_named_country %in% c(eu_real, eu_synth)) |>
  dplyr::mutate(
    Treated = EU,
    event_time = Year - REFERENCE_YEAR,
    EU_group = dplyr::if_else(EU == 1L, "Real_EU", "Synthetic_EU"),
    treated_factor = factor(dplyr::if_else(
      EU == 1L,
      Domicile_named_country,
      "CONTROL"
    ))
  ) |>
  dplyr::left_join(regression_weights, by = "Domicile_named_country")

if (anyNA(panel_main$regression_weight)) {
  stop("At least one baseline country has no regression weight.", call. = FALSE)
}

missing_treated <- setdiff(eu_real, unique(panel_main$Domicile_named_country))
if (length(missing_treated) > 0L) {
  stop(
    "The analysis panel is missing treated countries\n",
    paste(missing_treated, collapse = ", "),
    call. = FALSE
  )
}

# World Bank controls are cached after the first successful pull.
iso_custom_match <- c(
  "Cyprus (European Union)" = "CY",
  "United States of America" = "US",
  "Korea, Republic of" = "KR",
  "Tanzania, United Republic of" = "TZ",
  "Vietnam [Viet Nam]" = "VN",
  "Russian Federation" = "RU"
)

iso_map <- data.frame(
  Domicile_named_country = unique(panel_main$Domicile_named_country),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    iso2c = countrycode::countrycode(
      Domicile_named_country,
      origin = "country.name",
      destination = "iso2c",
      custom_match = iso_custom_match,
      warn = TRUE
    ),
    iso3c = countrycode::countrycode(iso2c, "iso2c", "iso3c", warn = FALSE)
  )

unmatched_iso <- iso_map |>
  dplyr::filter(is.na(iso2c))

if (nrow(unmatched_iso) > 0L) {
  warning(
    "No World Bank code for ",
    paste(unmatched_iso$Domicile_named_country, collapse = ", "),
    call. = FALSE
  )
}

wb_panel_cache <- file.path(CHECKPOINT_DIR, "world_bank_panel_raw.rds")
wb_codes <- c(
  GDPpc = "NY.GDP.PCAP.PP.KD",
  Pop15_65 = "SP.POP.1564.TO",
  YouthUnempl = "SL.UEM.1524.ZS",
  EduSpendGDP = "SE.XPD.TOTL.GD.ZS"
)

if (file.exists(wb_panel_cache)) {
  wb_raw <- readRDS(wb_panel_cache)
} else {
  wb_raw <- WDI::WDI(
    country = stats::na.omit(iso_map$iso2c),
    indicator = wb_codes,
    start = ANALYSIS_START_YEAR,
    end = 2024L
  )
}

wb_raw <- normalize_wdi_columns(
  wb_raw,
  indicators = wb_codes,
  required_indicators = c("GDPpc", "Pop15_65", "YouthUnempl"),
  required_keys = c("iso2c", "country", "year"),
  source_label = "World Bank panel data"
)
saveRDS(wb_raw, wb_panel_cache)

lebanon_2024 <- wb_raw$country == "Lebanon" & wb_raw$year == 2024L
lebanon_history <- wb_raw$country == "Lebanon" & wb_raw$year %in% c(2022L, 2023L)

if (any(lebanon_2024)) {
  wb_raw$GDPpc[lebanon_2024] <- mean(wb_raw$GDPpc[lebanon_history], na.rm = TRUE)
  wb_raw$YouthUnempl[lebanon_2024] <- mean(
    wb_raw$YouthUnempl[lebanon_history],
    na.rm = TRUE
  )
}

wb_clean <- wb_raw |>
  dplyr::transmute(
    iso2c,
    Year = as.integer(year),
    GDPpc,
    Pop15_65,
    YouthUnempl,
    EduSpendGDP
  )

panel_wb <- panel_main |>
  dplyr::mutate(did = EU * Post) |>
  dplyr::left_join(
    iso_map[c("Domicile_named_country", "iso2c")],
    by = "Domicile_named_country"
  ) |>
  dplyr::left_join(wb_clean, by = c("iso2c", "Year")) |>
  dplyr::filter(
    !is.na(GDPpc),
    !is.na(Pop15_65),
    !is.na(YouthUnempl)
  ) |>
  dplyr::mutate(
    c_log_GDPpc = log(GDPpc) - mean(log(GDPpc), na.rm = TRUE),
    c_log_Pop15_65 = log(Pop15_65) - mean(log(Pop15_65), na.rm = TRUE),
    c_YouthUnempl = YouthUnempl - mean(YouthUnempl, na.rm = TRUE),
    treated_factor = factor(dplyr::if_else(
      EU == 1L,
      Domicile_named_country,
      "CONTROL"
    ))
  )

eu27_names <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia",
  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece",
  "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg",
  "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia",
  "Slovenia", "Spain", "Sweden"
)

wages <- wages_raw |>
  dplyr::filter(Country %in% eu27_names)

missing_wage_countries <- setdiff(eu27_names, wages$Country)
if (length(missing_wage_countries) > 0L) {
  warning(
    "wages.xlsx is missing ",
    paste(missing_wage_countries, collapse = ", "),
    call. = FALSE
  )
}

eu_iso2 <- c(
  "AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE",
  "GR", "HU", "IE", "IT", "LV", "LT", "LU", "MT", "NL", "PL", "PT",
  "RO", "SK", "SI", "ES", "SE"
)

wb_eu_cache <- file.path(CHECKPOINT_DIR, "world_bank_eu_2022_raw.rds")
wb_eu_codes <- c(
  gdp_ppp_pc = "NY.GDP.PCAP.PP.KD",
  gini = "SI.POV.GINI",
  youth_unemp = "SL.UEM.1524.ZS",
  educ_spend_gdp = "SE.XPD.TOTL.GD.ZS",
  tertiary_enroll = "SE.TER.ENRR",
  net_migration = "SM.POP.NETM"
)

if (file.exists(wb_eu_cache)) {
  wb_eu_2022 <- readRDS(wb_eu_cache)
} else {
  wb_eu_2022 <- WDI::WDI(
    country = eu_iso2,
    indicator = wb_eu_codes,
    start = 2022L,
    end = 2022L
  )
}

wb_eu_2022 <- normalize_wdi_columns(
  wb_eu_2022,
  indicators = wb_eu_codes,
  required_indicators = "gdp_ppp_pc",
  required_keys = "iso2c",
  source_label = "World Bank EU 2022 data"
) |>
  dplyr::select(iso2c, dplyr::all_of(names(wb_eu_codes)))
saveRDS(wb_eu_2022, wb_eu_cache)

openxlsx::write.xlsx(panel_wb, file.path(OUTPUT_DIR, "panel_wb.xlsx"), overwrite = TRUE)

cleaning_checkpoint <- list(
  data_source_roles = data_source_roles,
  pre_treatment_years = pre_treatment_years,
  pre_treatment_window_audit = pre_treatment_window_audit,
  trajectory_eligibility = trajectory_eligibility,
  apps_clean = apps_clean,
  panel_main = panel_main,
  panel_wb = panel_wb,
  donor_pool = donor_pool,
  weights = weights,
  eu_synth = eu_synth,
  regression_weights = regression_weights,
  iso_map = iso_map,
  wages = wages,
  wb_eu_2022 = wb_eu_2022
)

saveRDS(cleaning_checkpoint, file.path(CHECKPOINT_DIR, "02_cleaning.rds"))

message("02_cleaning.R complete")
