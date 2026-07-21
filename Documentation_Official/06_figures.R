# 06_figures.R
# Build and save every figure from stored model and result objects.

if (!exists("robustness_checkpoint", inherits = FALSE)) {
  stop("Run scripts 01 through 05 first.", call. = FALSE)
}

theme_eer <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "grey15"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = "grey35"),
      plot.caption = ggplot2::element_text(
        size = base_size - 3,
        colour = "grey40",
        hjust = 0,
        margin = ggplot2::margin(t = 8)
      ),
      axis.line = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.4, "lines")
    )
}

accent <- "#1F4E79"

p_event <- ggplot2::ggplot(
  event_coefficients,
  ggplot2::aes(x = event_time, y = estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    width = 0.15,
    colour = accent
  ) +
  ggplot2::geom_line(colour = accent) +
  ggplot2::geom_point(colour = accent, size = 2.2) +
  ggplot2::labs(
    x = "Years relative to 2020",
    y = "Log difference in applications",
    title = "Event study for EU applicants and the synthetic control"
  ) +
  theme_eer()

leads_coefficients <- broom::tidy(m_leads_only, conf.int = TRUE) |>
  dplyr::filter(grepl("^event_time::", term)) |>
  dplyr::mutate(
    event_time = as.integer(sub(
      "^event_time::(-?[0-9]+):EU$",
      "\\1",
      term
    ))
  ) |>
  dplyr::arrange(event_time)

p_leads_only <- ggplot2::ggplot(
  leads_coefficients,
  ggplot2::aes(x = event_time, y = estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    width = 0.15,
    colour = accent
  ) +
  ggplot2::geom_line(colour = accent) +
  ggplot2::geom_point(colour = accent, size = 2.2) +
  ggplot2::labs(
    x = "Years relative to 2020",
    y = "Log difference in applications",
    title = "Leads-only test without anticipation effects"
  ) +
  theme_eer()

p_saturated_event <- ggplot2::ggplot(
  saturated_event_coefficients,
  ggplot2::aes(x = event_time, y = estimate)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    width = 0.15,
    colour = accent
  ) +
  ggplot2::geom_line(colour = accent) +
  ggplot2::geom_point(colour = accent, size = 1.6) +
  ggplot2::facet_wrap(~ Country_UCAS, scales = "free_y") +
  ggplot2::labs(
    x = "Years relative to 2020",
    y = "Estimated log effect",
    title = "Country-saturated event study"
  ) +
  theme_eer(base_size = 10)

fee_pre <- 9250
fee_post <- 22200
log_fee_change <- log(fee_post) - log(fee_pre)
arithmetic_fee_change <- (fee_post - fee_pre) / fee_pre

fee_mechanical <- panel_main |>
  dplyr::filter(Domicile_named_country %in% eu_real) |>
  dplyr::group_by(Domicile_named_country) |>
  dplyr::summarise(
    app_pre = mean(Applicants[Year == 2020L], na.rm = TRUE),
    app_post = mean(Applicants[Year == 2021L], na.rm = TRUE),
    eta_mechanical = ((app_post - app_pre) / app_pre) / arithmetic_fee_change,
    .groups = "drop"
  ) |>
  dplyr::rename(Country_UCAS = Domicile_named_country)

decomposition_data <- baseline_country_effects |>
  dplyr::mutate(
    eta_log = Estimate_TWFE / log_fee_change,
    eta_log_low = conf_low / log_fee_change,
    eta_log_high = conf_high / log_fee_change,
    eta_arc = (exp(Estimate_TWFE) - 1) / arithmetic_fee_change
  ) |>
  dplyr::left_join(fee_mechanical, by = "Country_UCAS") |>
  dplyr::arrange(eta_log) |>
  dplyr::mutate(Country_UCAS = factor(Country_UCAS, levels = Country_UCAS))

p_elasticity_decomposition <- ggplot2::ggplot(
  decomposition_data,
  ggplot2::aes(y = Country_UCAS)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = eta_log_low, xmax = eta_log_high),
    height = 0.15,
    colour = "grey45"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = eta_log, shape = "Structural log elasticity"),
    size = 2.6,
    colour = accent
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = eta_arc, shape = "DiD arc elasticity"),
    size = 2.4,
    colour = "#6BAED6"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = eta_mechanical, shape = "Mechanical fee elasticity"),
    size = 2.4,
    colour = "black"
  ) +
  ggplot2::scale_shape_manual(values = c(16, 17, 1), name = NULL) +
  ggplot2::labs(
    x = "Demand elasticity",
    y = NULL,
    title = "Post-Brexit demand elasticity measures"
  ) +
  theme_eer()

covid_wide <- covid_coefficients |>
  dplyr::transmute(
    Country_UCAS,
    period,
    estimate = estimate / log_fee_change
  ) |>
  tidyr::pivot_wider(names_from = period, values_from = estimate) |>
  dplyr::arrange(dplyr::desc(`Post-COVID`)) |>
  dplyr::mutate(Country_UCAS = factor(Country_UCAS, levels = Country_UCAS))

covid_long <- covid_wide |>
  tidyr::pivot_longer(
    cols = c(`Post-Brexit`, `Post-COVID`),
    names_to = "period",
    values_to = "estimate"
  )

p_covid_split <- ggplot2::ggplot(
  covid_long,
  ggplot2::aes(y = Country_UCAS)
) +
  ggplot2::geom_segment(
    data = covid_wide,
    ggplot2::aes(
      x = `Post-Brexit`,
      xend = `Post-COVID`,
      yend = Country_UCAS
    ),
    colour = "grey75"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = estimate, colour = period),
    size = 2.5
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  ggplot2::scale_colour_manual(
    values = c("Post-Brexit" = "#4C78A8", "Post-COVID" = "#E45756"),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Demand elasticity",
    y = NULL,
    title = "Immediate and persistent post-Brexit responses"
  ) +
  theme_eer()

basic_models <- list()
for (sample_name in names(pooled_robustness_models)) {
  for (weighting in names(pooled_robustness_models[[sample_name]])) {
    label <- paste(tools::toTitleCase(sample_name), tools::toTitleCase(weighting))
    basic_models[[label]] <- pooled_robustness_models[[sample_name]][[weighting]]
  }
}
basic_models[["Ireland only"]] <- m_ireland_only

basic_robustness_data <- dplyr::bind_rows(lapply(names(basic_models), function(name) {
  broom::tidy(basic_models[[name]], conf.int = TRUE) |>
    dplyr::filter(term == "EU:Post") |>
    dplyr::mutate(model = name)
})) |>
  dplyr::mutate(model = factor(model, levels = rev(names(basic_models))))

p_basic_robustness <- ggplot2::ggplot(
  basic_robustness_data,
  ggplot2::aes(x = estimate, y = model)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    height = 0.17,
    colour = accent
  ) +
  ggplot2::geom_point(colour = accent, size = 2.5) +
  ggplot2::labs(
    x = "EU × Post effect on log applicants",
    y = NULL,
    title = "Robustness across donor pools and weighting choices"
  ) +
  theme_eer()

saturated_models_flat <- list()
for (sample_name in names(saturated_robustness_models)) {
  for (weighting in names(saturated_robustness_models[[sample_name]])) {
    label <- paste(tools::toTitleCase(sample_name), substr(weighting, 1L, 1L))
    saturated_models_flat[[label]] <- saturated_robustness_models[[sample_name]][[weighting]]
  }
}
saturated_models_flat[["Ireland"]] <- sat_ireland_only

saturated_robustness_data <- dplyr::bind_rows(lapply(
  names(saturated_models_flat),
  function(model_name) {
    broom::tidy(saturated_models_flat[[model_name]], conf.int = TRUE) |>
      dplyr::filter(grepl("treated_factor::.*:Post$", term)) |>
      dplyr::mutate(
        model = model_name,
        Country_UCAS = term |>
          stringr::str_remove("^treated_factor::") |>
          stringr::str_remove(":Post$")
      )
  }
)) |>
  dplyr::mutate(model = factor(model, levels = names(saturated_models_flat)))

saturated_rank <- saturated_robustness_data |>
  dplyr::group_by(Country_UCAS) |>
  dplyr::summarise(average = mean(estimate), .groups = "drop") |>
  dplyr::arrange(average)

highlight_countries <- c(
  saturated_rank$Country_UCAS[1L],
  saturated_rank$Country_UCAS[nrow(saturated_rank)]
)

saturated_average <- saturated_robustness_data |>
  dplyr::group_by(model) |>
  dplyr::summarise(estimate = mean(estimate), .groups = "drop")

p_saturated_robustness <- ggplot2::ggplot(
  saturated_robustness_data,
  ggplot2::aes(x = model, y = estimate, group = Country_UCAS)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_line(colour = "grey75", alpha = 0.55) +
  ggplot2::geom_line(
    data = saturated_robustness_data |>
      dplyr::filter(Country_UCAS %in% highlight_countries),
    ggplot2::aes(colour = Country_UCAS),
    linewidth = 1
  ) +
  ggplot2::geom_line(
    data = saturated_average,
    ggplot2::aes(x = model, y = estimate, group = 1),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.2
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Effect on log applicants",
    colour = NULL,
    title = "Country-saturated robustness across specifications"
  ) +
  theme_eer(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

elasticity_control_data <- dplyr::bind_rows(lapply(
  names(saturated_control_models),
  function(model_name) {
    extract_saturated(
      saturated_control_models[[model_name]],
      "effect"
    ) |>
      dplyr::mutate(model = model_name)
  }
)) |>
  dplyr::mutate(
    elasticity = effect / log_fee_change,
    elasticity_low = conf_low / log_fee_change,
    elasticity_high = conf_high / log_fee_change,
    model = factor(model, levels = names(saturated_control_models))
  )

p_controlled_elasticities <- ggplot2::ggplot(
  elasticity_control_data,
  ggplot2::aes(x = elasticity, y = Country_UCAS, colour = model)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = elasticity_low, xmax = elasticity_high),
    height = 0.12,
    alpha = 0.7
  ) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::labs(
    x = "Implied demand elasticity",
    y = NULL,
    colour = NULL,
    title = "Country elasticities across control specifications"
  ) +
  theme_eer()

p_sdid_comparison <- sdid_comparison |>
  tidyr::pivot_longer(
    cols = c(Estimate_TWFE, Estimate_synthetic_DiD),
    names_to = "estimator",
    values_to = "estimate"
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = reorder(Country_UCAS, estimate))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65") +
  ggplot2::geom_point(ggplot2::aes(colour = estimator), size = 2.2) +
  ggplot2::scale_colour_manual(
    values = c(
      Estimate_TWFE = accent,
      Estimate_synthetic_DiD = "#D95F02"
    ),
    labels = c("Saturated TWFE", "Synthetic DiD"),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Estimated effect",
    y = NULL,
    title = "Saturated TWFE and synthetic DiD estimates"
  ) +
  theme_eer()

PV_grid <- seq(500, 9250, by = 50)
pv_curve_data <- data.frame(
  PV = PV_grid,
  elasticity = abs(beta_pooled) / (log(Ppost) - log(PV_grid))
)

p_pv_curve <- ggplot2::ggplot(
  pv_curve_data,
  ggplot2::aes(x = PV, y = elasticity)
) +
  ggplot2::geom_line(colour = accent, linewidth = 1) +
  ggplot2::geom_vline(xintercept = L, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(
    data = aggregate_pv_table |>
      dplyr::filter(Starting_Salary %in% c(20000, 26000, 32000, 42000)),
    ggplot2::aes(x = PV_Loan, y = Elasticity_magnitude),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.8
  ) +
  ggplot2::geom_text(
    data = aggregate_pv_table |>
      dplyr::filter(Starting_Salary %in% c(20000, 26000, 32000, 42000)),
    ggplot2::aes(
      x = PV_Loan,
      y = Elasticity_magnitude,
      label = paste0("£", Starting_Salary)
    ),
    inherit.aes = FALSE,
    nudge_y = 0.06,
    size = 3.2
  ) +
  ggplot2::labs(
    x = "Present value of the UK student loan in GBP",
    y = "Implied elasticity magnitude",
    title = "Loan value and the implied demand elasticity"
  ) +
  theme_eer()

p_calibration <- ggplot2::ggplot(
  sweep_all,
  ggplot2::aes(x = value, y = p)
) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#8B0000") +
  ggplot2::geom_line(colour = accent) +
  ggplot2::geom_point(colour = accent, size = 2) +
  ggplot2::facet_wrap(~ param, scales = "free_x") +
  ggplot2::labs(
    x = "Parameter value",
    y = "p-value of the income gradient",
    title = "Calibration sensitivity of the corrected income gradient"
  ) +
  theme_eer()

p_placebo_denominator <- placebo_summary |>
  dplyr::mutate(
    denominator = factor(
      denominator,
      levels = c(
        "Equal k (common price change)",
        "Real k_i (loan-adjusted)",
        "Placebo k (income only)"
      )
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = slope, y = denominator)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = conf.low, xmax = conf.high, colour = significant_5pct),
    height = 0.15
  ) +
  ggplot2::geom_point(ggplot2::aes(colour = significant_5pct), size = 2.4) +
  ggplot2::scale_colour_manual(
    values = c(`TRUE` = "#8B0000", `FALSE` = "grey55"),
    guide = "none"
  ) +
  ggplot2::labs(
    x = "Income-gradient slope",
    y = NULL,
    title = "Real and placebo elasticity denominators"
  ) +
  theme_eer()

p_standardized_correlations <- standardized_correlations |>
  dplyr::mutate(variable = reorder(variable, estimate)) |>
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = variable)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = conf.low,
      xmax = conf.high,
      colour = significant_5pct
    ),
    height = 0.15
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = significant_5pct),
    size = 2.3
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("R² = %.2f", r2)),
    hjust = -0.2,
    vjust = -0.8,
    size = 3,
    colour = "grey35"
  ) +
  ggplot2::scale_colour_manual(
    values = c(`TRUE` = "#8B0000", `FALSE` = "grey60"),
    guide = "none"
  ) +
  ggplot2::labs(
    x = "Standardized association with corrected elasticity",
    y = NULL,
    title = "Country correlates of the loan-corrected elasticity"
  ) +
  theme_eer()

shift_data <- elas_df |>
  dplyr::select(
    Country,
    log_gdp,
    equal_k = elasticity_TWFE,
    corrected = elasticity_new_TWFE
  ) |>
  dplyr::arrange(log_gdp) |>
  dplyr::mutate(position = dplyr::row_number())

p_gradient_shift <- ggplot2::ggplot(
  shift_data,
  ggplot2::aes(x = position)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      xend = position,
      y = equal_k,
      yend = corrected,
      colour = log_gdp
    ),
    arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed")
  ) +
  ggplot2::geom_point(ggplot2::aes(y = equal_k), colour = "grey55") +
  ggplot2::geom_point(
    ggplot2::aes(y = corrected, colour = log_gdp),
    size = 2.3
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(y = equal_k),
    method = "lm",
    se = FALSE,
    colour = "grey55",
    linetype = "dashed"
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(y = corrected),
    method = "lm",
    se = FALSE,
    colour = accent
  ) +
  ggplot2::scale_x_continuous(
    breaks = shift_data$position,
    labels = shift_data$Country
  ) +
  ggplot2::scale_colour_gradient(low = "#8B0000", high = "#4B0082") +
  ggplot2::labs(
    x = NULL,
    y = "Demand elasticity",
    colour = "Log GDP per capita",
    title = "Country corrections and the income gradient"
  ) +
  theme_eer() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))

gradient_plot_data <- gradient_coefficients |>
  dplyr::mutate(model = factor(model, levels = rev(model)))

p_gradient_coefficients <- ggplot2::ggplot(
  gradient_plot_data,
  ggplot2::aes(x = estimate, y = model)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    height = 0.15,
    colour = accent
  ) +
  ggplot2::geom_point(colour = accent, size = 2.3) +
  ggplot2::labs(
    x = "Slope on log GDP per capita",
    y = NULL,
    title = "Income-gradient slope across elasticity specifications"
  ) +
  theme_eer()

k_data <- elas_df |>
  dplyr::select(Country, k_i) |>
  dplyr::arrange(k_i) |>
  dplyr::mutate(Country = factor(Country, levels = Country))

p_k_compression <- ggplot2::ggplot(
  k_data,
  ggplot2::aes(x = k_i, y = Country)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = k_common, xend = k_i, yend = Country),
    colour = "grey75"
  ) +
  ggplot2::geom_vline(xintercept = k_common, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(colour = accent, size = 2.3) +
  ggplot2::labs(
    x = "Loan-adjusted log price change",
    y = NULL,
    title = "Unequal compression of the effective price increase"
  ) +
  theme_eer()

figures <- list(
  event_study = p_event,
  leads_only = p_leads_only,
  saturated_event_study = p_saturated_event,
  elasticity_decomposition = p_elasticity_decomposition,
  covid_split = p_covid_split,
  pooled_robustness = p_basic_robustness,
  saturated_robustness = p_saturated_robustness,
  controlled_elasticities = p_controlled_elasticities,
  synthetic_did_comparison = p_sdid_comparison,
  present_value_curve = p_pv_curve,
  calibration_sensitivity = p_calibration,
  placebo_denominator = p_placebo_denominator,
  standardized_correlations = p_standardized_correlations,
  gradient_shift = p_gradient_shift,
  gradient_coefficients = p_gradient_coefficients,
  price_compression = p_k_compression
)

figure_dimensions <- data.frame(
  name = names(figures),
  width = c(7, 7, 11, 8, 8, 8, 10, 9, 8, 7, 8, 8, 8, 9, 8, 7),
  height = c(4.5, 4.5, 9, 7, 7, 5, 6, 7, 7, 4.5, 4.5, 4, 5, 5, 4.5, 6),
  stringsAsFactors = FALSE
)

figure_manifest <- lapply(seq_len(nrow(figure_dimensions)), function(index) {
  figure_name <- figure_dimensions$name[index]
  figure_file <- file.path(FIGURE_DIR, paste0(figure_name, ".png"))

  ggplot2::ggsave(
    filename = figure_file,
    plot = figures[[figure_name]],
    width = figure_dimensions$width[index],
    height = figure_dimensions$height[index],
    units = "in",
    dpi = 320,
    bg = "white"
  )

  data.frame(
    figure = figure_name,
    file = figure_file,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

saveRDS(
  list(figures = figures, manifest = figure_manifest),
  file.path(CHECKPOINT_DIR, "06_figures.rds")
)

utils::write.csv(
  figure_manifest,
  file.path(OUTPUT_DIR, "06_figure_manifest.csv"),
  row.names = FALSE
)

utils::capture.output(
  utils::sessionInfo(),
  file = file.path(OUTPUT_DIR, "session_info.txt")
)

print(figure_manifest)
message("06_figures.R complete")
