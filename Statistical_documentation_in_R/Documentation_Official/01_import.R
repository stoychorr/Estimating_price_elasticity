# 01_import.R
# Import the two local source workbooks and define shared pipeline settings.

script_directory <- function() {
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

PROJECT_DIR    <- script_directory()
CHECKPOINT_DIR <- file.path(PROJECT_DIR, "checkpoints")
OUTPUT_DIR     <- file.path(PROJECT_DIR, "output")
FIGURE_DIR     <- file.path(PROJECT_DIR, "figures")

set.seed(20260721)

invisible(lapply(
  c(CHECKPOINT_DIR, OUTPUT_DIR, FIGURE_DIR),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_packages <- c(
  "broom", "countrycode", "dplyr", "fixest", "ggplot2", "ggrepel",
  "HonestDiD", "knitr", "lmtest", "openxlsx", "psych", "readxl", "Rsolnp",
  "sandwich", "stringr", "synthdid", "tidyr", "WDI", "writexl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the missing packages before running the pipeline\n",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

input_files <- c(
  applicants_historical = file.path(PROJECT_DIR, "data_filtered_2014_2024.xlsx"),
  applicants_main       = file.path(PROJECT_DIR, "data_filtered.xlsx"),
  wages                 = file.path(PROJECT_DIR, "wages.xlsx")
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Place the following input files beside 01_import.R\n",
    paste(basename(missing_inputs), collapse = ", "),
    call. = FALSE
  )
}

input_checksums <- data.frame(
  file = basename(input_files),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  input_checksums,
  file.path(OUTPUT_DIR, "01_input_checksums.csv"),
  row.names = FALSE
)

apps_historical_input <- readxl::read_excel(
  input_files[["applicants_historical"]],
  sheet = "full_data"
)
apps_main_input <- readxl::read_excel(
  input_files[["applicants_main"]],
  sheet = "full_data"
)
wages_raw <- readxl::read_excel(input_files[["wages"]], sheet = "main_df")

required_app_columns <- c("Domicile_named_country", "Year", "Applicants")
required_wage_columns <- c("Country", "Wage_EURO")
missing_wage_columns <- setdiff(required_wage_columns, names(wages_raw))
if (length(missing_wage_columns) > 0L) {
  stop(
    "The main_df sheet in wages.xlsx is missing columns\n",
    paste(missing_wage_columns, collapse = ", "),
    call. = FALSE
  )
}

prepare_applicant_source <- function(data, source_name) {
  missing_columns <- setdiff(required_app_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      source_name,
      " full_data is missing columns\n",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data <- data |>
    dplyr::mutate(
      Domicile_named_country = as.character(Domicile_named_country),
      Year = as.integer(Year),
      Applicants = as.numeric(Applicants)
    )

  if (anyNA(data[c("Domicile_named_country", "Year", "Applicants")])) {
    stop(
      source_name,
      " contains missing country, year, or applicant values.",
      call. = FALSE
    )
  }

  if (any(data$Applicants <= 0)) {
    stop(
      source_name,
      " contains nonpositive applicant values, which cannot be logged.",
      call. = FALSE
    )
  }

  country_year_totals <- data |>
    dplyr::group_by(Domicile_named_country, Year) |>
    dplyr::summarise(
      source_rows = dplyr::n(),
      Applicants = sum(Applicants),
      .groups = "drop"
    )

  aggregation_summary <- data.frame(
    source = source_name,
    source_rows = nrow(data),
    collapsed_country_years = nrow(country_year_totals),
    country_years_with_multiple_rows = sum(country_year_totals$source_rows > 1L),
    applicants_before = sum(data$Applicants),
    applicants_after = sum(country_year_totals$Applicants),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(all.equal(
    aggregation_summary$applicants_before,
    aggregation_summary$applicants_after,
    tolerance = 1e-10
  ))) {
    stop(
      source_name,
      " country-year aggregation changed the total applicant count.",
      call. = FALSE
    )
  }

  list(
    data = country_year_totals |>
      dplyr::select(Domicile_named_country, Year, Applicants),
    summary = aggregation_summary,
    duplicates = country_year_totals |>
      dplyr::filter(source_rows > 1L) |>
      dplyr::mutate(source = source_name)
  )
}

historical_import <- prepare_applicant_source(
  apps_historical_input,
  "data_filtered_2014_2024.xlsx"
)
main_import <- prepare_applicant_source(
  apps_main_input,
  "data_filtered.xlsx"
)

apps_historical_raw <- historical_import$data
apps_raw <- main_import$data
aggregation_summary <- dplyr::bind_rows(
  historical_import$summary,
  main_import$summary
)

wages_raw <- wages_raw |>
  dplyr::mutate(
    Country = as.character(Country),
    Wage_EURO = as.numeric(Wage_EURO)
  )

utils::write.csv(
  dplyr::bind_rows(
    historical_import$duplicates,
    main_import$duplicates
  ),
  file.path(OUTPUT_DIR, "01_duplicate_country_years.csv"),
  row.names = FALSE
)

utils::write.csv(
  aggregation_summary,
  file.path(OUTPUT_DIR, "01_country_year_aggregation_summary.csv"),
  row.names = FALSE
)

message(
  "Imported historical applicant years ",
  min(apps_historical_raw$Year),
  "--",
  max(apps_historical_raw$Year),
  " and main applicant years ",
  min(apps_raw$Year),
  "--",
  max(apps_raw$Year),
  "."
)

eu_real <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus (European Union)",
  "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany",
  "Greece", "Hungary", "Italy", "Latvia", "Lithuania", "Luxembourg",
  "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia",
  "Slovenia", "Spain", "Sweden"
)

non_countries <- c(
  "England", "Scotland", "Wales", "Northern Ireland", "Alderney",
  "Guernsey", "Jersey", "Isle of Man", "Anguilla", "Bermuda",
  "British Antarctic Territory", "British Indian Ocean Territory",
  "British Virgin Islands", "Cayman Islands", "Falkland Islands (Malvinas)",
  "Gibraltar", "Montserrat", "St Helena, Ascension & Tristan da Cunha",
  "Turks and Caicos Islands", "Aland Islands", "Aruba", "Canary Islands",
  "Christmas Island", "Curacao", "Faroe Islands", "French Guiana",
  "French Polynesia", "Greenland", "Guadeloupe", "Guam", "Hong Kong",
  "Macao", "Martinique", "Netherlands Antilles", "New Caledonia",
  "Norfolk Island", "Northern Mariana Islands", "Puerto Rico", "Reunion",
  "Sint Maarten (Dutch Part)", "St Barthelemy", "St Martin (French part)",
  "Virgin Islands (US)", "Cyprus (Non-European Union)",
  "Cyprus (Not otherwise specified)", "All",
  "No information provided (Overseas)", "Not Known", "Stateless",
  "USSR (Not In Use)", "French West Indies (Not In Use)"
)

ANALYSIS_START_YEAR <- 2016L
TREATMENT_YEAR      <- 2021L
REFERENCE_YEAR      <- 2020L
N_NEAREST_DONORS    <- 3L

import_checkpoint <- list(
  apps_historical_raw = apps_historical_raw,
  apps_raw = apps_raw,
  aggregation_summary = aggregation_summary,
  wages_raw = wages_raw,
  input_checksums = input_checksums,
  eu_real = eu_real,
  non_countries = non_countries,
  settings = list(
    seed = 20260721L,
    analysis_start_year = ANALYSIS_START_YEAR,
    treatment_year = TREATMENT_YEAR,
    reference_year = REFERENCE_YEAR,
    n_nearest_donors = N_NEAREST_DONORS
  )
)

saveRDS(import_checkpoint, file.path(CHECKPOINT_DIR, "01_import.rds"))

message("01_import.R complete")
