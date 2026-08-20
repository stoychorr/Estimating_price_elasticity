# Reproducible analysis pipeline

Place `data_filtered_2014_2024.xlsx`, `data_filtered.xlsx`, and `wages.xlsx` in this folder. Both applicant workbooks must contain a `full_data` sheet with `Domicile_named_country`, `Year`, and `Applicants`. The wage workbook must contain a `main_df` sheet with `Country` and `Wage_EURO`.

If either `full_data` sheet contains several disjoint rows for the same country and year, `01_import.R` sums them into the intended country-year `ALL` series before any model runs. The script writes the affected cells to `output/01_duplicate_country_years.csv` and verifies totals separately for each workbook.

Donor construction uses only 2014--2020 from `data_filtered_2014_2024.xlsx`. The main econometric panel uses only `data_filtered.xlsx` from 2016 onward. The selected donor window and country-level eligibility checks are written to `output/02_pre_treatment_window.csv` and `output/02_trajectory_eligibility.csv`.

The separation between the two applicant sources is recorded in `output/02_data_source_roles.csv` so the historical extract cannot silently enter the main regressions.

The scripts stop before reading the data if an R package is missing and print the missing package names. `HonestDiD` may be installed from its source repository when it is unavailable from the configured package repository.

```r
install.packages("remotes")
remotes::install_github("asheshrambachan/HonestDiD")
```

Run the numbered scripts in one fresh R session in this order.

```r
source("01_import.R")
source("02_cleaning.R")
source("03_descriptives.R")
source("04_baseline.R")
source("05_robustness.R")
source("06_figures.R")
rmarkdown::render("paper.Rmd")
```

The same sequence can be run with one command.

```r
source("run_all.R")
```

After scripts 01 through 05 have completed, numerical results can be printed directly in the R console without regenerating any graphs.

```r
source("Outputs_raw.R")
```

The script prints labelled figure inputs, numbered-table inputs, model-scale coefficients, complete model summaries, and statistical tests. It creates no files or graphs. The same material remains available in the `Outputs_raw` list for interactive inspection after the script finishes.

The first successful run downloads and caches the required World Bank series. Later runs reuse the cached raw responses. Delete only the two `world_bank_*_raw.rds` files inside `checkpoints` when a deliberate data refresh is required.

World Bank columns are normalized after both downloads and cache reads. This supports WDI installations that return friendly aliases as well as installations that return raw indicator codes, and it uses `iso2c` consistently for joins.

Generated workbooks are stored in `output`, plots in `figures`, and stage checkpoints in `checkpoints`. The source workbooks are never overwritten.

The uploaded notebooks were mapped in their supplied order. The control-group notebook feeds `01_import.R` and `02_cleaning.R`. The main documentation notebook is distributed across `02_cleaning.R` through `06_figures.R`. The present-value notebook feeds the loan and country-gradient sections of `05_robustness.R` and `06_figures.R`. The manually maintained beta workbook is no longer an input because the pipeline derives it from the fitted models and exports an auditable copy to `output/betas_saturated_model.xlsx`.

`paper.Rmd` contains the supplied manuscript text. Its five numbered tables and thirteen numbered figures are generated from pipeline checkpoints.

| Manuscript item | Generated object or file |
| --- | --- |
| Table 1 | `03_descriptives.rds` |
| Table 2 | `02_cleaning.rds` donor weights |
| Table 3 | `05_robustness.rds` control models |
| Table 4 | `05_robustness.rds` alternative samples |
| Table 5 | `05_robustness.rds` estimator comparison |
| Figure 1 | `figures/event_study.png` |
| Figure 2 | `figures/leads_only.png` |
| Figure 3 | `figures/saturated_event_study.png` |
| Figure 4 | `figures/elasticity_decomposition.png` |
| Figure 5 | `figures/gradient_coefficients.png` |
| Figure 6 | `figures/covid_split.png` |
| Figure 7 | `figures/controlled_elasticities.png` |
| Figure 8 | `figures/pooled_robustness.png` |
| Figure 9 | `figures/saturated_robustness.png` |
| Figure 10 | `figures/calibration_sensitivity.png` |
| Figure 11 | `figures/present_value_curve.png` |
| Figure 12 | `figures/gradient_shift.png` |
| Figure 13 | `figures/standardized_correlations.png` |
