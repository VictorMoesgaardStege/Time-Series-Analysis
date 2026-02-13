Project structure for Assignment 1

.
├─ data/
│  └─ DST_BIL54.csv
├─ R/
│  ├─ read_data.R          # functions: read + split
│  ├─ helpers.R            # small helpers (time variable, plotting utilities)
│  ├─ ols_trend.R          # OLS fit + SEs + forecasts + residual checks
│  ├─ wls_trend.R          # WLS (lambda weights) + forecasts
│  └─ rls.R                # recursive least squares + forgetting + RMSE grid
├─ scripts/
│  ├─ 01_plot_data.R
│  ├─ 02_matrix_form.R
│  ├─ 03_ols.R
│  ├─ 04_wls.R
│  └─ 05_rls.R
├─ report/
│  ├─ assignment1.Rmd      # generates the single PDF you hand in
│  └─ figures/             # saved plots (optional)
├─ output/
│  ├─ tables/              # forecast tables (csv) (optional)
│  └─ models/              # saved RDS objects (optional)
├─ assignment1.pdf
├─ README.md
└─ renv/ + renv.lock       # (optional but recommended) reproducible packages
