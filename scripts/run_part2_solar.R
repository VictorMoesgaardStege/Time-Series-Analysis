# Run script for Assignment 3, Part 2.
#
# Expected project structure (from README.md):
# .
# ├─ data/
# │  └─ datasolar.csv
# ├─ R/
# │  └─ solar_part2_functions.R
# ├─ scripts/
# │  └─ run_part2_solar.R
# ├─ report/figures/
# └─ output/{tables,models}/
#
# Run from the project root:
#   Rscript scripts/run_part2_solar.R

source("R/solar_part2_functions.R")

results <- run_solar_part2(
  data_path = "data/datasolar.csv",
  phi1 = -0.38,
  Phi1 = -0.94,
  mu = 5.72,
  sigma2_eps = 0.222,
  h = 12,
  figure_dir = "report/figures",
  table_dir = "output/tables",
  model_dir = "output/models"
)

cat("\nResidual diagnostic summary:\n")
print(results$diagnostics)

cat("\n12-month forecast table:\n")
print(results$forecast[, c("horizon", "date", "forecast", "lower", "upper")])

cat("\nFiles written to:\n")
cat("  - report/figures/part2_observed_solar_series.png\n")
cat("  - report/figures/part2_residual_diagnostics.png\n")
cat("  - report/figures/part2_forecast_with_intervals.png\n")
cat("  - output/tables/part2_solar_forecast_table.csv\n")
cat("  - output/tables/part2_residual_diagnostic_tests.csv\n")
cat("  - output/models/part2_solar_results.rds\n")
