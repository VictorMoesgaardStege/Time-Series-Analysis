# Assignment 3 - Question 1: Stability

rm(list = ls())
source("R/stability_helpers.R")

seed <- 20260411
n <- 200
n_rep <- 5
lag_max <- 30
burn_in <- 300

cases <- data.frame(
  question = c("1.1-1.2", "1.3", "1.4", "1.5", "1.6"),
  phi1 = c(-0.6, -0.6, 0.6, -0.7, -0.75),
  phi2 = c( 0.5, -0.3, -0.3, -0.3, -0.3)
)

for (dir_path in c("report/figures", "output/tables", "output/models")) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

all_summaries <- vector("list", nrow(cases))
names(all_summaries) <- cases$question

for (i in seq_len(nrow(cases))) {
  phi1 <- cases$phi1[i]
  phi2 <- cases$phi2[i]
  qlab <- cases$question[i]

  res <- summarize_case(
    phi1 = phi1,
    phi2 = phi2,
    n = n,
    n_rep = n_rep,
    lag_max = lag_max,
    seed = seed + i - 1,
    burn_in = burn_in
  )

  all_summaries[[i]] <- res

  base_name <- paste0(
    "q1_",
    gsub("\\.", "", qlab),
    "_phi1_", gsub("-", "m", as.character(phi1)),
    "_phi2_", gsub("-", "m", as.character(phi2))
  )

  plot_simulations(
    res$simulations,
    phi1 = phi1,
    phi2 = phi2,
    file = file.path("report/figures", paste0(base_name, "_simulations.png"))
  )

  plot_acfs(
    res$simulations,
    phi1 = phi1,
    phi2 = phi2,
    lag_max = lag_max,
    file = file.path("report/figures", paste0(base_name, "_acf.png"))
  )

  root_tbl <- data.frame(
    question = qlab,
    phi1 = phi1,
    phi2 = phi2,
    stationary = res$stationary,
    root = res$roots$root,
    real = res$roots$real,
    imaginary = res$roots$imaginary,
    modulus = res$roots$modulus
  )

  write.csv(
    root_tbl,
    file = file.path("output/tables", paste0(base_name, "_roots.csv")),
    row.names = FALSE
  )
}

root_summary <- do.call(
  rbind,
  lapply(seq_along(all_summaries), function(i) {
    res <- all_summaries[[i]]
    data.frame(
      question = names(all_summaries)[i],
      phi1 = res$phi1,
      phi2 = res$phi2,
      stationary = res$stationary,
      min_root_modulus = min(res$roots$modulus),
      max_root_modulus = max(res$roots$modulus)
    )
  })
)

write.csv(root_summary, "output/tables/q1_stationarity_summary.csv", row.names = FALSE)
saveRDS(all_summaries, "output/models/q1_stability_results.rds")

print(root_summary)

# -----------------------------
# Suggested concise comments for the report
# -----------------------------
cat("\nSuggested interpretation for the report:\n")
cat("1.1-1.2: phi1 = -0.6, phi2 = 0.5 gives complex roots with modulus < 1, so the process is stationary.\n")
cat("         The realizations fluctuate around zero and the ACF shows a damped oscillating pattern.\n")
cat("1.3: phi1 = -0.6, phi2 = -0.3 is stationary because both roots are inside the unit circle.\n")
cat("     The ACF decays without the same oscillatory behavior as in 1.1-1.2.\n")
cat("1.4: phi1 = 0.6, phi2 = -0.3 is also stationary, but the sign change of phi1 changes the correlation pattern.\n")
cat("     Expect more alternating behavior in the ACF and a different short-run dynamic in the paths.\n")
cat("1.5: phi1 = -0.7, phi2 = -0.3 is on the stationarity boundary because one root equals 1.\n")
cat("     The process is not stationary; realizations can look persistent and the sample ACF decays very slowly.\n")
cat("1.6: phi1 = -0.75, phi2 = -0.3 is non-stationary because one root has modulus > 1.\n")
cat("     Realizations typically show stronger drift/persistence and the empirical ACF becomes misleading.\n")
cat("1.7: Use both time-series plots and the ACF. The ACF alone can hide level changes, bursts, and near-nonstationary behavior.\n")
