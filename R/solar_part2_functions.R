# Assignment 3, Part 2: Predicting monthly solar power
# Helper functions following the project structure in README.md.
# The code uses the fixed model given in the assignment:
# (1 + phi1 B)(1 + Phi1 B^12)(log(Y_t) - mu) = eps_t

# -----------------------------
# Utility helpers
# -----------------------------

ensure_dirs <- function(paths) {
  for (p in paths) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

lag_vec <- function(x, lag = 1L) {
  lag <- as.integer(lag)
  n <- length(x)
  if (lag < 0L) stop("lag must be non-negative")
  if (lag == 0L) return(x)
  c(rep(NA_real_, lag), x[seq_len(n - lag)])
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

# -----------------------------
# Data loading
# -----------------------------

read_solar_data <- function(path = "data/datasolar.csv") {
  if (!file.exists(path)) {
    stop(
      paste0(
        "Could not find '", path, "'. Place datasolar.csv in data/ as described in README.md."
      )
    )
  }

  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names_lower <- tolower(names(dat))

  exact_priority <- c("power", "generation", "gen", "mwh", "yt", "y")
  y_col <- NULL

  for (nm in exact_priority) {
    idx <- which(names_lower == nm)
    if (length(idx) > 0 && is.numeric(dat[[idx[1]]])) {
      y_col <- idx[1]
      break
    }
  }

  if (is.null(y_col)) {
    candidate_idx <- grep("(power|generation|gen|mwh|solar)", names_lower)
    candidate_idx <- candidate_idx[vapply(dat[candidate_idx], is.numeric, logical(1))]
    if (length(candidate_idx) > 0) {
      y_col <- candidate_idx[1]
    }
  }

  if (is.null(y_col)) {
    numeric_cols <- which(vapply(dat, is.numeric, logical(1)))
    exclude <- which(names_lower %in% c("year", "month"))
    numeric_cols <- setdiff(numeric_cols, exclude)

    if (length(numeric_cols) == 0L) {
      stop("No suitable numeric response column found in datasolar.csv.")
    }
    y_col <- numeric_cols[1]
  }

  y <- as.numeric(dat[[y_col]])
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("The response variable must contain strictly positive finite values, since log(Y_t) is used.")
  }

  date <- NULL

  if ("date" %in% names_lower) {
    idx <- which(names_lower == "date")[1]
    try_formats <- c("%Y-%m-%d", "%d-%m-%Y", "%m/%d/%Y", "%Y/%m/%d", "%Y-%m")
    for (fmt in try_formats) {
      tmp <- as.Date(dat[[idx]], format = fmt)
      if (sum(!is.na(tmp)) == length(tmp)) {
        date <- tmp
        break
      }
    }
  }

  if (is.null(date) && all(c("year", "month") %in% names_lower)) {
    year_idx <- which(names_lower == "year")[1]
    month_idx <- which(names_lower == "month")[1]
    date <- as.Date(sprintf("%04d-%02d-01", dat[[year_idx]], dat[[month_idx]]))
  }

  if (is.null(date)) {
    date <- seq.Date(from = as.Date("2008-01-01"), by = "month", length.out = length(y))
  }

  data.frame(
    t = seq_along(y),
    date = date,
    Y = y
  )
}

# -----------------------------
# Model formulas
# -----------------------------

compute_transformed_series <- function(y, mu) {
  log(y) - mu
}

compute_one_step_residuals <- function(X, phi1, Phi1) {
  n <- length(X)
  eps_hat <- rep(NA_real_, n)
  xhat_1step <- rep(NA_real_, n)

  # X_t + phi1 X_{t-1} + Phi1 X_{t-12} + phi1*Phi1 X_{t-13} = eps_t
  for (tt in 14:n) {
    xhat_1step[tt] <- -phi1 * X[tt - 1] - Phi1 * X[tt - 12] - phi1 * Phi1 * X[tt - 13]
    eps_hat[tt] <- X[tt] - xhat_1step[tt]
  }

  data.frame(
    t = seq_len(n),
    xhat_1step = xhat_1step,
    eps_hat = eps_hat
  )
}

forecast_seasonal_ar <- function(X, h, phi1, Phi1) {
  n <- length(X)
  x_all <- c(X, rep(NA_real_, h))
  x_fc <- rep(NA_real_, h)

  for (k in seq_len(h)) {
    idx <- n + k
    x_all[idx] <- -phi1 * x_all[idx - 1] - Phi1 * x_all[idx - 12] - phi1 * Phi1 * x_all[idx - 13]
    x_fc[k] <- x_all[idx]
  }

  x_fc
}

compute_ar1_prediction_intervals <- function(x_forecast, mu, sigma2_eps, phi1, level = 0.95) {
  h <- length(x_forecast)

  # Ignore seasonal part as requested in 2.3
  # (1 + phi1 B)X_t = eps_t  <=>  X_t = a X_{t-1} + eps_t, where a = -phi1
  a <- -phi1
  z <- qnorm(1 - (1 - level) / 2)

  pred_var <- vapply(
    seq_len(h),
    function(k) sigma2_eps * sum(a^(2 * (0:(k - 1)))),
    numeric(1)
  )
  pred_se <- sqrt(pred_var)

  lower_x <- x_forecast - z * pred_se
  upper_x <- x_forecast + z * pred_se

  data.frame(
    horizon = seq_len(h),
    forecast_log = mu + x_forecast,
    lower_log = mu + lower_x,
    upper_log = mu + upper_x,
    pred_var = pred_var,
    pred_se = pred_se,
    forecast = exp(mu + x_forecast),
    lower = exp(mu + lower_x),
    upper = exp(mu + upper_x)
  )
}

# -----------------------------
# Diagnostics and tests
# -----------------------------

residual_test_table <- function(residuals) {
  e <- residuals[is.finite(residuals)]
  n <- length(e)
  if (n < 8L) stop("Too few residuals available for diagnostics.")

  lb_lag_6 <- min(6L, max(1L, floor(n / 3)))
  lb_lag_12 <- min(12L, max(1L, floor(n / 2)))

  data.frame(
    check = c(
      "Residual mean",
      "Residual SD",
      sprintf("Ljung-Box p-value (lag %d)", lb_lag_6),
      sprintf("Ljung-Box p-value (lag %d)", lb_lag_12),
      "Shapiro-Wilk p-value"
    ),
    value = c(
      mean(e),
      sd(e),
      unname(Box.test(e, lag = lb_lag_6, type = "Ljung-Box")$p.value),
      unname(Box.test(e, lag = lb_lag_12, type = "Ljung-Box")$p.value),
      unname(shapiro.test(e)$p.value)
    )
  )
}

# -----------------------------
# Plotting
# -----------------------------

save_observed_plot <- function(solar_df, file) {
  png(file, width = 1200, height = 700, res = 130)
  par(mar = c(4.2, 4.5, 3.0, 1.2), mgp = c(2.5, 0.8, 0))
  plot(
    solar_df$date, solar_df$Y,
    type = "o", pch = 1,
    xlab = "Time",
    ylab = "Generation (MWh)",
    main = "Monthly solar power generation",
    sub = "Observed electricity generation used for seasonal forecasting"
  )
  grid()
  dev.off()
}

save_residual_diagnostics_plot <- function(solar_df, residual_df, file) {
  e <- residual_df$eps_hat
  keep <- is.finite(e)

  png(file, width = 1400, height = 1000, res = 130)
  par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.0, 1.2), mgp = c(2.5, 0.8, 0))

  plot(
    solar_df$date[keep], e[keep],
    type = "o", pch = 16,
    xlab = "Time",
    ylab = expression(hat(epsilon)[t*"|"*t-1]),
    main = "Residuals over time"
  )
  abline(h = 0, lty = 2)
  grid()

  hist(
    e[keep],
    main = "Histogram of residuals",
    xlab = expression(hat(epsilon)[t*"|"*t-1]),
    border = "white"
  )

  qqnorm(e[keep], main = "Normal Q-Q plot of residuals", pch = 16)
  qqline(e[keep], lty = 2)

  acf(e[keep], lag.max = 12, main = "Residual ACF")

  dev.off()
}

save_forecast_only_plot <- function(solar_df, forecast_table, file) {
  last_date <- tail(solar_df$date, 1)
  future_dates <- forecast_table$date
  ylim <- range(c(solar_df$Y, forecast_table$forecast), na.rm = TRUE)

  png(file, width = 1300, height = 760, res = 130)
  par(mar = c(4.2, 4.8, 3.2, 1.2), mgp = c(2.7, 0.8, 0))

  plot(
    solar_df$date, solar_df$Y,
    type = "o", pch = 1,
    ylim = ylim,
    xlim = range(c(solar_df$date, future_dates)),
    xlab = "Time",
    ylab = "Generation (MWh)",
    main = "Observed series and 12-month ahead forecast",
    sub = "Point forecasts from the specified seasonal AR model"
  )
  grid()

  lines(solar_df$date, solar_df$Y, type = "o", pch = 1)
  lines(future_dates, forecast_table$forecast, type = "o", pch = 16, lty = 2, lwd = 2)
  abline(v = last_date, lty = 2)

  legend(
    "topright",
    legend = c("Observed", "Forecast"),
    lty = c(1, 2),
    pch = c(1, 16),
    bty = "n"
  )

  dev.off()
}

save_forecast_interval_plot <- function(solar_df, forecast_table, file) {
  last_date <- tail(solar_df$date, 1)
  future_dates <- forecast_table$date
  ylim <- range(c(solar_df$Y, forecast_table$lower, forecast_table$upper), na.rm = TRUE)

  png(file, width = 1300, height = 760, res = 130)
  par(mar = c(4.2, 4.8, 3.2, 1.2), mgp = c(2.7, 0.8, 0))

  plot(
    solar_df$date, solar_df$Y,
    type = "o", pch = 1,
    ylim = ylim,
    xlim = range(c(solar_df$date, future_dates)),
    xlab = "Time",
    ylab = "Generation (MWh)",
    main = "Observed series, forecasts, and 95% prediction intervals",
    sub = "12-month ahead forecasts from the specified seasonal AR model"
  )
  grid()

  polygon(
    x = c(future_dates, rev(future_dates)),
    y = c(forecast_table$lower, rev(forecast_table$upper)),
    border = NA,
    col = gray(0.85)
  )

  lines(solar_df$date, solar_df$Y, type = "o", pch = 1)
  lines(future_dates, forecast_table$forecast, type = "o", pch = 16, lty = 2, lwd = 2)
  lines(future_dates, forecast_table$lower, lty = 3)
  lines(future_dates, forecast_table$upper, lty = 3)
  abline(v = last_date, lty = 2)

  legend(
    "topright",
    legend = c("Observed", "Forecast", "95% prediction interval"),
    lty = c(1, 2, 3),
    pch = c(1, 16, NA),
    bty = "n"
  )

  dev.off()
}

# -----------------------------
# Main analysis wrapper
# -----------------------------

run_solar_part2 <- function(
  data_path = "data/datasolar.csv",
  phi1 = -0.38,
  Phi1 = -0.94,
  mu = 5.72,
  sigma2_eps = 0.222,
  h = 12,
  figure_dir = "report/figures",
  table_dir = "output/tables",
  model_dir = "output/models"
) {
  ensure_dirs(c(figure_dir, table_dir, model_dir))

  solar_df <- read_solar_data(data_path)
  X <- compute_transformed_series(solar_df$Y, mu = mu)

  residual_df <- compute_one_step_residuals(X, phi1 = phi1, Phi1 = Phi1)
  forecast_x <- forecast_seasonal_ar(X, h = h, phi1 = phi1, Phi1 = Phi1)

  forecast_tab <- compute_ar1_prediction_intervals(
    x_forecast = forecast_x,
    mu = mu,
    sigma2_eps = sigma2_eps,
    phi1 = phi1,
    level = 0.95
  )

  future_dates <- seq.Date(from = tail(solar_df$date, 1), by = "month", length.out = h + 1L)[-1]
  forecast_tab$date <- future_dates
  forecast_tab$forecast <- round(forecast_tab$forecast, 2)
  forecast_tab$lower <- round(forecast_tab$lower, 2)
  forecast_tab$upper <- round(forecast_tab$upper, 2)

  diag_tab <- residual_test_table(residual_df$eps_hat)

  # Separate tables for 2.2 and 2.3
  forecast_only_tab <- forecast_tab[, c("horizon", "date", "forecast")]
  forecast_pi_tab <- forecast_tab[, c("horizon", "date", "forecast", "lower", "upper", "pred_se")]

  # Save outputs
  save_observed_plot(
    solar_df,
    file.path(figure_dir, "part2_observed_solar_series.png")
  )
  save_residual_diagnostics_plot(
    solar_df,
    residual_df,
    file.path(figure_dir, "part2_residual_diagnostics.png")
  )
  save_forecast_only_plot(
    solar_df,
    forecast_tab,
    file.path(figure_dir, "part2_forecast_only.png")
  )
  save_forecast_interval_plot(
    solar_df,
    forecast_tab,
    file.path(figure_dir, "part2_forecast_with_intervals.png")
  )

  write.csv(
    forecast_only_tab,
    file.path(table_dir, "part2_solar_forecast_table.csv"),
    row.names = FALSE
  )

  write.csv(
    forecast_pi_tab,
    file.path(table_dir, "part2_solar_forecast_intervals_table.csv"),
    row.names = FALSE
  )

  write.csv(
    diag_tab,
    file.path(table_dir, "part2_residual_diagnostic_tests.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      solar_data = solar_df,
      X = X,
      residuals = residual_df,
      forecast = forecast_tab,
      forecast_only = forecast_only_tab,
      forecast_intervals = forecast_pi_tab,
      diagnostics = diag_tab,
      parameters = list(phi1 = phi1, Phi1 = Phi1, mu = mu, sigma2_eps = sigma2_eps)
    ),
    file.path(model_dir, "part2_solar_results.rds")
  )

  invisible(
    list(
      solar_data = solar_df,
      X = X,
      residuals = residual_df,
      forecast = forecast_tab,
      forecast_only = forecast_only_tab,
      forecast_intervals = forecast_pi_tab,
      diagnostics = diag_tab
    )
  )
}