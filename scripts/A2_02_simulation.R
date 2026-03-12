# ============================================================
# Module: simulate_and_plot_seasonal_models.R
# Purpose:
#   Simulate the six seasonal models from the assignment and
#   plot each model as:
#     1) time series
#     2) ACF
#     3) PACF
#   in one large matrix plot.
#
# Notes:
#   arima.sim() has no direct seasonal argument, so the seasonal
#   AR and MA terms are inserted manually at lags 12, 13, etc.
# ============================================================

simulate_and_plot_seasonal_models <- function(
    n = 300,
    burnin = 300,
    seasonal_period = 12,
    seed = 123,
    max_lag = 50
) {
  set.seed(seed)

  s <- seasonal_period

  # ----------------------------------------------------------
  # Helper: simulate one model
  # ----------------------------------------------------------
  simulate_model <- function(name, ar = NULL, ma = NULL, n = 300, burnin = 300) {
    x <- arima.sim(
      n = n,
      model = list(ar = ar, ma = ma),
      n.start = burnin
    )
    list(name = name, x = as.numeric(x), ar = ar, ma = ma)
  }

  # ----------------------------------------------------------
  # 2.1  (1,0,0) x (0,0,0)_12
  # phi1 = 0.6
  #
  # (1 - 0.6 B) Y_t = e_t
  # ----------------------------------------------------------
  m1 <- simulate_model(
    name = "2.1  (1,0,0)x(0,0,0)[12]\nphi1 = 0.6",
    ar = c(0.6),
    n = n,
    burnin = burnin
  )

  # ----------------------------------------------------------
  # 2.2  (0,0,0) x (1,0,0)_12
  # Phi1 = -0.9
  #
  # Assignment notation:
  #   Phi(B^12) = 1 + Phi1 B^12 = 1 - 0.9 B^12
  # Recursive form gives lag-12 AR coefficient = 0.9
  # ----------------------------------------------------------
  ar_12 <- rep(0, s)
  ar_12[s] <- 0.9

  m2 <- simulate_model(
    name = "2.2  (0,0,0)x(1,0,0)[12]\nPhi1 = -0.9",
    ar = ar_12,
    n = n,
    burnin = burnin
  )

  # ----------------------------------------------------------
  # 2.3  (1,0,0) x (0,0,1)_12
  # phi1 = 0.9, Theta1 = -0.7
  #
  # (1 - 0.9B)Y_t = (1 - 0.7B^12)e_t
  # AR lag 1 = 0.9
  # MA lag 12 = -0.7
  # ----------------------------------------------------------
  ma_12_a <- rep(0, s)
  ma_12_a[s] <- -0.7

  m3 <- simulate_model(
    name = "2.3  (1,0,0)x(0,0,1)[12]\nphi1 = 0.9, Theta1 = -0.7",
    ar = c(0.9),
    ma = ma_12_a,
    n = n,
    burnin = burnin
  )

  # ----------------------------------------------------------
  # 2.4  (1,0,0) x (1,0,0)_12
  # phi1 = -0.6, Phi1 = -0.8
  #
  # (1 - 0.6B)(1 - 0.8B^12)
  # = 1 - 0.6B - 0.8B^12 + 0.48B^13
  #
  # Recursive AR coefficients:
  # lag 1  = 0.6
  # lag 12 = 0.8
  # lag 13 = -0.48
  # ----------------------------------------------------------
  ar_13 <- rep(0, s + 1)
  ar_13[1]  <- 0.6
  ar_13[12] <- 0.8
  ar_13[13] <- -0.48

  m4 <- simulate_model(
    name = "2.4  (1,0,0)x(1,0,0)[12]\nphi1 = -0.6, Phi1 = -0.8",
    ar = ar_13,
    n = n,
    burnin = burnin
  )

  # ----------------------------------------------------------
  # 2.5  (0,0,1) x (0,0,1)_12
  # theta1 = 0.4, Theta1 = -0.8
  #
  # (1 + 0.4B)(1 - 0.8B^12)
  # = 1 + 0.4B - 0.8B^12 - 0.32B^13
  #
  # MA coefficients:
  # lag 1  = 0.4
  # lag 12 = -0.8
  # lag 13 = -0.32
  # ----------------------------------------------------------
  ma_13 <- rep(0, s + 1)
  ma_13[1]  <- 0.4
  ma_13[12] <- -0.8
  ma_13[13] <- -0.32

  m5 <- simulate_model(
    name = "2.5  (0,0,1)x(0,0,1)[12]\ntheta1 = 0.4, Theta1 = -0.8",
    ma = ma_13,
    n = n,
    burnin = burnin
  )

  # ----------------------------------------------------------
  # 2.6  (0,0,1) x (1,0,0)_12
  # theta1 = -0.4, Phi1 = 0.7
  #
  # MA part: 1 - 0.4B
  # Seasonal AR part: 1 + 0.7B^12
  #
  # Recursive AR lag 12 = -0.7
  # MA lag 1 = -0.4
  # ----------------------------------------------------------
  ar_12_b <- rep(0, s)
  ar_12_b[s] <- -0.7

  m6 <- simulate_model(
    name = "2.6  (0,0,1)x(1,0,0)[12]\ntheta1 = -0.4, Phi1 = 0.7",
    ar = ar_12_b,
    ma = c(-0.4),
    n = n,
    burnin = burnin
  )

  models <- list(m1, m2, m3, m4, m5, m6)

  # ----------------------------------------------------------
  # Plot all models in one large matrix:
  # each row = one model
  # columns = series, ACF, PACF
  # ----------------------------------------------------------
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(
    mfrow = c(length(models), 3),
    mar = c(3, 3, 3, 1),
    oma = c(2, 2, 4, 1)
  )

  for (m in models) {
    # Time series
    plot(
      m$x,
      type = "l",
      main = paste(m$name, "\nTime series"),
      xlab = "Time",
      ylab = "Y_t"
    )

    # ACF
    acf(
      m$x,
      lag.max = max_lag,
      main = paste(m$name, "\nACF")
    )

    # PACF
    pacf(
      m$x,
      lag.max = max_lag,
      main = paste(m$name, "\nPACF")
    )
  }

  mtext(
    "Simulated seasonal models: time series, ACF and PACF",
    outer = TRUE,
    cex = 1.4,
    font = 2
  )

  invisible(models)
}

# ============================================================
# Example usage
# ============================================================
models <- simulate_and_plot_seasonal_models()