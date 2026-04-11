# Helper functions for Assignment 3, Question 1 (Stability)
# Course notation:
#   X_t + phi1 X_{t-1} + phi2 X_{t-2} = eps_t
#
# In recursive form:
#   X_t = -phi1 X_{t-1} - phi2 X_{t-2} + eps_t
#
# This manual simulator works for both stationary and non-stationary cases.

simulate_ar2_single <- function(phi1, phi2, n = 200, burn_in = 300,
                                sigma_eps = 1, x_init = c(0, 0)) {
  total_n <- n + burn_in
  x <- numeric(total_n)
  eps <- rnorm(total_n, mean = 0, sd = sigma_eps)

  # initial values
  x[1] <- x_init[1]
  x[2] <- x_init[2]

  # recursion:
  # X_t + phi1 X_{t-1} + phi2 X_{t-2} = eps_t
  # => X_t = -phi1 X_{t-1} - phi2 X_{t-2} + eps_t
  for (t in 3:total_n) {
    x[t] <- -phi1 * x[t - 1] - phi2 * x[t - 2] + eps[t]
  }

  x[(burn_in + 1):total_n]
}

simulate_ar2_reps <- function(phi1, phi2, n = 200, n_rep = 5, seed = 20260411,
                              burn_in = 300, sigma_eps = 1) {
  set.seed(seed)
  sims <- replicate(
    n_rep,
    simulate_ar2_single(
      phi1 = phi1,
      phi2 = phi2,
      n = n,
      burn_in = burn_in,
      sigma_eps = sigma_eps
    )
  )
  sims <- as.matrix(sims)
  colnames(sims) <- paste0("sim_", seq_len(n_rep))
  sims
}

get_stationarity_roots <- function(phi1, phi2) {
  # AR polynomial in z:
  # z^2 + phi1 z + phi2 = 0
  roots <- polyroot(c(phi2, phi1, 1))
  data.frame(
    root = c("z1", "z2"),
    real = Re(roots),
    imaginary = Im(roots),
    modulus = Mod(roots)
  )
}

is_stationary_ar2 <- function(phi1, phi2, tol = 1e-8) {
  roots <- polyroot(c(phi2, phi1, 1))
  all(Mod(roots) < 1 - tol)
}

get_theoretical_acf <- function(phi1, phi2, lag_max = 30) {
  if (!is_stationary_ar2(phi1, phi2)) {
    return(rep(NA_real_, lag_max + 1))
  }
  as.numeric(ARMAacf(ar = -c(phi1, phi2), lag.max = lag_max))
}

get_empirical_acfs <- function(sims, lag_max = 30, plot = FALSE) {
  apply(sims, 2, function(x) {
    acf(x, lag.max = lag_max, plot = plot, na.action = na.pass)$acf[, 1, 1]
  })
}

plot_simulations <- function(sims, phi1, phi2, file = NULL) {
  if (!is.null(file)) {
    png(file, width = 1200, height = 700, res = 120)
    on.exit(dev.off(), add = TRUE)
  }

  matplot(
    sims,
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "t",
    ylab = expression(X[t]),
    main = bquote(
      "Five simulated realizations of AR(2): " ~
      X[t] + .(phi1) * X[t-1] + .(phi2) * X[t-2] == epsilon[t]
    )
  )
  grid()
  legend(
    "topright",
    legend = colnames(sims),
    col = seq_len(ncol(sims)),
    lty = 1,
    lwd = 2,
    bty = "n"
  )
}

plot_acfs <- function(sims, phi1, phi2, lag_max = 30, file = NULL) {
  emp_acfs <- get_empirical_acfs(sims, lag_max = lag_max)
  theo_acf <- get_theoretical_acf(phi1, phi2, lag_max = lag_max)
  lags <- 0:lag_max
  stationary <- is_stationary_ar2(phi1, phi2)

  if (!is.null(file)) {
    png(file, width = 1200, height = 700, res = 120)
    on.exit(dev.off(), add = TRUE)
  }

  ylim <- range(emp_acfs, na.rm = TRUE)
  if (stationary) {
    ylim <- range(c(emp_acfs, theo_acf), na.rm = TRUE)
  }

  plot(
    lags, emp_acfs[, 1],
    type = "l",
    lwd = 1.5,
    col = 1,
    xlab = "Lag k",
    ylab = "ACF",
    ylim = ylim,
    main = if (stationary) {
      bquote(
        "Empirical ACFs and theoretical " * rho(k) *
        " for AR(2): " ~ X[t] + .(phi1) * X[t-1] + .(phi2) * X[t-2] == epsilon[t]
      )
    } else {
      bquote(
        "Empirical ACFs for non-stationary AR(2): " ~
        X[t] + .(phi1) * X[t-1] + .(phi2) * X[t-2] == epsilon[t]
      )
    }
  )
  grid()

  if (ncol(emp_acfs) > 1) {
    for (j in 2:ncol(emp_acfs)) {
      lines(lags, emp_acfs[, j], lwd = 1.5, col = j)
    }
  }

  if (stationary) {
    lines(lags, theo_acf, lwd = 3, col = "black")
    legend(
      "topright",
      legend = c(colnames(sims), "theoretical rho(k)"),
      col = c(seq_len(ncol(sims)), "black"),
      lty = 1,
      lwd = c(rep(1.5, ncol(sims)), 3),
      bty = "n"
    )
  } else {
    legend(
      "topright",
      legend = c(colnames(sims), "no theoretical ACF (non-stationary)"),
      col = c(seq_len(ncol(sims)), NA),
      lty = c(rep(1, ncol(sims)), NA),
      lwd = c(rep(1.5, ncol(sims)), NA),
      bty = "n"
    )
  }
}

summarize_case <- function(phi1, phi2, n = 200, n_rep = 5, lag_max = 30,
                           seed = 20260411, burn_in = 300) {
  sims <- simulate_ar2_reps(
    phi1 = phi1,
    phi2 = phi2,
    n = n,
    n_rep = n_rep,
    seed = seed,
    burn_in = burn_in
  )

  emp_acfs <- get_empirical_acfs(sims, lag_max = lag_max)
  theo_acf <- get_theoretical_acf(phi1, phi2, lag_max = lag_max)
  roots <- get_stationarity_roots(phi1, phi2)
  stationary <- is_stationary_ar2(phi1, phi2)

  list(
    phi1 = phi1,
    phi2 = phi2,
    stationary = stationary,
    roots = roots,
    simulations = sims,
    empirical_acfs = emp_acfs,
    theoretical_acf = theo_acf
  )
}

