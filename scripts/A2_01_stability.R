# ============================================================
# AR(2) stability, invertibility, and autocorrelation analysis
# Model:
#   X_t + phi1 * X_{t-1} + phi2 * X_{t-2} = epsilon_t
# with epsilon_t ~ WN(0, 1)
# ============================================================

# ----------------------------
# 1. Define coefficients
# ----------------------------
phi1 <- -0.7
phi2 <- -0.2

cat("AR(2) model:\n")
cat("X_t +", phi1, "* X_{t-1} +", phi2, "* X_{t-2} = epsilon_t\n\n")

# In recursive form:
# X_t = -phi1 * X_{t-1} - phi2 * X_{t-2} + epsilon_t
a1 <- -phi1
a2 <- -phi2

cat("Recursive form:\n")
cat("X_t =", a1, "* X_{t-1} +", a2, "* X_{t-2} + epsilon_t\n\n")


# ----------------------------
# 2. Question 1.1:
#    Check stationarity using roots
# ----------------------------
# The AR polynomial is:
#   1 + phi1*z + phi2*z^2 = 0
# The process is stationary if all roots satisfy |z| > 1

coef_poly <- c(phi2, phi1, 1)   # corresponds to phi2*z^2 + phi1*z + 1 = 0
roots <- polyroot(coef_poly)

cat("Roots of the characteristic equation:\n")
print(roots)
cat("\nModuli of roots:\n")
print(Mod(roots))
cat("\n")

if (all(Mod(roots) > 1)) {
  cat("Conclusion: The process is stationary, because all roots lie outside the unit circle.\n\n")
} else {
  cat("Conclusion: The process is NOT stationary, because at least one root lies inside the unit circle.\n\n")
}


# ----------------------------
# 3. Question 1.2:
#    Is the process invertible?
# ----------------------------
# Invertibility is formally an MA-property, not an AR-property.
# However, a stationary AR process has a convergent MA(infinity) representation.
# So in this sense, the AR(2) can be inverted to an MA(infinity) form.

cat("Invertibility:\n")
cat("A pure AR(2) model is not usually discussed in terms of MA invertibility.\n")
cat("However, since the process is stationary, it admits a convergent MA(infinity) representation.\n")
cat("Thus, in that sense, it is invertible.\n\n")


# ----------------------------
# 4. Question 1.3:
#    Compute rho(k) as a function of phi1 and phi2
# ----------------------------
# For
#   X_t + phi1*X_{t-1} + phi2*X_{t-2} = epsilon_t
# the ACF satisfies:
#   rho(0) = 1
#   rho(1) = -phi1 / (1 + phi2)
#   rho(k) = -phi1*rho(k-1) - phi2*rho(k-2),  k >= 2

rho1 <- -phi1 / (1 + phi2)

cat("Autocorrelation formula:\n")
cat("rho(0) = 1\n")
cat("rho(1) = -phi1 / (1 + phi2) =", rho1, "\n")
cat("rho(k) = -phi1 * rho(k-1) - phi2 * rho(k-2), for k >= 2\n\n")


# ----------------------------
# 5. Question 1.4:
#    Compute and plot rho(k) up to lag 30
# ----------------------------
n_lag <- 30
rho <- numeric(n_lag + 1)   # store rho(0), rho(1), ..., rho(30)

rho[1] <- 1                 # rho(0)
rho[2] <- rho1              # rho(1)

for (k in 3:(n_lag + 1)) {
  rho[k] <- -phi1 * rho[k - 1] - phi2 * rho[k - 2]
}

lags <- 0:n_lag

cat("Autocorrelation values up to lag 30:\n")
acf_table <- data.frame(
  lag = lags,
  rho = round(rho, 6)
)
print(acf_table)
cat("\n")


# ----------------------------
# 6. Plot the ACF
# ----------------------------
plot(
  lags, rho,
  type = "h",
  lwd = 2,
  xlab = "Lag k",
  ylab = expression(rho(k)),
  main = "Autocorrelation Function for AR(2)",
  ylim = c(min(0, rho), 1)
)
points(lags, rho, pch = 16)
abline(h = 0, col = "red", lty = 2)


# ----------------------------
# 7. Optional: compare with built-in ARMAacf
# ----------------------------
# R uses the convention:
#   X_t = ar[1] X_{t-1} + ar[2] X_{t-2} + epsilon_t
# so we must pass ar = c(-phi1, -phi2)

rho_builtin <- ARMAacf(ar = c(-phi1, -phi2), lag.max = n_lag)

cat("Check against built-in ARMAacf:\n")
comparison <- data.frame(
  lag = lags,
  recursive = round(rho, 6),
  builtin   = round(as.numeric(rho_builtin), 6)
)
print(comparison)