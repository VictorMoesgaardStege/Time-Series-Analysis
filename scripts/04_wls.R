# scripts/04_wls.R
# Step 4: WLS local linear trend model (forgetting weights)

source("R/read_data.R")

# ---------- Load data ----------
dat <- load_bil54(file = file.path("data", "DST_BIL54.csv"))
Dtrain <- dat$train
Dtest  <- dat$test

# ---------- Helpers ----------
make_x <- function(time_vec) {
  lt <- as.POSIXlt(time_vec)
  (lt$year + 1900) + (lt$mon) / 12
}

dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("report", "figures"), recursive = TRUE, showWarnings = FALSE)

x_train <- make_x(Dtrain$time)
y_train <- as.numeric(Dtrain$total)

x_test  <- make_x(Dtest$time)
y_test  <- as.numeric(Dtest$total)

# Design matrix: [1, x]
X <- cbind(1, x_train)
y <- as.matrix(y_train)

N <- nrow(X)
p <- ncol(X)

# ---------- Core WLS function ----------
fit_wls_local <- function(x_train, y_train, lambda) {
  X <- cbind(1, x_train)
  y <- as.matrix(y_train)
  N <- nrow(X); p <- ncol(X)

  # Local variance-cov matrix Σ is diagonal.
  # We choose Σ_ii = 1 / lambda^(N - i)  =>  W = Σ^{-1} has diag weights = lambda^(N - i)
  w <- lambda^(N - (1:N))  # i=N => lambda^0 = 1 (highest at latest timepoint)
  W <- diag(w) # W = Σ^{-1} is diagonal with entries w_i

  XtWX <- t(X) %*% W %*% X # X' W X
  XtWy <- t(X) %*% W %*% y # X' W y

  theta_hat <- solve(XtWX) %*% XtWy # (X' W X)^{-1} X' W y
  y_hat <- as.numeric(X %*% theta_hat) # fitted values on training set. Computes X theta_hat
  e <- y_train - y_hat # residuals on training set

  # Weighted sigma^2 estimate (common pragmatic choice)
  sigma2_hat <- as.numeric(t(e) %*% W %*% e) / (N - p) # computes e' W e / (N - p)
  sigma_hat  <- sqrt(sigma2_hat) # weighted residual standard error

  Vtheta <- sigma2_hat * solve(XtWX) # Var(theta_hat) = sigma^2 * (X' W X)^{-1}
  se_theta <- sqrt(diag(Vtheta)) # standard errors for theta_hat

  list(
    lambda = lambda,
    theta_hat = theta_hat,
    se_theta = se_theta,
    sigma2_hat = sigma2_hat,
    sigma_hat = sigma_hat,
    w = w,
    X = X,
    W = W,
    XtWX_inv = solve(XtWX),
    fitted = y_hat,
    resid = e
  )
}

# ---------- OLS (for comparison in 4.5 plots) ----------
fit_ols <- lm(y_train ~ x_train)
ols_pred <- predict(fit_ols,
                    newdata = data.frame(x_train = x_test),
                    interval = "prediction",
                    level = 0.95)

# ---------- Part 4 with lambda = 0.9 ----------
lambda <- 0.9
wls <- fit_wls_local(x_train, y_train, lambda)

theta1_hat <- as.numeric(wls$theta_hat[1])
theta2_hat <- as.numeric(wls$theta_hat[2])
se1 <- as.numeric(wls$se_theta[1])
se2 <- as.numeric(wls$se_theta[2])

cat("\n=== WLS local trend (lambda =", lambda, ") ===\n")
cat(sprintf("theta1_hat = %.6f (SE %.6f)\n", theta1_hat, se1))
cat(sprintf("theta2_hat = %.6f (SE %.6f)\n", theta2_hat, se2))
cat(sprintf("sigma_hat  = %.6f\n", wls$sigma_hat))

# ---------- 4.1: Describe Σ (save relevant parts) ----------
# Σ is diagonal with entries Σ_ii = 1 / w_i, where w_i = lambda^(N-i).
# For OLS/global model, Σ = I (i.e. equal variances and no weighting).
Sigma_diag <- 1 / wls$w  # diag entries of Σ
Sigma_snip <- data.frame(
  i = c(1:5, (N-4):N),
  time = format(Dtrain$time[c(1:5, (N-4):N)], "%Y-%m"),
  Sigma_ii = Sigma_diag[c(1:5, (N-4):N)],
  weight_wi = wls$w[c(1:5, (N-4):N)]
)

write.csv(Sigma_snip,
          file = file.path("output", "tables", sprintf("04_wls_sigma_snip_lambda_%s.csv", lambda)),
          row.names = FALSE)

# ---------- 4.2: Plot weights vs time ----------

# Note: we can plot weights w_t vs time (calendar time) or vs index t in training set.

# png(file.path("report", "figures", sprintf("04_wls_weights_lambda_%s.png", lambda)),
#     width = 1200, height = 700, res = 150)

# t_idx <- 1:N
# plot(t_idx, wls$w, type = "l", lwd = 2,
#      xlab = "t (index in training set)",
#      ylab = "Weight w_t = λ^(N-t)",
#      main = sprintf("WLS local weights vs time (λ = %.2f)", lambda))
# points(t_idx, wls$w, pch = 16)
# grid()
# dev.off()

# cat("\nHighest weight is at latest training timepoint t = N (", format(Dtrain$time[N], "%Y-%m"), ")\n", sep="")

# ---------- 4.2: Plot weights vs time (calendar time) ----------
png(file.path("report", "figures", 
              sprintf("04_wls_weights_lambda_%s.png", lambda)),
    width = 1200, height = 700, res = 150)

plot(Dtrain$time, wls$w,
     type = "l",
     lwd = 2,
     xlab = "Time",
     ylab = expression(w[t] == lambda^(N - t)),
     main = sprintf("WLS local weights over time (λ = %.2f)", lambda))

points(Dtrain$time, wls$w, pch = 16)
grid()
dev.off()

cat("\nHighest weight is at latest training observation:",
    format(Dtrain$time[N], "%Y-%m"), "\n")

# ---------- 4.3: Sum of weights ----------
sum_w <- sum(wls$w)
cat(sprintf("\nSum of WLS weights (lambda=%.2f): %.6f\n", lambda, sum_w))
cat(sprintf("Sum of weights in OLS (all weights = 1): %d\n", N))

write.csv(
  data.frame(lambda = lambda, N = N, sum_weights = sum_w, sum_weights_ols = N),
  file = file.path("output", "tables", sprintf("04_wls_sum_weights_lambda_%s.csv", lambda)),
  row.names = FALSE
)

# ---------- 4.4: Save parameter table ----------
param_tbl <- data.frame(
  model = "WLS_local",
  lambda = lambda,
  theta1_hat = theta1_hat,
  se_theta1 = se1,
  theta2_hat = theta2_hat,
  se_theta2 = se2,
  sigma_hat = wls$sigma_hat
)

write.csv(param_tbl,
          file = file.path("output", "tables", sprintf("04_wls_params_lambda_%s.csv", lambda)),
          row.names = FALSE)

# Also save a "fit on training"
png(file.path("report", "figures", sprintf("04_wls_fit_train_lambda_%s.png", lambda)),
    width = 1200, height = 700, res = 150)

plot(x_train, y_train,
     pch = 16,
     xlab = "x (year + month/12)",
     ylab = "Total vehicles (millions)",
     main = sprintf("WLS fit on training set (λ = %.2f): observations and estimated mean", lambda))
ord <- order(x_train)
lines(x_train[ord], wls$fitted[ord], lwd = 2)
dev.off()

# ---------- 4.5: Forecast next 12 months (WLS) ----------
# Predictions for test x:
Xnew <- cbind(1, x_test)
y_pred_wls <- as.numeric(Xnew %*% wls$theta_hat)

# Optional: prediction intervals for WLS (approx, using weighted sigma_hat)
# Var(prediction at x0) = sigma^2 * (1 + x0' (X'WX)^{-1} x0)
tcrit <- qt(0.975, df = N - p)
pred_var <- apply(Xnew, 1, function(x0) {
  x0 <- as.matrix(x0, ncol = 1)
  as.numeric(wls$sigma2_hat * (1 + t(x0) %*% wls$XtWX_inv %*% x0))
})
pred_se <- sqrt(pred_var)
pi_lower_wls <- y_pred_wls - tcrit * pred_se
pi_upper_wls <- y_pred_wls + tcrit * pred_se

forecast_tbl <- data.frame(
  time = format(Dtest$time, "%Y-%m"),
  x = x_test,
  wls_pred = y_pred_wls,
  wls_pi_lower = pi_lower_wls,
  wls_pi_upper = pi_upper_wls,
  ols_pred = as.numeric(ols_pred[, "fit"]),
  ols_pi_lower = as.numeric(ols_pred[, "lwr"]),
  ols_pi_upper = as.numeric(ols_pred[, "upr"])
)

write.csv(forecast_tbl,
          file = file.path("output", "tables", sprintf("04_wls_forecast_2024_lambda_%s.csv", lambda)),
          row.names = FALSE)

cat("\nSaved forecast table to:",
    file.path("output", "tables", sprintf("04_wls_forecast_2024_lambda_%s.csv", lambda)), "\n")

# Plot training obs + OLS and WLS predictions for test (with PI segments)
png(file.path("report", "figures", sprintf("04_ols_vs_wls_forecast_lambda_%s.png", lambda)),
    width = 1200, height = 800, res = 150)

x_all <- c(x_train, x_test)
y_all <- c(y_train, y_test)
plot(x_train, y_train,
     pch = 16,
     xlab = "x (year + month/12)",
     ylab = "Total vehicles (millions)",
     main = sprintf("Forecast comparison for 2024: OLS vs WLS (λ = %.2f)", lambda),
     xlim = range(x_all),
     ylim = (range(y_all)) + 0.12)

# OLS forecast means + PI
points(x_test, forecast_tbl$ols_pred, pch = 17)
segments(x0 = x_test, y0 = forecast_tbl$ols_pi_lower,
         x1 = x_test, y1 = forecast_tbl$ols_pi_upper)

# WLS forecast means + PI
points(x_test, forecast_tbl$wls_pred, pch = 15)
segments(x0 = x_test, y0 = forecast_tbl$wls_pi_lower,
         x1 = x_test, y1 = forecast_tbl$wls_pi_upper)

legend("topleft",
       legend = c("Training obs", "OLS forecast mean", "OLS 95% PI", "WLS forecast mean", "WLS 95% PI"),
       pch = c(16, 17, NA, 15, NA),
       lty = c(NA, NA, 1, NA, 1),
       bty = "n")
grid()
dev.off()

# ---------- 4.6 Optional: repeat for multiple lambdas ----------
lambdas <- c(0.99, 0.9, 0.8, 0.7, 0.6)
multi_params <- data.frame()
multi_forecasts <- data.frame()

for (lam in lambdas) {
  wfit <- fit_wls_local(x_train, y_train, lam)

  multi_params <- rbind(multi_params, data.frame(
    lambda = lam,
    theta1_hat = as.numeric(wfit$theta_hat[1]),
    theta2_hat = as.numeric(wfit$theta_hat[2]),
    se_theta1 = as.numeric(wfit$se_theta[1]),
    se_theta2 = as.numeric(wfit$se_theta[2]),
    sigma_hat = wfit$sigma_hat
  ))

  ypred <- as.numeric(cbind(1, x_test) %*% wfit$theta_hat)
  multi_forecasts <- rbind(multi_forecasts, data.frame(
    lambda = lam,
    time = format(Dtest$time, "%Y-%m"),
    x = x_test,
    y_pred = ypred
  ))
}

write.csv(multi_params,
          file = file.path("output", "tables", "04_wls_params_multiple_lambdas.csv"),
          row.names = FALSE)
write.csv(multi_forecasts,
          file = file.path("output", "tables", "04_wls_forecast_multiple_lambdas_2024.csv"),
          row.names = FALSE)

# Plot how forecasts change with lambda (means only)
png(file.path("report", "figures", "04_wls_forecast_multiple_lambdas.png"),
    width = 1200, height = 700, res = 150)

plot(x_train, y_train, pch = 16,
     xlab = "x (year + month/12)",
     ylab = "Total vehicles (millions)",
     main = "WLS forecasts for 2024 for multiple λ (means only)",
     xlim = range(c(x_train, x_test)))

# Add OLS mean line over test for reference
lines(x_test[order(x_test)], forecast_tbl$ols_pred[order(x_test)], lwd = 2)

# Add each lambda's forecast as a line
for (lam in lambdas) {
  sub <- multi_forecasts[multi_forecasts$lambda == lam, ]
  ord <- order(sub$x)
  lines(sub$x[ord], sub$y_pred[ord], lwd = 2)
}

legend("topleft",
       legend = c("Training obs", "OLS mean (2024)", paste0("WLS mean, λ=", lambdas)),
       pch = c(16, NA, rep(NA, length(lambdas))),
       lty = c(NA, 1, rep(1, length(lambdas))),
       bty = "n")
grid()
dev.off()

cat("\nSaved multi-lambda parameter+forecast tables and plots.\n")