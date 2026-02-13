# scripts/03_ols.R
# Step 3: OLS global linear trend model

source("R/read_data.R")

# ---------- Load data ----------
dat <- load_bil54(file = file.path("data", "DST_BIL54.csv"))
Dtrain <- dat$train
Dtest  <- dat$test

# ---------- Construct x_t ----------
# x_t = year + month/12 (Jan -> year + 0/12, Feb -> year + 1/12, ...)
make_x <- function(time_vec) {
  lt <- as.POSIXlt(time_vec)
  (lt$year + 1900) + (lt$mon) / 12
}

x_train <- make_x(Dtrain$time)
y_train <- Dtrain$total

x_test <- make_x(Dtest$time)
y_test <- Dtest$total  # not used for fitting, but useful later

# ---------- OLS via matrix formula ----------
# Model: y = X theta + e, X = [1, x]
X <- cbind(1, x_train)
y <- as.matrix(y_train)

theta_hat <- solve(t(X) %*% X) %*% (t(X) %*% y)  # (X'X)^{-1} X'y
theta1_hat <- as.numeric(theta_hat[1])
theta2_hat <- as.numeric(theta_hat[2])

# Residuals + sigma^2 estimate
y_hat <- as.numeric(X %*% theta_hat)
resid <- y_train - y_hat

N <- nrow(X)
p <- ncol(X) # 2 params
sigma2_hat <- sum(resid^2) / (N - p)
sigma_hat  <- sqrt(sigma2_hat)

# Var(theta_hat) = sigma^2 * (X'X)^{-1}
XtX_inv <- solve(t(X) %*% X)
se_theta <- sqrt(diag(sigma2_hat * XtX_inv))
se_theta1 <- as.numeric(se_theta[1])
se_theta2 <- as.numeric(se_theta[2])

cat("\n=== OLS (matrix) estimates on training set ===\n")
cat(sprintf("theta1_hat = %.6f (SE %.6f)\n", theta1_hat, se_theta1))
cat(sprintf("theta2_hat = %.6f (SE %.6f)\n", theta2_hat, se_theta2))
cat(sprintf("sigma_hat  = %.6f\n", sigma_hat))

# ---------- Cross-check with lm ----------
fit <- lm(total ~ x, data = data.frame(total = y_train, x = x_train))
cat("\n=== lm() cross-check ===\n")
print(summary(fit)$coefficients)

# ---------- Create folders ----------
dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("report", "figures"), recursive = TRUE, showWarnings = FALSE)

# ---------- 3.2 Plot estimated mean line with observations ----------
png(file.path("report", "figures", "03_ols_fit_train.png"), width = 1200, height = 700, res = 150)
plot(x_train, y_train,
     pch = 16,
     xlab = "x (year + month/12)",
     ylab = "Total vehicles (millions)",
     main = "OLS fit on training set: observations and estimated mean")
# Sort by x to draw a nice line
ord <- order(x_train)
lines(x_train[ord], y_hat[ord], lwd = 2)
dev.off()

# ---------- 3.3 Forecast test set with prediction intervals ----------
# Use predict.lm to get standard prediction intervals
pred <- predict(fit,
                newdata = data.frame(x = x_test),
                interval = "prediction",
                level = 0.95)

# pred is a matrix with columns: fit, lwr, upr
forecast_tbl <- data.frame(
  time = format(Dtest$time, "%Y-%m"),
  x = x_test,
  y_pred = pred[, "fit"],
  pi_lower = pred[, "lwr"],
  pi_upper = pred[, "upr"]
)

# Save table
write.csv(forecast_tbl,
          file = file.path("output", "tables", "03_ols_forecast_2024.csv"),
          row.names = FALSE)

cat("\nSaved forecast table to output/tables/03_ols_forecast_2024.csv\n")
cat("\nForecasts (first rows):\n")
print(head(forecast_tbl, 12))

# ---------- 3.4 Plot fitted model + training data + forecasts + PI ----------
# For a continuous fitted line over training+test x-range:
x_all <- c(x_train, x_test)
x_grid <- seq(min(x_all), max(x_all), length.out = 300)
y_grid <- predict(fit, newdata = data.frame(x = x_grid))

png(file.path("report", "figures", "03_ols_fit_and_forecast.png"), width = 1200, height = 700, res = 150)

plot(x_train, y_train,
     pch = 16,
     xlab = "x (year + month/12)",
     ylab = "Total vehicles (millions)",
     main = "OLS: training fit and 2024 forecasts with 95% prediction intervals",
     xlim = range(x_all))

# fitted mean line
lines(x_grid, y_grid, lwd = 2)

# forecast points
points(x_test, forecast_tbl$y_pred, pch = 17_
