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
     xlab = "(Year)",
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
# ---------- 3.3 Plot forecasts with labels ----------
# Build forecast table
forecast_tbl <- data.frame(
  time       = format(Dtest$time, "%Y-%m"),
  year       = x_test,
  predicted  = pred[, "fit"],
  lower_95   = pred[, "lwr"],
  upper_95   = pred[, "upr"]
)


forecast_tbl[, 3:5] <- round(forecast_tbl[, 3:5], 4)

# Save table
write.csv(forecast_tbl,
          file = file.path("output", "tables", "03_ols_forecast_2024.csv"),
          row.names = FALSE)

cat("\nSaved forecast table to output/tables/03_ols_forecast_2024.csv\n")

# Print first rows nicely
print(forecast_tbl)


png(file.path("report", "figures", "03_ols_forecast_test.png"),
    width = 1200, height = 700, res = 150)

x_all <- c(x_train, x_test)
y_all <- c(y_train, pred[, "lwr"], pred[, "upr"])

plot(x_train, y_train,
     pch = 16,
     col = "black",
     xlim = range(x_all),
     ylim = range(y_all),
     xlab = "Year",
     ylab = "Total vehicles (millions)",
     main = "OLS Forecast with 95% Prediction Interval")

# Training fitted line
ord_train <- order(x_train)
lines(x_train[ord_train], y_hat[ord_train],
      col = "blue", lwd = 2)

# Forecast line
ord_test <- order(x_test)
lines(x_test[ord_test], pred[ord_test, "fit"],
      col = "red", lwd = 2)

# Prediction intervals
lines(x_test[ord_test], pred[ord_test, "lwr"],
      col = "red", lty = 2)
lines(x_test[ord_test], pred[ord_test, "upr"],
      col = "red", lty = 2)

# Vertical line separating train/test
abline(v = min(x_test), lty = 3, col = "darkgray")

# Legend
legend("topleft",
       legend = c("Observed (training)",
                  "Fitted mean (training)",
                  "Forecast (test)",
                  "95% prediction interval",
                  "Train/Test split"),
       col = c("black", "blue", "red", "red", "darkgray"),
       pch = c(16, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 2, 3),
       lwd = c(NA, 2, 2, 1, 1),
       bty = "n")

dev.off()



# ---------- 3.3 Plot forecasts + prediction intervals + actual test data ----------
png(file.path("report", "figures", "03_ols_forecast_test_vs_actual.png"),
    width = 1200, height = 700, res = 150)

# Ensure ordering for nice lines
ord_train <- order(x_train)
ord_test  <- order(x_test)

# Axis ranges must include training, forecast intervals, AND actual test values
x_all <- c(x_train, x_test)
y_all <- c(y_train,
           pred[, "lwr"], pred[, "upr"],
           y_test)

plot(x_train, y_train,
     pch = 16,
     col = "black",
     xlim = range(x_all),
     ylim = range(y_all),
     xlab = "Year",
     ylab = "Total vehicles (millions)",
     main = "OLS Forecast vs Actual Test Data (95% Prediction Interval)")

# Training fitted line (estimated mean on training set)
lines(x_train[ord_train], y_hat[ord_train], col = "blue", lwd = 2)

# Train/Test split marker
abline(v = min(x_test), lty = 3, col = "darkgray")

# Forecast mean line (test period)
lines(x_test[ord_test], pred[ord_test, "fit"], col = "red", lwd = 2)

# Prediction interval (test period)
lines(x_test[ord_test], pred[ord_test, "lwr"], col = "red", lty = 2)
lines(x_test[ord_test], pred[ord_test, "upr"], col = "red", lty = 2)

# Actual test observations (held out)
points(x_test, y_test, pch = 17, col = "darkgreen")  # triangle markers

legend("topleft",
       legend = c("Observed (training)",
                  "Fitted mean (training)",
                  "Train/Test split",
                  "Forecast mean (test)",
                  "95% prediction interval",
                  "Actual (test)"),
       col = c("black", "blue", "darkgray", "red", "red", "darkgreen"),
       pch = c(16, NA, NA, NA, NA, 17),
       lty = c(NA, 1, 3, 1, 2, NA),
       lwd = c(NA, 2, 1, 2, 1, NA),
       bty = "n")

dev.off()

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
points(x_test, forecast_tbl$y_pred, pch = 17)

# prediction interval "band" for test points (as vertical segments)
segments(x0 = x_test, y0 = forecast_tbl$pi_lower,
         x1 = x_test, y1 = forecast_tbl$pi_upper)

# optional: connect forecast means
ord2 <- order(x_test)
lines(x_test[ord2], forecast_tbl$y_pred[ord2])

legend("topleft",
       legend = c("Training obs", "Fitted mean", "Forecast mean (2024)", "95% PI (2024)"),
       pch = c(16, NA, 17, NA),
       lty = c(NA, 1, 1, NA),
       bty = "n")
dev.off()

# ---------- 3.6 Residual diagnostics ----------
# ---------- Studentized Residual Plot ----------
dir.create(file.path("report", "figures"), recursive = TRUE, showWarnings = FALSE)

png(file.path("report", "figures", "03_studentized_residuals.png"),
    width = 1200, height = 700, res = 150)

stud_res <- rstudent(fit)

plot(x_train, stud_res,
     type = "b",
     pch = 16,
     col = "black",
     xlab = "Year",
     ylab = "Studentized Residuals",
     main = "Studentized Residuals over Time")

abline(h = 0, col = "red")
abline(h = c(-2, 2), col = "blue", lty = 2)

dev.off()

# Residuals vs fitted
png(file.path("report", "figures", "03_resid_vs_fitted.png"), width = 1200, height = 700, res = 150)
plot(fitted(fit), e,
     pch = 16,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs fitted")
abline(h = 0)
dev.off()

# Histogram
png(file.path("report", "figures", "03_resid_hist.png"), width = 1200, height = 700, res = 150)
hist(e, main = "Residual histogram", xlab = "Residuals")
dev.off()

# QQ plot
png(file.path("report", "figures", "03_resid_qq.png"), width = 1200, height = 700, res = 150)
qqnorm(e, main = "Normal Q-Q plot of residuals")
qqline(e)
dev.off()

# ACF of residuals
png(file.path("report", "figures", "03_resid_acf.png"), width = 1200, height = 700, res = 150)
acf(e, main = "ACF of residuals")
dev.off()

# Ljung-Box test for autocorrelation (choose a lag; 12 is common for monthly)
lb <- Box.test(e, lag = 12, type = "Ljung-Box")
cat("\nLjung-Box test (lag 12):\n")
print(lb)

# Shapiro-Wilk normality test (note: can be sensitive)
sh <- shapiro.test(e)
cat("\nShapiro-Wilk test:\n")
print(sh)

cat("\nResidual SD (training):", sd(e), "\n")
