# Assignment 3 - Question 3: ARX model for the heating of a box
#
# Run this script from the repo root, or open it in RStudio while the working
# directory is set to the repo root (Session > Set Working Directory > To Project).
# The script auto-detects its location when opened in RStudio.

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Working directory + output folders
# ---------------------------------------------------------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  repo_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(repo_root)
}
cat("Working directory:", getwd(), "\n")

for (d in c("report/figures", "output/tables", "output/models")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# 3.1  Read data and plot the three non-lagged time series
# ---------------------------------------------------------------------------
box <- read.csv("data/box_data_60min.csv")
box$tdate <- as.POSIXct(box$tdate, tz = "UTC")

png("report/figures/q3_31_time_series.png", width = 900, height = 700)
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1), oma = c(3, 0, 0, 0))
plot(box$tdate, box$Ph,     type = "l", ylab = "Ph (W)",     xlab = "", main = "Heating power")
plot(box$tdate, box$Tdelta, type = "l", ylab = "Tdelta (°C)", xlab = "", main = "Temperature difference (inside - outside)")
plot(box$tdate, box$Gv,     type = "l", ylab = "Gv (W/m²)",  xlab = "", main = "Vertical solar radiation")
mtext("Time", side = 1, outer = TRUE, line = 1)
dev.off()

# Comment (3.1):
# Ph is ~50-60 W at night and drops when Gv is high (solar gain replaces electrical heating).
# Tdelta tracks inside-minus-outside temperature; higher when colder outside.
# Gv spikes during daytime only. Negative Ph-Gv relationship and positive Ph-Tdelta
# relationship are expected: more radiation or milder outside → less heating needed.

# ---------------------------------------------------------------------------
# 3.2  Train / test split  (thour 1..167 = training, 168..231 = test)
# ---------------------------------------------------------------------------
train <- box[box$thour <= 167, ]
test  <- box[box$thour >  167, ]

cat("Training rows:", nrow(train),
    "| Last training point:", format(max(train$tdate)), "\n")
cat("Test rows:", nrow(test),
    "| First test point:", format(min(test$tdate)), "\n")

# ---------------------------------------------------------------------------
# 3.3  Exploratory analysis on training data
# ---------------------------------------------------------------------------

# --- Scatter plots ---------------------------------------------------------
png("report/figures/q3_33_scatter.png", width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(train$Tdelta, train$Ph, pch = 16, cex = 0.5,
     xlab = "Tdelta (°C)", ylab = "Ph (W)", main = "Ph vs Tdelta")
plot(train$Gv,     train$Ph, pch = 16, cex = 0.5,
     xlab = "Gv (W/m²)",   ylab = "Ph (W)", main = "Ph vs Gv")
dev.off()

# --- ACF / PACF of Ph ------------------------------------------------------
png("report/figures/q3_33_acf_ph.png", width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
acf(train$Ph,  lag.max = 24, main = "ACF of Ph")
pacf(train$Ph, lag.max = 24, main = "PACF of Ph")
dev.off()

# --- Raw CCF: qualitative only (see 3.4 for prewhitened version) -----------
# NOTE (Textbook §8.7, p.269): Raw CCF is confounded by autocorrelation in the inputs.
# Use for qualitative direction only; prewhitening (3.4) gives the correct impulse response.
png("report/figures/q3_33_ccf_raw.png", width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
ccf(train$Tdelta, train$Ph, lag.max = 24, main = "CCF (raw): Tdelta -> Ph")
ccf(train$Gv,     train$Ph, lag.max = 24, main = "CCF (raw): Gv -> Ph")
dev.off()

# Comment (3.3):
# ACF of Ph decays slowly → strong AR dynamics needed. PACF cuts off after 1-2 lags.
# Raw CCF: Tdelta → Ph positive and persistent; Gv → Ph strongly negative.
# Scatter shows approximate linearity for Tdelta; Gv has a floor effect at Gv = 0
# (night-time). Because both inputs are autocorrelated, the raw CCF is biased — the
# prewhitened CCF in 3.4 gives the unbiased impulse response estimate.

# ---------------------------------------------------------------------------
# 3.4  Impulse response via prewhitening (Textbook §8.7, p.269)
# ---------------------------------------------------------------------------
max_lag_ir <- 10

# Helper: apply an AR filter (ar_fit object) to a series x
apply_ar_filter <- function(x, ar_fit) {
  p <- ar_fit$order
  if (p == 0) return(x)                 # AR(0) = identity filter
  n <- length(x)
  out <- x[(p + 1):n]
  for (k in seq_len(p)) {
    out <- out - ar_fit$ar[k] * x[(p + 1 - k):(n - k)]
  }
  out
}

# Step 1: fit AR to each input (order selected by AIC)
ar_td <- ar(train$Tdelta, order.max = 15, aic = TRUE, method = "ols")
ar_gv <- ar(train$Gv,     order.max = 15, aic = TRUE, method = "ols")
cat("AR order for Tdelta:", ar_td$order, "\n")
cat("AR order for Gv    :", ar_gv$order, "\n")

# Step 2: prewhiten inputs; apply same filter to Ph
alpha_td <- apply_ar_filter(train$Tdelta, ar_td)   # whitened Tdelta
beta_td  <- apply_ar_filter(train$Ph,     ar_td)   # Ph filtered by Tdelta AR

alpha_gv <- apply_ar_filter(train$Gv, ar_gv)       # whitened Gv
beta_gv  <- apply_ar_filter(train$Ph, ar_gv)       # Ph filtered by Gv AR

# Step 3: CCF of prewhitened pairs → estimated impulse response
png("report/figures/q3_34_impulse_response.png", width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
ccf(alpha_td, beta_td, lag.max = max_lag_ir,
    main = "Impulse response (prewhitened): Tdelta -> Ph", ylab = "CCF")
ccf(alpha_gv, beta_gv, lag.max = max_lag_ir,
    main = "Impulse response (prewhitened): Gv -> Ph", ylab = "CCF")
dev.off()

# Comment (3.4):
# Prewhitening removes input autocorrelation so the CCF estimates the true impulse response.
# Tdelta → Ph: positive response at lag 0 persisting several hours (thermal inertia of box).
# Gv → Ph: negative at lag 0 (solar gain reduces heating); effect decays within ~2 hours.
# The lag structure informs which lags to include in the ARX model.

# ---------------------------------------------------------------------------
# 3.5  Linear regression (no AR): Ph = ω1·Tdelta + ω2·Gv + ε
# ---------------------------------------------------------------------------
fit_lm <- lm(Ph ~ Tdelta + Gv, data = train)
cat("\n--- 3.5 Linear regression ---\n")
print(summary(fit_lm))

png("report/figures/q3_35_lm_diagnostics.png", width = 1100, height = 700)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot(train$tdate, train$Ph, type = "l", col = "black",
     ylab = "Ph (W)", xlab = "Time", main = "3.5 OLS: observed vs fitted")
lines(train$tdate, fitted(fit_lm), col = "red", lty = 2)
legend("topright", c("Observed", "Fitted"), col = c("black", "red"), lty = 1:2, cex = 0.8)
plot(train$tdate, residuals(fit_lm), type = "l",
     ylab = "Residual (W)", xlab = "Time", main = "Residuals over time")
abline(h = 0, col = "gray")
acf(residuals(fit_lm),  lag.max = 24, main = "ACF of residuals")
pacf(residuals(fit_lm), lag.max = 24, main = "PACF of residuals")
ccf(train$Tdelta, residuals(fit_lm), lag.max = 24, main = "CCF: Tdelta vs resid")
ccf(train$Gv,     residuals(fit_lm), lag.max = 24, main = "CCF: Gv vs resid")
dev.off()

write.csv(as.data.frame(summary(fit_lm)$coefficients),
          "output/tables/q3_35_lm_coefficients.csv")

# Comment (3.5):
# Significant autocorrelation in residuals: ACF/PACF show clear AR pattern. CCF of residuals
# with inputs also shows remaining structure. OLS i.i.d. assumption is violated → need ARX.

# ---------------------------------------------------------------------------
# 3.6  First-order ARX: Ph = −ϕ1·Ph_{t-1} + ω1·Tdelta + ω2·Gv + ε
# ---------------------------------------------------------------------------
fit_arx1 <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train)
cat("\n--- 3.6 ARX(1) ---\n")
print(summary(fit_arx1))

png("report/figures/q3_36_arx1_diagnostics.png", width = 1100, height = 700)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot(train$tdate, train$Ph, type = "l", col = "black",
     ylab = "Ph (W)", xlab = "Time", main = "3.6 ARX(1): observed vs fitted")
lines(train$tdate, fitted(fit_arx1), col = "red", lty = 2)
legend("topright", c("Observed", "Fitted"), col = c("black", "red"), lty = 1:2, cex = 0.8)
plot(train$tdate, residuals(fit_arx1), type = "l",
     ylab = "Residual (W)", xlab = "Time", main = "Residuals over time")
abline(h = 0, col = "gray")
acf(residuals(fit_arx1),  lag.max = 24, main = "ACF of residuals")
pacf(residuals(fit_arx1), lag.max = 24, main = "PACF of residuals")
ccf(train$Tdelta, residuals(fit_arx1), lag.max = 24, main = "CCF: Tdelta vs resid")
ccf(train$Gv,     residuals(fit_arx1), lag.max = 24, main = "CCF: Gv vs resid")
dev.off()

write.csv(as.data.frame(summary(fit_arx1)$coefficients),
          "output/tables/q3_36_arx1_coefficients.csv")

# Comment (3.6):
# Adding Ph.l1 substantially reduces residual autocorrelation; fit improves markedly.
# Some structure may remain at longer lags → higher orders evaluated in 3.7.

# ---------------------------------------------------------------------------
# 3.7  ARX order selection: AIC and BIC vs. order p = 1..max_order
#      Order p means p AR lags and lags 0..p for each exogenous input.
# ---------------------------------------------------------------------------
max_order <- 10

aic_vec <- numeric(max_order)
bic_vec <- numeric(max_order)

for (p in seq_len(max_order)) {
  regressors <- c(paste0("Ph.l",     1:p),
                  paste0("Tdelta.l", 0:p),
                  paste0("Gv.l",     0:p))
  # Keep only columns that actually exist (caps at the lag depth in the data)
  regressors <- intersect(regressors, names(box))
  fit_p      <- lm(as.formula(paste("Ph ~", paste(regressors, collapse = " + "))),
                   data = train)
  aic_vec[p] <- AIC(fit_p)
  bic_vec[p] <- BIC(fit_p)
}

ic_df <- data.frame(order = 1:max_order, AIC = aic_vec, BIC = bic_vec)
write.csv(ic_df, "output/tables/q3_37_aic_bic.csv", row.names = FALSE)

png("report/figures/q3_37_aic_bic.png", width = 700, height = 480)
par(mar = c(4, 4, 2, 1))
ylims <- range(c(aic_vec, bic_vec))
plot(1:max_order, aic_vec, type = "b", pch = 16, col = "blue",
     ylim = ylims, xlab = "Model order (p)",
     ylab = "Information criterion", main = "AIC and BIC vs. ARX order")
lines(1:max_order, bic_vec, type = "b", pch = 16, col = "red", lty = 2)
abline(v = which.min(aic_vec), col = "blue", lty = 3)
abline(v = which.min(bic_vec), col = "red",  lty = 3)
legend("topright", c("AIC", "BIC"), col = c("blue", "red"), lty = 1:2, pch = 16)
dev.off()

cat("AIC minimum at order:", which.min(aic_vec), "\n")
cat("BIC minimum at order:", which.min(bic_vec), "\n")

# Comment (3.7):
# BIC penalises extra parameters more heavily (k·ln(N) vs 2k for AIC) and therefore
# selects a more parsimonious model. AIC may favour a slightly higher order. The chosen
# order balances in-sample fit with out-of-sample parsimony.

# ---------------------------------------------------------------------------
# 3.8  One-step RMSE on the test set vs. model order
#      RMSE = sqrt( (1/64) * sum_{t=168}^{231} ε²_{t|t-1} )
# ---------------------------------------------------------------------------
rmse_vec <- numeric(max_order)

for (p in seq_len(max_order)) {
  regressors <- intersect(
    c(paste0("Ph.l", 1:p), paste0("Tdelta.l", 0:p), paste0("Gv.l", 0:p)),
    names(box)
  )
  fit_p      <- lm(as.formula(paste("Ph ~", paste(regressors, collapse = " + "))),
                   data = train)
  pred_test  <- predict(fit_p, newdata = test)   # one-step: observed lags in test data
  rmse_vec[p] <- sqrt(mean((test$Ph - pred_test)^2))
}

rmse_df <- data.frame(order = 1:max_order, RMSE = rmse_vec)
write.csv(rmse_df, "output/tables/q3_38_rmse.csv", row.names = FALSE)

png("report/figures/q3_38_rmse.png", width = 700, height = 480)
par(mar = c(4, 4, 2, 1))
plot(1:max_order, rmse_vec, type = "b", pch = 16, col = "darkgreen",
     xlab = "Model order (p)", ylab = "RMSE (W)",
     main = "One-step RMSE on test set vs. ARX order")
abline(v = which.min(rmse_vec), col = "darkgreen", lty = 2)
dev.off()

cat("RMSE minimum at order:", which.min(rmse_vec), "\n")

# Comment (3.8):
# Test-set RMSE measures actual predictive performance on unseen data. It may point to a
# different order than AIC/BIC (often lower, avoiding overfitting).

# ---------------------------------------------------------------------------
# 3.9  Multi-step simulation over the full period
#      AR lags are computed iteratively from the simulated series;
#      input lags (Tdelta, Gv) use the observed values throughout.
# ---------------------------------------------------------------------------
# Choose model order (update p_best manually after inspecting 3.7 and 3.8)
p_best <- which.min(bic_vec)   # <-- change this if AIC / RMSE suggest a different order
cat("\nUsing ARX(", p_best, ") for multi-step simulation.\n")

regressors_best <- intersect(
  c(paste0("Ph.l", 1:p_best), paste0("Tdelta.l", 0:p_best), paste0("Gv.l", 0:p_best)),
  names(box)
)
fit_best <- lm(as.formula(paste("Ph ~", paste(regressors_best, collapse = " + "))),
               data = train)
saveRDS(fit_best, paste0("output/models/q3_arx_order", p_best, ".rds"))
write.csv(as.data.frame(summary(fit_best)$coefficients),
          paste0("output/tables/q3_39_arx", p_best, "_coefficients.csv"))

# Iterative simulation
Ph_sim <- box$Ph                       # initialise with observed values
coefs  <- coef(fit_best)               # named coefficient vector including "(Intercept)"
ar_lags    <- grep("^Ph\\.l",     regressors_best, value = TRUE)
input_lags <- grep("^(Tdelta|Gv)", regressors_best, value = TRUE)

for (t in seq(p_best + 1, nrow(box))) {
  ar_part    <- sum(coefs[ar_lags]    * sapply(ar_lags,    function(nm) {
    lag_k <- as.integer(sub("Ph.l", "", nm))
    Ph_sim[t - lag_k]          # use the iterated (simulated) value
  }))
  input_part <- sum(coefs[input_lags] * as.numeric(box[t, input_lags]))
  Ph_sim[t]  <- coefs["(Intercept)"] + ar_part + input_part
}

sim_df <- data.frame(tdate  = box$tdate,
                     Ph_obs = box$Ph,
                     Ph_sim = Ph_sim,
                     set    = ifelse(box$thour <= 167, "train", "test"))
write.csv(sim_df, "output/tables/q3_39_simulation.csv", row.names = FALSE)

png("report/figures/q3_39_multistep_simulation.png", width = 1000, height = 480)
par(mar = c(4, 4, 2, 1))
plot(sim_df$tdate, sim_df$Ph_obs, type = "l", col = "black",
     ylab = "Ph (W)", xlab = "Time",
     main = paste0("Multi-step simulation — ARX(", p_best, ")"))
lines(sim_df$tdate, sim_df$Ph_sim, col = "red", lty = 2)
abline(v = as.numeric(max(train$tdate)), col = "gray60", lty = 3, lwd = 1.5)
legend("topright",
       c("Observed", "Simulated", "Train / test split"),
       col = c("black", "red", "gray60"), lty = c(1, 2, 3), lwd = c(1, 1, 1.5), cex = 0.8)
dev.off()

# Comment (3.9):
# Multi-step simulation propagates AR lags from the simulated (not observed) series, so
# errors accumulate over time. The model may still capture the broad daily pattern while
# showing larger deviations than one-step predictions. In a real-time setting the observed
# inputs (Tdelta, Gv) would also need to be forecast, adding further uncertainty — unless
# they can be obtained from a weather forecast service.

# ---------------------------------------------------------------------------
# 3.10  Summary
# ---------------------------------------------------------------------------
cat("\n=== 3.10 Summary ===\n")
cat("Best order by AIC :", which.min(aic_vec), "\n")
cat("Best order by BIC :", which.min(bic_vec), "\n")
cat("Best order by RMSE:", which.min(rmse_vec), "\n")
cat("Model used for simulation: ARX(", p_best, ")\n\n")
cat("Coefficient table for chosen model:\n")
print(summary(fit_best)$coefficients)

# Conclusion (3.10):
# The ARX model successfully captures the transfer-function dynamics between weather
# variables (Tdelta, Gv) and heating demand (Ph). Including AR lags removes the residual
# autocorrelation present in the plain regression. BIC selects a parsimonious order while
# the test-set RMSE confirms predictive accuracy. Multi-step simulation shows the model
# tracks the daily heating pattern well in-sample and reasonably in the test period,
# though accuracy degrades as the simulation horizon grows.
