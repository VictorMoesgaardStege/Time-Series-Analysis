# Assignment 3 - Question 3: ARX model for the heating of a box
#
# Open in RStudio from the repo root, or run with working directory = repo root.
# All figures → report/figures/*.pdf   (include in Overleaf with \includegraphics)
# All tables  → report/tables/*.tex    (include in Overleaf with \input)

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Working directory, packages, output folders
# ---------------------------------------------------------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
}
cat("Working directory:", getwd(), "\n")

if (!requireNamespace("xtable", quietly = TRUE)) install.packages("xtable")
library(xtable)

for (d in c("report/figures", "report/tables", "output/tables", "output/models")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Plot style — mirrors assignment3_q1_stability / stability_helpers.R
# ---------------------------------------------------------------------------
FW   <- 6.2   # one-column Overleaf text width (inches)
FH_S <- 3.2   # single-row / side-by-side panels
FH_D <- 4.8   # 2×3 diagnostic grid
FH_3 <- 7.0   # three stacked panels
FH_2 <- 4.2   # two stacked panels

COL_OBS  <- "black"
COL_FIT  <- "red"
COL_AIC  <- "steelblue"
COL_BIC  <- "firebrick"
COL_RMSE <- "darkgreen"
LWD      <- 2      # default line width throughout

# Results collector — appended throughout; written to markdown at the end
md <- character(0)
md_add <- function(...) md <<- c(md, paste0(...))

# Convenience wrapper: save a named coefficient table as a .tex file
save_coef_tex <- function(fit, file, caption, label) {
  tbl <- as.data.frame(summary(fit)$coefficients)
  colnames(tbl) <- c("Estimate", "Std.\\ Error", "$t$-value", "$p$-value")
  xt <- xtable(tbl, caption = caption, label = label, digits = c(0, 4, 4, 3, 4))
  print(xt,
        booktabs               = TRUE,
        sanitize.text.function = identity,
        file                   = file,
        caption.placement      = "top")
}

# ---------------------------------------------------------------------------
# 3.1  Read data and plot the three non-lagged time series
# ---------------------------------------------------------------------------
box <- read.csv("data/box_data_60min.csv")
box$tdate <- as.POSIXct(box$tdate, tz = "UTC")

pdf("report/figures/q3_31_time_series.pdf", width = FW, height = FH_3)
par(mfrow = c(3, 1), mar = c(2, 4.5, 1.8, 1), oma = c(3, 0, 0, 0))

plot(box$tdate, box$Ph, type = "l", lwd = LWD,
     ylab = expression(P[h]~(W)), xlab = "",
     main = expression("Heating power  " * P[h]))
grid()

plot(box$tdate, box$Tdelta, type = "l", lwd = LWD,
     ylab = expression(T[Delta]~(degree*C)), xlab = "",
     main = expression("Temperature difference  " * T[Delta]))
grid()

plot(box$tdate, box$Gv, type = "l", lwd = LWD,
     ylab = expression(G[v]~(W/m^2)), xlab = "",
     main = expression("Vertical solar radiation  " * G[v]))
grid()

mtext("Time", side = 1, outer = TRUE, line = 1.5)
dev.off()

md_add("# Question 3 Results Summary\n")
md_add("## 3.1 Time Series")
md_add("**Figure:** `report/figures/q3_31_time_series.pdf`")
md_add("**Ph range:** ", round(min(box$Ph),1), " – ", round(max(box$Ph),1), " W")
md_add("**Tdelta range:** ", round(min(box$Tdelta),1), " – ", round(max(box$Tdelta),1), " °C")
md_add("**Gv range:** ", round(min(box$Gv),1), " – ", round(max(box$Gv),1), " W/m²")
md_add("**N total:** ", nrow(box), " hourly observations\n")

# ---------------------------------------------------------------------------
# 3.2  Train / test split
# ---------------------------------------------------------------------------
# thour starts at 19 (first 18 rows lost to lag computation), so the assignment's
# "2013-02-06 00:00" corresponds to thour = 185, which is row 167 of the data.
train <- box[1:167, ]
test  <- box[168:nrow(box), ]
cat("Training rows:", nrow(train), "| Last:", format(max(train$tdate)), "\n")
cat("Test rows    :", nrow(test),  "| First:", format(min(test$tdate)), "\n")

md_add("## 3.2 Train/Test Split")
md_add("**Training set:** ", nrow(train), " obs (rows 1–167, thour up to 185), up to ", format(max(train$tdate)))
md_add("**Test set:** ", nrow(test), " obs, thour 168–231, from ", format(min(test$tdate)), "\n")

# ---------------------------------------------------------------------------
# 3.3  Exploratory analysis (training set)
# ---------------------------------------------------------------------------

# Figure A: scatter plots
pdf("report/figures/q3_33_scatter.pdf", width = FW, height = FH_S)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

plot(train$Tdelta, train$Ph, pch = 16, cex = 0.5,
     xlab = expression(T[Delta]~(degree*C)),
     ylab = expression(P[h]~(W)),
     main = expression(P[h]~vs~T[Delta]))
grid()

plot(train$Gv, train$Ph, pch = 16, cex = 0.5,
     xlab = expression(G[v]~(W/m^2)),
     ylab = expression(P[h]~(W)),
     main = expression(P[h]~vs~G[v]))
grid()
dev.off()

# Figure B: ACF/PACF of Ph and raw CCF (2×2)
pdf("report/figures/q3_33_acf_ccf.pdf", width = FW, height = FH_D)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2.2, 1))

acf(train$Ph,  lag.max = 24, main = expression("ACF of  " * P[h]))
grid()
pacf(train$Ph, lag.max = 24, main = expression("PACF of  " * P[h]))
grid()
ccf(train$Tdelta, train$Ph, lag.max = 24,
    main = expression("CCF (raw):  " * T[Delta] %->% P[h]), ylab = "CCF")
grid()
ccf(train$Gv, train$Ph, lag.max = 24,
    main = expression("CCF (raw):  " * G[v] %->% P[h]), ylab = "CCF")
grid()
dev.off()

md_add("## 3.3 Exploratory Analysis")
md_add("**Figure A (scatter):** `report/figures/q3_33_scatter.pdf`")
md_add("**Figure B (ACF/PACF/CCF):** `report/figures/q3_33_acf_ccf.pdf`")
md_add("**ACF Ph:** slow exponential decay — strong AR dynamics needed")
md_add("**PACF Ph:** cuts off sharply after lag 1–2 — suggests low-order AR")
md_add("**Raw CCF Tdelta→Ph:** positive at lag 0, persists several lags")
md_add("**Raw CCF Gv→Ph:** strongly negative at lag 0 (solar gain replaces heating)")
md_add("**Note:** raw CCF is biased when input is autocorrelated — see prewhitened version in 3.4\n")

# ---------------------------------------------------------------------------
# 3.4  Impulse response via joint FIR regression (lags 0..10 for both inputs)
# ---------------------------------------------------------------------------
max_lag_ir <- 10

fir_dat <- data.frame(
  Ph = train$Ph
)

for (k in 0:max_lag_ir) {
  fir_dat[[paste0("Tdelta_l", k)]] <- lag_vec(train$Tdelta, k)
  fir_dat[[paste0("Gv_l", k)]]     <- lag_vec(train$Gv, k)
}

fir_dat <- na.omit(fir_dat)

fit_fir <- lm(Ph ~ . - 1, data = fir_dat)
cat("\n--- 3.4 Joint FIR regression ---\n")
print(summary(fit_fir))

save_coef_tex(fit_fir,
              file    = "report/tables/q3_34_fir_coef.tex",
              caption = "Estimated FIR coefficients for the joint impulse-response model in 3.4.",
              label   = "tab:q3_fir_coef")

coef_hat <- coef(fit_fir)
h_T <- coef_hat[paste0("Tdelta_l", 0:max_lag_ir)]
h_G <- coef_hat[paste0("Gv_l", 0:max_lag_ir)]

irf_df <- data.frame(
  lag      = 0:max_lag_ir,
  TdeltaIR = as.numeric(h_T),
  GvIR     = as.numeric(h_G)
)

write.csv(irf_df, "output/tables/q3_34_impulse_response.csv", row.names = FALSE)

pdf("report/figures/q3_34_impulse_response.pdf", width = FW, height = FH_S)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.2, 1))

plot(irf_df$lag, irf_df$TdeltaIR,
     type = "h", lwd = 3,
     xlab = "Lag (hours)",
     ylab = expression(h[Tdelta](k)),
     main = expression("Impulse response:  " * T[Delta] %->% P[h]))
abline(h = 0, lty = 2)
grid()

plot(irf_df$lag, irf_df$GvIR,
     type = "h", lwd = 3,
     xlab = "Lag (hours)",
     ylab = expression(h[Gv](k)),
     main = expression("Impulse response:  " * G[v] %->% P[h]))
abline(h = 0, lty = 2)
grid()
dev.off()

peak_T_lag <- irf_df$lag[which.max(abs(irf_df$TdeltaIR))]
peak_G_lag <- irf_df$lag[which.max(abs(irf_df$GvIR))]

md_add("## 3.4 Impulse Response (Joint FIR Regression)")
md_add("**Figure:** `report/figures/q3_34_impulse_response.pdf`")
md_add("**Table:** `report/tables/q3_34_fir_coef.tex`")
md_add("**CSV output:** `output/tables/q3_34_impulse_response.csv`")
md_add("**Method:** joint FIR regression with lags 0–10 for both Tdelta and Gv")
md_add("**Largest |Tdelta impulse-response coefficient| at lag:** ", peak_T_lag, " h")
md_add("**Largest |Gv impulse-response coefficient| at lag:** ", peak_G_lag, " h")
md_add("**Comment:** pre-whitened CCF is a standard textbook identification idea for single-input transfer-function models, but with two potentially correlated inputs it is less reliable. Therefore the impulse responses are estimated jointly here through finite distributed lags.\n")
# ---------------------------------------------------------------------------
# Helper: 2×3 diagnostic panel used for both 3.5 and 3.6
# ---------------------------------------------------------------------------
plot_diagnostics <- function(fit, tdate, Tdelta, Gv, title_prefix) {
  res  <- residuals(fit)
  fitt <- fitted(fit)
  obs  <- fitt + res          # same as the response

  # Panel 1: observed vs fitted
  plot(tdate, obs, type = "l", lwd = LWD, col = COL_OBS,
       ylab = expression(P[h]~(W)), xlab = "Time",
       main = bquote(.(title_prefix) * ": observed vs fitted"))
  lines(tdate, fitt, lwd = LWD, col = COL_FIT, lty = 2)
  legend("topright", c("Observed", "Fitted"),
         col = c(COL_OBS, COL_FIT), lty = 1:2, lwd = LWD, bty = "n")
  grid()

  # Panel 2: residuals over time
  plot(tdate, res, type = "l", lwd = LWD,
       ylab = "Residual (W)", xlab = "Time",
       main = "Residuals over time")
  abline(h = 0, col = "gray50", lwd = 1)
  grid()

  # Panels 3–6: ACF, PACF, CCF
  acf(res,  lag.max = 24, main = "ACF of residuals")
  grid()
  pacf(res, lag.max = 24, main = "PACF of residuals")
  grid()
  ccf(Tdelta, res, lag.max = 24,
      main = expression("CCF:  " * T[Delta]~vs~residuals), ylab = "CCF")
  grid()
  ccf(Gv, res, lag.max = 24,
      main = expression("CCF:  " * G[v]~vs~residuals), ylab = "CCF")
  grid()
}

# ---------------------------------------------------------------------------
# 3.5  Linear regression: Ph = omega1*Tdelta + omega2*Gv + epsilon
# ---------------------------------------------------------------------------
fit_lm <- lm(Ph ~ Tdelta + Gv, data = train)
cat("\n--- 3.5 Linear regression ---\n"); print(summary(fit_lm))

save_coef_tex(fit_lm,
              file    = "report/tables/q3_35_lm_coef.tex",
              caption = "OLS coefficient estimates for the linear regression model (3.5).",
              label   = "tab:lm_coef")

pdf("report/figures/q3_35_lm_diagnostics.pdf", width = FW, height = FH_D)
par(mfrow = c(2, 3), mar = c(4, 4.5, 2.2, 1))
plot_diagnostics(fit_lm, train$tdate, train$Tdelta, train$Gv,
                 title_prefix = "OLS")
dev.off()

lm_s    <- summary(fit_lm)
lm_coef <- coef(lm_s)
md_add("## 3.5 Linear Regression (no AR)")
md_add("**Figure:** `report/figures/q3_35_lm_diagnostics.pdf`")
md_add("**Table:** `report/tables/q3_35_lm_coef.tex`")
md_add("**R²:** ", round(lm_s$r.squared, 4), " | **Adj R²:** ", round(lm_s$adj.r.squared, 4))
md_add("**Residual std error:** ", round(lm_s$sigma, 3), " W")
md_add("**omega1 (Tdelta):** ", round(lm_coef["Tdelta","Estimate"],4),
       " (p=", format(lm_coef["Tdelta","Pr(>|t|)"], digits=3), ")")
md_add("**omega2 (Gv):** ", round(lm_coef["Gv","Estimate"],4),
       " (p=", format(lm_coef["Gv","Pr(>|t|)"], digits=3), ")")
md_add("**Residual ACF:** significant autocorrelation remains → AR terms needed")
md_add("**Residual CCF:** structure visible vs both inputs → transfer function model needed\n")

# ---------------------------------------------------------------------------
# 3.6  ARX(1): Ph = -phi1*Ph_{t-1} + omega1*Tdelta + omega2*Gv + epsilon
# ---------------------------------------------------------------------------
fit_arx1 <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train)
cat("\n--- 3.6 ARX(1) ---\n"); print(summary(fit_arx1))

save_coef_tex(fit_arx1,
              file    = "report/tables/q3_36_arx1_coef.tex",
              caption = "OLS coefficient estimates for the ARX(1) model (3.6).",
              label   = "tab:arx1_coef")

pdf("report/figures/q3_36_arx1_diagnostics.pdf", width = FW, height = FH_D)
par(mfrow = c(2, 3), mar = c(4, 4.5, 2.2, 1))
plot_diagnostics(fit_arx1, train$tdate, train$Tdelta, train$Gv,
                 title_prefix = "ARX(1)")
dev.off()

arx1_s    <- summary(fit_arx1)
arx1_coef <- coef(arx1_s)
md_add("## 3.6 ARX(1)")
md_add("**Figure:** `report/figures/q3_36_arx1_diagnostics.pdf`")
md_add("**Table:** `report/tables/q3_36_arx1_coef.tex`")
md_add("**R²:** ", round(arx1_s$r.squared, 4), " | **Adj R²:** ", round(arx1_s$adj.r.squared, 4))
md_add("**Residual std error:** ", round(arx1_s$sigma, 3), " W")
md_add("**phi1 (Ph.l1):** ", round(arx1_coef["Ph.l1","Estimate"],4),
       " (p=", format(arx1_coef["Ph.l1","Pr(>|t|)"], digits=3), ")")
md_add("**omega1 (Tdelta):** ", round(arx1_coef["Tdelta","Estimate"],4),
       " (p=", format(arx1_coef["Tdelta","Pr(>|t|)"], digits=3), ")")
md_add("**omega2 (Gv):** ", round(arx1_coef["Gv","Estimate"],4),
       " (p=", format(arx1_coef["Gv","Pr(>|t|)"], digits=3), ")")
md_add("**Improvement over OLS:** residual autocorrelation substantially reduced\n")

# ---------------------------------------------------------------------------
# 3.7  Model order selection: AIC and BIC vs. p = 1..10
# ---------------------------------------------------------------------------
max_order <- 10
aic_vec   <- numeric(max_order)
bic_vec   <- numeric(max_order)

for (p in seq_len(max_order)) {
  regs      <- intersect(c(paste0("Ph.l", 1:p),
                           paste0("Tdelta.l", 0:(p-1)),
                           paste0("Gv.l",     0:(p-1))), names(box))
  fit_p     <- lm(as.formula(paste("Ph ~", paste(regs, collapse = " + "))), data = train)
  aic_vec[p] <- AIC(fit_p)
  bic_vec[p] <- BIC(fit_p)
}

ic_df <- data.frame(Order = 1:max_order, AIC = round(aic_vec, 2), BIC = round(bic_vec, 2))
write.csv(ic_df, "output/tables/q3_37_aic_bic.csv", row.names = FALSE)
print(xtable(ic_df,
             caption = "AIC and BIC for ARX models of order $p = 1, \\ldots, 10$ (3.7).",
             label   = "tab:ic", digits = c(0, 0, 2, 2)),
      booktabs = TRUE, include.rownames = FALSE,
      sanitize.text.function = identity,
      file = "report/tables/q3_37_ic.tex", caption.placement = "top")

pdf("report/figures/q3_37_aic_bic.pdf", width = FW, height = FH_S)
par(mar = c(4.5, 4.5, 2, 1))
ylims <- range(c(aic_vec, bic_vec))
plot(1:max_order, aic_vec, type = "b", lwd = LWD, pch = 16, col = COL_AIC,
     ylim = ylims, xlab = "Model order p",
     ylab = "Information criterion",
     main = "AIC and BIC vs. ARX model order")
lines(1:max_order, bic_vec, type = "b", lwd = LWD, pch = 16, col = COL_BIC, lty = 2)
abline(v = which.min(aic_vec), col = COL_AIC, lty = 3, lwd = 1.5)
abline(v = which.min(bic_vec), col = COL_BIC,  lty = 3, lwd = 1.5)
legend("topright", c("AIC", "BIC"), col = c(COL_AIC, COL_BIC),
       lty = 1:2, lwd = LWD, pch = 16, bty = "n")
grid()
dev.off()

cat("AIC minimum at order:", which.min(aic_vec), "\n")
cat("BIC minimum at order:", which.min(bic_vec), "\n")

md_add("## 3.7 Model Order Selection (AIC / BIC)")
md_add("**Figure:** `report/figures/q3_37_aic_bic.pdf`")
md_add("**Table:** `report/tables/q3_37_ic.tex`")
md_add("**Best order by AIC:** ", which.min(aic_vec), " (AIC = ", round(min(aic_vec),2), ")")
md_add("**Best order by BIC:** ", which.min(bic_vec), " (BIC = ", round(min(bic_vec),2), ")")
md_add("**AIC values (p=1..10):** ", paste(round(aic_vec,1), collapse=", "))
md_add("**BIC values (p=1..10):** ", paste(round(bic_vec,1), collapse=", "), "\n")

# ---------------------------------------------------------------------------
# 3.8  One-step RMSE on test set vs. model order
# ---------------------------------------------------------------------------
rmse_vec <- numeric(max_order)
for (p in seq_len(max_order)) {
  regs      <- intersect(c(paste0("Ph.l", 1:p),
                           paste0("Tdelta.l", 0:(p-1)),
                           paste0("Gv.l",     0:(p-1))), names(box))
  fit_p     <- lm(as.formula(paste("Ph ~", paste(regs, collapse = " + "))), data = train)
  pred_test <- predict(fit_p, newdata = test)
  resid_test  <- test$Ph - pred_test
  stopifnot(length(resid_test) == 64)            # guard: assignment specifies 1/64
  rmse_vec[p] <- sqrt(sum(resid_test^2) / 64)
}

rmse_df <- data.frame(Order = 1:max_order, RMSE = round(rmse_vec, 4))
write.csv(rmse_df, "output/tables/q3_38_rmse.csv", row.names = FALSE)
print(xtable(rmse_df,
             caption = "One-step RMSE (W) on the test set for ARX models of order $p = 1, \\ldots, 10$ (3.8).",
             label   = "tab:rmse", digits = c(0, 0, 4)),
      booktabs = TRUE, include.rownames = FALSE,
      sanitize.text.function = identity,
      file = "report/tables/q3_38_rmse.tex", caption.placement = "top")

pdf("report/figures/q3_38_rmse.pdf", width = FW, height = FH_S)
par(mar = c(4.5, 4.5, 2, 1))
plot(1:max_order, rmse_vec, type = "b", lwd = LWD, pch = 16, col = COL_RMSE,
     xlab = "Model order p", ylab = expression("RMSE  "(W)),
     main = "One-step test-set RMSE vs. ARX model order")
abline(v = which.min(rmse_vec), col = COL_RMSE, lty = 3, lwd = 1.5)
legend("topright", paste0("Min at p = ", which.min(rmse_vec)),
       col = COL_RMSE, lty = 3, lwd = 1.5, bty = "n")
grid()
dev.off()

cat("RMSE minimum at order:", which.min(rmse_vec), "\n")

md_add("## 3.8 Test-Set RMSE vs. Model Order")
md_add("**Figure:** `report/figures/q3_38_rmse.pdf`")
md_add("**Table:** `report/tables/q3_38_rmse.tex`")
md_add("**Best order by RMSE:** ", which.min(rmse_vec), " (RMSE = ", round(min(rmse_vec),4), " W)")
md_add("**RMSE values (p=1..10):** ", paste(round(rmse_vec,3), collapse=", "), "\n")

# ---------------------------------------------------------------------------
# 3.9  Multi-step simulation over the full period
# ---------------------------------------------------------------------------
# p_best defaults to BIC-selected order — change manually if 3.7/3.8 suggest otherwise
p_best <- which.min(bic_vec)
cat("\nUsing ARX(", p_best, ") for multi-step simulation.\n")

regs_best <- intersect(c(paste0("Ph.l", 1:p_best),
                         paste0("Tdelta.l", 0:(p_best-1)),
                         paste0("Gv.l",     0:(p_best-1))), names(box))
fit_best  <- lm(as.formula(paste("Ph ~", paste(regs_best, collapse = " + "))), data = train)
saveRDS(fit_best, paste0("output/models/q3_arx_order", p_best, ".rds"))

save_coef_tex(fit_best,
              file    = paste0("report/tables/q3_39_arx", p_best, "_coef.tex"),
              caption = paste0("Coefficient estimates for the chosen ARX(", p_best, ") model."),
              label   = paste0("tab:arx", p_best, "_coef"))

# Iterative simulation: AR lags from simulated Ph, input lags from observed data
coefs      <- coef(fit_best)
ar_lags    <- grep("^Ph\\.l",      regs_best, value = TRUE)
input_lags <- grep("^(Tdelta|Gv)", regs_best, value = TRUE)
Ph_sim     <- box$Ph

for (t in seq(p_best + 1, nrow(box))) {
  ar_part    <- sum(coefs[ar_lags] * sapply(ar_lags, function(nm)
                  Ph_sim[t - as.integer(sub("Ph.l", "", nm))]))
  input_part <- sum(coefs[input_lags] * as.numeric(box[t, input_lags]))
  Ph_sim[t]  <- coefs["(Intercept)"] + ar_part + input_part
}

sim_df <- data.frame(tdate = box$tdate, Ph_obs = box$Ph, Ph_sim = Ph_sim)
write.csv(sim_df, "output/tables/q3_39_simulation.csv", row.names = FALSE)

split_t <- as.numeric(max(train$tdate))

pdf("report/figures/q3_39_multistep_simulation.pdf", width = FW, height = FH_S)
par(mar = c(4.5, 4.5, 2, 1))
plot(sim_df$tdate, sim_df$Ph_obs, type = "l", lwd = LWD, col = COL_OBS,
     ylab = expression(P[h]~(W)), xlab = "Time",
     main = bquote("Multi-step simulation  \u2014  ARX(" * .(p_best) * ")"))
lines(sim_df$tdate, sim_df$Ph_sim, lwd = LWD, col = COL_FIT, lty = 2)
abline(v = split_t, col = "gray40", lty = 3, lwd = 1.5)
legend("topright",
       c("Observed", "Simulated", "Train / test split"),
       col = c(COL_OBS, COL_FIT, "gray40"),
       lty = c(1, 2, 3), lwd = c(LWD, LWD, 1.5), bty = "n")
grid()
dev.off()

# ---------------------------------------------------------------------------
# 3.10  Print summary to console + finish markdown
# ---------------------------------------------------------------------------
sim_rmse_train <- sqrt(mean((sim_df$Ph_obs[1:167] -
                              sim_df$Ph_sim[1:167])^2, na.rm = TRUE))
sim_rmse_test  <- sqrt(mean((sim_df$Ph_obs[168:nrow(box)] -
                              sim_df$Ph_sim[168:nrow(box)])^2, na.rm = TRUE))

cat("\n=== 3.10 Summary ===\n")
cat("Best order — AIC:", which.min(aic_vec),
    "| BIC:", which.min(bic_vec),
    "| RMSE:", which.min(rmse_vec), "\n")
cat("Model used for simulation: ARX(", p_best, ")\n\n")
print(summary(fit_best)$coefficients)

md_add("## 3.9 Multi-Step Simulation")
md_add("**Figure:** `report/figures/q3_39_multistep_simulation.pdf`")
md_add("**Table (coefficients):** `report/tables/q3_39_arx", p_best, "_coef.tex`")
md_add("**Model used:** ARX(", p_best, ")")
md_add("**Simulation RMSE — training period:** ", round(sim_rmse_train, 3), " W")
md_add("**Simulation RMSE — test period:** ",     round(sim_rmse_test,  3), " W")
md_add("**Coefficients:**")
for (nm in names(coef(fit_best))) md_add("  - ", nm, ": ", round(coef(fit_best)[nm], 4))

md_add("\n## 3.10 Summary")
md_add("| Criterion | Best order |")
md_add("|-----------|------------|")
md_add("| AIC       | ", which.min(aic_vec), " |")
md_add("| BIC       | ", which.min(bic_vec), " |")
md_add("| Test RMSE | ", which.min(rmse_vec), " |")
md_add("| Chosen    | ", p_best, " (BIC) |\n")

writeLines(md, "report/q3_results_summary.md")

cat("\nAll PDFs  → report/figures/\n")
cat("All .tex  → report/tables/\n")
cat("Summary   → report/q3_results_summary.md\n")
