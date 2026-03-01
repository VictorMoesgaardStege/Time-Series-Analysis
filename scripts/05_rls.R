# scripts/05_rls.R
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

options(scipen = 999)
##### 5.1 - 5.2


# RLS: theta_t = theta_{t-1} + R_t^(-1) * X_t * [Y_t - X_t^T * theta_{t-1}]
# RLS: R_t = R_{t-1} + X_t * X_t^T

rls_estimation <- function(x_data, y_data, R0, theta0) {
  n <- length(x_data)
  
  # Initialize
  R <- R0
  theta <- theta0
  
  # Storage
  theta_history <- matrix(NA, nrow = n, ncol = 2)
  colnames(theta_history) <- c("theta1_intercept", "theta2_slope")
  
  predictions <- numeric(n)      # One-step ahead predictions
  fitted_values <- numeric(n)    # Fitted values (after update)
  residuals <- numeric(n)        # One-step ahead residuals
  
  # RLS Loop
  for (t in 1:n) {
    
    # Design vector for current observation
    X_t <- matrix(c(1, x_data[t]), nrow = 2)
    Y_t <- y_data[t]
    
    # One-step ahead prediction using theta from PREVIOUS time step
    # y_hat_{t|t-1} = X_t^T * theta_{t-1}
    predictions[t] <- as.numeric(t(X_t) %*% theta)
    
    # One-step ahead residual
    residuals[t] <- Y_t - predictions[t]
    
    # Update R
    R <- R + X_t %*% t(X_t)
    
    # Update theta
    theta <- theta + solve(R) %*% X_t %*% residuals[t]
    
    # Store theta for this time step
    theta_history[t, ] <- as.numeric(theta)
    
    # Fitted value (after update)
    fitted_values[t] <- as.numeric(t(X_t) %*% theta)
  }
  
  # Return results
  return(list(
    theta_final = as.numeric(theta),
    theta_history = theta_history,
    R_final = R,
    predictions = predictions,
    fitted_values = fitted_values,
    residuals = residuals
  ))
}
# Run RLS for first 3 observations
results_5_2 <- rls_estimation(
  x_data = x_train[1:3],
  y_data = y_train[1:3],
  R0 = diag(0.1, 2),
  theta0 = c(0, 0)
)

print("Theta at t=3:")
print(results_5_2$theta_final)

print("Theta history:")
print(results_5_2$theta_history)


print("Predictions for t1, t2 and t3:")
print(as.numeric(results_5_2$theta_history[1,1] + 2018 * results_5_2$theta_history[1,2]))
print(as.numeric(results_5_2$theta_history[2,1] + (2018+(1/12)) * results_5_2$theta_history[2,2]))
print(as.numeric(results_5_2$theta_history[3,1] + (2018+(2/12)) * results_5_2$theta_history[3,2]))


##### 5.3

# OLS estimate
model_ols <- lm(total ~ year, data = Dtrain)

# RLS with given R0
results_5_3 <- rls_estimation(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.1, 2),  # Small R0 for convergence
  theta0 = c(0, 0)
)

# Compare
print("OLS coefficients:")
print(coef(model_ols))

print("RLS theta_N:")
print(results_5_3$theta_final)

print("Difference:")
print(results_5_3$theta_final - coef(model_ols))

# RLS with much higher uncertainty
results_5_3 <- rls_estimation(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.00000001, 2),  # Small R0 for convergence
  theta0 = c(0, 0)
)

# Compare
print("OLS coefficients:")
print(coef(model_ols))

print("RLS theta_N:")
print(results_5_3$theta_final)

print("Difference:")
print(results_5_3$theta_final - coef(model_ols))

##### 5.4
rls_estimation_forgetting <- function(x_data, y_data, R0, theta0, lambda) {
  n <- length(x_data)
  
  # Initialize
  R <- R0
  theta <- theta0
  
  # Storage
  theta_history <- matrix(NA, nrow = n, ncol = 2)
  colnames(theta_history) <- c("theta1_intercept", "theta2_slope")
  
  predictions <- numeric(n)      # One-step ahead predictions
  fitted_values <- numeric(n)    # Fitted values (after update)
  residuals <- numeric(n)        # One-step ahead residuals
  
  # RLS Loop
  for (t in 1:n) {
    
    # Design vector for current observation
    X_t <- matrix(c(1, x_data[t]), nrow = 2)
    Y_t <- y_data[t]
    
    # One-step ahead prediction using theta from PREVIOUS time step
    # y_hat_{t|t-1} = X_t^T * theta_{t-1}
    predictions[t] <- as.numeric(t(X_t) %*% theta)
    
    # One-step ahead residual
    residuals[t] <- Y_t - predictions[t]
    
    # Update R
    R <- lambda * R + X_t %*% t(X_t)
    
    # Update theta
    theta <- theta + solve(R) %*% X_t %*% residuals[t]
    
    # Store theta for this time step
    theta_history[t, ] <- as.numeric(theta)
    
    # Fitted value (after update)
    fitted_values[t] <- as.numeric(t(X_t) %*% theta)
  }
  
  # Return results
  return(list(
    theta_final = as.numeric(theta),
    theta_history = theta_history,
    R_final = R,
    predictions = predictions,
    fitted_values = fitted_values,
    residuals = residuals
  ))
}

# Run RLS with forgetting for lambda = 0.7 and lambda = 0.99
results_5_4_07 <- rls_estimation_forgetting(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.1, 2),
  theta0 = c(0, 0),
  lambda = 0.7
)

results_5_4_099 <- rls_estimation_forgetting(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.1, 2),
  theta0 = c(0, 0),
  lambda = 0.99
)
# Remove burn-in period
burnin <- 5
t_plot <- burnin:length(x_train)

# Create 2 plots (one for each parameter)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Plot theta_1 (intercept)
plot(t_plot, results_5_4_07$theta_history[t_plot, 1], 
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "theta_1 (intercept)",
     main = "Evolution of theta_1 (intercept)")
lines(t_plot, results_5_4_099$theta_history[t_plot, 1], 
      col = "red", lwd = 2)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99"), 
       col = c("blue", "red"), lwd = 2)
grid()

# Plot theta_2 (slope)
plot(t_plot, results_5_4_07$theta_history[t_plot, 2], 
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "theta_2 (slope)",
     main = "Evolution of theta_2 (slope)")
lines(t_plot, results_5_4_099$theta_history[t_plot, 2], 
      col = "red", lwd = 2)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99"), 
       col = c("blue", "red"), lwd = 2)
grid()



# Run RLS with forgetting for lambda = 0.7 and lambda = 0.99
results_5_4_07 <- rls_estimation_forgetting(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.00000001, 2),
  theta0 = c(0, 0),
  lambda = 0.7
)

results_5_4_099 <- rls_estimation_forgetting(
  x_data = x_train,
  y_data = y_train,
  R0 = diag(0.00000001, 2),
  theta0 = c(0, 0),
  lambda = 0.99
)
# Remove burn-in period
burnin <- 5
t_plot <- burnin:length(x_train)

# Create 2 plots (one for each parameter)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Plot theta_1 (intercept)
plot(t_plot, results_5_4_07$theta_history[t_plot, 1], 
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "theta_1 (intercept)",
     main = "Evolution of theta_1 (intercept)")
lines(t_plot, results_5_4_099$theta_history[t_plot, 1], 
      col = "red", lwd = 2)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99"), 
       col = c("blue", "red"), lwd = 2)
grid()

# Plot theta_2 (slope)
plot(t_plot, results_5_4_07$theta_history[t_plot, 2], 
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "theta_2 (slope)",
     main = "Evolution of theta_2 (slope)")
lines(t_plot, results_5_4_099$theta_history[t_plot, 2], 
      col = "red", lwd = 2)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99"), 
       col = c("blue", "red"), lwd = 2)
grid()

print(results_5_4_07$theta_final)
print(results_5_4_099$theta_final)

##### Q 5.5
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
burnin <- 4
pred_plot <- burnin:length(x_train)
residuals_07 <- results_5_4_07$predictions[pred_plot]-y_train[pred_plot]
residuals_099 <- results_5_4_099$predictions[pred_plot]-y_train[pred_plot]

plot(pred_plot, results_5_4_099$predictions[pred_plot], 
     type = "l", col = "red", lwd = 2,
     xlab = "Time", ylab = "prediction in millions",
     main = "One-step prediction of both lambda values")
lines(pred_plot, results_5_4_07$predictions[pred_plot],
      col = "blue", lwd = 2)
lines(pred_plot, y_train[pred_plot],
      col = "black", lwd = 1)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99", "y_t"), 
       col = c("blue", "red", "black"), lwd = 2)
grid()

plot(pred_plot, residuals_099, 
     type = "l", col = "red", lwd = 2,
     xlab = "Time", ylab = "residual in millions",
     main = "Residuals of both lambda values")
lines(pred_plot, residuals_07, 
      col = "blue", lwd = 2)
legend("bottomright", 
       legend = c("lambda = 0.7", "lambda = 0.99"), 
       col = c("blue", "red"), lwd = 2)
grid()
abline(h = 0, lty = 2, col = "black", lwd = 1)

##### Q 5.6
rls_estimation_forgetting_k_horizon <- function(x_data, y_data, R0, theta0, lambda, max_horizon) {
  n <- length(x_data)
  
  R <- R0
  theta <- theta0
  
  theta_history <- matrix(NA, nrow = n, ncol = 2)
  colnames(theta_history) <- c("theta1_intercept", "theta2_slope")
  
  predictions <- numeric(n)
  fitted_values <- numeric(n)
  residuals <- numeric(n)
  
  horizon_residuals <- matrix(NA, nrow = n, ncol = max_horizon)
  
  for (t in 1:n) {
    X_t <- matrix(c(1, x_data[t]), nrow = 2)
    Y_t <- y_data[t]
    
    predictions[t] <- as.numeric(t(X_t) %*% theta)
    residuals[t] <- Y_t - predictions[t]
    
    R <- lambda * R + X_t %*% t(X_t)
    theta <- theta + solve(R) %*% X_t %*% residuals[t]
    
    theta_history[t, ] <- as.numeric(theta)
    fitted_values[t] <- as.numeric(t(X_t) %*% theta)
    
    for (h in 1:max_horizon) {
      future_t <- t + h
      if (future_t > n) break
      
      X_future <- matrix(c(1, x_data[future_t]), nrow = 2)
      y_hat <- as.numeric(t(X_future) %*% theta)
      horizon_residuals[t, h] <- y_data[future_t] - y_hat
    }
  }
  
  rmse_k <- numeric(max_horizon)
  for (h in 1:max_horizon) {
    burn_in <- 4
    valid_residuals <- horizon_residuals[(burn_in + 1):(n - h), h]
    rmse_k[h] <- sqrt(mean(valid_residuals^2, na.rm = TRUE))
  }
  
  return(list(
    theta_final = as.numeric(theta),
    theta_history = theta_history,
    R_final = R,
    predictions = predictions,
    fitted_values = fitted_values,
    residuals = residuals,
    horizon_residuals = horizon_residuals,
    rmse_k = rmse_k
  ))
}

# Loop - save into matrix
lambdas <- seq(0.5, 0.99, by = 0.01)
rmse_results <- matrix(NA, nrow = 12, ncol = length(lambdas))

for (i in seq_along(lambdas)) {
  result <- rls_estimation_forgetting_k_horizon(
    x_data = x_train,
    y_data = y_train,
    R0 = diag(0.00000001, 2),
    theta0 = c(0, 0),
    lambda = lambdas[i],
    max_horizon = 12
  )
  rmse_results[, i] <- result$rmse_k
}

library(ggplot2)
library(reshape2)

# Convert matrix to long format
colnames(rmse_results) <- round(lambdas, 2)
rownames(rmse_results) <- 1:12

rmse_long <- melt(rmse_results, varnames = c("horizon", "lambda"), value.name = "rmse")

ggplot(rmse_long, aes(x = lambda, y = rmse, group = horizon, color = factor(horizon))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = colorRampPalette(c("darkblue", "lightblue"))(12),
    name = "Horizon (k)"
  ) +
  labs(
    title = "RMSE by lambda and forecast horizon",
    x = expression(lambda),
    y = expression(RMSE[k])
  ) +
  theme_bw()
lambdas <- seq(0.5, 0.99, by = 0.01)

for (k in 1:12) {
  valid <- which(lambdas <= 0.9)
  min_idx <- valid[which.min(rmse_results[k, valid])]
  cat("k =", k, "| optimal lambda =", lambdas[min_idx], "| RMSE =", round(rmse_results[k, min_idx], 6), "\n")
}


##### Q5.7
optimal_lambdas <- c(0.50, 0.50, 0.50, 0.55, 0.61, 0.64, 0.66, 0.67, 0.68, 0.68, 0.66, 0.65)

predictions_test <- numeric(12)

for (k in 1:12) {
  result <- rls_estimation_forgetting_k_horizon(
    x_data = x_train,
    y_data = y_train,
    R0 = diag(0.00000001, 2),
    theta0 = c(0, 0),
    lambda = optimal_lambdas[k],
    max_horizon = k  # only need up to horizon k
  )
  
  # theta at t=N is the last row of theta_history
  theta_N <- result$theta_history[nrow(result$theta_history), ]
  predictions_test[k] <- theta_N[1] + theta_N[2] * x_test[k]
}

ols_pred <- read.csv("output/tables/03_ols_forecast_2024.csv")
wls_pred <- read.csv("output/tables/04_wls_forecast_2024_lambda_0.9.csv")

str(ols_pred)
str(wls_pred)

# Plot
par(mar = c(7, 4, 4, 2))
plot(x_test, y_test,
     type = "l", lwd = 2, col = "black",
     xlab = "Time", ylab = "Total vehicles",
     main = "5.7 - Test set predictions: OLS vs WLS vs RLS",
     ylim = range(c(y_test, predictions_test, ols_pred$predicted, wls_pred$wls_pred)))
points(x_test, predictions_test, col = "purple", pch = 17)
points(x_test, ols_pred$predicted, col = "blue", pch = 15)
points(x_test, wls_pred$wls_pred, col = "red", pch = 16)
par(xpd = TRUE)
legend("bottom", inset = c(0, -0.35),
       legend = c("Observed", "RLS optimal λ per k", "OLS", "WLS (λ=0.9)"),
       col = c("black", "purple", "blue", "red"),
       pch = c(NA, 17, 15, 16),
       lwd = c(2, NA, NA, NA),
       horiz = TRUE, bty = "n", cex = 0.85)
par(xpd = FALSE)
grid()
