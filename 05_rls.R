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
