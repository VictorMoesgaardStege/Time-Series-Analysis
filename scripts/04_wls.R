source("R/read_data.R")

# ---------- Load data ----------
dat <- load_bil54(file = file.path("data", "/DST_BIL54.csv"))
Dtrain <- dat$train
Dtest  <- dat$test


# Q1
N = nrow(Dtrain)
lambda = 0.9
SIGMA <- diag(N)
for (i in 1:N) {
  SIGMA[i,i] <- 1/lambda^(N-i)
}

SIGMA_inv = solve(SIGMA)

#Q2
weights <- diag(SIGMA_inv)
t <- 1:N

plot(t, weights, type="l", lwd=2, col="darkblue",
     xlab="Time", ylab="Weight (λ^k)",
     main=paste("WLS Weights vs Time (λ =", lambda, ")"))
points(t, weights, pch=16, col="darkblue")
grid()

#Q3
sum(weights)


#Q4