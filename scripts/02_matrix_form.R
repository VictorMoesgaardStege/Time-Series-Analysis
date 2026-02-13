# scripts/02_matrix_form.R
source("R/read_data.R")

dat <- load_bil54(file = file.path("data", "DST_BIL54.csv"))
Dtrain <- dat$train

head(Dtrain, 3)

# First 3 time points
D3 <- Dtrain[1:3, ]

# Define x_t exactly like in section 1 (year + month/12)
yr  <- as.POSIXlt(D3$time)$year + 1900
mon <- as.POSIXlt(D3$time)$mon
x3  <- yr + mon / 12

# Output vector y (3x1)
y3 <- as.matrix(D3$total)

# Design matrix X (3x2): [1, x_t]
X3 <- cbind(1, x3)

# Print "actual values" with max 3 digits (decimals)
cat("\nFirst 3 time points (raw rows):\n")
print(D3[, c("time", "total")])

cat("\nVector y (3x1), rounded to 3 decimals:\n")
print(round(y3, 3))

cat("\nDesign matrix X (3x2) = [1, x_t], rounded to 3 decimals:\n")
print(round(X3, 3))

# Optional: label rows/cols nicely for screenshots
colnames(X3) <- c("1", "x_t")
rownames(X3) <- paste0("t=", 1:3)
rownames(y3) <- paste0("t=", 1:3)

cat("\nLabeled (nice for report screenshot):\n")
cat("\ny:\n"); print(round(y3, 3))
cat("\nX:\n"); print(round(X3, 3))
