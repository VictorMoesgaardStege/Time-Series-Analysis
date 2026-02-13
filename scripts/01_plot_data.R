# scripts/01_plot_data.R
# Section 1: Plot data (training set only)

source("R/read_data.R")

# Load data
dat <- load_bil54(file = file.path("data", "DST_BIL54.csv"))
Dtrain <- dat$train
Dtest  <- dat$test

# --- 1.1 Make time variable x ---
# Requirement: 2018-Jan -> 2018, 2018-Feb -> 2018 + 1/12, etc.
# We compute x from the POSIXct time column
train_year  <- as.POSIXlt(Dtrain$time)$year + 1900
train_month <- as.POSIXlt(Dtrain$time)$mon  # 0..11
x <- train_year + train_month / 12

# Sanity check: first few x values should be 2018, 2018+1/12, 2018+2/12, ...
print(head(x, 6))

# --- Plot: training data vs x ---
# Create output folder for figures (optional but convenient for the report)
if (!dir.exists(file.path("report", "figures"))) {
  dir.create(file.path("report", "figures"), recursive = TRUE)
}

png(filename = file.path("report", "figures", "01_training_vs_x.png"),
    width = 1200, height = 700, res = 150)

plot(
  x, Dtrain$total,
  pch = 16,
  xlab = "x (year + month/12)",
  ylab = "Total vehicles (millions)",
  main = "Training data: total vehicles vs x"
)
lines(x, Dtrain$total)

dev.off()

# Also show the plot in the interactive device
plot(
  x, Dtrain$total,
  pch = 16,
  xlab = "x (year + month/12)",
  ylab = "Total vehicles (millions)",
  main = "Training data: total vehicles vs x"
)
lines(x, Dtrain$total)

# --- Optional extras for writing 1.2 (describe the time series) ---

# A) Plot against calendar time (sometimes easier to interpret than x)
png(filename = file.path("report", "figures", "01_training_vs_time.png"),
    width = 1200, height = 700, res = 150)

plot(
  Dtrain$time, Dtrain$total,
  type = "l",
  xlab = "Time",
  ylab = "Total vehicles (millions)",
  main = "Training data: total vehicles over time"
)
points(Dtrain$time, Dtrain$total, pch = 16)

dev.off()

# B) Monthly increments (to see growth/seasonality more clearly)
# Note: diff() loses first observation
delta_total <- diff(Dtrain$total)

png(filename = file.path("report", "figures", "01_training_monthly_change.png"),
    width = 1200, height = 700, res = 150)

plot(
  Dtrain$time[-1], delta_total,
  type = "h",
  xlab = "Time",
  ylab = "Monthly change in total (millions)",
  main = "Training data: month-to-month change"
)
abline(h = 0)

dev.off()

# Console summary (useful for the description in 1.2)
cat("\nTraining period:\n")
cat("Start:", format(min(Dtrain$time)), "\n")
cat("End  :", format(max(Dtrain$time)), "\n")
cat("N    :", nrow(Dtrain), "\n")

cat("\nTotal vehicles (millions) summary (training):\n")
print(summary(Dtrain$total))
