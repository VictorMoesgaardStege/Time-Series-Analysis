# Read the data
df <- read.csv("box_data_60min.csv", stringsAsFactors = FALSE)

# Convert timestamp column to datetime
df$tdate <- as.POSIXct(df$tdate, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Create hourly time series from Ph
ph_ts <- ts(df$Ph, frequency = 24)

# Plot ACF and PACF
par(mfrow = c(1, 2))
acf(ph_ts, lag.max = 30, main = "ACF of Ph")
pacf(ph_ts, lag.max = 30, main = "PACF of Ph")
par(mfrow = c(1, 1))