# R/read_data.R

load_bil54 <- function(file = NULL, teststart = "2024-01-01", tz = "UTC") {
  # If no file is provided, assume standard repo layout: data/DST_BIL54.csv
  if (is.null(file)) {
    if (!requireNamespace("here", quietly = TRUE)) {
      stop("Package 'here' is required. Install it with install.packages('here').")
    }
    file <- here::here("data", "DST_BIL54.csv")
  }

  D <- read.csv(file)

  # Parse monthly time stamps
  D$time <- as.POSIXct(paste0(D$time, "-01"), "%Y-%m-%d", tz = tz)

  # Time covariate used in assignment: year + month/12
  D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

  # Convert to millions
  D$total <- as.numeric(D$total) / 1e6

  # Split train/test
  teststart <- as.POSIXct(teststart, tz = tz)
  Dtrain <- D[D$time < teststart, ]
  Dtest  <- D[D$time >= teststart, ]

  list(D = D, train = Dtrain, test = Dtest)
}
