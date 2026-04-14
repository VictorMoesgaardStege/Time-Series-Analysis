solar_df <- read_solar_data("data/datasolar.csv")
head(solar_df)
summary(solar_df$Y)

dat <- read.csv("data/datasolar.csv", check.names = FALSE)
str(dat)
head(dat)
