

### Test to see of OLE works with non-cal dates

set.seed(123)

## Load libraries and external functions
library(nimbleCarbon)
source('Functions.R')

start_date <- 1000
end_date <- 1900
sample_size <- 80
k <- 10

n <- end_date - start_date

## Generate probability
x <- seq(0, 100, length.out = n)
K <- 1     # carrying capacity (max value)
r <- 0.1   # growth rate
x0 <- 50   # midpoint

y <- K / (1 + exp(-r * (x - x0)))

plot(y, type = "l")

dates <- round(sample(seq(start_date,end_date, length.out = n), sample_size, prob = y))
sorted_dates <- sort(dates, decreasing = F)
sorted_dates <- sorted_dates[1:k]

OLE_res <- OLE.test(dd = sorted_dates, alpha = 0.05)
sorted_dates[1]
OLE_res




