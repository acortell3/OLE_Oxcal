####### SCRIPT TO SIMULATE DATES USING A LOGISTIC FUNCTION

set.seed(123)

## Create vector and logistic function
tr <- seq(6500,5500) ## Time vector

## Logistic growth
K <- 1000 ## Carrying capacity
r <- 0.01 ## Growth rate
ti <- seq(1,length(tr)) ## Time

pop <- K / (1 + exp(-r*(ti-(max(ti)/2))))

prob_pop <- pop / max(pop)

## Sample dates
ndates <- 80 ## I do 80 dates because OLE function doesn't work with more than 97 dates
dates <- sample(tr,ndates,prob_pop,replace = TRUE)

####### PLOT THE THEORETICAL RESULTS #######

# Plot the theoretical logistic curve (probabilities)
plot(tr, prob_pop, type = "l", lwd = 2, col = "blue",
     ylab = "Probability / Frequency", xlab = "Time",
     main = "Logistic Sampling of Dates")

# Add histogram of sampled dates (scaled to match probabilities)
hist(dates, breaks = 30, freq = FALSE, add = TRUE,
     col = rgb(1, 0, 0, 0.3), border = NA)

# Add a legend
legend("topleft", legend = c("Logistic probability", "Sampled dates"),
       col = c("blue", rgb(1, 0, 0, 0.3)), lwd = c(2, 10),
       pch = c(NA, 15), pt.cex = 2, bty = "n")

####### SPD AND CALIBRATION
library(rcarbon)

## Produce errors
er <- sample(seq(20,80,20),ndates,replace = TRUE)

## Uncalibrate
uncal_dates <- uncalibrate(dates)

## Back-calibrate
cal_dates <- calibrate(uncal_dates$ccCRA, errors = er)
spdates <- spd(cal_dates, timeRange = c(6700,5600))
plot(spdates)

## Save info and dates for pop models
YearsBP <- data.frame("Year" = uncal_dates$ccCRA,
		      "Sd" = er)

saveRDS(YearsBP,"../Simu_data/Uncal_YearsBP_100dates.rds")

######### APPROACH WITH DTRAPEZOID
library(nimbleCarbon)

modelPlot(dLogisticGrowth,a=6500,b=5000,params=c(k=0.1,r=0.01),alpha=1,col='black') #the sampling distribution where k is the proportion of K (carrying capacity).
x <- replicate(80,rLogisticGrowth(a=6500,b=5000,k=0.1,r=0.01)) #sample 80 dates from the distribution above in calendar time (cannot use n=80 as nimble is unhappy with that)
dates <- uncalibrate(x,CRAerrors=20)[,c('rCRA','rError')]

saveRDS(dates,"../Simu_data/Uncal_YearsBP_100dates.rds") ## It purposedly overwrites the above





