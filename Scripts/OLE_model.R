
########## SCRIPT TO DO AN OPTIMAL LINEAR ESTIMATION, BASED ON THE SIMULATED DATES FROM Simu_dates.R

## Load necessary packages
library(rcarbon)
library(sExtinct)
library(weibulltools)

## Load dates
dates <- readRDS("../Simu_data/Uncal_YearsBP_100dates.rds")
colnames(dates) <- c("Year", "Sd") 


## Since the dates have been uncalibrated in the simulation, I need to calibrate them back here
cal_dates <- calibrate(dates$Year, dates$Sd)
spdates <- spd(cal_dates, timeRange = c(6700,5600))
plot(spdates) ## Everything is cool

### Apply OLE
source("functions.R")

### First, OLE for medians
## Function for extracting the median of a calibrated date
dates_median <- rep(NA,nrow(dates))

for (i in 1:nrow(dates)){
	dates_median[i] <- medCal(calibrate(dates[i,1],dates[i,2]))
}

OLE_medians <- OLE.test(dd = dates_median, alpha = 0.05)
write.csv(OLE_medians,"../Results/OLE_medians.csv")


### Second, OLE with resampling
resamp <- 80 ## Number of resamples
OLE_resamp <- data.frame("Estimate" = rep(NA,resamp),
			 "upperCI" = rep(NA,resamp),
			 "lowerCI" = rep(NA,resamp),
			 "WeibullShape" = rep(NA,resamp))

dates_for_OLE <- matrix(nrow = 80, ncol = 80)

for (j in 1:resamp){
	date_vec <- rep(NA,nrow(dates))
	
	for (i in 1:nrow(dates)){
		date_vec[i] <- sample(cal_dates$grids[[i]][,1],1,cal_dates$grids[[i]][,2], replace = T)
	}
	
	dates_for_OLE[j,] <- date_vec
	OLE_res <- OLE.test.v(dd = date_vec,alpha = 0.05)

	OLE_resamp[j,] <- OLE_res

}

W_scale <- function(x,k){
	l <- x/(gamma(1+(1/k)))
	return(l)
}

# With medians
W_scale <- function(x,k){
	l <- (log(2)^(1/k)) / x
	return(l)
}

OLE_resamp$WeibullScale <- W_scale(x = OLE_resamp$Estimate, k = OLE_resamp$WeibullShape)

#wx <- seq(0, 30000, length.out = 1000)
curve(dweibull(x, shape = OLE_resamp$WeibullShape[1], scale = OLE_resamp$WeibullScale[1]))

for (i in 2:resamp){
	curve(dweibull(x, shape = OLE_resamp$WeibullShape[i], scale = OLE_resamp$WeibullScale[i]), add = TRUE)
}

weigh_WB <- c("Shape" = mean(OLE_resamp$WeibullShape), 
	      "Scale" = mean(OLE_resamp$WeibullScale))
mean_WB <- weigh_WB[2]*gamma(1+(1/weigh_WB[1]))
names(mean_WB) <- "E(X)"
qweibull(c(0.025, 0.975), shape = weigh_WB[1], scale = weigh_WB[2]) ## Confidence intervals are huuuuuge, so this approach doesn't work

SU <- (-log(0.05)/80)^-weigh_WB[1] 
upperCI <- max(dates_for_OLE) + ((max(dates_for_OLE)-min(dates_for_OLE))/(SU-1))

lowerCI <- min(dates_for_OLE) + ((min(dates_for_OLE)-max(dates_for_OLE))/(SU-1))

mean_WB
upperCI
lowerCI

OLE_resamp_res <- c("mean" = mean_WB, "upperCI" = upperCI, "lowerCI" = lowerCI)
write.csv(OLE_resamp_res,"../Results/OLE_resamp_res.csv")
