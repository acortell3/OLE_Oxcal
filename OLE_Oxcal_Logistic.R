

################################################################################
################################################################################

### Comparison of Oxcal/OLE using a logistic model

################################################################################
################################################################################


## Load packages
library(rcarbon)
library(nimbleCarbon)

## Functions
######## OLE calculation function taken from Key et al. 2024 since sExtinct has been retracted from CRAN
OLE.test <- function(dd, alpha){ ## dd is the dates and alpha is the confidence interval
  # records are sorted in a reverse order, as required by OLE method
  sights <- rev(sort(dd))
  # calculation of k, v, e, lambda and other values
  k <- length(sights)
  v <- (1/(k-1)) * sum(log((sights[1] - sights[k])/(sights[1] - sights[2:(k-1)])))
  e <- matrix(rep(1,k), ncol=1)
  SU<-(-log(alpha)/length(sights))^-v
  myfun <- function(i,j,v){(gamma(2*v+i)*gamma(v+j))/(gamma(v+i)*gamma(j))}
  lambda <- outer(1:k, 1:k, myfun, v=v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda)) 
  a <- as.vector(solve(t(e)%*%solve(lambda)%*%e)) * solve(lambda)%*%e
  # calculation of CI ("upperCI") and extinction time ("extest")
  upperCI<-max(sights) + ((max(sights)-min(sights))/(SU-1))
  extest<-sum(t(a)%*%sights)
  # return of results produced by the function
  res<-data.frame(Estimate=extest, upperCI=upperCI)
  return(res)	
}


################################################################################
######################     STEP 1. GENERATE THE MODEL     ######################
################################################################################

## Set parameters (these will be randomised and simulated n number of times)
r <- 0.01 ## Rate
s <- 6500 ## Starting date (range for simulation will be 10000-3000)
e <- s-500 ## End date. I can also randomise this with range 300-1000
n <- 10 ## Number of dates (range for simulation will be 15-50)
k <- 0.00001 ## Initial proportion of carrying capacity
alpha <- 0.05 ## Significance OLE

## Generate the model 
dates <- replicate(n,rLogisticGrowth(a = s, b = e, k = k, r = r)) ## Logistic model
errors <- sample(seq(20,100,20), n, replace = TRUE)

## Uncalibrate
uncal_dates <- uncalibrate(dates)

## Back-calibrate
cal_dates <- calibrate(uncal_dates$ccCRA, errors = errors)

## SPD
spd_dates <- spd(cal_dates, timeRange = c(s,e))

par(mfrow = c(2,1))
hist(dates, breaks = 20)
plot(spd_dates)

################################################################################
######################         STEP 2. OLE median         ######################
################################################################################
med_dates <- rep(NA,n)

## Extract medians per date
for (i in 1:n){
  med_dates[i] <- round(median(cal_dates[i]$grids[[1]][,1]))
}  

## OLE on medians
OLE_med <- OLE.test(med_dates, alpha = alpha)
OLE_med


################################################################################
######################        STEP 3. OLE resampling      ######################
################################################################################
res_dates <- rep(NA,n)

## Sample date based on date probability density
for (i in 1:n){
  res_dates[i] <- sample(cal_dates[i]$grids[[1]][,1], 1, prob = cal_dates[i]$grids[[1]][,2])
}

## OLE on resampled
OLE_res <- OLE.test(res_dates, alpha = alpha)
OLE_res


################################################################################
######################           STEP 4. OxCal            ######################
################################################################################







