

### GENERATION PROCESSES

## Load libraries
library(rcarbon)
library(nimbleCarbon) ## for dLogisticGrowth
library(mc2d) ## For the PERT distribution

set.seed(1)

## Two samples sizes. It was said 50 and 100. In reality, sample size must be computed regarding the length of the process
## since it is not the same 50 dates for 1000 than for 10000 years. This sample size is compted ratio r = s/l

### Function 1
## Function 1. Generation of population dynamics
#' @title n_dates
#' Returns the number of dates based on the ratio. If n_dates = TRUE, it returns the ratio
#' @param a: Oldest date for the timeframe considered
#' @param b: Youngest date for the timeframe considered
#' @param r: If n_dates = FALSE (default) this must be specified and the function will return
#' a number of dates based on the ratio of dates per total length of the process
#' @param n_dates: Default is FALSE. It can be set to numeric, in which case, the function returns the 
#' ratio of dates per length of the process. If specified it must fulfill n_dates < a-b
#' @returns: Returns the number of dates based on the ratio (length*ratio) or the ratio (n_dates/l) if the 
#' dates are specified (n_dates != FALSE)
#' @export

n_dates <- function(a, b, r, n_dates = FALSE){
  l <- a-b
  res <- ifelse(isFALSE(n_dates), l*r, n_dates/l)
  ifelse(res >= 1, names(res) <- "N_dates", names(res) <- "ratio")
  return(res)
}

## We consider the following common data
start <- 7000 ## Start of the process
end <- 4000 ## End of the process

ts <- c(7000:4000) # Time span

## Ratios
rsm <- 0.02 ## Small ratio
rmd <- 0.04 ## Medium ratio
rlg <- 0.06 ## Large ratio
elg <- 0.08 ## Extra large ratio

## Sample sizes
ss_sm <- n_dates(a = start, b = end, r = rsm) ## Small
ss_md <- n_dates(a = start, b = end, r = rmd) ## Medium
ss_lg <- n_dates(a = start, b = end, r = rlg) ## Large
ss_el <- n_dates(a = start, b = end, r = elg) ## Extra-Large


################################################################################
#####################         Exponential        ###############################
################################################################################

## Generate calendar dates

## Parameter for the exponential distribution
lambda1 <- 0.05
names(lambda1) <- "Lambda"
#lambda2 <- 1
#lambda3 <- 2

pm_exp_sm <- sample(sort(ts),ss_sm,prob=dexp(seq(0,100,length.out = length(ts)),lambda1), replace = TRUE)
pm_exp_md <- sample(sort(ts),ss_md,prob=dexp(seq(0,100,length.out = length(ts)),lambda1), replace = TRUE)
pm_exp_lg <- sample(sort(ts),ss_lg,prob=dexp(seq(0,100,length.out = length(ts)),lambda1), replace = TRUE)
pm_exp_el <- sample(sort(ts),ss_el,prob=dexp(seq(0,100,length.out = length(ts)),lambda1), replace = TRUE)

## Uncalibrate
pm_exp_sm_unc <- uncalibrate(pm_exp_sm)
pm_exp_md_unc <- uncalibrate(pm_exp_md)
pm_exp_lg_unc <- uncalibrate(pm_exp_lg)
pm_exp_el_unc <- uncalibrate(pm_exp_el)

## Back-calibrate
pm_exp_sm_cal <- calibrate(pm_exp_sm_unc$ccCRA,errors = sample(seq(20,100,20),ss_sm, replace = TRUE))
pm_exp_md_cal <- calibrate(pm_exp_md_unc$ccCRA,errors = sample(seq(20,100,20),ss_md, replace = TRUE))
pm_exp_lg_cal <- calibrate(pm_exp_lg_unc$ccCRA,errors = sample(seq(20,100,20),ss_lg, replace = TRUE))
pm_exp_el_cal <- calibrate(pm_exp_el_unc$ccCRA,errors = sample(seq(20,100,20),ss_el, replace = TRUE))

## Produce SPDs
pm_exp_sm_spd <- spd(pm_exp_sm_cal, timeRange = c(start,end))
pm_exp_md_spd <- spd(pm_exp_md_cal, timeRange = c(start,end))
pm_exp_lg_spd <- spd(pm_exp_lg_cal, timeRange = c(start,end))
pm_exp_el_spd <- spd(pm_exp_el_cal, timeRange = c(start,end))

## Plot SPDs
par(mfrow = c(2,2))
plot(pm_exp_sm_spd, main = paste0("Sample size = ", ss_sm))
plot(pm_exp_md_spd, main = paste0("Sample size = ", ss_md))
plot(pm_exp_lg_spd, main = paste0("Sample size = ", ss_lg))
plot(pm_exp_el_spd, main = paste0("Sample size = ", ss_el))

Exp_pm <- list("Back-cal.dates" = list("Small.sample.size" = pm_exp_sm_cal,
                                       "Medium.sample.size" = pm_exp_md_cal,
                                       "Large.sample.size" = pm_exp_lg_cal,
                                       "Extralarge.sample.size" = pm_exp_el_cal),
               "SPDs" = list("Small.sample.size" = pm_exp_sm_spd,
                             "Medium.sample.size" = pm_exp_md_spd,
                             "Large.sample.size" = pm_exp_lg_spd,
                             "Extralarge.sample.size" = pm_exp_el_spd),
               "Parameters" = lambda1)


################################################################################
#####################     Negative exponential   ###############################
################################################################################


## Generate calendar dates

## Parameter is the same lambda as the one for the exponential distribution. Re-parameterisation below
 
pm_nexp_sm <- sample(sort(ts, decreasing = TRUE),ss_sm,prob=dexp(seq(0,100,length.out = length(ts)),1/lambda1), replace = TRUE)
pm_nexp_md <- sample(sort(ts, decreasing = TRUE),ss_md,prob=dexp(seq(0,100,length.out = length(ts)),1/lambda1), replace = TRUE)
pm_nexp_lg <- sample(sort(ts, decreasing = TRUE),ss_lg,prob=dexp(seq(0,100,length.out = length(ts)),1/lambda1), replace = TRUE)
pm_nexp_el <- sample(sort(ts, decreasing = TRUE),ss_el,prob=dexp(seq(0,100,length.out = length(ts)),1/lambda1), replace = TRUE)

## Uncalibrate
pm_nexp_sm_unc <- uncalibrate(pm_nexp_sm)
pm_nexp_md_unc <- uncalibrate(pm_nexp_md)
pm_nexp_lg_unc <- uncalibrate(pm_nexp_lg)
pm_nexp_el_unc <- uncalibrate(pm_nexp_el)

## Back-calibrate
pm_nexp_sm_cal <- calibrate(pm_nexp_sm_unc$ccCRA,errors = sample(seq(20,100,20),ss_sm, replace = TRUE))
pm_nexp_md_cal <- calibrate(pm_nexp_md_unc$ccCRA,errors = sample(seq(20,100,20),ss_md, replace = TRUE))
pm_nexp_lg_cal <- calibrate(pm_nexp_lg_unc$ccCRA,errors = sample(seq(20,100,20),ss_lg, replace = TRUE))
pm_nexp_el_cal <- calibrate(pm_nexp_el_unc$ccCRA,errors = sample(seq(20,100,20),ss_el, replace = TRUE))

## Produce SPDs
pm_nexp_sm_spd <- spd(pm_nexp_sm_cal, timeRange = c(start,end))
pm_nexp_md_spd <- spd(pm_nexp_md_cal, timeRange = c(start,end))
pm_nexp_lg_spd <- spd(pm_nexp_lg_cal, timeRange = c(start,end))
pm_nexp_el_spd <- spd(pm_nexp_el_cal, timeRange = c(start,end))

## Plot SPDs
par(mfrow = c(2,2))
plot(pm_nexp_sm_spd, main = paste0("Sample size = ", ss_sm))
plot(pm_nexp_md_spd, main = paste0("Sample size = ", ss_md))
plot(pm_nexp_lg_spd, main = paste0("Sample size = ", ss_lg))
plot(pm_nexp_el_spd, main = paste0("Sample size = ", ss_el))

Nexp_pm <- list("Back-cal.dates" = list("Small.sample.size" = pm_nexp_sm_cal,
                                       "Medium.sample.size" = pm_nexp_md_cal,
                                       "Large.sample.size" = pm_nexp_lg_cal,
                                       "Extralarge.sample.size" = pm_nexp_el_cal),
               "SPDs" = list("Small.sample.size" = pm_nexp_sm_spd,
                             "Medium.sample.size" = pm_nexp_md_spd,
                             "Large.sample.size" = pm_nexp_lg_spd,
                             "Extralarge.sample.size" = pm_nexp_el_spd),
               "Parameters" = lambda1)


################################################################################
#####################           Logistic         ###############################
################################################################################

## Generate calendar dates

## Parameters for the logistic growth
k <- 0.001
r <- 0.005
pars_log <- c(k,r)
names(pars_log) <- c("k","r")

pm_log_sm <- sample(sort(ts, decreasing = TRUE),ss_sm,
                    prob=dLogisticGrowth(x = c(start:end), a = start, b = end, k = k, r = r, log = FALSE), 
                    replace = TRUE)
pm_log_md <- sample(sort(ts, decreasing = TRUE),ss_md,
                    prob=dLogisticGrowth(x = c(start:end), a = start, b = end, k = k, r = r, log = FALSE),
                    replace = TRUE)
pm_log_lg <- sample(sort(ts, decreasing = TRUE),ss_lg,
                    prob=dLogisticGrowth(x = c(start:end), a = start, b = end, k = k, r = r, log = FALSE),
                    replace = TRUE)
pm_log_el <- sample(sort(ts, decreasing = TRUE),ss_el,
                    prob=dLogisticGrowth(x = c(start:end), a = start, b = end, k = k, r = r, log = FALSE),
                    replace = TRUE)

## Uncalibrate
pm_log_sm_unc <- uncalibrate(pm_log_sm)
pm_log_md_unc <- uncalibrate(pm_log_md)
pm_log_lg_unc <- uncalibrate(pm_log_lg)
pm_log_el_unc <- uncalibrate(pm_log_el)

## Back-calibrate
pm_log_sm_cal <- calibrate(pm_log_sm_unc$ccCRA,errors = sample(seq(20,100,20),ss_sm, replace = TRUE))
pm_log_md_cal <- calibrate(pm_log_md_unc$ccCRA,errors = sample(seq(20,100,20),ss_md, replace = TRUE))
pm_log_lg_cal <- calibrate(pm_log_lg_unc$ccCRA,errors = sample(seq(20,100,20),ss_lg, replace = TRUE))
pm_log_el_cal <- calibrate(pm_log_el_unc$ccCRA,errors = sample(seq(20,100,20),ss_el, replace = TRUE))

## Produce SPDs
pm_log_sm_spd <- spd(pm_log_sm_cal, timeRange = c(start,end))
pm_log_md_spd <- spd(pm_log_md_cal, timeRange = c(start,end))
pm_log_lg_spd <- spd(pm_log_lg_cal, timeRange = c(start,end))
pm_log_el_spd <- spd(pm_log_el_cal, timeRange = c(start,end))

## Plot SPDs
par(mfrow = c(2,2))
plot(pm_log_sm_spd, main = paste0("Sample size = ", ss_sm))
plot(pm_log_md_spd, main = paste0("Sample size = ", ss_md))
plot(pm_log_lg_spd, main = paste0("Sample size = ", ss_lg))
plot(pm_log_el_spd, main = paste0("Sample size = ", ss_el))

Log_pm <- list("Back-cal.dates" = list("Small.sample.size" = pm_log_sm_cal,
                                       "Medium.sample.size" = pm_log_md_cal,
                                       "Large.sample.size" = pm_log_lg_cal,
                                       "Extralarge.sample.size" = pm_log_el_cal),
               "SPDs" = list("Small.sample.size" = pm_log_sm_spd,
                             "Medium.sample.size" = pm_log_md_spd,
                             "Large.sample.size" = pm_log_lg_spd,
                             "Extralarge.sample.size" = pm_log_el_spd),
               "Parameters" = pars_log)


################################################################################
#####################            Normal          ###############################
################################################################################

## Generate calendar dates

## Parameters for the normal distribution
m <- round(mean(seq(0,100,length.out = length(ts)))) ## Mean of the vector for sampling
sd <- round(sd(seq(0,100,length.out = length(ts)))) ## Standard deviation of the vector for sampling
pars_nor <- c(m,sd)
names(pars_nor) <- c("Mean", "Sd")

pm_nor_sm <- sample(sort(ts),ss_sm,prob=dnorm(seq(0,100,length.out = length(ts)), mean = m, sd = sd), replace = TRUE)
pm_nor_md <- sample(sort(ts),ss_md,prob=dnorm(seq(0,100,length.out = length(ts)), mean = m, sd = sd), replace = TRUE)
pm_nor_lg <- sample(sort(ts),ss_lg,prob=dnorm(seq(0,100,length.out = length(ts)), mean = m, sd = sd), replace = TRUE)
pm_nor_el <- sample(sort(ts),ss_el,prob=dnorm(seq(0,100,length.out = length(ts)), mean = m, sd = sd), replace = TRUE)

## Uncalibrate
pm_nor_sm_unc <- uncalibrate(pm_nor_sm)
pm_nor_md_unc <- uncalibrate(pm_nor_md)
pm_nor_lg_unc <- uncalibrate(pm_nor_lg)
pm_nor_el_unc <- uncalibrate(pm_nor_el)

## Back-calibrate
pm_nor_sm_cal <- calibrate(pm_nor_sm_unc$ccCRA,errors = sample(seq(20,100,20),ss_sm, replace = TRUE))
pm_nor_md_cal <- calibrate(pm_nor_md_unc$ccCRA,errors = sample(seq(20,100,20),ss_md, replace = TRUE))
pm_nor_lg_cal <- calibrate(pm_nor_lg_unc$ccCRA,errors = sample(seq(20,100,20),ss_lg, replace = TRUE))
pm_nor_el_cal <- calibrate(pm_nor_el_unc$ccCRA,errors = sample(seq(20,100,20),ss_el, replace = TRUE))

## Produce SPDs
pm_nor_sm_spd <- spd(pm_nor_sm_cal, timeRange = c(start,end))
pm_nor_md_spd <- spd(pm_nor_md_cal, timeRange = c(start,end))
pm_nor_lg_spd <- spd(pm_nor_lg_cal, timeRange = c(start,end))
pm_nor_el_spd <- spd(pm_nor_el_cal, timeRange = c(start,end))

## Plot SPDs
par(mfrow = c(2,2))
plot(pm_nor_sm_spd, main = paste0("Sample size = ", ss_sm))
plot(pm_nor_md_spd, main = paste0("Sample size = ", ss_md))
plot(pm_nor_lg_spd, main = paste0("Sample size = ", ss_lg))
plot(pm_nor_el_spd, main = paste0("Sample size = ", ss_el))

Nor_pm <- list("Back-cal.dates" = list("Small.sample.size" = pm_nor_sm_cal,
                                       "Medium.sample.size" = pm_nor_md_cal,
                                       "Large.sample.size" = pm_nor_lg_cal,
                                       "Extralarge.sample.size" = pm_nor_el_cal),
               "SPDs" = list("Small.sample.size" = pm_nor_sm_spd,
                             "Medium.sample.size" = pm_nor_md_spd,
                             "Large.sample.size" = pm_nor_lg_spd,
                             "Extralarge.sample.size" = pm_nor_el_spd),
               "Parameters" = pars_nor)


################################################################################
#####################            PERT            ###############################
################################################################################

## Generate calendar dates

## Parameters for the PERT distribution
a <- min(seq(0,100,length.out = length(ts)))
b <- max(seq(0,100,length.out = length(ts)))
mo <- quantile(seq(0,100,length.out = length(ts)))[4]
pars_per <- c(a,b,mo)
names(pars_per) <-c("Lowest", "Highest", "Mode")

pm_per_sm <- sample(sort(ts, decreasing = TRUE),ss_sm,
                    prob=dpert(seq(0,100,length.out = length(ts)),min = a, max = b, mode = mo), replace = TRUE)
pm_per_md <- sample(sort(ts, decreasing = TRUE),ss_md,
                    prob=dpert(seq(0,100,length.out = length(ts)),min = a, max = b, mode = mo), replace = TRUE)
pm_per_lg <- sample(sort(ts, decreasing = TRUE),ss_lg,
                    prob=dpert(seq(0,100,length.out = length(ts)),min = a, max = b, mode = mo), replace = TRUE)
pm_per_el <- sample(sort(ts, decreasing = TRUE),ss_el,
                    prob=dpert(seq(0,100,length.out = length(ts)),min = a, max = b, mode = mo), replace = TRUE)

## Uncalibrate
pm_per_sm_unc <- uncalibrate(pm_per_sm)
pm_per_md_unc <- uncalibrate(pm_per_md)
pm_per_lg_unc <- uncalibrate(pm_per_lg)
pm_per_el_unc <- uncalibrate(pm_per_el)

## Back-calibrate
pm_per_sm_cal <- calibrate(pm_per_sm_unc$ccCRA,errors = sample(seq(20,100,20),ss_sm, replace = TRUE))
pm_per_md_cal <- calibrate(pm_per_md_unc$ccCRA,errors = sample(seq(20,100,20),ss_md, replace = TRUE))
pm_per_lg_cal <- calibrate(pm_per_lg_unc$ccCRA,errors = sample(seq(20,100,20),ss_lg, replace = TRUE))
pm_per_el_cal <- calibrate(pm_per_el_unc$ccCRA,errors = sample(seq(20,100,20),ss_el, replace = TRUE))

## Produce SPDs
pm_per_sm_spd <- spd(pm_per_sm_cal, timeRange = c(start,end))
pm_per_md_spd <- spd(pm_per_md_cal, timeRange = c(start,end))
pm_per_lg_spd <- spd(pm_per_lg_cal, timeRange = c(start,end))
pm_per_el_spd <- spd(pm_per_el_cal, timeRange = c(start,end))

## Plot SPDs
par(mfrow = c(2,2))
plot(pm_per_sm_spd, main = paste0("Sample size = ", ss_sm))
plot(pm_per_md_spd, main = paste0("Sample size = ", ss_md))
plot(pm_per_lg_spd, main = paste0("Sample size = ", ss_lg))
plot(pm_per_el_spd, main = paste0("Sample size = ", ss_el))

Per_pm <- list("Back-cal.dates" = list("Small.sample.size" = pm_per_sm_cal,
                                       "Medium.sample.size" = pm_per_md_cal,
                                       "Large.sample.size" = pm_per_lg_cal,
                                       "Extralarge.sample.size" = pm_per_el_cal),
               "SPDs" = list("Small.sample.size" = pm_per_sm_spd,
                             "Medium.sample.size" = pm_per_md_spd,
                             "Large.sample.size" = pm_per_lg_spd,
                             "Extralarge.sample.size" = pm_per_el_spd),
               "Parameters" = pars_per)

Population_models <- list("Exponential" = Exp_pm,
                          "Neg.exponential" = Nexp_pm,
                          "Logistic" = Log_pm,
                          "Normal" = Nor_pm,
                          "PERT" = Per_pm)

saveRDS(Population_models, "Population_models.rds")






