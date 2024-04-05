

################################################################################
################################################################################

### Comparison of Oxcal/OLE using a logistic model

################################################################################
################################################################################


## Load packages
library(rcarbon)
library(nimbleCarbon)
library(oxcAAR)
library(foreach)
library(doParallel)

## Settings
n_sim <- 5
cores <- 10
set.seed(123)

## Load functions
source("Functions.R")

### Start parallel
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- foreach (j = 1:n_sim, .packages = c('rcarbon', 'nimbleCarbon','oxcAAR')) %dopar%{
  
  
  ################################################################################
  ######################     STEP 1. GENERATE THE MODEL     ######################
  ################################################################################
  
  
  ## Set parameters (these will be randomised and simulated n number of times)
  r <- 0.01 ## Rate
  s <- sample(seq(3000,10000),1) ## Starting date 
  e <- s-sample(seq(300,1000),1) ## End date. 
  n <- sample(seq(15,50),1) ## Number of dates 
  k <- 0.00001 ## Initial proportion of carrying capacity
  alpha <- 0.05 ## Significance OLE
  
  ## Generate the model 
  dates <- replicate(n,rLogisticGrowth(a = s, b = e, k = k, r = r)) ## Logistic model
  errors <- sample(seq(20,100,20), n, replace = TRUE)
  
  ## Calibrate
  cal_dates <- calibrate(dates, errors = errors)
  
  
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
  
  
  ################################################################################
  ######################       STEP 4. OxCal Trapezoid      ######################
  ################################################################################
  
  quickSetupOxcal(path = getwd())
  
  ## Produce file for OxCal script
  oxcalScriptGen(c14age=dates,errors=errors,fn=paste0('trap_set.oxcal'),model='trapezoid')
  
  ## Read oxcal script
  oxScript_trap  <- readLines(paste0('trap_set.oxcal'))
  
  ## Execute oxcal script
  res_file_trap <- executeOxcalScript(oxScript_trap)
  
  ## Store result
  res_trap <- readOxcalOutput(res_file_trap)
  
  ## Human readable
  res_trap <- parseFullOxcalOutput(res_trap)
  
  
  ################################################################################
  ######################        STEP 5. OxCal Uniform       ######################
  ################################################################################
  
  ## Produce file for OxCal script
  oxcalScriptGen(c14age=dates,errors=errors,fn=paste0('uni_set.oxcal'),model='uniform')
  
  ## Read oxcal script
  oxScript_uni  <- readLines(paste0('uni_set.oxcal'))
  
  ## Execute oxcal script
  res_file_uni <- executeOxcalScript(oxScript_uni)
  
  ## Store result
  res_uni <- readOxcalOutput(res_file_uni)
  
  ## Human readable
  res_uni <- parseFullOxcalOutput(res_uni)
  
  
  ################################################################################
  ######################       STEP 6. OxCal Gaussian       ######################
  ################################################################################
  
  ## Produce file for OxCal script
  oxcalScriptGen(c14age=dates,errors=errors,fn=paste0('gau_set.oxcal'),model='gaussian')
  
  ## Read oxcal script
  oxScript_gau  <- readLines(paste0('gau_set.oxcal'))
  
  ## Execute oxcal script
  res_file_gau <- executeOxcalScript(oxScript_gau)
  
  ## Store result
  res_gau <- readOxcalOutput(res_file_gau)
  
  ## Human readable
  res_gau <- parseFullOxcalOutput(res_gau)
  
  
  ## Produce results
  results <- list(OLE_med, OLE_res, res_trap, res_uni, res_gau)
  
}

stopCluster(cl)
saveRDS(results,"results.rds")


