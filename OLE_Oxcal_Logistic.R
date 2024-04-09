

################################################################################
################################################################################

### Comparison of Oxcal/OLE using a logistic model

################################################################################
################################################################################


## Load packages
library(rcarbon)
library(nimbleCarbon)
library(oxcAAR)
library(parallel)

## Settings
n_sim <- 5
cores <- 5
set.seed(123)

## Load functions
source("Functions.R")

################################################################################
######################     STEP 1. GENERATE THE MODEL     ######################
################################################################################

## Parameters
r <- rep(0.01, n_sim) ## Rate
s <- sample(seq(3000,10000),n_sim, replace = TRUE) ## Starting date 
e <- s-sample(seq(300,1000),n_sim, replace = TRUE) ## End date. 
n <- sample(seq(15,50),n_sim, replace = TRUE) ## Number of dates 
k <- rep(0.00001, n_sim) ## Initial proportion of carrying capacity
alpha <- rep(0.05, n_sim) ## Significance OLE

pars <- data.frame("r" = r,
                   "s" = s, 
                   "e" = e, 
                   "n" = n,
                   "k" = k,
                   "alpha" = alpha)

## Generate the model
chrono <- vector("list", length = n_sim)

for(i in 1:n_sim){
  dates <- replicate(pars$n[i], rLogisticGrowth(a = pars$s[i], b = pars$e[i], k = pars$k[i],
                                                r = pars$r[i])) ## Logistic model
  errors <- sample(seq(20,100,20), pars$n[i], replace = TRUE)
  
  dates_df <- data.frame("dates" = dates, "errors" = errors)
  chrono[[i]] <- dates_df
}

## Calibrate dates for OLE
cal_dates <- vector("list", length = n_sim)

for (i in 1:n_sim){
  cal_dates[[i]] <- calibrate(chrono[[i]]$dates, chrono[[i]]$errors)
}


################################################################################
######################         STEP 2. OLE median         ######################
################################################################################

## Extract medians per date
medians_list <- vector("list", length = n_sim)

for (i in 1:n_sim){
  
  ## Set appropriate dataset
  dataset <- cal_dates[[i]]
  n_samples <- pars$n[i]
  med_dates <- rep(NA,n_samples)
  
  ## Extract median per date
  for (j in 1:n_samples){
    med_dates[j] <- round(median(dataset$grids[[j]][,1]))
  }
  
  ## Assign to dataset
  medians_list[[i]] <- med_dates 
}

## OLE on medians
OLE_med <- data.frame("Estimate" = rep(NA,n_sim),
                      "upperCI" = rep(NA,n_sim),
                      "lowerCI" = rep(NA,n_sim))

for (i in 1:n_sim){
  OLE_med[i,] <- OLE.test(medians_list[[i]], alpha = alpha)
}


################################################################################
######################        STEP 3. OLE resampling      ######################
################################################################################

## Extract medians per date
random_list <- vector("list", length = n_sim)

for (i in 1:n_sim){
  
  ## Set appropriate dataset
  dataset <- cal_dates[[i]]
  n_samples <- pars$n[i]
  ran_dates <- rep(NA,n_samples)
  
  ## Extract median per date
  for (j in 1:n_samples){
    ran_dates[j] <- sample(dataset$grids[[j]][,1], 1, prob = dataset$grids[[j]][,2])
  }
  
  ## Assign to dataset
  random_list[[i]] <- ran_dates 
}

## OLE on resampled
OLE_res <- data.frame("Estimate" = rep(NA,n_sim),
                      "upperCI" = rep(NA,n_sim),
                      "lowerCI" = rep(NA,n_sim))

for (i in 1:n_sim){
  OLE_res[i,] <- OLE.test(random_list[[i]], alpha = alpha)
}


################################################################################
######################            STEP 4. OxCal           ######################
################################################################################

## Save Oxcal scripts
for (i in 1:n_sim){
  oxcalScriptGen(c14age=chrono[[i]]$dates,errors=chrono[[i]]$errors,
                 fn=paste0('OxScripts/trap_set_',i,'.oxcal'),model='trapezoid')
  oxcalScriptGen(c14age=chrono[[i]]$dates,errors=chrono[[i]]$errors,
                 fn=paste0('OxScripts/uni_set_',i,'.oxcal'),model='uniform')
  oxcalScriptGen(c14age=chrono[[i]]$dates,errors=chrono[[i]]$errors,
                 fn=paste0('OxScripts/gau_set_',i,'.oxcal'),model='gaussian')
}

## Fit depending on prior
cl  <- makeCluster(cores)
x  <- 1:n_sim

oxcal_trap  <- parLapply(cl=cl,X=x,fun=logi_trap) ## Trapezoid
oxcal_uni  <- parLapply(cl=cl,X=x,fun=logi_uni) ## Uniform
oxcal_gau  <- parLapply(cl=cl,X=x,fun=logi_gau) ## Gaussian

################################################################################
######################        STEP 5. Store results       ######################
################################################################################


results <- list(OLE_med, OLE_res, oxcal_trap, oxcal_uni, oxcal_gau)
saveRDS(results, "results.rds")



