
######### SCRIPT COMBINING ALL THE METHODS PROPOSED

## Load libraries
library(rcarbon)
library(oxcAAR) 
library(stringr)
library(Rextinct)
library(parallel)

## Oxcal setup
quickSetupOxcal() #this should be replaced by actual reference to the file path where oxcal is downloaded (otherwise would install it every time)

## Load required functions
source('functions.R')

set.seed(123)

## Load data and utilities
dates <- readRDS("../Simu_data/simuls_up.rds")
ndates <- 80
sims <- 1000

## Arrange for sample sizes in OLE
dates_ss <- c(5,10,ndates/2,ndates)

########### OLE MEDIANS

## Prepare to store results
OLE_medians <- data.frame("Estimate" = numeric(),
			  "upperCI" = numeric(),
			  "lowerCI" = numeric(),
			  "start_date" = numeric(),
			  "r" = numeric(),
			  "sd" = numeric(),
			  "sampled_dates" = numeric())
## Measure time
time_OLE_medians <- system.time({

## Subset dates from df
for (z in 1:sims){
	dates_subset <- dates[(ndates*z-79):(ndates*z),]
	dates_subset <- dates_subset[order(dates_subset$Year, decreasing = T),] ## Need this for different samples sizes

	############ OLE WITH MEDIANS

	dates_median <- rep(NA,nrow(dates_subset))
	for (i in 1:nrow(dates_subset)){
		dates_median[i] <- medCal(calibrate(dates_subset[i,1],dates_subset[i,2]))
	}
	
	## If the two first dates are the same, it gives NA. Sneaky way to sort that
	if (dates_median[1] == dates_median[2]){
		dates_median[1] <- dates_median[1]+1
	}
		
	for (i in 1:length(dates_ss)){
		dates_median_ss <- dates_median[1:dates_ss[i]]
		OLE_med_res <- OLE.test(dd = dates_median_ss, alpha = 0.05)
		OLE_medians[nrow(OLE_medians)+1,] <- c(OLE_med_res,dates_subset$start_date[1],dates_subset$r[1],dates_subset$Sd[1],dates_ss[i])
	}
}
})

saveRDS(OLE_medians,"../Results/OLE_medians_up.rds")
saveRDS(time_OLE_medians,"../Results/time_OLE_medians_up.rds")

## to save space
rm(OLE_medians)
rm(time_OLE_medians)

######### OLE RESAMPLE CALDATE (Our proposal)

rsmp <- 1000 ## Number of resamples. Shared between the two OLE resampling methods

## Prepare to store results
OLE_resamp_caldate <- data.frame("Estimate" = numeric(),
	        	 	 "upperCI" = numeric(),
			 	 "lowerCI" = numeric(),
			 	 "start_date" = numeric(),
			 	 "r" = numeric(),
			 	 "sd" = numeric(),
			 	 "sampled_dates" = numeric())


time_OLE_resamp_caldate <- system.time({
for (z in 1:(sims*4)){
	## Temporal OLE results
	temp_OLE <- data.frame("Estimate" = numeric(),
			       "upperCI" = numeric(),
		       	       "lowerCI" = numeric())
	
	dates_subset <- dates[(ndates*z-79):(ndates*z),]
	
	# Calibrate dates
	cal_dates <- calibrate(dates_subset[,1],dates_subset[,2])
	for (i in 1:rsmp){
		valid <- F ## I'll use this to check whether the matrix is singluar
		
		while(!valid){
			# Resample
			resamp_dates <- rep(NA,ndates)
			for (k in 1:ndates){
				resamp_dates[k] <- sample(cal_dates$grids[[k]][,1],size = 1,cal_dates$grids[[k]][,2], replace = T)
			}
			resamp_dates <- resamp_dates[order(resamp_dates, decreasing = T)]
			
			## If the two first dates are the same, it gives NA. Sneaky way to sort that
			if (resamp_dates[1] == resamp_dates[2]){
				resamp_dates[1] <- resamp_dates[1]+1
			}

			for (j in 1:length(dates_ss)){
				resamp_dates <- resamp_dates[1:dates_ss[j]]
			
				## Check that matrix is singlular and produce results only if it is
				OLE_res <- OLE.test(dd = resamp_dates, alpha = 0.05)
		
				if (!any(is.na(OLE_res))){
					temp_OLE[i,] <- OLE_res
					valid <- TRUE
				}
			}
		}
	}
	
	OLE_resamp_caldate[z,] <- c(apply(temp_OLE,2,mean),dates_subset$start_date[1],dates_subset$r[1],dates_subset$Sd[1],dates_ss[z])
}
})

saveRDS(OLE_resamp_caldate,"../Results/OLE_resamp_caldate_up.rds")
saveRDS(time_OLE_resamp_caldate,"../Results/time_OLE_resamp_caldate_up.rds")

## to save space
rm(OLE_resamp_caldate)
rm(time_OLE_resamp_caldate)

######### OLE RESAMPLE NORM DIST (as in Key et al., 2021 and Bebber and Key 2022)

## Prepare to store results
OLE_resamp_norm <- data.frame("Estimate" = numeric(),
	                      "upperCI" = numeric(),
			      "lowerCI" = numeric(),
			      "start_date" = numeric(),
			      "r" = numeric(),
			      "sd" = numeric(),
			      "sampled_dates" = numeric())

time_OLE_resamp_norm <- system.time({
for (z in 1:(sims*4)){
	## Temporal OLE results
	temp_OLE <- data.frame("Estimate" = numeric(),
			       "upperCI" = numeric(),
		       	       "lowerCI" = numeric())
	
	dates_subset <- dates[(ndates*z-79):(ndates*z),]
	
	# Calibrate dates
	cal_dates <- calibrate(dates_subset[,1],dates_subset[,2])
	
	for (i in 1:rsmp){
		valid <- F ## I'll use this to check whether the matrix is singluar
		
		while(!valid){
			# Resample
			resamp_dates <- rep(NA,ndates)
			for (k in 1:ndates){
				dat_mean <- mean(cal_dates$grids[[k]][,1])
				dat_sd <- (max(cal_dates$grids[[k]][,1])-dat_mean)/2
				resamp_dates[k] <- rnorm(1,dat_mean,dat_sd)
			}
			resamp_dates <- resamp_dates[order(resamp_dates, decreasing = T)]
			## If the two first dates are the same, it gives NA. Sneaky way to sort that
			if (resamp_dates[1] == resamp_dates[2]){
				resamp_dates[1] <- resamp_dates[1]+1
			}
			
			for (j in 1:length(dates_ss)){
				resamp_dates <- resamp_dates[1:dates_ss[j]]
			
				## Check that matrix is singlular and produce results only if it is
				OLE_res <- OLE.test(dd = resamp_dates, alpha = 0.05)
		
				if (!any(is.na(OLE_res))){
					temp_OLE[i,] <- OLE_res
					valid <- TRUE
				}
			}
		}
	}
	
	OLE_resamp_norm[z,] <- c(apply(temp_OLE,2,mean),dates_subset$start_date[1],dates_subset$r[1],dates_subset$Sd[1],dates_ss[z])
}
})

saveRDS(OLE_resamp_norm,"../Results/OLE_resamp_norm_up.rds")
saveRDS(time_OLE_resamp_norm,"../Results/time_OLE_resamp_norm_up.rds")

## to save space
rm(OLE_resamp_norm)
rm(time_OLE_resamp_norm)

######### OLE RESAMPLE UNI DIST (as in Djakovic et al. 2022)

## Prepare to store results
OLE_resamp_uni <- data.frame("Estimate" = numeric(),
	                     "upperCI" = numeric(),
			     "lowerCI" = numeric(),
			     "start_date" = numeric(),
			     "r" = numeric(),
			     "sd" = numeric(),
			     "sampled_dates" = numeric())

time_OLE_resamp_uni <- system.time({
for (z in 1:(sims*4)){
	## Temporal OLE results
	temp_OLE <- data.frame("Estimate" = numeric(),
			       "upperCI" = numeric(),
		       	       "lowerCI" = numeric())
	
	dates_subset <- dates[(ndates*z-79):(ndates*z),]
	
	# Calibrate dates
	cal_dates <- calibrate(dates_subset[,1],dates_subset[,2])
	
	for (i in 1:rsmp){
		valid <- F ## I'll use this to check whether the matrix is singluar
		
		while(!valid){
			# Resample
			resamp_dates <- rep(NA,ndates)
			for (k in 1:ndates){
				resamp_dates[k] <- runif(1,min(cal_dates$grids[[k]][,1]),max(cal_dates$grids[[k]][,1]))
			}
			resamp_dates <- resamp_dates[order(resamp_dates, decreasing = T)]
			## If the two first dates are the same, it gives NA. Sneaky way to sort that
			if (resamp_dates[1] == resamp_dates[2]){
				resamp_dates[1] <- resamp_dates[1]+1
			}
			
			for (j in 1:length(dates_ss)){
				resamp_dates <- resamp_dates[1:dates_ss[j]]
			
				## Check that matrix is singlular and produce results only if it is
				OLE_res <- OLE.test(dd = resamp_dates, alpha = 0.05)
		
				if (!any(is.na(OLE_res))){
					temp_OLE[i,] <- OLE_res
					valid <- TRUE
				}
			}
		}
	}
	
	OLE_resamp_uni[z,] <- c(apply(temp_OLE,2,mean),dates_subset$start_date[1],dates_subset$r[1],dates_subset$Sd[1],dates_ss[z])
}
})

saveRDS(OLE_resamp_uni,"../Results/OLE_resamp_uni_up.rds")
saveRDS(time_OLE_resamp_uni,"../Results/time_OLE_resamp_uni_up.rds")

## to save space
rm(OLE_resamp_uni)
rm(time_OLE_resamp_uni)

######### Oxcal

# Cores
ncores <- 10

time_oxcal <- system.time({
oxcal_res <- mclapply(1:sims, function(z) {
			      dates_subset <- dates[(ndates*z-79):(ndates*z), ]
			      list(unif  = oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "uniform"),
				   trap  = oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "trapezoid")
			      )
			       }, mc.cores = ncores)
})

saveRDS(oxcal_res,"../Results/oxcal_up.rds")
saveRDS(time_oxcal,"../Results/time_oxcal_up.rds")

## to save space
rm(oxcal_res)
rm(time_oxcal)

######### CRIWM

sim <- sims
C14_errors <- c(20,50)
r_vals <- c(0.01,0.02,0.03,0.04)

# Combinations of (i, r, sd)
comb <- expand.grid(sim = 1:sim,r = r_vals,sd = C14_errors)

ncores <- 10

time_criwm <- system.time({

  criwm_res <- mclapply(1:nrow(comb), function(idx) {
    row <- comb[idx, ]

    # Build file name
    filepath <- paste0("../Simu_data/Uncal_YearsBP_sim_up_",row$sim, "_r_",row$r,"_error_",row$sd,".txt")

    # Run criwm
    res <- criwm(filepath)

    # Extract values
    data.frame(Estimate = res$criwm[2,2],upperCI = res$criwm[2,3],lowerCI = res$criwm[2,1],sim = row$sim,r = row$r,sd = row$sd)
  }, mc.cores = ncores)
})

# Single data frame
CRIWM <- do.call(rbind, criwm_res)

saveRDS(CRIWM,"../Results/CRIWM_up.rds")
saveRDS(time_criwm,"../Results/time_criwm_up.rds")

## to save space
rm(CRIWM)
rm(time_criwm)
