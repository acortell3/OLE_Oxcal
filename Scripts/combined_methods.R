
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
source('oxcal_functions.R')
source('functions.R')

## Load data and utilities
dates <- readRDS("../Simu_data/simuls.rds")
ndates <- 80

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
system.time({

## Subset dates from df
for (z in 1:10){
	dates_subset <- dates[(ndates*z-79):(ndates*z),]
	dates_subset <- dates_subset[order(dates_subset$Year, decreasing = T),] ## Need this for different samples sizes

	############ OLE WITH MEDIANS

	dates_median <- rep(NA,nrow(dates_subset))
	for (i in 1:nrow(dates_subset)){
		dates_median[i] <- medCal(calibrate(dates_subset[i,1],dates_subset[i,2]))
	}

	for (i in 1:length(dates_ss)){
		dates_median_ss <- dates_median[1:dates_ss[i]]
		OLE_med_res <- OLE.test(dd = dates_median_ss, alpha = 0.05)
		OLE_medians[nrow(OLE_medians)+1,] <- c(OLE_med_res,dates_subset$start_date[1],dates_subset$r[1],dates_subset$Sd[1],dates_ss[i])
	}
}
})


######### OLE RESAMPLE CALDATE (Our proposal)
## Ok, so what they do in (Key et al. 2021, Bebber & Key, 2022 or Djakovic et al., 2022) is basically resample from the mean of the radiocarbon date with half standard deviation. They sample 10,000 times and from the results, the just get the mean of the estimated means and CIs. I'll do the last bit (although it is dodgy), but I'm not resampling from the radiocarbon date because this does not take into account the calibration error. I'll sample from the probability distribution of the calibrated date. 

rsmp <- 1000 ## Number of resamples. Shared between the two OLE resampling methods

## Prepare to store results
OLE_resamp_caldate <- data.frame("Estimate" = numeric(),
	        	 	 "upperCI" = numeric(),
			 	 "lowerCI" = numeric(),
			 	 "start_date" = numeric(),
			 	 "r" = numeric(),
			 	 "sd" = numeric(),
			 	 "sampled_dates" = numeric())


system.time({
for (z in 1:40){
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

######### OLE RESAMPLE NORM DIST (as in Key et al., 2021 and Bebber and Key 2022)
## Ok, so what they do in (Key et al. 2021, Bebber & Key, 2022 or Djakovic et al., 2022) is basically resample from the mean of the radiocarbon date with half standard deviation. They sample 10,000 times and from the results, the just get the mean of the estimated means and CIs. I'll do the last bit (although it is dodgy), but I'm not resampling from the radiocarbon date because this does not take into account the calibration error. I'll sample from the probability distribution of the calibrated date. 


## Prepare to store results
OLE_resamp_norm <- data.frame("Estimate" = numeric(),
	                      "upperCI" = numeric(),
			      "lowerCI" = numeric(),
			      "start_date" = numeric(),
			      "r" = numeric(),
			      "sd" = numeric(),
			      "sampled_dates" = numeric())

system.time({
for (z in 1:40){
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

######### OLE RESAMPLE UNI DIST (as in Djakovic et al. 2022)
## Ok, so what they do in (Key et al. 2021, Bebber & Key, 2022 or Djakovic et al., 2022) is basically resample from the mean of the radiocarbon date with half standard deviation. They sample 10,000 times and from the results, the just get the mean of the estimated means and CIs. I'll do the last bit (although it is dodgy), but I'm not resampling from the radiocarbon date because this does not take into account the calibration error. I'll sample from the probability distribution of the calibrated date. 


## Prepare to store results
OLE_resamp_uni <- data.frame("Estimate" = numeric(),
	                     "upperCI" = numeric(),
			     "lowerCI" = numeric(),
			     "start_date" = numeric(),
			     "r" = numeric(),
			     "sd" = numeric(),
			     "sampled_dates" = numeric())

system.time({
for (z in 1:40){
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

######### Oxcal

# Cores
ncores <- 10

system.time({
oxcal_res <- mclapply(1:10, function(z) {
			      dates_subset <- dates[(ndates*z-79):(ndates*z), ]
			      list(unif  = oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "uniform"),
				   trap  = oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "trapezoid")
			      )
			       }, mc.cores = ncores)
})

######### CRIWM

## Need this to recover files by name. They match the ones in Simu_dates.R
sim <- 1
C14_errors <- c(20,50)
r <- c(0.01,0.02,0.03,0.04)
index <- 1

## Prepare to store results
CRIWM <- data.frame("Estimate" = numeric(),
	            "upperCI" = numeric(),
		    "lowerCI" = numeric(),
		    "sim" = numeric(),
		    "r" = numeric(),
		    "sd" = numeric())
system.time({
for (i in 1:sim){
	for (k in 1:length(r)){
		for (j in 1:length(C14_errors)){
			creiwm_res_temp <- criwm(paste0("../Simu_data/Uncal_YearsBP_sim_",i,"_r_",r[k],"_error_",C14_errors[j],".txt"))
			CRIWM[index,] <- c(creiwm_res_temp$criwm[2,2],
					   creiwm_res_temp$criwm[2,3],
					   creiwm_res_temp$criwm[2,1],
					   sim[i],r[k],C14_errors[j])
			index <- index + 1
		}
	}
}
})
