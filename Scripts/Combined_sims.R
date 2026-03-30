
######### SCRIPT COMBINING ALL THE METHODS PROPOSED

#### Authors: Alfredo Cortell-Nicolau & Enrico R. Crema

#################################################
########### HOLOCENE DATES
#################################################


## Load libraries
library(rcarbon)
library(oxcAAR) 
library(stringr)
library(Rextinct)
library(parallel)

## Oxcal setup
quickSetupOxcal() #this should be replaced by actual reference to the file path where oxcal is downloaded (otherwise would install it every time)

## Load required functions
source('Functions.R')

## Load data and utilities
dates <- readRDS("../Simu_data/simuls_hol.rds")
dates_seq <- 130
sims <- nrow(dates) / dates_seq

C14_errors <- c(20,50)
r_vals <- c(0.01,0.03,0.06)
sample_size <- c(10,40,80)

ncores <- 35
dates_OLE <- c(5,10,40,80)

########## OLE MEDIANS
## Measure time
time_OLE_medians <- system.time({

OLE_medians_list <- mclapply(
			     seq_len(sims),
			     function(z){
				     set.seed(z)
				     pre_dates_subset <- dates[(dates_seq*z - (dates_seq-1)):(dates_seq*z),]
				     res_list <- list()
				     row_id <- 1
				     for (h in seq_along(sample_size)){
					     dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
					     ndates <- sample_size[h]
					     st <- which(pre_dates_subset$sample_size == ndates)[1]
					     dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
					     
					     ## OLE WITH MEDIANS
					     dates_median <- numeric(nrow(dates_subset))
					     
					     for (i in seq_len(nrow(dates_subset))){
						     dates_median[i] <- medCal(calibrate(dates_subset[i,1], dates_subset[i,2]))}
					     dates_median <- sort(dates_median, decreasing = TRUE)
					     
					     for (i in seq_along(dates_ss)){
						     dates_median_ss <- dates_median[1:dates_ss[i]]
						     
						     ## Avoid ties
						     for (b in 2:length(dates_median_ss)){
							     if (dates_median_ss[b] == dates_median_ss[b-1]){
								     dates_median_ss[b] <- dates_median_ss[b] + 1
							     }
						     }
						     
						     ## Run OLE
						     OLE_med_res <- OLE.test(dd = dates_median_ss, alpha = 0.05)
						     res_list[[row_id]] <- data.frame(Estimate = OLE_med_res[1],upperCI = OLE_med_res[2],lowerCI = OLE_med_res[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[i],sample_size = sample_size[h],simID_seed = z)
						     row_id <- row_id + 1
					     }
				     }
				     do.call(rbind, res_list)
			     },
			     mc.cores = ncores
)
OLE_medians <- do.call(rbind, OLE_medians_list)
})

saveRDS(OLE_medians,"../Results/OLE_medians_hol.rds")
saveRDS(time_OLE_medians,"../Results/time_OLE_medians_hol.rds")

## to save space
rm(OLE_medians)
rm(time_OLE_medians)

######### OLE RESAMPLE CALDATE (Our proposal)

rsmp <- 1000

time_OLE_resamp_caldate <- system.time({

OLE_resamp_list <- mclapply(
			    seq_len(sims),
			    function(z){
				    
				    set.seed(z)
				    
				    pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
				    
				    res_list <- list()
				    row_id <- 1
				    
				    for (h in seq_along(sample_size)){
					    dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
					    ndates <- sample_size[h]
					    st <- which(pre_dates_subset$sample_size == ndates)[1]
					    dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
					    
					    # Calibrate dates
					    cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
					    for (j in seq_along(dates_ss)){
						    temp_OLE <- matrix(NA, rsmp, 3)
						    
						    for (i in seq_len(rsmp)){
							    valid <- FALSE
							    while(!valid){
								    # Resample calibrated dates
								    resamp_dates <- vapply(seq_len(ndates),function(k){sample(cal_dates$grids[[k]][,1],size = 1,prob = cal_dates$grids[[k]][,2],replace = TRUE)},numeric(1))
								    resamp_dates <- sort(resamp_dates, decreasing = TRUE)
								    resamp_dates_temp <- resamp_dates[1:dates_ss[j]]
								    # Break ties
								    for (b in 2:length(resamp_dates_temp)){
									    if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
										    resamp_dates_temp[b] <- resamp_dates_temp[b]+1
									    }
								    }
								    OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
								    if (!any(is.na(OLE_res))){
									    temp_OLE[i,] <- OLE_res
									    valid <- TRUE
								    }
							    }
						    }
						    means <- colMeans(temp_OLE)
						    
						    res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
						    row_id <- row_id + 1
					    }
				    }
				    
				    do.call(rbind, res_list)
			    },
			    mc.cores = ncores
)
OLE_resamp_caldate <- do.call(rbind, OLE_resamp_list)
})

saveRDS(OLE_resamp_caldate,"../Results/OLE_resamp_caldate_hol.rds")
saveRDS(time_OLE_resamp_caldate,"../Results/time_OLE_resamp_caldate_hol.rds")

## to save space
rm(OLE_resamp_caldate)
rm(time_OLE_resamp_caldate)

######### OLE RESAMPLE NORM DIST (as in Key et al., 2021 and Bebber and Key 2022)

time_OLE_resamp_norm <- system.time({

OLE_resamp_norm_list <- mclapply(
				 seq_len(sims),
				 function(z){
					 
					 set.seed(z)
					 
					 pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
					 
					 res_list <- list()
					 row_id <- 1
					 
					 for (h in seq_along(sample_size)){
						 dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
						 ndates <- sample_size[h]
						 st <- which(pre_dates_subset$sample_size == ndates)[1]
						 dates_subset <- pre_dates_subset[st:(st + ndates - 1),]

						 # Calibrate dates
						 cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
						 
						 for (j in seq_along(dates_ss)){						 
							 temp_OLE <- matrix(NA, rsmp, 3)
							 
							 for (i in seq_len(rsmp)){
								 valid <- FALSE
								 while(!valid){
									 # Resample from normal distributions
									 resamp_dates <- vapply(seq_len(ndates),function(k){dat_mean <- mean(cal_dates$grids[[k]][,1])
												dat_sd <- (max(cal_dates$grids[[k]][,1]) - dat_mean)/2
												rnorm(1, dat_mean, dat_sd)},numeric(1))
									 
									 resamp_dates <- sort(resamp_dates, decreasing = TRUE)
									 resamp_dates_temp <- resamp_dates[1:dates_ss[j]]							 
									 # Break ties
									 for (b in 2:length(resamp_dates_temp)){
										 if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
											 resamp_dates_temp[b] <- resamp_dates_temp[b]+1
										 }
									 }
									 OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
									 if (!any(is.na(OLE_res))){
										 temp_OLE[i,] <- OLE_res
										 valid <- TRUE
									 }
								 }
							 }
							 
							 means <- colMeans(temp_OLE)
							 res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
							 row_id <- row_id + 1
						 }
					 }
					 
					 do.call(rbind, res_list)
				 },
				 mc.cores = ncores
)
OLE_resamp_norm <- do.call(rbind, OLE_resamp_norm_list)
})

saveRDS(OLE_resamp_norm,"../Results/OLE_resamp_norm_hol.rds")
saveRDS(time_OLE_resamp_norm,"../Results/time_OLE_resamp_norm_hol.rds")

## to save space
rm(OLE_resamp_norm)
rm(time_OLE_resamp_norm)

######### OLE RESAMPLE UNI DIST (Djakovic et al. 2022)

time_OLE_resamp_unif <- system.time({

OLE_resamp_unif_list <- mclapply(
				 seq_len(sims),function(z){
					 set.seed(z)
					 pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
					 
					 res_list <- list()
					 row_id <- 1
					 
					 for (h in seq_along(sample_size)){
						 dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
						 ndates <- sample_size[h]
						 
						 st <- which(pre_dates_subset$sample_size == ndates)[1]
						 dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
						 
						 # Calibrate dates
						 cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
						 
						 for (j in seq_along(dates_ss)){
							 temp_OLE <- matrix(NA, rsmp, 3)
							 
							 for (i in seq_len(rsmp)){
								 valid <- FALSE
								 
								 while(!valid){
									 # Resample from uniform distributions
									 resamp_dates <- vapply(seq_len(ndates),function(k){
													runif(1,min(cal_dates$grids[[k]][,1]),max(cal_dates$grids[[k]][,1]))},numeric(1))
									 resamp_dates <- sort(resamp_dates, decreasing = TRUE)
									 resamp_dates_temp <- resamp_dates[1:dates_ss[j]]
									 
									 # Break ties
									 for (b in 2:length(resamp_dates_temp)){
										 if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
											 resamp_dates_temp[b] <- resamp_dates_temp[b]+1}
									 }
									 OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
									 
									 if (!any(is.na(OLE_res))){
										 temp_OLE[i,] <- OLE_res
										 valid <- TRUE
									 }
								 }
							 }
							 
							 means <- colMeans(temp_OLE)
							 res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
							 row_id <- row_id + 1
						 }
					 }
					 
					 do.call(rbind, res_list)
				 },
				 mc.cores = ncores
)
OLE_resamp_unif <- do.call(rbind, OLE_resamp_unif_list)
})

saveRDS(OLE_resamp_unif,"../Results/OLE_resamp_unif_hol.rds")
saveRDS(time_OLE_resamp_unif,"../Results/time_OLE_resamp_unif_hol.rds")

## to save space
rm(OLE_resamp_unif)
rm(time_OLE_resamp_unif)

######### Oxcal
## Arrange for sample sizes in Oxcal and CRIWM
sim_index <- expand.grid(dts = 1:(nrow(dates)/dates_seq), ESS = sample_size) ## For indexing sims


## Oxcal uniform
time_oxcal_unif <- system.time({
  oxcal_unif <- mclapply(
	1:nrow(sim_index),function(z) {
		dts_id <- sim_index$dts[z]
		ESS <- sim_index$ESS[z]
		
		pre_dates_subset <- dates[(dates_seq*dts_id - (dates_seq-1)):(dates_seq*dts_id),]
		dates_subset <- pre_dates_subset[pre_dates_subset$sample_size == ESS,]
		
		est <- oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "uniform")
		est$Sd <- unique(dates_subset$Sd)
	        est$start_date <- unique(dates_subset$start_date)
	        est$r <- unique(dates_subset$r)
	        est$sample_size <- nrow(dates_subset)
	        est$ESS <- ESS
	
	        return(est)
	},
	mc.cores = ncores
  )
})

saveRDS(oxcal_unif,"../Results/oxcal_unif_hol.rds")
saveRDS(time_oxcal_unif,"../Results/time_oxcal_unif_hol.rds")

## to save space
rm(oxcal_unif)
rm(time_oxcal_unif)

## Oxcal trapezoid
time_oxcal_trap <- system.time({
  oxcal_trap <- mclapply(
	1:nrow(sim_index),function(z) {
		dts_id <- sim_index$dts[z]
		ESS <- sim_index$ESS[z]
		
		pre_dates_subset <- dates[(dates_seq*dts_id - (dates_seq-1)):(dates_seq*dts_id),]
		dates_subset <- pre_dates_subset[pre_dates_subset$sample_size == ESS,]
		
		est <- oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "trapezoid")
		est$Sd <- unique(dates_subset$Sd)
	        est$start_date <- unique(dates_subset$start_date)
	        est$r <- unique(dates_subset$r)
	        est$sample_size <- nrow(dates_subset)
	        est$ESS <- ESS
	
	        return(est)
	},
	mc.cores = ncores
  )
})

saveRDS(oxcal_trap,"../Results/oxcal_trap_hol.rds")
saveRDS(time_oxcal_trap,"../Results/time_oxcal_trap_hol.rds")

## to save space
rm(oxcal_trap)
rm(time_oxcal_trap)


######### CRIWM
sim <- 1000

# Combinations of (i, r, sd)
comb <- expand.grid(sim = 1:sim,r = r_vals,Sd = C14_errors, ESS = sample_size)

time_criwm <- system.time({

  criwm_res <- mclapply(1:nrow(comb), function(idx) {
    set.seed(idx)
    row <- comb[idx, ]
    
    # Build file name
    filepath <- paste0("../Simu_data/Uncal_YearsBP_sim_hol_",row$sim, "_r_",row$r,"_error_",row$Sd,"_sample_size_",row$ESS,".txt")

    # Run criwm
    res <- criwm(filepath, signor_lipps="arr")

    # Extract values
    r_v <- row$r
    Sd_v <- row$Sd
    ESS_v <- row$ESS
    data.frame(Estimate = res$criwm[2,2],upperCI = res$criwm[2,3],lowerCI = res$criwm[2,1],start_date = unique(dates$start_date)[comb$sim[idx]],r = unique(r_v),Sd = unique(Sd_v),ESS = unique(ESS_v),sample_size = unique(ESS_v),simID_seed = idx)
  }, mc.cores = ncores)
})

# Single data frame
CRIWM <- do.call(rbind, criwm_res)

saveRDS(CRIWM,"../Results/CRIWM_hol.rds")
saveRDS(time_criwm,"../Results/time_criwm_hol.rds")

## to save space
rm(CRIWM)
rm(time_criwm)


#################################################
########### UPPER PALAEOLITHIC DATES
#################################################

## Load data and utilities
dates <- readRDS("../Simu_data/simuls_upl.rds")

## Utilities picked up from the beginning of the script

########## OLE MEDIANS
## Measure time
time_OLE_medians <- system.time({

OLE_medians_list <- mclapply(
			     seq_len(sims),
			     function(z){
				     set.seed(z)
				     pre_dates_subset <- dates[(dates_seq*z - (dates_seq-1)):(dates_seq*z),]
				     res_list <- list()
				     row_id <- 1
				     for (h in seq_along(sample_size)){
					     dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
					     ndates <- sample_size[h]
					     st <- which(pre_dates_subset$sample_size == ndates)[1]
					     dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
					     
					     ## OLE WITH MEDIANS
					     dates_median <- numeric(nrow(dates_subset))
					     
					     for (i in seq_len(nrow(dates_subset))){
						     dates_median[i] <- medCal(calibrate(dates_subset[i,1], dates_subset[i,2]))}
					     dates_median <- sort(dates_median, decreasing = TRUE)
					     
					     for (i in seq_along(dates_ss)){
						     dates_median_ss <- dates_median[1:dates_ss[i]]
						     
						     ## Avoid ties
						     for (b in 2:length(dates_median_ss)){
							     if (dates_median_ss[b] == dates_median_ss[b-1]){
								     dates_median_ss[b] <- dates_median_ss[b]+1
							     }
						     }
						     
						     ## Run OLE
						     OLE_med_res <- OLE.test(dd = dates_median_ss, alpha = 0.05)
						     res_list[[row_id]] <- data.frame(Estimate = OLE_med_res[1],upperCI = OLE_med_res[2],lowerCI = OLE_med_res[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[i],sample_size = sample_size[h],simID_seed = z)
						     row_id <- row_id + 1
					     }
				     }
				     do.call(rbind, res_list)
			     },
			     mc.cores = ncores
)
OLE_medians <- do.call(rbind, OLE_medians_list)
})

saveRDS(OLE_medians,"../Results/OLE_medians_upl.rds")
saveRDS(time_OLE_medians,"../Results/time_OLE_medians_upl.rds")

## to save space
rm(OLE_medians)
rm(time_OLE_medians)

######## OLE RESAMPLE CALDATE (Our proposal)

rsmp <- 1000

time_OLE_resamp_caldate <- system.time({

OLE_resamp_list <- mclapply(
			    seq_len(sims),
			    function(z){
				    
				    set.seed(z)
				    
				    pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
				    
				    res_list <- list()
				    row_id <- 1
				    
				    for (h in seq_along(sample_size)){
					    dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
					    ndates <- sample_size[h]
					    st <- which(pre_dates_subset$sample_size == ndates)[1]
					    dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
					    
					    # Calibrate dates
					    cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
					    for (j in seq_along(dates_ss)){
						    temp_OLE <- matrix(NA, rsmp, 3)
						    
						    for (i in seq_len(rsmp)){
							    valid <- FALSE
							    while(!valid){
								    # Resample calibrated dates
								    resamp_dates <- vapply(seq_len(ndates),function(k){sample(cal_dates$grids[[k]][,1],size = 1,prob = cal_dates$grids[[k]][,2],replace = TRUE)},numeric(1))
								    resamp_dates <- sort(resamp_dates, decreasing = TRUE)
								    resamp_dates_temp <- resamp_dates[1:dates_ss[j]]
								    # Break ties
								    for (b in 2:length(resamp_dates_temp)){
									    if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
										    resamp_dates_temp[b] <- resamp_dates_temp[b]+1
									    }
								    }
								    OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
								    if (!any(is.na(OLE_res))){
									    temp_OLE[i,] <- OLE_res
									    valid <- TRUE
								    }
							    }
						    }
						    means <- colMeans(temp_OLE)
						    
						    res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
						    row_id <- row_id + 1
					    }
				    }
				    
				    do.call(rbind, res_list)
			    },
			    mc.cores = ncores
)
OLE_resamp_caldate <- do.call(rbind, OLE_resamp_list)
})

saveRDS(OLE_resamp_caldate,"../Results/OLE_resamp_caldate_upl.rds")
saveRDS(time_OLE_resamp_caldate,"../Results/time_OLE_resamp_caldate_upl.rds")

## to save space
rm(OLE_resamp_caldate)
rm(time_OLE_resamp_caldate)

######### OLE RESAMPLE NORM DIST (as in Key et al., 2021 and Bebber and Key 2022)

time_OLE_resamp_norm <- system.time({

OLE_resamp_norm_list <- mclapply(
				 seq_len(sims),
				 function(z){
					 
					 set.seed(z)
					 
					 pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
					 
					 res_list <- list()
					 row_id <- 1
					 
					 for (h in seq_along(sample_size)){
						 dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
						 ndates <- sample_size[h]
						 st <- which(pre_dates_subset$sample_size == ndates)[1]
						 dates_subset <- pre_dates_subset[st:(st + ndates - 1),]

						 # Calibrate dates
						 cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
						 
						 for (j in seq_along(dates_ss)){						 
							 temp_OLE <- matrix(NA, rsmp, 3)
							 
							 for (i in seq_len(rsmp)){
								 valid <- FALSE
								 while(!valid){
									 # Resample from normal distributions
									 resamp_dates <- vapply(seq_len(ndates),function(k){dat_mean <- mean(cal_dates$grids[[k]][,1])
												dat_sd <- (max(cal_dates$grids[[k]][,1]) - dat_mean)/2
												rnorm(1, dat_mean, dat_sd)},numeric(1))
									 
									 resamp_dates <- sort(resamp_dates, decreasing = TRUE)
									 resamp_dates_temp <- resamp_dates[1:dates_ss[j]]							 
									 # Break ties
									 for (b in 2:length(resamp_dates_temp)){
										 if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
											 resamp_dates_temp[b] <- resamp_dates_temp[b]+1
										 }
									 }
									 OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
									 if (!any(is.na(OLE_res))){
										 temp_OLE[i,] <- OLE_res
										 valid <- TRUE
									 }
								 }
							 }
							 
							 means <- colMeans(temp_OLE)
							 res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
							 row_id <- row_id + 1
						 }
					 }
					 
					 do.call(rbind, res_list)
				 },
				 mc.cores = ncores
)
OLE_resamp_norm <- do.call(rbind, OLE_resamp_norm_list)
})

saveRDS(OLE_resamp_norm,"../Results/OLE_resamp_norm_upl.rds")
saveRDS(time_OLE_resamp_norm,"../Results/time_OLE_resamp_norm_upl.rds")

## to save space
rm(OLE_resamp_norm)
rm(time_OLE_resamp_norm)

######### OLE RESAMPLE UNI DIST (Djakovic et al. 2022)

time_OLE_resamp_unif <- system.time({

OLE_resamp_unif_list <- mclapply(
				 seq_len(sims),function(z){
					 set.seed(z)
					 pre_dates_subset <- dates[(dates_seq*z-(dates_seq-1)):(dates_seq*z),]
					 
					 res_list <- list()
					 row_id <- 1
					 
					 for (h in seq_along(sample_size)){
						 dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
						 ndates <- sample_size[h]
						 
						 st <- which(pre_dates_subset$sample_size == ndates)[1]
						 dates_subset <- pre_dates_subset[st:(st + ndates - 1),]
						 
						 # Calibrate dates
						 cal_dates <- calibrate(dates_subset[,1], dates_subset[,2])
						 
						 for (j in seq_along(dates_ss)){
							 temp_OLE <- matrix(NA, rsmp, 3)
							 
							 for (i in seq_len(rsmp)){
								 valid <- FALSE
								 
								 while(!valid){
									 # Resample from uniform distributions
									 resamp_dates <- vapply(seq_len(ndates),function(k){
													runif(1,min(cal_dates$grids[[k]][,1]),max(cal_dates$grids[[k]][,1]))},numeric(1))
									 resamp_dates <- sort(resamp_dates, decreasing = TRUE)
									 resamp_dates_temp <- resamp_dates[1:dates_ss[j]]
									 
									 # Break ties
									 for (b in 2:length(resamp_dates_temp)){
										 if (resamp_dates_temp[b] == resamp_dates_temp[b-1]){
											 resamp_dates_temp[b] <- resamp_dates_temp[b]+1}
									 }
									 OLE_res <- as.numeric(OLE.test(dd = resamp_dates_temp, alpha = 0.05))
									 
									 if (!any(is.na(OLE_res))){
										 temp_OLE[i,] <- OLE_res
										 valid <- TRUE
									 }
								 }
							 }
							 
							 means <- colMeans(temp_OLE)
							 res_list[[row_id]] <- data.frame(Estimate = means[1],upperCI = means[2],lowerCI = means[3],start_date = dates_subset$start_date[1],r = dates_subset$r[1],Sd = dates_subset$Sd[1],ESS = dates_ss[j],sample_size = sample_size[h],simID_seed = z)
							 row_id <- row_id + 1
						 }
					 }
					 
					 do.call(rbind, res_list)
				 },
				 mc.cores = ncores
)
OLE_resamp_unif <- do.call(rbind, OLE_resamp_unif_list)
})

saveRDS(OLE_resamp_unif,"../Results/OLE_resamp_unif_upl.rds")
saveRDS(time_OLE_resamp_unif,"../Results/time_OLE_resamp_unif_upl.rds")

## to save space
rm(OLE_resamp_unif)
rm(time_OLE_resamp_unif)

######### Oxcal
## Arrange for sample sizes in Oxcal and CRIWM
sim_index <- expand.grid(dts = 1:(nrow(dates)/dates_seq), ESS = sample_size) ## For indexing sims


## Oxcal uniform
time_oxcal_unif <- system.time({
  oxcal_unif <- mclapply(
	1:nrow(sim_index),function(z) {
		dts_id <- sim_index$dts[z]
		ESS <- sim_index$ESS[z]
		
		pre_dates_subset <- dates[(dates_seq*dts_id - (dates_seq-1)):(dates_seq*dts_id),]
		dates_subset <- pre_dates_subset[pre_dates_subset$sample_size == ESS,]
		
		est <- oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "uniform")
		est$Sd <- unique(dates_subset$Sd)
	        est$start_date <- unique(dates_subset$start_date)
	        est$r <- unique(dates_subset$r)
	        est$sample_size <- nrow(dates_subset)
	        est$ESS <- ESS
	
	        return(est)
	},
	mc.cores = ncores
  )
})

saveRDS(oxcal_unif,"../Results/oxcal_unif_upl.rds")
saveRDS(time_oxcal_unif,"../Results/time_oxcal_unif_upl.rds")

## to save space
rm(oxcal_unif)
rm(time_oxcal_unif)

## Oxcal trapezoid
time_oxcal_trap <- system.time({
  oxcal_trap <- mclapply(
	1:nrow(sim_index),function(z) {
		dts_id <- sim_index$dts[z]
		ESS <- sim_index$ESS[z]
		
		pre_dates_subset <- dates[(dates_seq*dts_id - (dates_seq-1)):(dates_seq*dts_id),]
		dates_subset <- pre_dates_subset[pre_dates_subset$sample_size == ESS,]
		
		est <- oxcalRunner(c14age = round(dates_subset$Year),error  = dates_subset$Sd,model  = "trapezoid")
		est$Sd <- unique(dates_subset$Sd)
	        est$start_date <- unique(dates_subset$start_date)
	        est$r <- unique(dates_subset$r)
	        est$sample_size <- nrow(dates_subset)
	        est$ESS <- ESS
	
	        return(est)
	},
	mc.cores = ncores
  )
})

saveRDS(oxcal_trap,"../Results/oxcal_trap_upl.rds")
saveRDS(time_oxcal_trap,"../Results/time_oxcal_trap_upl.rds")

## to save space
rm(oxcal_trap)
rm(time_oxcal_trap)


######### CRIWM

# Combinations of (i, r, sd)
comb <- expand.grid(sim = 1:sim,r = r_vals,Sd = C14_errors, ESS = sample_size)

time_criwm <- system.time({

  criwm_res <- mclapply(1:nrow(comb), function(idx) {
    set.seed(idx)
    row <- comb[idx, ]
    
    # Build file name
    filepath <- paste0("../Simu_data/Uncal_YearsBP_sim_upl_",row$sim, "_r_",row$r,"_error_",row$Sd,"_sample_size_",row$ESS,".txt")

    # Run criwm
    res <- criwm(filepath, signor_lipps="arr")

    # Extract values
    r_v <- row$r
    Sd_v <- row$Sd
    ESS_v <- row$ESS
    data.frame(Estimate = res$criwm[2,2],upperCI = res$criwm[2,3],lowerCI = res$criwm[2,1],start_date = unique(dates$start_date)[comb$sim[idx]],r = unique(r_v),Sd = unique(Sd_v),ESS = unique(ESS_v),sample_size = unique(ESS_v),simID_seed = idx)
  }, mc.cores = ncores)
})

# Single data frame
CRIWM <- do.call(rbind, criwm_res)

saveRDS(CRIWM,"../Results/CRIWM_upl.rds")
saveRDS(time_criwm,"../Results/time_criwm_upl.rds")

## to save space
rm(CRIWM)
rm(time_criwm)
