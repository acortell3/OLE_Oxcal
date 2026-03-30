

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

ncores <- 15
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
					     # Select ESSs for OLE. Cannot be higher than sample size.
					     dates_ss <- dates_OLE[1:which(dates_OLE == sample_size[h])]
					     ## Sample size
					     ndates <- sample_size[h]
					     ## Index of the first dates that equals sample size
					     st <- which(pre_dates_subset$sample_size == ndates)[1]
					     ## Dates subset (it's correct)
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
								     dates_median_ss[b] <- dates_median_ss[b] +1 
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
