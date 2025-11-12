
########## SCRIPT COMBINING ALL THE METHODS PROPOSED

## Load libraries
library(rcarbon)
library(oxcAAR) 
library(stringr)
library(Rextinct)

## Load required functions
source('oxcal_functions.R')
source('functions.R')

## Load data
nfiles <- 1600/2 ## I know I have 1600 usable files in my Simu_data folder, 800 rds and 800 txt

## Need this to recover files by name. They match the ones in Simu_dates.R
sim <- 100
C14_errors <- c(20,50)
r <- c(0.01,0.02,0.03,0.04)

for (h in 1:nfiles){
	for (i in 1:sim){
		for (k in 1:length(r)){
			for (j in 1:length(C14_errors)){
				dates  <- readRDS(paste0('../Simu_data/Uncal_YearsBP_sim_',i,'_r_',r[k],'_error_',C14_errors[j],'.rds'))
				colnames(dates) <- c("Year","Sd")
				dates_text  <- read.table(paste0('../Simu_data/Uncal_YearsBP_sim_',i,'_r_',r[k],'_error_',C14_errors[j],'.txt'))
			}
		}
	}
}


## Back-calibrate dates (for OLE)
cal_dates <- calibrate(dates$Year, dates$Sd)

## Arrange for sample sizes in OLE
dates_ss <- c(5,10,round(nrow(dates)/2),nrow(dates))

OLE_medians <- list("d5" = NA, # To store results on medians
		    "d10" = NA,
		    "half" = NA,
		    "full" = NA)

OLE_resample <- list("d5" = NA, # To store results on resample
		     "d10" = NA,
		     "half" = NA,
		     "full" = NA)
## OLE medians
dates_median <- rep(NA,nrow(dates))
for (i in 1:nrow(dates)){
	dates_median[i] <- medCal(calibrate(dates[i,1],dates[i,2]))
}

## sample dates acccording to sample size
for (i in 1:length(dates_ss)){
	dates_median_temp <- sample(dates_median,dates_ss[i])
	OLE_medians[[i]] <- OLE.test(dd = dates_median_temp, alpha = 0.05)
}

## OLE resample
rsmp <- 100 ## Number of resamples

for (j in 1:length(dates_ss)){
	## Create object to store with given sample size
	temp_object <- data.frame("Estimate" = rep(NA,rsmp),
				  "upperCI" = rep(NA,rsmp),
				  "lowerCI" = rep(NA,rsmp))
	for (k in 1:rsmp){
		## Resample dates
		date_vec <- rep(NA,nrow(dates))
		for (i in 1:nrow(dates)){
			date_vec[i] <- sample(cal_dates$grids[[i]][,1],1,cal_dates$grids[[i]][,2], replace = T)
		}
	
		## Subset according to sample sizes
		sampled_dates <- sample(date_vec,dates_ss[j])
		temp_object[k,] <- OLE.test(dd = sampled_dates, alpha = 0.05)
	}

	OLE_resample[[j]] <- temp_object
}

## CRIWM
criwm_res <- criwm("dates_text.txt")

## Oxcal
quickSetupOxcal() #this should be replaced by actual reference to the file path where oxcal is downloaded (otherwise would install it every time)
unif.res  <- oxcalRunner(c14age=round(x$Year),error=x$Sd,model='uniform')
trap.res  <- oxcalRunner(c14age=round(x$Year),error=x$Sd,model='trapezoid')
