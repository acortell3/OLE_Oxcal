####### SCRIPT TO SIMULATE DATES USING A LOGISTIC FUNCTION

set.seed(123)

## Load libraries
library(nimbleCarbon)
library(rcarbon)

## Set pars for different simulations
C14_errors <- c(20,50)
r <- c(0.01,0.03,0.06)

########## Younger dates HOLOCENE
start_date <- sample(1000:11000, 1000, replace = F)
ndates <- 80
i_row <- 1 # To update df row index
ESS_criwm <- c(10,ndates/2,ndates)


## Create data frame to store simulations
simuls <- data.frame("YearsBP" = NA,
		     "Sd" = NA,
		     "start_date" = NA,
		     "r" = NA)
## Generate dates
for (i in 1:length(start_date)){
	for (j in 1:length(r)){
		for (k in 1:length(C14_errors)){
			modelPlot(dLogisticGrowth,a=start_date[i],b=start_date[i]-1000,params=c(k=0.1,r=r[j]),alpha=1,col='black') #the sampling distribution where k is the proportion of K (carrying capacity).
			x <- replicate(ndates,rLogisticGrowth(a=start_date[i],b=start_date[i]-1000,k=0.1,r=r[j])) #sample 80 dates from the distribution above in calendar time (cannot use n=80 as nimble is unhappy with that)
			YearsBP <- uncalibrate(x,CRAerrors=C14_errors[k])[,c('rCRA','rError')]
			simuls[c(i_row:(i_row+ndates-1)),] <- c(YearsBP[1],YearsBP[2],as.data.frame(rep(start_date[i],ndates)),as.data.frame(rep(r[j],ndates)))
			i_row <- 1+nrow(simuls) ## Update the row to store the information
    			# Need to save it like this for CWIRM
			for (h in 1:length(ESS_criwm)){
				YearsBP_ESS <- YearsBP[sample(nrow(YearsBP),ESS_criwm[h]),]
				write.table(YearsBP_ESS,file = paste0("../Simu_data/Uncal_YearsBP_sim_hol_",i,"_r_",r[j],"_error_",C14_errors[k],"_ESS_",ESS_criwm[h],".txt"), sep = "\t", row.names = F) ## It purposedly overwrites the above
			}
		}
	}
}

saveRDS(simuls,"../Simu_data/simuls_hol.rds")

##### OLDER UPPER PALAEOLITHIC

## Set pars for different simulations
start_date <- sample(30000:40000, 1000, replace = F)
i_row <- 1 # To update df row index

## Create data frame to store simulations
simuls <- data.frame("YearsBP" = NA,
		     "Sd" = NA,
		     "start_date" = NA,
		     "r" = NA)
## Generate dates
for (i in 1:length(start_date)){
	for (j in 1:length(r)){
		for (k in 1:length(C14_errors)){
			modelPlot(dLogisticGrowth,a=start_date[i],b=start_date[i]-1000,params=c(k=0.1,r=r[j]),alpha=1,col='black') #the sampling distribution where k is the proportion of K (carrying capacity).
			x <- replicate(ndates,rLogisticGrowth(a=start_date[i],b=start_date[i]-1000,k=0.1,r=r[j])) #sample 80 dates from the distribution above in calendar time (cannot use n=80 as nimble is unhappy with that)
			YearsBP <- uncalibrate(x,CRAerrors=C14_errors[k])[,c('rCRA','rError')]
			simuls[c(i_row:(i_row+ndates-1)),] <- c(YearsBP[1],YearsBP[2],as.data.frame(rep(start_date[i],ndates)),as.data.frame(rep(r[j],ndates)))
			i_row <- 1+nrow(simuls) ## Update the row to store the information
    			# Need to save it like this for CWIRM
			for (h in 1:length(ESS_criwm)){
				YearsBP_ESS <- YearsBP[sample(nrow(YearsBP),ESS_criwm[h]),]
				write.table(YearsBP_ESS,file = paste0("../Simu_data/Uncal_YearsBP_sim_upl_",i,"_r_",r[j],"_error_",C14_errors[k],"_ESS_",ESS_criwm[h],".txt"), sep = "\t", row.names = F) ## It purposedly overwrites the above
			}
		}
	}
}

saveRDS(simuls,"../Simu_data/simuls_upl.rds")


