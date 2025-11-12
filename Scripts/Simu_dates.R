####### SCRIPT TO SIMULATE DATES USING A LOGISTIC FUNCTION

set.seed(123)

## Load libraries
library(nimbleCarbon)
library(rcarbon)

## Set pars for different simulations
C14_errors <- c(20,50)
r <- c(0.01,0.02,0.03,0.04)
start_date <- sample(1000:45000, 100 , replace = F)

for (i in 1:length(start_date)){
	for (j in 1:length(r)){
		for (k in 1:length(C14_errors)){
			modelPlot(dLogisticGrowth,a=start_date[i],b=start_date[i]-1000,params=c(k=0.1,r=r[j]),alpha=1,col='black') #the sampling distribution where k is the proportion of K (carrying capacity).
			x <- replicate(80,rLogisticGrowth(a=start_date[i],b=start_date[i]-1000,k=0.1,r=r[j])) #sample 80 dates from the distribution above in calendar time (cannot use n=80 as nimble is unhappy with that)
        		YearsBP <- uncalibrate(x,CRAerrors=C14_errors[k])[,c('rCRA','rError')]
			colnames(YearsBP) <- c("Year", "Sd")
    			saveRDS(YearsBP,paste0("../Simu_data/Uncal_YearsBP_sim_",i,"_r_",r[j],"_error_",C14_errors[k],".rds")) ## It purposedly overwrites the above
    			write.table(YearsBP,file = paste0("../Simu_data/Uncal_YearsBP_sim_",i,"_r_",r[j],"_error_",C14_errors[k],".txt"), sep = "\t", row.names = F) ## It purposedly overwrites the above
		}
	}
}

## Save start dates for later comparison
saveRDS(start_date,"../Simu_data/start_dates.rds")

