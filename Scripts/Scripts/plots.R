##### Just some preliminary plots

## Load librares
library(Rmisc) ## for CI
library(rcarbon)
library(coda)

## Load and prepare data
NC_trapezoid <- readRDS("../Results/NC_trapezoid.rds")
NC_uniform <- readRDS("../Results/NC_uniform.rds")
OLE_medians <- read.csv("../Results/OLE_medians.csv")[,-1]
OLE_resamp <- t(read.csv("../Results/OLE_resamp_res.csv"))
OLE_resamp <- OLE_resamp[-1,]
OLE_resamp <- as.numeric(OLE_resamp)
names(OLE_resamp) <- names(OLE_medians)

mean_NC_tr <- mean(NC_trapezoid$samples$chain1[,'a'])
NC_trapezoid <- HPDinterval(NC_trapezoid$samples$chain1[,'a'])

mean_NC_un <- mean(NC_uniform$samples$chain1[,'alpha[2]'])
NC_uniform <- HPDinterval(NC_uniform$samples$chain1[,'alpha[2]'])

## Get CIs for Oxcal (from mean and sd extracted from Oxcal results)
Oxcal_uniform <- c("mean" = 4117+1950, "upper" = 4209+1950, "lower" = 4072+1950)
Oxcal_trapezoid <- c("mean" = 4074+1950, "upper" = 4288+1950, "lower" = 3957+1950)

## Merge everything into a single dataset
results <- data.frame("mean" = c(as.integer(OLE_medians[1]),as.integer(OLE_resamp[1]),mean_NC_un, mean_NC_tr,Oxcal_uniform[1],Oxcal_trapezoid[1]),
		      "upperCI" = c(as.integer(OLE_medians[2]), as.integer(OLE_resamp[2]), NC_uniform[2], NC_trapezoid[2], Oxcal_uniform[2], Oxcal_trapezoid[2]), 
		      "lowerCI" = c(as.integer(OLE_medians[3]), as.integer(OLE_resamp[3]), NC_uniform[1], NC_trapezoid[1], Oxcal_uniform[3], Oxcal_trapezoid[3]))

rownames(results) <- c("OLE_medians", "OLE_resamp", "NC_uniform", "NC_trapezoid", "Oxcal_uniform", "Oxcal_trapezoid")

## Target dates
# Uncalibrated dates
effective_oldest_date <- c(6500,40) ## Sd is arbitrary
sampled_oldest_date <- c(6306,60) ## As extracted from the script Simu_dates.R

# back-calibration
#uncalibrate
uncal_eod <- uncalibrate(effective_oldest_date[1])
uncal_sod <- uncalibrate(sampled_oldest_date[1])

#back-calibrate
cal_eod <- calibrate(uncal_eod$ccCRA, errors = effective_oldest_date[2])
cal_sod <- calibrate(uncal_sod$ccCRA, error = sampled_oldest_date[2])

### Generate means and CIs for the observed dates
cal_eod_stats <- c("median" = unlist(summary(cal_eod)[2]), "upper" = unlist(as.numeric(sub(" .*","",summary(cal_eod)[5]))), "lower" = unlist(as.numeric(sub(".*to","",summary(cal_eod)[6])))) 

cal_sod_stats <- c("median" = unlist(summary(cal_sod)[2]), "upper" = unlist(as.numeric(sub(" .*","",summary(cal_sod)[5]))), "lower" = unlist(as.numeric(sub(".*to","",summary(cal_sod)[6])))) 

#### Let's go with the plots
## prepare the colors for the distributions
light_cols <- c("gold", "brown1","cadetblue1","deepskyblue","seagreen1","indianred1")
dark_cols <- c("gold4","brown4","cadetblue4","deepskyblue4","seagreen4","indianred4")

plot(x = c(0,7), y = c(cal_eod_stats[1],cal_eod_stats[1]), type = "l", col = "red", ylim = c(4500,7000), ylab = "Age", xlab = "Method", xaxt = "n")
axis(1, at = 1:6, labels = c("OLE med","OLE res",
			     "NC uni","NC tra",
			     "Oxc uni","Oxc tra"))
legend("bottomright", legend = c("real oldest date","oldest C14 sample"), lty = c(1,1), col = c("red","blue"))
lines(x = c(0,7), y = c(cal_eod_stats[2],cal_eod_stats[2]), lty = 2, col = "pink")
lines(x = c(0,7), y = c(cal_eod_stats[3],cal_eod_stats[3]), lty = 2, col = "pink")

lines(x = c(0,7), y = c(cal_sod_stats[1],cal_sod_stats[1]), lty = 1, col = "blue")
lines(x = c(0,7), y = c(cal_sod_stats[2],cal_sod_stats[2]), lty = 2, col = "lightblue")
lines(x = c(0,7), y = c(cal_sod_stats[3],cal_sod_stats[3]), lty = 2, col = "lightblue")

for (i in 1:6){
	polygon(x = c(i-0.05,i+0.05,i+0.05,i-0.05), y = c(results[i,2],results[i,2],results[i,3], results[i,3]), col = light_cols[i])
	lines(x = c(i-0.1, i+0.1), y = c(results[i,1], results[i,1]), lwd = 4, col = dark_cols[i])
}

