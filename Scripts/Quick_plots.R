
############### PLOTS

## Load packages
library(stringr)
library(rcarbon)

## Load data
sim_data <- readRDS("../sim_data.rds")
results <- readRDS("../results.rds")

## Run necessary functions
## Generate how many times the actual start is captured
is_correct <- function(x,up,lo){
  res <- as.factor(ifelse(x <= up & x >= lo, "hit", "miss"))
  return(res)
}

## Extracts the values from the OxCal output
extract_vals <- function(x){
  x <- x[3]
  upperCI <- as.numeric(str_extract_all(x,"\\d+")[[1]][3])
  lowerCI <- as.numeric(str_extract_all(x,"\\d+")[[1]][6])
  return(data.frame("upperCI" = upperCI, "lowerCI" = lowerCI))
}

## Assign values and objects
starts <- sim_data[[4]]$s
ndates <- sim_data[[4]]$n

OLE_med_res <- results[[1]]
OLE_ran_res <- results[[2]]
OxCal_tra <- lapply(results[[3]],extract_vals)
OxCal_tra <- do.call(rbind,OxCal_tra)
OxCal_uni <- lapply(results[[4]],extract_vals)
OxCal_uni <- do.call(rbind,OxCal_uni)
OxCal_gau <- lapply(results[[5]],extract_vals)
OxCal_gau <- do.call(rbind,OxCal_gau)

## Cal start for oxcal
cal_start <- calibrate(starts[1], error = 20) ## Select 20 sd
cal_start_med <- round(median(cal_start$grids$`1`$calBP)) ## median
cal_start_mea <- round(mean(cal_start$grids$`1`$calBP)) ## mean
cal_start_upCI <- hpdi(cal_start)[[1]][1,1]
cal_start_loCI <- hpdi(cal_start)[[1]][3,2]

## Compute how many times does it get it right
OLE_med_res$Correct <- is_correct(cal_start_med,OLE_med_res$upperCI,OLE_med_res$lowerCI)
OLE_ran_res$Correct <- is_correct(cal_start_med,OLE_ran_res$upperCI,OLE_ran_res$lowerCI)
OxCal_tra$Correct <- is_correct(cal_start_med,OxCal_tra$upperCI,OxCal_tra$lowerCI)
OxCal_uni$Correct <- is_correct(cal_start_med,OxCal_uni$upperCI,OxCal_uni$lowerCI)
OxCal_gau$Correct <- is_correct(cal_start_med,OxCal_gau$upperCI,OxCal_gau$lowerCI)

### Plot accuracy
### OLE
par(mfrow = c(2,1))
plot(ndates,rep(cal_start_med,length(ndates)), ylim = c(4000,10000), type = "l", col = "black", 
     main = "Accuracy OLE Medians", ylab = "Years BP", xlab = "Number of Dates")
polygon(x = c(rep(min(ndates),2),rep(max(ndates),2)), 
        y = c(cal_start_loCI,cal_start_upCI,cal_start_upCI,cal_start_loCI),
        border = FALSE, col = adjustcolor("gray", alpha.f = 0.3))
lines(ndates,rep(cal_start_upCI,length(ndates)), col = "gray")
lines(ndates,rep(cal_start_loCI,length(ndates)), col = "gray")
points(ndates,OLE_med_res$upperCI, pch = 16, col = adjustcolor("lightblue", alpha.f = 0.3))
points(ndates,OLE_med_res$lowerCI, pch = 16, col = adjustcolor("darkblue", alpha.f = 0.3))
legend("topright", legend = c("Upper confidence interval","Lower confidence interval"),
       col = c("lightblue", "darkblue"), pch = 16)

plot(ndates,rep(cal_start_med,length(ndates)), ylim = c(4000,10000), type = "l", col = "black", 
     main = "Accuracy OLE Random", ylab = "Years BP", xlab = "Number of Dates")
polygon(x = c(rep(min(ndates),2),rep(max(ndates),2)), 
        y = c(cal_start_loCI,cal_start_upCI,cal_start_upCI,cal_start_loCI),
        border = FALSE, col = adjustcolor("gray", alpha.f = 0.3))
lines(ndates,rep(cal_start_upCI,length(ndates)), col = "gray")
lines(ndates,rep(cal_start_loCI,length(ndates)), col = "gray")
points(ndates,OLE_ran_res$upperCI, pch = 16, col = adjustcolor("lightblue", alpha.f = 0.3))
points(ndates,OLE_ran_res$lowerCI, pch = 16, col = adjustcolor("darkblue", alpha.f = 0.3))
legend("topright", legend = c("Upper confidence interval","Lower confidence interval"),
       col = c("lightblue", "darkblue"), pch = 16)

### OxCal
par(mfrow = c(3,1))
plot(ndates,rep(cal_start_med,length(ndates)), ylim = c(4500,8000), type = "l", col = "black", 
     main = "Accuracy Oxcal Trapezoid", ylab = "Years BP", xlab = "Number of Dates")
polygon(x = c(rep(min(ndates),2),rep(max(ndates),2)), 
        y = c(cal_start_loCI,cal_start_upCI,cal_start_upCI,cal_start_loCI),
        border = FALSE, col = adjustcolor("gray", alpha.f = 0.3))
lines(ndates,rep(cal_start_upCI,length(ndates)), col = "gray")
lines(ndates,rep(cal_start_loCI,length(ndates)), col = "gray")
points(ndates,OxCal_tra$upperCI, pch = 16, col = adjustcolor("lightblue", alpha.f = 0.3))
points(ndates,OxCal_tra$lowerCI, pch = 16, col = adjustcolor("darkblue", alpha.f = 0.3))
legend("topright", legend = c("Upper confidence interval","Lower confidence interval"),
       col = c("lightblue", "darkblue"), pch = 16)

plot(ndates,rep(cal_start_med,length(ndates)), ylim = c(4500,8000), type = "l", col = "black", 
     main = "Accuracy Oxcal uniform", ylab = "Years BP", xlab = "Number of Dates")
polygon(x = c(rep(min(ndates),2),rep(max(ndates),2)), 
        y = c(cal_start_loCI,cal_start_upCI,cal_start_upCI,cal_start_loCI),
        border = FALSE, col = adjustcolor("gray", alpha.f = 0.3))
lines(ndates,rep(cal_start_upCI,length(ndates)), col = "gray")
lines(ndates,rep(cal_start_loCI,length(ndates)), col = "gray")
points(ndates,OxCal_uni$upperCI, pch = 16, col = adjustcolor("lightblue", alpha.f = 0.3))
points(ndates,OxCal_uni$lowerCI, pch = 16, col = adjustcolor("darkblue", alpha.f = 0.3))
legend("topright", legend = c("Upper confidence interval","Lower confidence interval"),
       col = c("lightblue", "darkblue"), pch = 16)

plot(ndates,rep(cal_start_med,length(ndates)), ylim = c(4500,8000), type = "l", col = "black", 
     main = "Accuracy Oxcal Gaussian", ylab = "Years BP", xlab = "Number of Dates")
polygon(x = c(rep(min(ndates),2),rep(max(ndates),2)), 
        y = c(cal_start_loCI,cal_start_upCI,cal_start_upCI,cal_start_loCI),
        border = FALSE, col = adjustcolor("gray", alpha.f = 0.3))
lines(ndates,rep(cal_start_upCI,length(ndates)), col = "gray")
lines(ndates,rep(cal_start_loCI,length(ndates)), col = "gray")
points(ndates,OxCal_gau$upperCI, pch = 16, col = adjustcolor("lightblue", alpha.f = 0.3))
points(ndates,OxCal_gau$lowerCI, pch = 16, col = adjustcolor("darkblue", alpha.f = 0.3))
legend("topright", legend = c("Upper confidence interval","Lower confidence interval"),
       col = c("lightblue", "darkblue"), pch = 16)


## Correctness plots
par(mfrow = c(2,1))
cdplot(ndates,OLE_med_res$Correct, col = c("lightblue", "darkblue"), main = "Correct OLE Median", 
       xlab = "Number of dates", ylab = "")
cdplot(ndates,OLE_ran_res$Correct, col = c("lightblue", "darkblue"), main = "Correct OLE Random", 
       xlab = "Number of dates", ylab = "")

par(mfrow = c(2,1))
cdplot(ndates,OxCal_tra$Correct, col = c("lightblue", "darkblue"), main = "Correct Oxcal Trapezoid", 
       xlab = "Number of dates", ylab = "")
cdplot(ndates,OxCal_uni$Correct, col = c("lightblue", "darkblue"), main = "Correct Oxcal Uniform", 
       xlab = "Number of dates", ylab = "")
#cdplot(ndates,OxCal_gau$Correct, col = c("lightblue", "darkblue"), main = "Correct Oxcal Gaussian", 
#       xlab = "Number of dates", ylab = "")


