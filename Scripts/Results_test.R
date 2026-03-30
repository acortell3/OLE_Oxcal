

####### Script to plot and visualise the results of the paper <whateva>

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")


### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
#a###### HOLOCENE
## Load data
hol_data_full <- full_data[full_data$Period == "Holocene",]
hol_data_full <- hol_data_full[complete.cases(hol_data_full),]
hol_data_full$accuracy <- ifelse(hol_data_full$upperCI>hol_data_full$start_date & hol_data_full$lowerCI<hol_data_full$start_date,TRUE,FALSE)
nrow(hol_data_full)
hol_data_full$precision <- abs(hol_data_full$upperCI - hol_data_full$lowerCI)
n <- 80
k <- 40



hol_10_10 <- hol_data_full[hol_data_full$sample_size == n & hol_data_full$ESS == k,]

prop.table(table(hol_10_10$accuracy))
mean(hol_10_10$precision)

## Individual checks

## Load libraries
library(rcarbon)
library(oxcAAR) 
library(stringr)
library(Rextinct)
library(parallel)

## Oxcal setup
#quickSetupOxcal() #this should be replaced by actual reference to the file path where oxcal is downloaded (otherwise would install it every time)

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

## Values for resampling
sd_val <- 5365
r_val <- 0.03
Sd_val <- 20
ESS_val <- 5
ss_val <- 80
simID_seed_val <- 1017



date_check <- dates[dates$start_date == sd_val & dates$r == r_val & dates$Sd == Sd_val & dates$sample_size == ss_val,]

date_check_medians <- medCal(calibrate(date_check[,1],date_check[,2]))

date_check_medians5 <- sort(date_check_medians, decreasing = TRUE)[1:5]
date_check_medians5[3] <- date_check_medians5[2]+1
date_check_medians5[5] <- date_check_medians5[4]+1

res <- OLE.test(dd = date_check_medians5, alpha = 0.05)





