

####### Script to plot and visualise the results of the paper <whateva>

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")


### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE
## Load data
hol_data_full <- full_data[full_data$Period == "Holocene",]
n <- 40
k <- 10


hol_10_10 <- hol_data_full[hol_data_full$sample_size == n & hol_data_full$ESS == k,]

round(prop.table(table(hol_10_10$start_date < hol_10_10$upperCI & hol_10_10$start_date > hol_10_10$lowerCI))[2],3)






