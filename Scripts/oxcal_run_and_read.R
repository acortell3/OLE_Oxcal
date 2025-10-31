library(oxcAAR) 
library(stringr)
source('oxcal_functions.R')
x  <- readRDS('../Simu_data/Uncal_YearsBP_100dates.rds')
quickSetupOxcal() #this should be replaced by actual reference to the file path where oxcal is downloaded (otherwise would install it every time)
unif.res  <- oxcalRunner(c14age=round(x$Year),error=x$Sd,model='uniform')
trap.res  <- oxcalRunner(c14age=round(x$Year),error=x$Sd,model='trapezoid')

