###### CREATE OXCAL SCRIPTS TO SEND TO OXCAL

## Load dates
dates <- readRDS("../Simu_data/Uncal_YearsBP_100dates.rds")
colnames(dates) <- c("Year", "Sd")

## Load functions
source("functions.R")

## Load packages

### Generate oxcal scripts
oxcalScriptGen(c14age = round(dates$Year), errors = dates$Sd, fn = "./Oxcal_scripts/uniform.oxcal", model = "uniform")

oxcalScriptGen(c14age = round(dates$Year), errors = dates$Sd, fn = "./Oxcal_scripts/trapezoid.oxcal", model = "trapezoid")








