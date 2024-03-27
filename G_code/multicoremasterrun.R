library(oxcAAR) #Load oxcAAR library
library(rcarbon) #to handle 14C dates
library(parallel)
library('sExtinct')
library('nimbleCarbon')
library('dplyr')

source('oxcalScript.R')
source('randomgen.R')
source('oxcalcleaner.R')
source('oxcalcleanerunif.R')
source('cleanthree.R')
source('cleanthreeunif.R')
source('cleanfour.R')
source('cleanfourunif.R')
source('cleanfive.R')
source('cleanfiveunif.R')
source('cleansix.R')
source('cleansixunif.R')
source('OLE_Median.R')
source('OLE_CI.R')

# Generate a list of n sample sets to run in parallel

n  <- 100 #loopruns
# datalist  <- vector('list',length=n)
# nsamples <- c(10,50,100)
nsamples  <- 10 #number of dates
true.start  <- 5000 
true.end  <- 3000
# rateofadoption <- c(0.01, 0.05, 0.1)
rateofadoption <- 0.01
error <- 20

quickSetupOxcal(path = "gugan/Documents/Diss_Code")

unif_df <- data.frame(Precision = numeric(),Accuracy = numeric(), Distance_from_start = numeric())
trap_df <- data.frame(Precision = numeric(),Accuracy = numeric(), Distance_from_start = numeric())
OLE.Median.results <- data.frame(Precision = numeric(),Accuracy = numeric(), Distance_from_start = numeric())
OLE.CI.results <- data.frame(Precision = numeric(),Accuracy = numeric(), Distance_from_start = numeric())

mainrun <- function(noofsamples,start,end,adoptionrate,errorrate){
  radiocarbon.ages<<- list()
  for (x in 1:n) {
    radiocarbon.ages<<-c(radiocarbon.ages,list(sort(sim_dates(noofsamples,start,end,adoptionrate,0.01))))
  }
  error.dates <<- rep(errorrate,nsamples) #20 years error for each date
}

genscript <- function (){
i <- 1
for (dates in radiocarbon.ages)
{
  # print(i)
  oxcalScriptGen(c14age=dates,errors=20,fn=paste0('trap_set_',i,'.oxcal'),model='trapezoid')
  oxcalScriptGen(c14age=dates,errors=20,fn=paste0('unif_set_',i,'.oxcal'),model='uniform')
  i <- i + 1
  }
}


runFun  <- function(x)
{
  library(oxcAAR)
  source('oxcalScript.R')
  quickSetupOxcal(path = "gugan/Documents/Diss_Code")
  oxScript.trap  <- readLines(paste0('trap_set_',x,'.oxcal'))
  res.file.trap <- executeOxcalScript(oxScript.trap)
  res.trap <- readOxcalOutput(res.file.trap)
  results  <- res.trap
  ii=grep("ocd\\[4\\].posterior.comment",results)
  rngtrap = results[ii]
}

runFun1  <- function(x)
{
  library(oxcAAR)
  source('oxcalScript.R')
  quickSetupOxcal(path = "gugan/Documents/Diss_Code")
  oxScript.unif  <- readLines(paste0('unif_set_',x,'.oxcal'))
  res.file.unif <- executeOxcalScript(oxScript.unif)
  res.unif <- readOxcalOutput(res.file.unif)
  results  <- res.unif
  ii=grep("ocd\\[3\\].posterior.comment",results)
  rng = results[ii]
}

##################
#Code Starts Here#
##################

mainrun(nsamples,true.start,true.end,rateofadoption,error)

genscript()

cl  <- makeCluster(8)
x  <- 1:n
out  <- parLapply(cl=cl,X=x,fun=runFun)
out1  <- parLapply(cl=cl,X=x,fun=runFun1)
stopCluster(cl)

# out is a list where in each slot you have the oxcal output of your model. For example you could extract from the first slot the start of the phase:

for (trap in out){
    oxcalcleaner(trap)
}
write.csv(trap_df, (paste0("sample_",(paste0(nsamples,(paste0("_rate_",(paste0(rateofadoption,"_Trap_OxCal.csv")))))))))

for (unif in out1){
  oxcalcleanerunif(unif)
}
write.csv(unif_df, (paste0("sample_",(paste0(nsamples,(paste0("_rate_",(paste0(rateofadoption,"_Unif_OxCal.csv")))))))))


for (dates in radiocarbon.ages)
{
  OLE.Median(dates, error)
  OLE.CI(dates)
}
write.csv(OLE.Median.results, (paste0("sample_",(paste0(nsamples,(paste0("_rate_",(paste0(rateofadoption,"_OLE.Median.results.csv")))))))))
write.csv(OLE.CI.results, (paste0("sample_",(paste0(nsamples,(paste0("_rate_",(paste0(rateofadoption,"_OLE.CI.results.csv")))))))))

# rm(unif_df)
# rm(trap_df)
# rm(OLE.true_estimate.median.results)
# rm(OLE.CI.results)


