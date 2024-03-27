OLE.CI = function(radioages_input) {

# radioages_input<-sort(replicate(100,rLogisticGrowth(a=5000, b=4000, k=0.01, r=0.01)), decreasing = FALSE)

radioages_input <- unique(radioages_input)
# length(radioages_input)
# Resampling Approach (accounts for chronological uncertainty but not for sampling error)

nsim  <- 1000

alpha <- 0.05

res2  <- numeric(length=nsim) #Place holder for results

radiocarbonages  <- calibrate(radioages_input,rep(1,length(radioages_input))) #Calibrating RC-Ages to Calender dates

grids <- radiocarbonages$grids #extracting grids of all dates' probabilities

distri_random <- vector('list',length=nsim) #place holder for looped results

for (i in 1:nsim)
{
  set.seed(i)
  distri_random[[i]]  <- lapply(grids,function(radioages_input){sample(radioages_input$calBP,size=1,prob=radioages_input$PrDens)}) |> unlist()
  distri_random[[i]]  <- distri_random[[i]] + runif(length(distri_random[[i]]))
  distri_random[[i]]  <- distri_random[[i]] - runif(length(distri_random[[i]]))
  distri_random[[i]]  <- sort(distri_random[[i]])
  distri_random[[i]] <- unique(distri_random[[i]])
  occurance <- rep(1, each = length(distri_random[[i]]))
  ole_df <- data.frame(distri_random[[i]], occurance)
  res2[i] <- OLE(ole_df,alpha)$Estimate
  
}

# hist(res2)

CIInterv <- quantile(res2,prob=c(0.025,0.975),na.rm = TRUE) |> round()

accuracy <- 0

# Define value
x1 <- 5000   

# Apply between function
range32 <- between(x1, CIInterv[1], CIInterv[2])

if (is.nan(range32) == "FALSE" & is.na(range32) == "FALSE"){

if (range32 == "TRUE")
  
{
  accuracy <- accuracy + 1
}

margin <-  CIInterv[2] - CIInterv[1]

distancefromstart <- 0

if (accuracy == 0){
  
  distancefromstart <- min(c(abs(5000 - CIInterv[2]),abs(5000 - CIInterv[1])))
  
}

new_row <- c(margin, accuracy, distancefromstart)

OLE.CI.results[nrow(OLE.CI.results) + 1, ] <<- new_row
 }
}
