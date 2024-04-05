
############# FUNCTIONS USED ####################

######## OLE calculation function taken from Key et al. 2024 since sExtinct has been retracted from CRAN
OLE.test <- function(dd, alpha){ ## dd is the dates and alpha is the confidence interval
  # records are sorted in a reverse order, as required by OLE method
  sights <- rev(sort(dd))
  # calculation of k, v, e, lambda and other values
  k <- length(sights)
  v <- (1/(k-1)) * sum(log((sights[1] - sights[k])/(sights[1] - sights[2:(k-1)])))
  e <- matrix(rep(1,k), ncol=1)
  SU<-(-log(alpha)/length(sights))^-v
  myfun <- function(i,j,v){(gamma(2*v+i)*gamma(v+j))/(gamma(v+i)*gamma(j))}
  lambda <- outer(1:k, 1:k, myfun, v=v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda)) 
  a <- as.vector(solve(t(e)%*%solve(lambda)%*%e)) * solve(lambda)%*%e
  # calculation of CI ("upperCI") and extinction time ("extest")
  upperCI<-max(sights) + ((max(sights)-min(sights))/(SU-1))
  lowerCI<-min(sights) + ((min(sights)-max(sights))/(SU-1))
  extest<-sum(t(a)%*%sights)
  # return of results produced by the function
  res<-data.frame(Estimate=extest, upperCI=upperCI, lowerCI=lowerCI)
  return(res)	
}

oxcalScriptGen = function(c14age,errors,fn,model=c("gaussian","uniform","trapezoid"))
{
  id  <- paste0('id',1:length(c14age))	
  export <- file(fn) #create export file
  cat("Options(){Resolution=5;SD1=FALSE;SD2=TRUE;BCAD=FALSE;};\n",file=fn,append=FALSE) #Start Sequence#
  cat("Plot(){\n",file=fn,append=TRUE) #Start Sequence#
  cat("Phase(){\n",file=fn,append=TRUE) #Start Sequence#
  
  # Phase Boundary Start####
  cat("Sequence(){\n", file = fn, append = TRUE) #Start Sequence#
  if (model=="gaussian")
  {
    cat(paste0('Sigma_Boundary("Start Phase");\n'), file = fn, append = TRUE)
  }
  if (model=="uniform")
  {
    cat(paste0('Boundary("Start Phase");\n'), file = fn, append = TRUE)
  }
  if (model=="trapezoid")
  {
    cat(paste0('Boundary("Start Phase"){\n'), file = fn, append = TRUE)
    cat(paste0('Start("Start of Start Phase");\n'), file = fn, append = TRUE)
    cat(paste0('Transition("Period of Start Phase");\n'), file = fn, append = TRUE)
    cat(paste0('End("End of Start Phase");\n'), file = fn, append = TRUE)
    cat('};\n', file = fn, append = TRUE)
  }
  
  
  # Actual Dates#####
  cat(paste0('Phase("Phase")\n'), file = fn, append = TRUE)
  cat('{\n', file = fn, append = TRUE)
  
  #start with dates to be combined
  
  for (i in 1:length(c14age))
  {
    cat(paste('R_Date(','\"',id[i],'\",',c14age[i],',',errors[i],');\n', sep = ""),file = fn, append = TRUE)
  } 
  
  cat('};\n',file=fn,append=TRUE)
  
  # Phase Boundary End####
  if (model=="gaussian")
  {
    cat(paste0('Sigma_Boundary("End Phase");\n'), file = fn, append = TRUE)
  }
  if (model=="uniform")
  {
    cat(paste0('Boundary("End Phase");\n'), file = fn, append = TRUE)
  }
  if (model=="trapezoid")
  {
    cat(paste0('Boundary("End Phase"){\n'), file = fn, append = TRUE)
    cat(paste0('Start("Start of End Phase");\n'), file = fn, append = TRUE)
    cat(paste0('Transition("Period of End Phase");\n'), file = fn, append = TRUE)
    cat(paste0('End("End of End Phase");\n'), file = fn, append = TRUE)
    cat('};\n', file = fn, append = TRUE)
  }
  cat('};\n', file = fn, append = TRUE)
  cat('};\n', file = fn, append = TRUE) 
  cat('};\n', file = fn, append = TRUE) 
  close(export)
}



