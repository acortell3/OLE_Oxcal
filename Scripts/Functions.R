
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

## OLE returning v Weibull parameter
OLE.test.v <- function(dd, alpha){ ## dd is the dates and alpha is the confidence interval
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
  res<-data.frame(Estimate=extest, upperCI=upperCI, lowerCI=lowerCI,W.shape = v)
  return(res)	
}

### Functions to fit oxcal
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

## Function to extract Oxcal results
oxcalReadjs <- function(x,model,path="")
{
  oxcalinput <- readLines(paste0(path,model,".oxcal"))
  oxcaloutput <- readLines(paste0(path,model,".js"))
  ovAg=as.numeric(strsplit(oxcaloutput[grep("model.overallAgreement",oxcaloutput)],"[=,;]")[[1]][2])
  
  print("Extracting data from oxcal output...")
  
  ii = which(is.na(x$grp)) #dates not part of a R_Combine Group
  
  pb <- txtProgressBar(min=1, max=length(ii), style=3)
  
  x$agreement = NA
  x$combine.agreement = NA
  x$pval = NA
    
  for (i in 1:length(ii))
  {
    setTxtProgressBar(pb, i)
    index=grep(paste0(x$id[ii[i]]," Posterior"),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".posterior.agreement"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    x$agreement[ii[i]]=as.numeric(strsplit(tmp2,"[=,;]")[[1]][2])	
  }
  close(pb)
  
  combineGroups = unique(x$grp)
  combineGroups = combineGroups[which(!is.na(combineGroups))]
  pb <- txtProgressBar(min=1, max=length(combineGroups), style=3)
  
  
  for (i in 1:length(unique(combineGroups)))
  {
    setTxtProgressBar(pb, i)
    #Extract Agreement
    index=grep(paste0('"Combination',combineGroups[i],' Posterior "'),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".posterior.agreement"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    x$combine.agreement[which(x$grp==combineGroups[i])] = as.numeric(strsplit(tmp2,"[=,;]")[[1]][2])	
    
    #Extract X-square Test
    index=grep(paste0('"Combination',combineGroups[i],'"'),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".likelihood.testText"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    df = as.numeric(gsub('^.*df=\\s*|\\s*T.*$', '', tmp2))
    t  = as.numeric(gsub('^.*T=\\s*|\\s*\\(.*$', '', tmp2))
    x$pval[which(x$grp==combineGroups[i])]=1-pchisq(t,df)
  }
  close(pb)
  
  
  print("Done.")
  return(list(ovAg=ovAg,df=x))
}


