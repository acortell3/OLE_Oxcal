
############# FUNCTIONS USED ####################

#a####### OLE calculation function taken from Key et al. 2024 since sExtinct has been retracted from CRAN
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
  
  ## Correction to avoid matrix singularity without stopping the loop
  if(any(!is.finite(lambda)) || qr(lambda)$rank < ncol(lambda)){
	  res <- NA
	  return(res)
  } else {
	  a <- as.vector(solve(t(e)%*%solve(lambda)%*%e)) * solve(lambda)%*%e
  	  # calculation of CI ("upperCI") and extinction time ("extest")
  	  upperCI<-max(sights) + ((max(sights)-min(sights))/(SU-1))
  	  lowerCI<-min(sights) + ((min(sights)-max(sights))/(SU-1))
  	  extest<-sum(t(a)%*%sights)
  	  # return of results produced by the function
  	  res<-data.frame(Estimate=extest, upperCI=upperCI, lowerCI=lowerCI)
  	  return(res)
  }
}

### OXCAL FUNCTIONS

oxcalRunner <- function(c14age,errors,fn=tempfile(fileext='oxcal'),model=c("uniform","trapezoid"))
{
	require(rcarbon)
	require(oxcAAR)

	# Routine for generating and storing oxcal script on fn
	id  <- paste0('id',1:length(c14age))	
	export <- file(fn) #create export file
	cat("Options(){Resolution=5;SD1=FALSE;SD2=TRUE;BCAD=FALSE;};\n",file=fn,append=FALSE) #Start Sequence#
	cat("Plot(){\n",file=fn,append=TRUE) #Start Sequence#
	cat("Phase(){\n",file=fn,append=TRUE) #Start Sequence#

	# Phase Boundary Start####
	cat("Sequence(){\n", file = fn, append = TRUE) #Start Sequence#
	if (model=="uniform")
	{
		cat(paste0('Boundary("Target");\n'), file = fn, append = TRUE)
	}
	if (model=="trapezoid")
	{
		cat(paste0('Boundary("Start Phase"){\n'), file = fn, append = TRUE)
		cat(paste0('Start("Target");\n'), file = fn, append = TRUE)
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
	modelscript <- readLines(fn)

	### Run OxCal ----
	na.check <- TRUE #OxCal for some reason does occasionaly yield '...' in the posterior range. This happens randomly, so the whle loop forces a re-run of the script to avoid reporting NA.
	na.rerun <- -1 #this counts how many times the while loop was used to rerun. If na.rerun is 0, it means that there were no issues.
	while(na.check)
	{
	na.rerun <- na.rerun + 1
	fit <- executeOxcalScript(modelscript) # executes the oxcal script, loc is where the output is stored
	oxcaloutput <- readLines(fit) #read output


	### Retreive Key Infor from output 
	# overall agreement index
	ovAg=as.numeric(strsplit(oxcaloutput[grep("model.overallAgreement",oxcaloutput)],"[=,;]")[[1]][2])

	# retrieve agreement index for each date
	res <- data.frame(CRA=c14age,CRAError=errors)
	res$agreement  <- NA
	    
	for (i in 1:nrow(res))
	{
		index=grep(paste0('id',i," Posterior"),oxcaloutput)
		tmp=oxcaloutput[index]	
		position=strsplit(tmp,"[.]")[[1]][1]
		position=gsub("\\[","\\\\[",position)
		position=gsub("\\]","\\\\]",position)
		index2=grep(paste0(position,".posterior.agreement"),oxcaloutput)
		tmp2=oxcaloutput[index2]
		res$agreement[i]=as.numeric(strsplit(tmp2,"[=,;]")[[1]][2])	
	}


	# retrieve 2-sigma posterior range for start date
	target <- grep('Target',oxcaloutput)
	# find out ocd number of target parameter (arrival date)
	ocd  <- sub("\\..*$","",oxcaloutput[target[1]])
	ocd  <- paste0("^",gsub("(\\[|\\])","\\\\\\1",ocd))
	posranges  <- grep(paste0(ocd,'.posterior.range\\[2\\]'),oxcaloutput)

	# Utility function for extracting key values
	extract_nums_after_equals <- function(x) {
		# find text inside [...] after the first '='
		m <- regexec("=\\s*\\[([^]]*)\\]", x)
		caps <- regmatches(x, m)

		lapply(caps, function(z) {
			       # if no match, return numeric(0)
			       if (length(z) < 2) return(numeric(0))

			       inside <- trimws(z[2])
			       if (inside == "") return(numeric(0))  # empty brackets

			       as.numeric(strsplit(inside, "\\s*,\\s*")[[1]])
})
	}

	post <- extract_nums_after_equals(oxcaloutput[posranges])
	post.df  <- do.call(rbind,lapply(post,function(y){if(length(y)==3) y else NULL})) |> as.data.frame()
	post.df[,1] <- BCADtoBP(post.df[,1])
	post.df[,2] <- BCADtoBP(post.df[,2])
	colnames(post.df)  <- c('Start','End','P')
	if (all(!is.na(post.df$Start))&all(!is.na(post.df$End))){na.check <- FALSE} #NA check
	}
	return(list(na.rerun=na.rerun,overallAgreement=ovAg,individualAgreement=res,posteriorRange=post.df))
}

