
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
	return(list(overallAgreement=ovAg,individualAgreement=res,posteriorRange=post.df))
}
