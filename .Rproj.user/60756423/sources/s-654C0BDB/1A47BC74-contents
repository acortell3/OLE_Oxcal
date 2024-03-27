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
