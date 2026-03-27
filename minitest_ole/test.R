library(nimbleCarbon)
library(rcarbon)
library(dplyr)
source('../Scripts/Functions.R')

nsim <- 500
r <- 0.01
cra.error  <- 50
n <- c(10,40,80)
k <- c(10,40,80)

res <- data.frame(sim=NA,start.date=NA,n=NA,k=NA,sort=NA,type=NA,est=NA,lo=NA,hi=NA,hit=NA)

for (i in 1:nsim)
{
	print(i)
	set.seed(i)
	#simulate date
	true.start <- sample(1000:11000,1,size=1)
	for (j in 1:length(n))
	{
		calendar.dates  <- replicate(n[j],rLogisticGrowth(a=true.start,b=true.start-1000,k=0.1,r=r))
		c14.dates  <- uncalibrate(calendar.dates)$rCRA
		med.dates <- medCal(calibrate(c14.dates,rep(cra.error,n[j]),verbose=F))
		calendar.dates  <- calendar.dates + runif(n[j],min=0.1,max=0.5)
		med.dates  <- med.dates + runif(n[j],min=0.1,max=0.5)
		calendar.dates  <- sort(calendar.dates,decreasing=T)
		med.dates  <- sort(med.dates,decreasing=T)

		for (z in 1:length(k))
		{
			if (n[j]>k[z])
			{
				#if n>k
				#earliest sample approach
				calendar.dates.early  <- calendar.dates[1:k[z]]
				med.dates.early  <- med.dates[1:k[z]]
				res.cal <- OLE.test(dd=calendar.dates.early,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort='yes',type='calendar',est=res.cal$Estimate,lo=res.cal$upperCI,hi=res.cal$lowerCI,hit=ifelse(res.cal$upperCI>true.start & res.cal$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
				res.c14 <- OLE.test(dd=med.dates.early,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort='yes',type='c14med',est=res.c14$Estimate,lo=res.c14$upperCI,hi=res.c14$lowerCI,hit=ifelse(res.c14$upperCI>true.start & res.c14$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
				#random sub-sample approach
				calendar.dates.rand  <- sample(calendar.dates,size=k[z],replace=F) |> sort(decreasing=TRUE)
				med.dates.rand  <- sample(calendar.dates,size=k[z],replace=F) |> sort(decreasing=TRUE)

				res.cal <- OLE.test(dd=calendar.dates.rand,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort='rand',type='calendar',est=res.cal$Estimate,lo=res.cal$upperCI,hi=res.cal$lowerCI,hit=ifelse(res.cal$upperCI>true.start & res.cal$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
				res.c14 <- OLE.test(dd=med.dates.rand,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort='rand',type='c14med',est=res.c14$Estimate,lo=res.c14$upperCI,hi=res.c14$lowerCI,hit=ifelse(res.c14$upperCI>true.start & res.c14$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
			}
			if (n[j]==k[z])
			{
				# if k=n
				res.cal <- OLE.test(dd=calendar.dates,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort=NA,type='calendar',est=res.cal$Estimate,lo=res.cal$upperCI,hi=res.cal$lowerCI,hit=ifelse(res.cal$upperCI>true.start & res.cal$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
				res.c14 <- OLE.test(dd=med.dates,alpha=0.05)
				tmp.res <- data.frame(sim=i,start.date=true.start,n=n[j],k=k[z],sort=NA,type='c14med',est=res.c14$Estimate,lo=res.c14$upperCI,hi=res.c14$lowerCI,hit=ifelse(res.c14$upperCI>true.start & res.c14$lowerCI < true.start,TRUE,FALSE))
				res <- rbind.data.frame(res,tmp.res)
			}
		}
	}

}
		

res <- res[-1,]
res$sort[is.na(res$sort)] <- 'na'

results.accuracy <- aggregate(hit~n+k+sort+type,data=res,FUN=mean)
results.precision <- aggregate(abs(lo-hi) ~ n + k + sort + type, data=res,FUN=mean)

















