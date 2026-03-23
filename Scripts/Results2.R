library(dplyr)
library(tidyr)
library(ggplot2)

res <- readRDS('../Results/full_output.rds')
res$precision  <- abs(res$upperCI-res$lowerCI)
res$accuracy <- ifelse(res$lowerCI < res$start_date & res$upperCI > res$start_date, 1, 0)
n.iter <- subset(res,Period=='Holocene' & r==0.01 & Sd==20 & ESS==10 &method=='CRIWM') |> nrow() 


res.holocene <- subset(res,Period=='Holocene')
res.pleistocene <- subset(res,Period=='Pleistocene')

#params
r  <- c(0.01,0.03,0.06)
Sd  <- c(20,50)
ESS  <- c(10,40,80)
method  <- c('OLE_medians','OLE_resamp_caldate','OLE_resamp_norm','OLE_resamp_unif','CRIWM','oxcal_trap','oxcal_unif') 


trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(r)){
	for (j in 1:length(Sd)){
	plot(NULL,xlim=c(1,26),ylim=c(0,1),axes=F,xlab='',ylab='Accuracy')
	title(main=paste0('r=',r[i],', Sd=',Sd[j]))
	counter  <- 1
		     for (k in 1:length(ESS)){
			     text(x=counter+2,y=1,paste('ESS=',ESS[k]),cex=1,pos=2)
				  for (m in 1:length(method))
				  {
					  ii  <- which(res.holocene$r==r[i]&res.holocene$Sd==Sd[j]&res.holocene$ESS==ESS[k]&res.holocene$method==method[m])
					  if(length(ii)>1)
					  {
					  tmp <- res.holocene[ii,]
					  rect(xleft=counter-trim,xright=counter+trim,ybottom=0,ytop=sum(tmp$accuracy)/nrow(tmp),border=NA,col='lightgrey')
axis(1,at=counter,label=method[m],las=2,cex=0.5)
counter  <- counter+1
				  }
				  }
	abline(v=counter-0.5,lty=2)
}
	axis(2,at=seq(0,1,0.25),las=2)
	abline(h=0.95,lty=3)
}}




trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(r)){
	for (j in 1:length(Sd)){
	plot(NULL,xlim=c(1,26),ylim=c(0,1200),axes=F,xlab='',ylab='Precision')
	title(main=paste0('r=',r[i],', Sd=',Sd[j]))
	counter  <- 1
		     for (k in 1:length(ESS)){
			     text(x=counter+2,y=1150,paste('ESS=',ESS[k]),cex=1,pos=2)
				  for (m in 1:length(method))
				  {
					  ii  <- which(res.holocene$r==r[i]&res.holocene$Sd==Sd[j]&res.holocene$ESS==ESS[k]&res.holocene$method==method[m])
					  if(length(ii)>1)
					  {
					  tmp <- res.holocene[ii,]
					  boxplot(tmp$precision,at=counter,add=T,outline=F,axes=F,lty=1)
axis(1,at=counter,label=method[m],las=2,cex=0.5)
counter  <- counter+1
				  }
				  }
	abline(v=counter-0.5,lty=2)
}
	axis(2,las=2)
	abline(h=0.95,lty=3)
}}
