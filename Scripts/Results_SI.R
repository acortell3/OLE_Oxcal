

###### Script to plot and visualise the results of the paper <whateva>

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")


### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE
## Load data
hol_data_full <- full_data[full_data$Period == "Holocene",]

## subset for OLE
hol_data <- hol_data_full[!(hol_data_full$method %in% unique(hol_data_full$method)[c(1,6,7)]),]

## indexes for looping
index_r <- unique(hol_data$r)
index_Sd <- unique(hol_data$Sd)
index_sample_size <- sort(unique(hol_data$sample_size))
index_method <- unique(hol_data$method)
pl_num <- length(unique(index_method))*4
axlabs <- c()
ylim_hol_lo <- -2000
ylim_hol_hi <- 750
index_ESS <- unique(hol_data$ESS)
index_method_OLE <- index_method[1:4]
index_SI <- c(1:4)
axlabs <- c()

## Shorter names for method
names_methods <- c("Median","True","Normal","Uniform")
names_accuracy <- c(paste0(names_methods,", k = 5"),paste0(names_methods,", k = 10"),paste0(names_methods,", k = 40"),paste0(names_methods,", k = 80"))


#### ACCURACY

png("../Figures/SI_1.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
ins <- c(rep("lightsteelblue1",4),rep("lightsteelblue2",4),rep("lightsteelblue3",4),rep("lightsteelblue4",4))
outs <- rep("gray97",length(ins))


par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		hol_data_ss <- hol_data_ss[order(hol_data_ss$ESS),]
		hol_data_ss$group <- paste("k =",hol_data_ss$ESS,"dates",hol_data_ss$method)
		hol_data_ss$group <- factor(hol_data_ss$group, levels = unique(hol_data_ss$group))

		inside <- hol_data_ss$start_date > hol_data_ss$lowerCI & hol_data_ss$start_date < hol_data_ss$upperCI
		res <- table(hol_data_ss$group, inside)
		mat_res <- t(res)
		mat_res <- rbind(mat_res[2,],mat_res[1,])
		colnames(mat_res) <- names_accuracy

		barplot(prop.table(mat_res,2), col = c("white"), las = 2, cex.names = 0.65, main = paste0("Holocene, r = ",index_r[h]," Sd = ",index_Sd[i]), border = NA)
			
		for (j in 1:ncol(mat_res)){
			xx = prop.table(mat_res,2)
			xx[,-j] = NA
			colnames(xx)[-j] = NA
			barplot(xx,col=c(ins[j],outs[j]), add=T, axes=F,xaxt = "n", border = NA) 
		}
		abline(h=0.95, lty = 2, lwd = 0.5)
	}
}
dev.off()

# Plot Precision
library(vioplot)
xlabs <- names_accuracy

png("../Figures/SI_3.png", res = 100, height = 1500, width = 1500)

hol_data$precision <- abs(hol_data$upperCI-hol_data$lowerCI)
colorinchis <- c("darkorange1","darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,9750),axes=F,xlab='',ylab='Precision')
		abline(h = seq(0,9000,100), lty = 2, col = adjustcolor("gray", alpha = 0.3))
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_ESS)){
			for (m in 1:length(index_method)){
				ii  <- which(hol_data$r==index_r[i]&hol_data$Sd==index_Sd[j]&hol_data$ESS==index_ESS[k]&hol_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- hol_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[k])
					counter  <- counter+1
				}
			}
			axis(1,at=c(1:16),label=xlabs,las=2,cex=0.5)
			abline(v=counter-0.5,lty=2)
		}
		axis(2,las=2)
		abline(h=0.95,lty=3)
	}
}

dev.off()

## Plots precision and accuracy
hol_data$accuracy <- ifelse(hol_data$start_date<hol_data$upperCI & hol_data$start_date >hol_data$lowerCI,TRUE,FALSE)

for (j in 1:length(index_ESS)){
	for (k in 1:length(index_method)){
		
		png(paste0("../Figures/SI_5_Holocene_k_",index_ESS[j],"_method_",index_method[k],".png"), res = 100, height = 1500, width = 1500)
		par(mfrow = c(3,2), mar=c(7,5,3,5))
		
		for (h in 1:length(index_r)){
			for (i in 1:length(index_Sd)){
				hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i] & hol_data$ESS == index_ESS[j] & hol_data$method == index_method[k],]
				lab <- paste0(paste0(0:11,'k-'),paste0(1:12,'k'))
				hol_data_ss$bracket  <- cut(hol_data_ss$precision,breaks=seq(0,12000,1000),labels = lab)
				pl_dt <- hol_data_ss |>
					filter(!is.na(bracket)) |>
					group_by(bracket) |>
					summarise(n = n(),accuracy = mean(accuracy == TRUE),.groups = "drop")
				if (length(lab) > nrow(pl_dt)){
					lab <- lab[1:nrow(pl_dt)]
				}
				barplot(height=pl_dt$n,space = 0,names.arg = lab,las=2,xlab='Precision',ylab='Number of Cases', col = "gold3", main = paste0("r = ",index_r[h]," Sd = ", index_Sd[i]))
				lines(x=(0:(length(pl_dt$accuracy)-1))+0.5,y=pl_dt$accuracy*(max(pl_dt$n)-20),ylim = c(0,max(pl_dt$n)),pch=21,type='b',lwd=2,cex=2, col = "darkolivegreen", bg = "darkolivegreen4")
				axis(side=4,at=c(0,max(pl_dt$n)*0.25,max(pl_dt$n)*0.5,max(pl_dt$n)*0.75,max(pl_dt$n)),labels=c(0,0.25,0.5,0.75,1.0),las=2)
				abline(h=max(pl_dt$n)*0.95,lty=2)
			}
		}
		dev.off()
	}
}

####################
######## PLEISTOCENE
ple_data_full <- full_data[full_data$Period == "Pleistocene",]

## subset for OLE
ple_data <- ple_data_full[!(ple_data_full$method %in% unique(ple_data_full$method)[c(1,6,7)]),]

## indexes for looping
index_r <- unique(ple_data$r)
index_Sd <- unique(ple_data$Sd)
index_sample_size <- sort(unique(ple_data$sample_size))
index_method <- unique(ple_data$method)
pl_num <- length(unique(index_method))*4
axlabs <- c()
ylim_ple_lo <- -3000
ylim_ple_hi <- 2000

#### ACCURACY

png("../Figures/SI_2.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
ins <- c(rep("lightsteelblue1",4),rep("lightsteelblue2",4),rep("lightsteelblue3",4),rep("lightsteelblue4",4))
outs <- rep("gray97",length(ins))


par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		ple_data_ss <- ple_data_ss[order(ple_data_ss$ESS),]
		ple_data_ss$group <- paste("k =",ple_data_ss$ESS,"dates",ple_data_ss$method)
		ple_data_ss$group <- factor(ple_data_ss$group, levels = unique(ple_data_ss$group))

		inside <- ple_data_ss$start_date > ple_data_ss$lowerCI & ple_data_ss$start_date < ple_data_ss$upperCI
		res <- table(ple_data_ss$group, inside)
		mat_res <- t(res)
		mat_res <- rbind(mat_res[2,],mat_res[1,])
		colnames(mat_res) <- names_accuracy

		barplot(prop.table(mat_res,2), col = c("white"), las = 2, cex.names = 0.65, main = paste0("Holocene, r = ",index_r[h]," Sd = ",index_Sd[i]), border = NA)
			
		for (j in 1:ncol(mat_res)){
			xx = prop.table(mat_res,2)
			xx[,-j] = NA
			colnames(xx)[-j] = NA
			barplot(xx,col=c(ins[j],outs[j]), add=T, axes=F,xaxt = "n", border = NA) 
		}
		abline(h=0.95, lty = 2, lwd = 0.5)
	}
}
dev.off()

png("../Figures/SI_4.png", res = 100, height = 1500, width = 1500)

ple_data$precision <- abs(ple_data$upperCI-ple_data$lowerCI)
colorinchis <- c("darkorange1","darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,9750),axes=F,xlab='',ylab='Precision')
		abline(h = seq(0,9000,100), lty = 2, col = adjustcolor("gray", alpha = 0.3))
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_ESS)){
			for (m in 1:length(index_method)){
				ii  <- which(ple_data$r==index_r[i]&ple_data$Sd==index_Sd[j]&ple_data$ESS==index_ESS[k]&ple_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- ple_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[k])
					counter  <- counter+1
				}
			}
			axis(1,at=c(1:16),label=xlabs,las=2,cex=0.5)
			abline(v=counter-0.5,lty=2)
		}
		axis(2,las=2)
		abline(h=0.95,lty=3)
	}
}

dev.off()

## Plots precision and accuracy
ple_data$accuracy <- ifelse(ple_data$start_date<ple_data$upperCI & ple_data$start_date >ple_data$lowerCI,TRUE,FALSE)

for (j in 1:length(index_ESS)){
	for (k in 1:length(index_method)){
		
		png(paste0("../Figures/SI_6_Pleistocene_k_",index_ESS[j],"_method_",index_method[k],".png"), res = 100, height = 1500, width = 1500)
		par(mfrow = c(3,2), mar=c(7,5,3,5))
		
		for (h in 1:length(index_r)){
			for (i in 1:length(index_Sd)){
				ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i] & ple_data$ESS == index_ESS[j] & ple_data$method == index_method[k],]
				lab <- paste0(paste0(0:11,'k-'),paste0(1:12,'k'))
				ple_data_ss$bracket  <- cut(ple_data_ss$precision,breaks=seq(0,12000,1000),labels = lab)
				pl_dt <- ple_data_ss |>
					filter(!is.na(bracket)) |>
					group_by(bracket) |>
					summarise(n = n(),accuracy = mean(accuracy == TRUE),.groups = "drop")
				if (length(lab) > nrow(pl_dt)){
					lab <- lab[1:nrow(pl_dt)]
				}
				barplot(height=pl_dt$n,space = 0,names.arg = lab,las=2,xlab='Precision',ylab='Number of Cases', col = "gold3", main = paste0("r = ",index_r[h]," Sd = ", index_Sd[i]))
				lines(x=(0:(length(pl_dt$accuracy)-1))+0.5,y=pl_dt$accuracy*(max(pl_dt$n)-20),ylim = c(0,max(pl_dt$n)),pch=21,type='b',lwd=2,cex=2, col = "darkolivegreen", bg = "darkolivegreen4")
				axis(side=4,at=c(0,max(pl_dt$n)*0.25,max(pl_dt$n)*0.5,max(pl_dt$n)*0.75,max(pl_dt$n)),labels=c(0,0.25,0.5,0.75,1.0),las=2)
				abline(h=max(pl_dt$n)*0.95,lty=2)
			}
		}
		dev.off()
	}
}
