

###### Script to plot and visualise the results of the paper <whateva>

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")


### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE
## Load data
hol_data_full <- full_data[full_data$Period == "Holocene",]

## subset for OLE
hol_data_OLE <- hol_data_full[!(hol_data_full$method %in% unique(hol_data_full$method)[c(1,6,7)]),]
## subset for rest
hol_data_rest <- hol_data_full[hol_data_full$method %in% unique(hol_data_full$method)[c(1,6,7)],]
hol_data_OLE_10 <- hol_data_OLE[hol_data_OLE$ESS == 10,]
hol_data <- rbind(hol_data_OLE_10,hol_data_rest)

## indexes for looping
index_r <- unique(hol_data$r)
index_Sd <- unique(hol_data$Sd)
index_sample_size <- sort(unique(hol_data$sample_size))
index_method <- unique(hol_data$method)
pl_num <- length(unique(index_method))*3
axlabs <- c()
ylim_hol_lo <- -2000
ylim_hol_hi <- 750

## Plot CIs

png("../Figures/Figure_Holocene.png", res = 160, height = 1500, width = 1500)

par(mfrow = c(3,2))

for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i],]
		
		plot(x = c(-1:(pl_num+2)), y = rep(0,(pl_num+4)), type = "l", lty = 2, xlim = c(0,pl_num+1), ylim = c(ylim_hol_lo,ylim_hol_hi), ylab = "centred age", xlab = "", xaxt = "n",  main = paste0("Holocene, r = ",index_r[h]," Sd = ",index_Sd[i]))
		polygon(x = c(0.5,7.5,7.5,0.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(7.5,14.5,14.5,7.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("red", alpha = 0.09), border = NA)
		polygon(x = c(14.5,21.5,21.5,14.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("red", alpha = 0.13), border = NA)
		text(x = c(3.5,10.5,17.5), y = rep(ylim_hol_lo+100,4), labels = c("10 dates", "40 dates", "80 dates"), cex  = 0.7)	
		index <- 1
		
		for (j in 1:length(index_sample_size)){
			for (k in 1:length(index_method)){
				sgl <- hol_data_ss[hol_data_ss$sample_size == index_sample_size[j] & hol_data_ss$method == index_method[k],]
				
				if (nrow(sgl)!= 0){
					sgl$upperCI_norm <- sgl$upperCI - sgl$start_date
					sgl$lowerCI_norm <- sgl$lowerCI - sgl$start_date
					
					lines(x = c(index,index), y = c(mean(sgl$upperCI_norm),mean(sgl$lowerCI_norm)), col = "lightblue3")
					lines(x = c(index-0.2,index+0.2), y = rep(mean(sgl$upperCI_norm),2), lwd = 2, col = "lightblue4")
					lines(x = c(index-0.2,index+0.2), y = rep(mean(sgl$lowerCI_norm),2), lwd = 2, col = "lightblue4")
					axlabs <- append(axlabs,index_method[k])	
					index <- index + 1
				}
			}
		}
		axis(side = 1, at = seq(1,pl_num), labels = axlabs, las = 2, cex.axis = 0.5)
		axlabs <- c()
	}
}

dev.off()

# Plot target within CI
names_accuracy <- c("OLE-median, n = 10", "OLE-true, n = 10","OLE-normal, n = 10","OLE-uniform, n = 10","CRIWM, n = 10","BPM-trapezoid, n = 10", "BPM-uniform, n = 10","OLE-median, n = 40", "OLE-true, n = 40","OLE-normal, n = 40","OLE-uniform, n = 40","CRIWM, n = 40","BPM-trapezoid, n = 40", "BPM-uniform, n = 40","OLE-median, n = 80", "OLE-true, n = 80","OLE-normal, n = 80","OLE-uniform, n = 80","CRIWM, n = 80","BPM-trapezoid, n = 80", "BPM-uniform, n = 80")

png("../Figures/Fig_1.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
ins <- c(rep("lightsteelblue1",7),rep("lightsteelblue3",7),rep("lightsteelblue4",7))
outs <- rep("gray97",length(ins))


par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		hol_data_ss <- hol_data_ss[order(hol_data_ss$sample_size),]
		hol_data_ss$group <- paste(hol_data_ss$sample_size,"dates",hol_data_ss$method)
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

library(vioplot)
# Plot Precision
xlabs <- c("OLE-median","OLE-true","OLE-normal","OLE-uniform","CRIWM","BPM-trapezoid","BPM-uniform")

png("../Figures/Fig_3.png", res = 100, height = 1500, width = 1500)

hol_data$precision <- abs(hol_data$upperCI-hol_data$lowerCI)
colorinchis <- c("darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,9750),axes=F,xlab='',ylab='Precision')
		abline(h = seq(0,9000,100), lty = 2, col = adjustcolor("gray", alpha = 0.3))
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_sample_size)){
			text(x=counter+4,y=9700,paste('sample size=',index_sample_size[k]),cex=1,pos=2)
			for (m in 1:length(index_method)){
				ii  <- which(hol_data$r==index_r[i]&hol_data$Sd==index_Sd[j]&hol_data$sample_size==index_sample_size[k]&hol_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- hol_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[k])
					#boxplot(tmp$precision,at=counter,add=T,outline=F,axes=F,lty=1,col = colorinchis[k])
					axis(1,at=counter,label=xlabs[m],las=2,cex=0.5)
					counter  <- counter+1
				}
			}
			abline(v=counter-0.5,lty=2)
		}
		axis(2,las=2)
		abline(h=0.95,lty=3)
	}
}

dev.off()

##### SUPPLEMENTARY
index_ESS <- unique(hol_data_OLE$ESS)
index_method_OLE <- index_method[1:4]
index_SI <- c(1:4)
axlabs <- c()

## Shorter names for method
method_OLE_names <- c("Medians", "rCaldate", "rNormal", "rUniform") 

for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_OLE_ss <- hol_data_OLE[hol_data_OLE$r == index_r[h] & hol_data_OLE$Sd == index_Sd[i],]
		png(paste0("../Figures/Figure_Holocene_OLE_SI_Sd_",index_Sd[i],"_r_",index_r[h],".png"), res = 160, height = 1500, width = 1500)
		par(mfrow=c(2,2))
		for (j in 1:length(index_method_OLE)){
			plot(x = c(-1:(pl_num+2)), y = rep(0,(pl_num+4)), type = "l", lty = 2, xlim = c(0,5), ylim = c(-1700,1700), ylab = "centred age", xlab = "", xaxt = "n",  main = paste0("r = ",index_r[h]," Sd = ",index_Sd[i]," method = ", method_OLE_names[j]))
			abline(h = seq(-1700,1700,100), lty = 2, col = adjustcolor("blue", alpha = 0.1))
			for (k in 1:length(index_ESS)){
			sgl <- hol_data_OLE_ss[hol_data_OLE_ss$sample_size == 80 & hol_data_OLE_ss$ESS == index_ESS[k] & hol_data_OLE_ss$method == index_method_OLE[j],]
			if (nrow(sgl)!= 0){
					sgl$upperCI_norm <- sgl$upperCI - sgl$start_date
					sgl$lowerCI_norm <- sgl$lowerCI - sgl$start_date
					lines(x = c(index_SI[k],index_SI[k]), y = c(mean(sgl$upperCI_norm),mean(sgl$lowerCI_norm)), col = "lightblue3")
					lines(x = c(index_SI[k]-0.2,index_SI[k]+0.2), y = rep(mean(sgl$upperCI_norm),2), lwd = 2, col = "lightblue4")
					lines(x = c(index_SI[k]-0.2,index_SI[k]+0.2), y = rep(mean(sgl$lowerCI_norm),2), lwd = 2, col = "lightblue4")
					axlabs <- append(axlabs,paste0("k = ",sgl$ESS[k]))	
				}	
			}
			axis(side = 1, at = seq(1,4), labels = axlabs, las = 2, cex.axis = 0.8)
			axlabs <- c()
			#dev.off()
		}
		dev.off()
	}
}

####################
######## PLEISTOCENE
ple_data_full <- full_data[full_data$Period == "Pleistocene",]

## subset for OLE
ple_data_OLE <- ple_data_full[!(ple_data_full$method %in% unique(ple_data_full$method)[c(1,6,7)]),]
## subset for rest
ple_data_rest <- ple_data_full[ple_data_full$method %in% unique(ple_data_full$method)[c(1,6,7)],]
ple_data_OLE_10 <- ple_data_OLE[ple_data_OLE$ESS == 10,]
ple_data <- rbind(ple_data_OLE_10,ple_data_rest)

## indexes for looping
index_r <- unique(ple_data$r)
index_Sd <- unique(ple_data$Sd)
index_sample_size <- sort(unique(ple_data$sample_size))
index_method <- unique(ple_data$method)
pl_num <- length(unique(index_method))*3
axlabs <- c()
ylim_ple_lo <- -3000
ylim_ple_hi <- 2000

png("../Figures/Figure_Pleistocene.png", res = 160, height = 1500, width = 1500)

par(mfrow = c(3,2))

for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		
		plot(x = c(-1:(pl_num+2)), y = rep(0,(pl_num+4)), type = "l", lty = 2, xlim = c(0,pl_num+1), ylim = c(ylim_ple_lo,ylim_ple_hi), ylab = "centred age", xlab = "", xaxt = "n",  main = paste0("Pleistocene, r = ",index_r[h]," Sd = ",index_Sd[i]))
		polygon(x = c(0.5,7.5,7.5,0.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(7.5,14.5,14.5,7.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("red", alpha = 0.09), border = NA)
		polygon(x = c(14.5,21.5,21.5,14.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("red", alpha = 0.13), border = NA)
		text(x = c(3.5,10.5,17.5), y = rep(ylim_ple_lo+100,4), labels = c("10 dates", "40 dates", "80 dates"), cex  = 0.7)	
		index <- 1
		
		for (j in 1:length(index_sample_size)){
			for (k in 1:length(index_method)){
				sgl <- ple_data_ss[ple_data_ss$sample_size == index_sample_size[j] & ple_data_ss$method == index_method[k],]
				
				if (nrow(sgl)!= 0){
					sgl$upperCI_norm <- sgl$upperCI - sgl$start_date
					sgl$lowerCI_norm <- sgl$lowerCI - sgl$start_date
					
					lines(x = c(index,index), y = c(mean(sgl$upperCI_norm),mean(sgl$lowerCI_norm)), col = "lightblue3")
					lines(x = c(index-0.2,index+0.2), y = rep(mean(sgl$upperCI_norm),2), lwd = 2, col = "lightblue4")
					lines(x = c(index-0.2,index+0.2), y = rep(mean(sgl$lowerCI_norm),2), lwd = 2, col = "lightblue4")
					axlabs <- append(axlabs,index_method[k])	
					index <- index + 1
				}
			}
		}
		axis(side = 1, at = seq(1,pl_num), labels = axlabs, las = 2, cex.axis = 0.5)
		axlabs <- c()
	}
}

dev.off()

# Plot target within CI
png("../Figures/Fig_2.png", res = 160, height = 1500, width = 1500)

par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		ple_data_ss <- ple_data_ss[order(ple_data_ss$sample_size),]
		ple_data_ss$group <- paste(ple_data_ss$sample_size,"dates",ple_data_ss$method)
		ple_data_ss$group <- factor(ple_data_ss$group, levels = unique(ple_data_ss$group))

		inside <- ple_data_ss$start_date > ple_data_ss$lowerCI & ple_data_ss$start_date < ple_data_ss$upperCI
		res <- table(ple_data_ss$group, inside)
		mat_res <- t(res)
		mat_res <- rbind(mat_res[2,],mat_res[1,])
		colnames(mat_res) <- names_accuracy

		barplot(prop.table(mat_res,2), col = c("white"), las = 2, cex.names = 0.65, main = paste0("Pleistocene, r = ",index_r[h]," Sd = ",index_Sd[i]), border = NA)
			
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
png("../Figures/Fig_4.png", res = 100, height = 1500, width = 1500)

ple_data$precision <- abs(ple_data$upperCI-ple_data$lowerCI)
colorinchis <- c("darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,11500),axes=F,xlab='',ylab='Precision')
		abline(h = seq(0,10000,100), lty = 2, col = adjustcolor("gray", alpha = 0.3))
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_sample_size)){
			text(x=counter+4,y=11450,paste('sample size=',index_sample_size[k]),cex=1,pos=2)
			for (m in 1:length(index_method)){
				ii  <- which(ple_data$r==index_r[i]&ple_data$Sd==index_Sd[j]&ple_data$sample_size==index_sample_size[k]&ple_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- ple_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[k])
					#boxplot(tmp$precision,at=counter,add=T,outline=F,axes=F,lty=1,col = colorinchis[k])
					axis(1,at=counter,label=xlabs[m],las=2,cex=0.5)
					counter  <- counter+1
				}
			}
			abline(v=counter-0.5,lty=2)
		}
		axis(2,las=2)
		abline(h=0.95,lty=3)
	}
}

dev.off()

##### SUPPLEMENTARY
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_OLE_ss <- ple_data_OLE[ple_data_OLE$r == index_r[h] & ple_data_OLE$Sd == index_Sd[i],]
		png(paste0("../Figures/Figure_Pleistocene_OLE_SI_Sd_",index_Sd[i],"_r_",index_r[h],".png"), res = 160, height = 1500, width = 1500)
		par(mfrow=c(2,2))
		for (j in 1:length(index_method_OLE)){
			plot(x = c(-1:(pl_num+2)), y = rep(0,(pl_num+4)), type = "l", lty = 2, xlim = c(0,5), ylim = c(-3000,3000), ylab = "centred age", xlab = "", xaxt = "n",  main = paste0("r = ",index_r[h]," Sd = ",index_Sd[i]," method = ", method_OLE_names[j]))
			abline(h = seq(-3000,3000,100), lty = 2, col = adjustcolor("blue", alpha = 0.1))
			for (k in 1:length(index_ESS)){
			sgl <- ple_data_OLE_ss[ple_data_OLE_ss$sample_size == 80 & ple_data_OLE_ss$ESS == index_ESS[k] & ple_data_OLE_ss$method == index_method_OLE[j],]
			if (nrow(sgl)!= 0){
					sgl$upperCI_norm <- sgl$upperCI - sgl$start_date
					sgl$lowerCI_norm <- sgl$lowerCI - sgl$start_date
					lines(x = c(index_SI[k],index_SI[k]), y = c(mean(sgl$upperCI_norm),mean(sgl$lowerCI_norm)), col = "lightblue3")
					lines(x = c(index_SI[k]-0.2,index_SI[k]+0.2), y = rep(mean(sgl$upperCI_norm),2), lwd = 2, col = "lightblue4")
					lines(x = c(index_SI[k]-0.2,index_SI[k]+0.2), y = rep(mean(sgl$lowerCI_norm),2), lwd = 2, col = "lightblue4")
					axlabs <- append(axlabs,paste0("k = ",sgl$ESS[k]))	
				}	
			}
			axis(side = 1, at = seq(1,4), labels = axlabs, las = 2, cex.axis = 0.8)
			axlabs <- c()
			#dev.off()
		}
		dev.off()
	}
}

####### Plots for calibration curve
## Utilities for for looping
## These are from before
#index_r <- unique(hol_data$r)
#index_Sd <- unique(hol_data$Sd)
#index_sample_size <- sort(unique(hol_data$sample_size))
#index_method <- unique(hol_data$method)
#pl_num <- length(unique(index_method))*3

## New
#index_period <- unique(full_data$Period)
cols <- c(rep("lightblue2",4),"lightblue3",rep("lightblue4",2)) 

for (h in 1:length(index_r)){
	for (k in 1:length(index_Sd)){
		png(paste0("../Figures/Figure_Cal_curves","_r_",index_r[h],"_Sd_",index_Sd[k],".png"), res = 100, height = 1500, width = 1500)
		par(mfcol = c(7,2))
		
		for (j in 1:length(index_method)){
			subset_plot <- full_data[full_data$Period == "Holocene"  & full_data$method == index_method[j] & full_data$r == index_r[h] & full_data$Sd == index_Sd[k] & full_data$sample_size == 80,]
			
			ylimvals <- c(-700,700)
			
			if (length(unique(subset_plot$ESS)) > 1){
				subset_plot <- subset_plot[subset_plot$ESS == 10,]
				ylimvals <- c(-5000,5000)
			}
			
			plot(x=c(11000:1000), y = rep(0,10001), col = "white", xlim = c(11000,1000), ylim = ylimvals, xlab = "Time", ylab = "Precision", main = paste0("Holocene ",index_method[j], " r = ",index_r[h], " Sd = ", index_Sd[k]))
			for (i in 1:nrow(subset_plot)){
				low <- subset_plot$upperCI[i] - subset_plot$start_date[i]
				high <- subset_plot$lowerCI[i] - subset_plot$start_date[i]
				lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = cols[j], lwd = 2) 
			}
		}
		
		for (j in 1:length(index_method)){
			subset_plot <- full_data[full_data$Period == "Pleistocene" & full_data$method == index_method[j] & full_data$r == index_r[h] & full_data$Sd == index_Sd[k] & full_data$sample_size == 80,]
			
			ylimvals <- c(-1500,1500)
			
			if (length(unique(subset_plot$ESS)) > 1){
				subset_plot <- subset_plot[subset_plot$ESS == 10,]
				ylimvals <- c(-6500,6500)
			}
			plot(x=c(40000:30000), y = rep(0,10001), col = "white", xlim = c(40000,30000), ylim = ylimvals, xlab = "Time", ylab = "Precision", main = paste0("Pleistocene ",index_method[j], " r = ", index_r[h], " Sd = ", index_Sd[k]))
			for (i in 1:nrow(subset_plot)){
				low <- subset_plot$upperCI[i] - subset_plot$start_date[i]
				high <- subset_plot$lowerCI[i] - subset_plot$start_date[i]
				lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = cols[j], lwd = 2) 
			}
		}
		dev.off()
	}
}

##### CHECK COMPUTATION TIMES
time_criwm_hol <- readRDS("../Results/time_criwm_hol.rds")
time_criwm_upl <- readRDS("../Results/time_criwm_upl.rds")
time_OLE_medians_hol <- readRDS("../Results/time_OLE_medians_hol.rds")
time_OLE_medians_upl <- readRDS("../Results/time_OLE_medians_upl.rds")
time_OLE_resamp_caldate_hol <- readRDS("../Results/time_OLE_resamp_caldate_hol.rds")
time_OLE_resamp_caldate_upl <- readRDS("../Results/time_OLE_resamp_caldate_upl.rds")
time_OLE_resamp_norm_hol <- readRDS("../Results/time_OLE_resamp_norm_hol.rds")
time_OLE_resamp_norm_upl <- readRDS("../Results/time_OLE_resamp_norm_upl.rds")
time_OLE_resamp_unif_hol <- readRDS("../Results/time_OLE_resamp_unif_hol.rds")
time_OLE_resamp_unif_upl <- readRDS("../Results/time_OLE_resamp_unif_upl.rds")
time_oxcal_unif_hol <- readRDS("../Results/time_oxcal_unif_hol.rds")
time_oxcal_unif_upl <- readRDS("../Results/time_oxcal_unif_upl.rds")
time_oxcal_trap_hol <- readRDS("../Results/time_oxcal_trap_hol.rds")
time_oxcal_trap_upl <- readRDS("../Results/time_oxcal_trap_upl.rds")

time_seconds <- data.frame("Holocene" = c(time_criwm_hol[3],time_OLE_medians_hol[3],time_OLE_resamp_unif_hol[3],time_OLE_resamp_norm_hol[3],time_OLE_resamp_caldate_hol[3],time_oxcal_unif_hol[3],time_oxcal_trap_hol[3]),
			  "Pleistocene"  = c(time_criwm_upl[3],time_OLE_medians_upl[3],time_OLE_resamp_unif_upl[3],time_OLE_resamp_norm_upl[3],time_OLE_resamp_caldate_upl[3],time_oxcal_unif_upl[3],time_oxcal_trap_upl[3]))


df_format <- as.data.frame(lapply(time_seconds, function(x) {
					  sprintf("%02.0f:%02.0f:%02.0f",x %/% 3600,(x %% 3600) %/% 60,x %% 60)}
			  )
)


