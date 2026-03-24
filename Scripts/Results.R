

####### Script to plot and visualise the results of the paper <whateva>

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
ylim_hol_lo <- -400
ylim_hol_hi <- 650

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
png("../Figures/Figure_Holocene_accuracy.png", res = 160, height = 1500, width = 1500)

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
png("../Figures/Figure_Holocene_Precision.png", res = 100, height = 1500, width = 1500)

hol_data$precision <- abs(hol_data$upperCI-hol_data$lowerCI)
colorinchis <- c("darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,1700),axes=F,xlab='',ylab='Precision')
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_sample_size)){
			text(x=counter+4,y=1650,paste('sample size=',index_sample_size[k]),cex=1,pos=2)
			for (m in 1:length(index_method)){
				ii  <- which(hol_data$r==index_r[i]&hol_data$Sd==index_Sd[j]&hol_data$sample_size==index_sample_size[k]&hol_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- hol_data[ii,]
					boxplot(tmp$precision,at=counter,add=T,outline=F,axes=F,lty=1,col = colorinchis[k])
					axis(1,at=counter,label=index_method[m],las=2,cex=0.5)
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
string_find <- unique(hol_data_OLE$method)
string_title <- c("medians","resamp cal curve","resamp normal","resamp uniform")

for (k in 1:length(string_find)){
	OLE_hol <- hol_data_OLE[hol_data_OLE$method == string_find[k],]
	
	cs <- list()
	index_nd <- unique(OLE_hol$ESS)
	index_r_si <- unique(OLE_hol$r)
	index_sd_si <- unique(OLE_hol$Sd)
	
	png(paste0("../Figures/Figure_SI_OLE_Hol_",k,".png"), res = 100, height = 1000, width = 1000)
	par(mfrow = c(4,2))
	for (j in 1:length(index_r_si)){
		for (h in 1:length(index_sd_si)){
			for (i in 1:length(index_nd)){
				cs[[i]] <- OLE_hol[OLE_hol$ESS == index_nd[i] & OLE_hol$r == index_r_si[j] & OLE_hol$Sd == index_sd_si[h],]
				cs[[i]]$Estimate_norm <- cs[[i]]$Estimate - cs[[i]]$start_date
				cs[[i]]$upperCI_norm <- cs[[i]]$upperCI - cs[[i]]$start_date
				cs[[i]]$lowerCI_norm <- cs[[i]]$lowerCI - cs[[i]]$start_date
			}
			
			plot(x=1, y = mean(cs[[1]]$Estimate_norm), pch = 16, col = "blue", ylim = c(-200,650), xlim = c(0,5), ylab = "Estimate", xlab = "ndates", xaxt = "n", main = paste0("Holocene OLE ",string_title[k],", sd = ",index_sd_si[h]," r = ",index_r_si[j]))
			axis(side = 1, at = seq(0,5), labels = c("",5,10,20,40,""))
			lines(x=rep(1,2), y = c(mean(cs[[1]]$upperCI_norm),mean(cs[[1]]$lowerCI_norm)), lwd = 2, col = "orange")
			abline(h = 0, lty = 2)
			for (i in 1:length(cs)){
				lines(x=rep(i,2), y = c(mean(cs[[i]]$upperCI_norm),mean(cs[[i]]$lowerCI_norm)), lwd = 2, col = "orange")
				points(x=i, y = mean(cs[[i]]$Estimate_norm), pch = 16, col = "blue")
			}
		}
	}
	dev.off()
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
ylim_ple_lo <- -500
ylim_ple_hi <- 1000

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
png("../Figures/Figure_Pleistocene_accuracy.png", res = 160, height = 1500, width = 1500)

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
png("../Figures/Figure_Pleistocene_Precision.png", res = 100, height = 1500, width = 1500)

ple_data$precision <- abs(ple_data$upperCI-ple_data$lowerCI)
colorinchis <- c("darkorange2","darkorange3","darkorange4")

#trim  <- 0.2
par(mfrow=c(3,2),mar=c(10,4,2,2))
for (i in 1:length(index_r)){
	for (j in 1:length(index_Sd)){
		plot(NULL,xlim=c(1,pl_num),ylim=c(0,1700),axes=F,xlab='',ylab='Precision')
		title(main=paste0('r=',index_r[i],', Sd=',index_Sd[j]))
		counter  <- 1
		for (k in 1:length(index_sample_size)){
			text(x=counter+4,y=1650,paste('sample size=',index_sample_size[k]),cex=1,pos=2)
			for (m in 1:length(index_method)){
				ii  <- which(ple_data$r==index_r[i]&ple_data$Sd==index_Sd[j]&ple_data$sample_size==index_sample_size[k]&ple_data$method==index_method[m])
				if(length(ii)>1){
					tmp <- ple_data[ii,]
					boxplot(tmp$precision,at=counter,add=T,outline=F,axes=F,lty=1,col = colorinchis[k])
					axis(1,at=counter,label=index_method[m],las=2,cex=0.5)
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
string_find <- unique(ple_data_OLE$method)
string_title <- c("medians","resamp cal curve","resamp normal","resamp uniform")

for (k in 1:length(string_find)){
	OLE_ple <- ple_data_OLE[ple_data_OLE$method == string_find[k],]
	
	cs <- list()
	index_nd <- unique(OLE_ple$ESS)
	index_r_si <- unique(OLE_ple$r)
	index_sd_si <- unique(OLE_ple$Sd)
	
	png(paste0("../Figures/Figure_SI_OLE_Ple_",k,".png"), res = 100, height = 1000, width = 1000)
	par(mfrow = c(4,2))
	for (j in 1:length(index_r_si)){
		for (h in 1:length(index_sd_si)){
			for (i in 1:length(index_nd)){
				cs[[i]] <- OLE_ple[OLE_ple$ESS == index_nd[i] & OLE_ple$r == index_r_si[j] & OLE_ple$Sd == index_sd_si[h],]
				cs[[i]]$Estimate_norm <- cs[[i]]$Estimate - cs[[i]]$start_date
				cs[[i]]$upperCI_norm <- cs[[i]]$upperCI - cs[[i]]$start_date
				cs[[i]]$lowerCI_norm <- cs[[i]]$lowerCI - cs[[i]]$start_date
			}
			
			plot(x=1, y = mean(cs[[1]]$Estimate_norm), pch = 16, col = "blue", ylim = c(-100,1000), xlim = c(0,5), ylab = "Estimate", xlab = "ndates", xaxt = "n", main = paste0("Pleistocene OLE ",string_title[k],", sd = ",index_sd_si[h]," r = ",index_r_si[j]))
			axis(side = 1, at = seq(0,5), labels = c("",5,10,20,40,""))
			lines(x=rep(1,2), y = c(mean(cs[[1]]$upperCI_norm),mean(cs[[1]]$lowerCI_norm)), lwd = 2, col = "orange")
			abline(h = 0, lty = 2)
			for (i in 1:length(cs)){
				lines(x=rep(i,2), y = c(mean(cs[[i]]$upperCI_norm),mean(cs[[i]]$lowerCI_norm)), lwd = 2, col = "orange")
				points(x=i, y = mean(cs[[i]]$Estimate_norm), pch = 16, col = "blue")
			}
		}
	}
	dev.off()
}
