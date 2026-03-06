

####### Script to plot and visualise the results of the paper <whateva>

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")


### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE
hol_data <- full_data[full_data$Period == "Holocene",]

## indexes for looping
index_r <- unique(hol_data$r)
index_Sd <- unique(hol_data$Sd)
index_ESS <- sort(unique(hol_data$ESS))
index_method <- unique(hol_data$method)
pl_num <- length(unique(index_method))*3 + 4
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
		polygon(x = c(0.5,4.5,4.5,0.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(4.5,11.5,11.5,4.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("blue", alpha = 0.05), border = NA)
		polygon(x = c(11.5,18.5,18.5,11.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(18.5,25.5,25.5,18.5), y = c(ylim_hol_lo,ylim_hol_lo,ylim_hol_hi,ylim_hol_hi), col = adjustcolor("blue", alpha = 0.05), border = NA)
		text(x = c(2.5,8.5,14.5,21.5), y = rep(ylim_hol_lo+100,4), labels = c("5 dates", "10 dates", "40 dates", "80 dates"), cex  = 0.7)	
		index <- 1
		
		for (j in 1:length(index_ESS)){
			for (k in 1:length(index_method)){
				sgl <- hol_data_ss[hol_data_ss$ESS == index_ESS[j] & hol_data_ss$method == index_method[k],]
				
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
png("../Figures/Figure_Holocene_hits.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
ins <- c(rep("lightsteelblue1",4),rep("lightsteelblue2",7),rep("lightsteelblue3",7),rep("lightsteelblue4",7))
outs <- rep("gray97",length(ins))


par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		hol_data_ss <- hol_data_ss[order(hol_data_ss$ESS),]
		hol_data_ss$group <- paste(hol_data_ss$ESS,"dates",hol_data_ss$method)
		hol_data_ss$group <- factor(hol_data_ss$group, levels = unique(hol_data_ss$group))

		inside <- hol_data_ss$start_date > hol_data_ss$lowerCI & hol_data_ss$start_date < hol_data_ss$upperCI
		res <- table(hol_data_ss$group, inside)
		mat_res <- t(res)
		mat_res <- rbind(mat_res[2,],mat_res[1,])


		#barplot(prop.table(mat_res,2), col = c("red","slateblue4","yellow3","black"), las = 2, cex.names = 0.65, main = paste0("Holocene, r = ",index_r[h]," Sd = ",index_Sd[i]))
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

####################
######## PLEISTOCENE
ple_data <- full_data[full_data$Period == "Pleistocene",]

## indexes for looping
index_r <- unique(ple_data$r)
index_Sd <- unique(ple_data$Sd)
index_ESS <- sort(unique(ple_data$ESS))
index_method <- unique(ple_data$method)
pl_num <- length(unique(index_method))*3 + 4
axlabs <- c()
ylim_ple_lo <- -500
ylim_ple_hi <- 1000

png("../Figures/Figure_Pleistocene.png", res = 160, height = 1500, width = 1500)

par(mfrow = c(3,2))

for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		
		plot(x = c(-1:(pl_num+2)), y = rep(0,(pl_num+4)), type = "l", lty = 2, xlim = c(0,pl_num+1), ylim = c(ylim_ple_lo,ylim_ple_hi), ylab = "centred age", xlab = "", xaxt = "n",  main = paste0("Pleistocene, r = ",index_r[h]," Sd = ",index_Sd[i]))
		polygon(x = c(0.5,4.5,4.5,0.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(4.5,11.5,11.5,4.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("blue", alpha = 0.05), border = NA)
		polygon(x = c(11.5,18.5,18.5,11.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("red", alpha = 0.05), border = NA)
		polygon(x = c(18.5,25.5,25.5,18.5), y = c(ylim_ple_lo,ylim_ple_lo,ylim_ple_hi,ylim_ple_hi), col = adjustcolor("blue", alpha = 0.05), border = NA)
		text(x = c(2.5,8.5,14.5,21.5), y = rep(ylim_ple_lo+100,4), labels = c("5 dates", "10 dates", "40 dates", "80 dates"), cex  = 0.7)
		index <- 1
		
		for (j in 1:length(index_ESS)){
			for (k in 1:length(index_method)){
				sgl <- ple_data_ss[ple_data_ss$ESS == index_ESS[j] & ple_data_ss$method == index_method[k],]
				
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
png("../Figures/Figure_Pleistocene_hits.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
ins <- c(rep("lightsteelblue1",4),rep("lightsteelblue2",7),rep("lightsteelblue3",7),rep("lightsteelblue4",7))
outs <- rep("gray97",length(ins))


par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		## Order to keep consistency with previous plot
		ple_data_ss <- ple_data_ss[order(ple_data_ss$ESS),]
		ple_data_ss$group <- paste(ple_data_ss$ESS,"dates",ple_data_ss$method)
		ple_data_ss$group <- factor(ple_data_ss$group, levels = unique(ple_data_ss$group))

		inside <- ple_data_ss$start_date > ple_data_ss$lowerCI & ple_data_ss$start_date < ple_data_ss$upperCI
		res <- table(ple_data_ss$group, inside)
		mat_res <- t(res)
		mat_res <- rbind(mat_res[2,],mat_res[1,])
    
		#barplot(prop.table(mat_res,2), col = c("red","slateblue4","yellow3","black"), las = 2, cex.names = 0.65, main = paste0("Holocene, r = ",index_r[h]," Sd = ",index_Sd[i]))
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


