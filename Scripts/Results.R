

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

#####################
####### PLEISTOCENE
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




