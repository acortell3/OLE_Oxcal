
######### CODE FOR "UNCERTAIN BEGINNINGS: A COMPARISON OF THE ACCURACY AND PRECISION OF METHODS ESTIMATING EXTREME CHRONOLOGICAL EVENTS", BY A.CORTELL-NICOLAU, E. R. CREMA, AND A. KEY


######### PLOT AND VISUALISE RESULTS

#### Authors: Alfredo Cortell-Nicolau & Enrico R. Crema

#### Load data, utilities and packages
full_data <- readRDS("../Results/full_output.rds")
full_data <- full_data[full_data$method != "OLE_resamp_norm" & full_data$method != "OLE_resamp_unif",]



### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE
## Load data
#hol_data_full <- full_data[full_data$Period == "Holocene",]
hol_data <- full_data[full_data$Period == "Holocene",]

## indexes for looping
index_r <- unique(hol_data$r)
index_Sd <- unique(hol_data$Sd)
index_sample_size <- sort(unique(hol_data$sample_size))
index_method <- unique(hol_data$method)
pl_num <- length(unique(index_method))*3
axlabs <- c()
ylim_hol_lo <- -2000
ylim_hol_hi <- 750

# Plot target within CI
names_accuracy <- c("CRIWM, n = 10","OLE-median, n = 10", "OLE-true, n = 10","BPM-trapezoid, n = 10", "BPM-uniform, n = 10","","CRIWM, n = 40","OLE-median, n = 40", "OLE-true, n = 40","BPM-trapezoid, n = 40", "BPM-uniform, n = 40","","CRIWM, n = 80","OLE-median, n = 80", "OLE-true, n = 80","BPM-trapezoid, n = 80", "BPM-uniform, n = 80")

png("../Figures/Fig_1.png", res = 160, height = 1500, width = 1500)

## Create colour palettes
col_group <- c(
 "#444444",  # dark grey
  "#E69F00",  # light orange
  "#D55E00",  # dark orange
  "#56B4E9",  # light blue
  "#0072B2"   # dark blue
)
ins <- c(col_group,"white",col_group,"white",col_group)
outs <- rep("gray97",length(ins))

par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		hol_data_ss <- hol_data[hol_data$r == index_r[h] & hol_data$Sd == index_Sd[i],]
		hol_data_ss <- hol_data_ss[hol_data_ss$sample_size == hol_data_ss$ESS,]
		## Order to keep consistency with previous plot
		hol_data_ss <- hol_data_ss[order(hol_data_ss$sample_size),]
		hol_data_ss$group <- paste(hol_data_ss$sample_size,"dates",hol_data_ss$method)
		hol_data_ss$group <- factor(hol_data_ss$group, levels = unique(hol_data_ss$group))

		inside <- hol_data_ss$start_date > hol_data_ss$lowerCI & hol_data_ss$start_date < hol_data_ss$upperCI
		idx <- which(hol_data_ss$nlegs > 1)
		inside[idx] <- sapply(idx, function(r) {
					      n <- hol_data_ss$nlegs[r]
					      start <- hol_data_ss$start_date[r]
					      upper <- as.numeric(hol_data_ss[r, 14 + 2 * (0:(n - 1))])
					      lower <- as.numeric(hol_data_ss[r, 15 + 2 * (0:(n - 1))])
					      any(start > lower & start < upper, na.rm = TRUE)})

		res <- table(hol_data_ss$group, inside)
		## Prepare blank rows
		res <- rbind(res[1:5,],c(1000,0),res[6:10,],c(1000,0),res[11:15,])
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
xlabs <- c("CRIWM","OLE-median","OLE-true","BPM-trapezoid","BPM-uniform")

png("../Figures/Fig_3.png", res = 100, height = 1500, width = 1500)

hol_data$precision <- abs(hol_data$upperCI-hol_data$lowerCI)
idx <- which(hol_data$nlegs > 1)
hol_data$precision[idx] <- sapply(idx, function(r){
					      n <- hol_data$nlegs[r]
					      uppers <- as.numeric(hol_data[r, 14 + 2 * (0:(n - 1))])
					      lowers <- as.numeric(hol_data[r, 15 + 2 * (0:(n - 1))])
					      return(sum(abs(uppers-lowers)))
					      })

colorinchis <- col_group

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
				ii  <- which(hol_data$r==index_r[i]&hol_data$Sd==index_Sd[j]&hol_data$sample_size==index_sample_size[k]&hol_data$method==index_method[m]&hol_data$ESS==index_sample_size[k])
				if(length(ii)>1){
					tmp <- hol_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[m])
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


####################
######## PLEISTOCENE
ple_data <- full_data[full_data$Period == "Pleistocene",]

## indexes for looping
index_r <- unique(ple_data$r)
index_Sd <- unique(ple_data$Sd)
index_sample_size <- sort(unique(ple_data$sample_size))
index_method <- unique(ple_data$method)
pl_num <- length(unique(index_method))*3
axlabs <- c()
ylim_ple_lo <- -3000
ylim_ple_hi <- 2000

# Plot target within CI
png("../Figures/Fig_2.png", res = 160, height = 1500, width = 1500)

par(mfrow = c(2,3), mar = c(8,4,4,2))
for (h in 1:length(index_r)){
	for (i in 1:length(index_Sd)){
		ple_data_ss <- ple_data[ple_data$r == index_r[h] & ple_data$Sd == index_Sd[i],]
		ple_data_ss <- ple_data_ss[ple_data_ss$sample_size == ple_data_ss$ESS,]
		## Order to keep consistency with previous plot
		ple_data_ss <- ple_data_ss[order(ple_data_ss$sample_size),]
		ple_data_ss$group <- paste(ple_data_ss$sample_size,"dates",ple_data_ss$method)
		ple_data_ss$group <- factor(ple_data_ss$group, levels = unique(ple_data_ss$group))

		inside <- ple_data_ss$start_date > ple_data_ss$lowerCI & ple_data_ss$start_date < ple_data_ss$upperCI	
		idx <- which(ple_data_ss$nlegs > 1)
		inside[idx] <- sapply(idx, function(r) {
					      n <- ple_data_ss$nlegs[r]
					      start <- ple_data_ss$start_date[r]
					      upper <- as.numeric(ple_data_ss[r, 14 + 2 * (0:(n - 1))])
					      lower <- as.numeric(ple_data_ss[r, 15 + 2 * (0:(n - 1))])
					      any(start > lower & start < upper, na.rm = TRUE)})
		
		res <- table(ple_data_ss$group, inside)
		## Prepare blank rows
		res <- rbind(res[1:5,],c(1000,0),res[6:10,],c(1000,0),res[11:15,])
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
idx <- which(hol_data$nlegs > 1)
ple_data$precision[idx] <- sapply(idx, function(r){
					      n <- ple_data$nlegs[r]
					      uppers <- as.numeric(hol_data[r, 14 + 2 * (0:(n - 1))])
					      lowers <- as.numeric(hol_data[r, 15 + 2 * (0:(n - 1))])
					      return(sum(abs(uppers-lowers)))
					      })
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
				ii  <- which(ple_data$r==index_r[i]&ple_data$Sd==index_Sd[j]&ple_data$sample_size==index_sample_size[k]&ple_data$method==index_method[m]&ple_data$ESS==index_sample_size[k])
				if(length(ii)>1){
					tmp <- ple_data[ii,]
					vioplot(tmp$precision, at = counter, add = T, outline = F, axes = F, lty =1 , col = colorinchis[m])
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


####### Plots for calibration curve
method_name <- c("CRIWM","OLE-medians","OLE-true","BPM-trapezoid","BPM-uniform")

for (h in 1:length(index_r)){
	for (k in 1:length(index_Sd)){
		png(paste0("../Figures/Figure_Cal_curves","_r_",index_r[h],"_Sd_",index_Sd[k],".png"), res = 100, height = 1500, width = 1500)
		par(mfcol = c(5,2))

		for (j in 1:length(index_method)){
			subset_plot <- full_data[full_data$Period == "Holocene"  & full_data$method == index_method[j] & full_data$r == index_r[h] & full_data$Sd == index_Sd[k] & full_data$sample_size == 80,]
			
			ylimvals <- c(-700,700)
			
			if (length(unique(subset_plot$ESS)) > 1){
				subset_plot <- subset_plot[subset_plot$ESS == 10,]
				ylimvals <- c(-5000,5000)
			}
			
			plot(x=c(11000:1000), y = rep(0,10001), col = "white", xlim = c(11000,1000), ylim = ylimvals, xlab = "Time", ylab = "Centered estimates", main = paste0("Holocene ",method_name[j], " r = ",index_r[h], " Sd = ", index_Sd[k]))
			for (i in 1:nrow(subset_plot)){
				if (!is.na(subset_plot$nleg[i]) && subset_plot$nleg[i]>1){
					n <- subset_plot$nleg[i]
					uppers <- as.numeric(subset_plot[i, 14 + 2 * (0:(n - 1))])
    					lowers <- as.numeric(subset_plot[i, 15 + 2 * (0:(n - 1))])
    
    					# Does 0 fall inside at least one of the legs?
    					inside <- any(0 > lowers & 0 < uppers, na.rm = TRUE)
   
				        leg_col <- ifelse(inside,"firebrick1","darkslategray4")	
    					#leg_col <- ifelse(inside, "darkslategray4", "firebrick1")
					for (x in 1:n){
						low <- subset_plot[i,13+x] - subset_plot$start_date[i]
						high <- subset_plot[i,14+x] - subset_plot$start_date[i]
						lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = leg_col, lwd = 2)
						abline(h=0,lty  = 2)}	
					} else {
						low <- subset_plot$upperCI[i] - subset_plot$start_date[i]
						high <- subset_plot$lowerCI[i] - subset_plot$start_date[i]
						lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = ifelse(0>high & 0 < low,"darkslategray4","firebrick1"), lwd = 2)
						abline(h=0,lty  = 2)
				}	
			}
		}
		
		for (j in 1:length(index_method)){
			subset_plot <- full_data[full_data$Period == "Pleistocene" & full_data$method == index_method[j] & full_data$r == index_r[h] & full_data$Sd == index_Sd[k] & full_data$sample_size == 80,]
			
			ylimvals <- c(-1500,1500)
			
			if (length(unique(subset_plot$ESS)) > 1){
				subset_plot <- subset_plot[subset_plot$ESS == 10,]
				ylimvals <- c(-6500,6500)
			}
			plot(x=c(40000:30000), y = rep(0,10001), col = "white", xlim = c(40000,30000), ylim = ylimvals, xlab = "Time", ylab = "Centered estimates", main = paste0("Pleistocene ",method_name[j], " r = ", index_r[h], " Sd = ", index_Sd[k]))
			
			for (i in 1:nrow(subset_plot)){
				if (!is.na(subset_plot$nleg[i]) && subset_plot$nleg[i]>1){
					n <- subset_plot$nleg[i]
					uppers <- as.numeric(subset_plot[i, 14 + 2 * (0:(n - 1))])
					lowers <- as.numeric(subset_plot[i, 15 + 2 * (0:(n - 1))])
					
					# Does 0 fall inside at least one of the legs?
    					inside <- any(0 > lowers & 0 < uppers, na.rm = TRUE)
   
				        leg_col <- ifelse(inside,"firebrick1","darkslategray4")	
    					#leg_col <- ifelse(inside, "darkslategray4", "firebrick1")
					for (x in 1:n){
						low <- subset_plot[i,13+x] - subset_plot$start_date[i]
						high <- subset_plot[i,14+x] - subset_plot$start_date[i]
						lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = leg_col, lwd = 2)
						abline(h=0,lty  = 2)}	
					} else {
						low <- subset_plot$upperCI[i] - subset_plot$start_date[i]
						high <- subset_plot$lowerCI[i] - subset_plot$start_date[i]
						lines(x = rep(subset_plot$start_date[i],2), y = c(low,high), col = ifelse(0>high & 0 < low,"darkslategray4","firebrick1"), lwd = 2)
						abline(h=0,lty  = 2)
				}	
			}
		}
		dev.off()
	}
}

##a### CHECK COMPUTATION TIMES
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





