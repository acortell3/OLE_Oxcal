

####### Script to plot and visualise the results of the paper <whateva>

### SELECT BEST OLE ACCORDING TU NUMBER OF DATES

#####################
####### HOLOCENE

string_find <- c("medians","resamp_caldate","resamp_norm","resamp_unif")
string_title <- c("medians","resamp cal curve","resamp normal","resamp uniform")

for (k in 1:length(string_find)){
	OLE_hol <- readRDS(paste0("../Results/OLE_",string_find[k],"_hol.rds"))
	
	cs <- list()
	index_nd <- unique(OLE_hol$sampled_dates)
	index_r <- unique(OLE_hol$r)
	index_sd <- unique(OLE_hol$Sd)
	
	png(paste0("../Figures/Figure_OLE_Hol_",k,".png"), res = 100, height = 1000, width = 1000)
	par(mfrow = c(4,2))
	for (j in 1:length(index_r)){
		for (h in 1:length(index_sd)){
			for (i in 1:length(index_nd)){
				cs[[i]] <- OLE_hol[OLE_hol$sampled_dates == index_nd[i] & OLE_hol$r == index_r[j] & OLE_hol$Sd == index_sd[h],]
				cs[[i]]$Estimate_norm <- cs[[i]]$start_date - cs[[i]]$Estimate
				cs[[i]]$upperCI_norm <- cs[[i]]$start_date - cs[[i]]$upperCI
				cs[[i]]$lowerCI_norm <- cs[[i]]$start_date - cs[[i]]$lowerCI
			}
			
			plot(x=1, y = mean(cs[[1]]$Estimate_norm), pch = 16, col = "blue", ylim = c(-600,100), xlim = c(0,5), ylab = "Estimate", xlab = "ndates", xaxt = "n", main = paste0("Holocene OLE ",string_title[k],", sd = ",index_sd[h]," r = ",index_r[j]))
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


#####################
####### UPPER PALAEOLITHIC

string_find <- c("medians","resamp_caldate","resamp_norm","resamp_unif")
string_title <- c("medians","resamp cal curve","resamp normal","resamp uniform")

for (k in 1:length(string_find)){
	OLE_up <- readRDS(paste0("../Results/OLE_",string_find[k],"_up.rds"))
	
	cs <- list()
	index_nd <- unique(OLE_up$sampled_dates)
	index_r <- unique(OLE_up$r)
	index_sd <- unique(OLE_up$Sd)
	
	png(paste0("../Figures/Figure_OLE_up_",k,".png"), res = 100, height = 1000, width = 1000)
	par(mfrow = c(4,2))
	for (j in 1:length(index_r)){
		for (h in 1:length(index_sd)){
			for (i in 1:length(index_nd)){
				cs[[i]] <- OLE_up[OLE_up$sampled_dates == index_nd[i] & OLE_up$r == index_r[j] & OLE_up$Sd == index_sd[h],]
				cs[[i]]$Estimate_norm <- cs[[i]]$start_date - cs[[i]]$Estimate
				cs[[i]]$upperCI_norm <- cs[[i]]$start_date - cs[[i]]$upperCI
				cs[[i]]$lowerCI_norm <- cs[[i]]$start_date - cs[[i]]$lowerCI
			}
			
			plot(x=1, y = mean(cs[[1]]$Estimate_norm), pch = 16, col = "blue", ylim = c(-900,100), xlim = c(0,5), ylab = "Estimate", xlab = "ndates", xaxt = "n", main = paste0("Upper Palaeolithic OLE ",string_title[k],", sd = ",index_sd[h]," r = ",index_r[j]))
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


######################## METHOD COMPARISON #############################
ox_trap <- readRDS("../Results/oxcal_trap_hol.rds")
ox_uni <- readRDS("../Results/oxcal_unif_hol.rds")
dates <- readRDS("../Simu_data/simuls_hol.rds")

single_val <- function(x){
	if (nrow(x) == 1){
		return(x)} else {
			x <- c(max(x[[1]], na.rm = T),
			       min(x[[2]], na.rm = T),
			       sum(x[[3]], na.rm = T))
			return(x)
	}
}

## Rearrange trapezoid Oxcal into df
red_ox_trap <- lapply(ox_trap,function(x){x[[3]] <- single_val(x[[3]])})
date_par_vals <- unique(dates[,-1])
trap_vals <- do.call(rbind,red_ox_trap)
oxcal_trap_vals <- cbind(trap_vals,date_par_vals)

oxcal_trap_vals$Start_norm <- oxcal_trap_vals$Start - oxcal_trap_vals$start_date
oxcal_trap_vals$End_norm <- oxcal_trap_vals$End - oxcal_trap_vals$start_date

## Rearrange uniform Oxcal into df
red_ox_uni <- lapply(ox_uni,function(x){x[[3]] <- single_val(x[[3]])})
date_par_vals <- unique(dates[,-1])
uni_vals <- do.call(rbind,red_ox_uni)
oxcal_uni_vals <- cbind(uni_vals,date_par_vals)

oxcal_uni_vals$Start_norm <- oxcal_uni_vals$Start - oxcal_uni_vals$start_date
oxcal_uni_vals$End_norm <- oxcal_uni_vals$End - oxcal_uni_vals$start_date

cs <- list(oxcal_trap_vals,oxcal_uni_vals)
#index_nd <- unique(OLE_up$sampled_dates)


for (k in 1:length(string_find)){
	OLE_hol <- readRDS(paste0("../Results/OLE_",string_find[k],"_hol.rds"))
	OLE_hol$Start_norm <- OLE_hol$upperCI - OLE_hol$start_date
	OLE_hol$Start_norm <- OLE_hol$lowerCI - OLE_hol$start_date
	OLE_hol <- OLE_hol[OLE_hol$sampled_dates == 10,]
	cs[[k+2]] <- OLE_hol
}


index_r <- unique(oxcal_trap_vals$r)
index_sd <- unique(oxcal_trap_vals$Sd)

plot(x=rep(1,2), y = rep(0,2), col = "white",xlim = c(0,8), ylim = c(-150,150),ylab = "Estimate", xlab = "ndates", xaxt = "n", main = paste0("Oxcal Trapezoid"))
axis(side = 1, at = seq(0,8), labels = c("","OxcTrap","OxcUnif","CRIWM","OLE_med","OLE_res","OLE_norm","OLE_unif",""))
abline(h = 0, lty = 2, col = "black")
for (i in 1:length(cs)){
	obj <- cs[[i]][cs[[i]]$r == 0.01 & cs[[i]]$Sd == 20,]
	lines(x=rep(i,2), y = c(mean(obj$Start_norm, na.rm = T),mean(obj$End_norm, na.rm = T)), lwd = 2, col = "dodgerblue")
	lines(x=c(i-0.1,i+0.1), y = c(mean(obj$Start_norm, na.rm = T),mean(obj$Start_norm, na.rm = T)), lwd = 2, col = "dodgerblue4")
	lines(x=c(i-0.1,i+0.1), y = c(mean(obj$End_norm, na.rm = T),mean(obj$End_norm, na.rm = T)), lwd = 2, col = "dodgerblue4")
}




## Explore structure from CRIWM
criwm <- readRDS("../Results/CRIWM_hol.rds")


