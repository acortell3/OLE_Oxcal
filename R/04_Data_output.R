
######### CODE FOR "UNCERTAIN BEGINNINGS: A COMPARISON OF THE ACCURACY AND PRECISION OF METHODS ESTIMATING EXTREME CHRONOLOGICAL EVENTS", BY A.CORTELL-NICOLAU, E. R. CREMA, AND A. KEY


######### SCRIPT COMBINING ALL THE METHODS PROPOSED

#### Authors: Alfredo Cortell-Nicolau & Enrico R. Crema

#######################################
######## This scripts homogeneises different data outputs into one single df
######################################

## Common utilities
path <- "../Results/"

### HOLOCENE
## Load rds

CRIWM_hol <- readRDS(paste0(path,"CRIWM_hol.rds"))
OLE_medians_hol <- readRDS(paste0(path,"OLE_medians_hol.rds"))
OLE_resamp_caldate_hol <- readRDS(paste0(path,"OLE_resamp_caldate_hol.rds"))
OLE_resamp_norm_hol <- readRDS(paste0(path,"OLE_resamp_norm_hol.rds"))
OLE_resamp_unif_hol <- readRDS(paste0(path,"OLE_resamp_unif_hol.rds"))
oxcal_trap_hol <- readRDS(paste0(path,"oxcal_trap_hol.rds"))
oxcal_unif_hol <- readRDS(paste0(path,"oxcal_unif_hol.rds"))

## Some OLEs are out of range. Discard those 
## Values to discard if error message
keep <- c("10","40","80")
OLE_medians_hol <- OLE_medians_hol[OLE_medians_hol$sample_size %in% keep,]
OLE_resamp_caldate_hol <- OLE_resamp_caldate_hol[OLE_resamp_caldate_hol$sample_size %in% keep,]
OLE_resamp_norm_hol <- OLE_resamp_norm_hol[OLE_resamp_norm_hol$sample_size %in% keep,]
OLE_resamp_unif_hol <- OLE_resamp_unif_hol[OLE_resamp_unif_hol$sample_size %in% keep,]

CRIWM_hol$overallAgreementIndex <- OLE_medians_hol$overallAgreementIndex <- OLE_resamp_caldate_hol$overallAgreementIndex <- OLE_resamp_norm_hol$overallAgreementIndex <- OLE_resamp_unif_hol$overallAgreementIndex <- NA

CRIWM_hol$nlegs <- OLE_medians_hol$nlegs <- OLE_resamp_caldate_hol$nlegs <- OLE_resamp_norm_hol$nlegs <- OLE_resamp_unif_hol$nlegs <- NA

## Builds oxcal dfs
## Some Oxcal results failed with "Error in BCAD(post.df[,2]: 0 BC/AD is not a valid year". Subset only the ones for which the time range was valid
oxcal_trap_hol <- oxcal_trap_hol[sapply(oxcal_trap_hol,length)>1]
oxcal_unif_hol <- oxcal_unif_hol[sapply(oxcal_unif_hol,length)>1]

## Oxca dfs
oxcal_trap_hol_df <- data.frame("Estimate" = rep(NA,length(oxcal_trap_hol)),
				"upperCI" = sapply(oxcal_trap_hol, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_trap_hol, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_trap_hol, function(x) x$start_date),
				"r" = sapply(oxcal_trap_hol, function(x) x$r),
				"Sd" = sapply(oxcal_trap_hol, function(x) x$Sd),
				"ESS" = sapply(oxcal_trap_hol, function(x) x$ESS),
				"sample_size" = sapply(oxcal_trap_hol, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_trap_hol)),
				"overallAgreementIndex" = sapply(oxcal_trap_hol, function(x) x$overallAgreement),
				"nlegs" = sapply(oxcal_trap_hol, function(x) nrow(x$posteriorRange)))

oxcal_unif_hol_df <- data.frame("Estimate" = rep(NA,length(oxcal_unif_hol)),
				"upperCI" = sapply(oxcal_unif_hol, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_unif_hol, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_unif_hol, function(x) x$start_date),
				"r" = sapply(oxcal_unif_hol, function(x) x$r),
				"Sd" = sapply(oxcal_unif_hol, function(x) x$Sd),
				"ESS" = sapply(oxcal_unif_hol, function(x) x$ESS),
				"sample_size" = sapply(oxcal_unif_hol, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_unif_hol)),
				"overallAgreementIndex" = sapply(oxcal_unif_hol, function(x) x$overallAgreement),
				"nlegs" = sapply(oxcal_unif_hol, function(x) nrow(x$posteriorRange)))

hol_df <- rbind(CRIWM_hol,OLE_medians_hol,OLE_resamp_caldate_hol,OLE_resamp_norm_hol,OLE_resamp_unif_hol,oxcal_trap_hol_df,oxcal_unif_hol_df)

hol_df$method <- c(rep("CRIWM",nrow(CRIWM_hol)),rep("OLE_medians",nrow(OLE_medians_hol)),rep("OLE_resamp_caldate",nrow(OLE_resamp_caldate_hol)),rep("OLE_resamp_norm",nrow(OLE_resamp_norm_hol)),rep("OLE_resamp_unif",nrow(OLE_resamp_unif_hol)),rep("oxcal_trap",nrow(oxcal_trap_hol_df)),rep("oxcal_unif",nrow(oxcal_unif_hol_df)))

hol_df$Period <- rep("Holocene",nrow(hol_df))

## The error had turned some columns to strings. Now we put all numerics back to numeric
hol_df[] <- lapply(hol_df, function(x){
			   y <- suppressWarnings(as.numeric(x))
			   if (all(is.na(x) == is.na(y))) y else x})

### PLEISTOCENE
## Load rds
CRIWM_upl <- readRDS(paste0(path,"CRIWM_upl.rds"))
OLE_medians_upl <- readRDS(paste0(path,"OLE_medians_upl.rds"))
OLE_resamp_caldate_upl <- readRDS(paste0(path,"OLE_resamp_caldate_upl.rds"))
OLE_resamp_norm_upl <- readRDS(paste0(path,"OLE_resamp_norm_upl.rds"))
OLE_resamp_unif_upl <- readRDS(paste0(path,"OLE_resamp_unif_upl.rds"))
oxcal_trap_upl <- readRDS(paste0(path,"oxcal_trap_upl.rds"))
oxcal_unif_upl <- readRDS(paste0(path,"oxcal_unif_upl.rds"))


CRIWM_upl$overallAgreementIndex <- OLE_medians_upl$overallAgreementIndex <- OLE_resamp_caldate_upl$overallAgreementIndex <- OLE_resamp_norm_upl$overallAgreementIndex <- OLE_resamp_unif_upl$overallAgreementIndex <- NA

CRIWM_upl$nlegs <- OLE_medians_upl$nlegs <- OLE_resamp_caldate_upl$nlegs <- OLE_resamp_norm_upl$nlegs <- OLE_resamp_unif_upl$nlegs <- NA

## Builds oxcal dfs
## Oxca dfs
oxcal_trap_upl_df <- data.frame("Estimate" = rep(NA,length(oxcal_trap_upl)),
				"upperCI" = sapply(oxcal_trap_upl, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_trap_upl, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_trap_upl, function(x) x$start_date),
				"r" = sapply(oxcal_trap_upl, function(x) x$r),
				"Sd" = sapply(oxcal_trap_upl, function(x) x$Sd),
				"ESS" = sapply(oxcal_trap_upl, function(x) x$ESS),
				"sample_size" = sapply(oxcal_trap_upl, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_trap_upl)),
				"overallAgreementIndex" = sapply(oxcal_trap_upl, function(x) x$overallAgreement),
				"nlegs" = sapply(oxcal_trap_upl, function(x) nrow(x$posteriorRange)))

oxcal_unif_upl_df <- data.frame("Estimate" = rep(NA,length(oxcal_unif_upl)),
				"upperCI" = sapply(oxcal_unif_upl, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_unif_upl, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_unif_upl, function(x) x$start_date),
				"r" = sapply(oxcal_unif_upl, function(x) x$r),
				"Sd" = sapply(oxcal_unif_upl, function(x) x$Sd),
				"ESS" = sapply(oxcal_unif_upl, function(x) x$ESS),
				"sample_size" = sapply(oxcal_unif_upl, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_unif_upl)),
				"overallAgreementIndex" = sapply(oxcal_unif_upl, function(x) x$overallAgreement),
				"nlegs" = sapply(oxcal_unif_upl, function(x) nrow(x$posteriorRange)))

upl_df <- rbind(CRIWM_upl,OLE_medians_upl,OLE_resamp_caldate_upl,OLE_resamp_norm_upl,OLE_resamp_unif_upl,oxcal_trap_upl_df,oxcal_unif_upl_df)

upl_df$method <- c(rep("CRIWM",nrow(CRIWM_upl)),rep("OLE_medians",nrow(OLE_medians_upl)),rep("OLE_resamp_caldate",nrow(OLE_resamp_caldate_upl)),rep("OLE_resamp_norm",nrow(OLE_resamp_norm_upl)),rep("OLE_resamp_unif",nrow(OLE_resamp_unif_upl)),rep("oxcal_trap",nrow(oxcal_trap_upl_df)),rep("oxcal_unif",nrow(oxcal_unif_upl_df)))

upl_df$Period <- rep("Pleistocene",nrow(upl_df))

## Add columns to account for oxcal's separated probability masses. The number 8/16 is obtained after unique(full_output$nlegs)

## HOLOCENE
ind_cis <- 1
for (i in ncol(hol_df):(ncol(hol_df)+7)){
	hol_df$newcol <- NA
	colnames(hol_df)[ncol(hol_df)] <- paste0("upperCI",ind_cis)
	hol_df$newcol <- NA
	colnames(hol_df)[ncol(hol_df)] <- paste0("lowerCI",ind_cis)
	ind_cis <- ind_cis + 1
}

## Remove everything with legs > 1
hol_df_1 <- hol_df[hol_df$nlegs == 1 | is.na(hol_df$nlegs),]
hol_df_mult <- hol_df[which(hol_df$nlegs > 1),]

## Subset only sims where posteriorRange > 1
oxcal_trap_hol_multiple <- Filter(function(x) nrow(x$posteriorRange) > 1,oxcal_trap_hol)
oxcal_unif_hol_multiple <- Filter(function(x) nrow(x$posteriorRange) > 1,oxcal_unif_hol)

for (i in 1:length(oxcal_trap_hol_multiple)){
	pr <- oxcal_trap_hol_multiple[[i]]$posteriorRange
	pr <- pr[,c(1,2)]
	values <- as.vector(t(pr))
	## Put max min CIs, although this is wrong. Just for the recrod
	hol_df_mult[i,2] <- max(pr)
	hol_df_mult[i,3] <- min(pr)
	hol_df_mult[i,c(14:(13+length(values)))] <- values
}

ind <- 1
for (i in (length(oxcal_trap_hol_multiple)+1):nrow(hol_df_mult)){
	pr <- oxcal_unif_hol_multiple[[ind]]$posteriorRange
	pr <- pr[,c(1,2)]
	values <- as.vector(t(pr))
	## Put max min CIs, although this is wrong. Just for the recrod
	hol_df_mult[i,2] <- max(pr)
	hol_df_mult[i,3] <- min(pr)
	hol_df_mult[i,c(14:(13+length(values)))] <- values
	ind <- ind + 1
}

hol_df <- rbind(hol_df_1,hol_df_mult)

## PLEISTOCENE
ind_cis <- 1
for (i in ncol(upl_df):(ncol(upl_df)+7)){
	upl_df$newcol <- NA
	colnames(upl_df)[ncol(upl_df)] <- paste0("upperCI",ind_cis)
	upl_df$newcol <- NA
	colnames(upl_df)[ncol(upl_df)] <- paste0("lowerCI",ind_cis)
	ind_cis <- ind_cis + 1
}

## Remove everything with legs > 1
upl_df_1 <- upl_df[upl_df$nlegs == 1 | is.na(upl_df$nlegs),]
upl_df_mult <- upl_df[which(upl_df$nlegs > 1),]

## Subset only sims where posteriorRange > 1
oxcal_trap_upl_multiple <- Filter(function(x) nrow(x$posteriorRange) > 1,oxcal_trap_upl)
oxcal_unif_upl_multiple <- Filter(function(x) nrow(x$posteriorRange) > 1,oxcal_unif_upl)

for (i in 1:length(oxcal_trap_upl_multiple)){
	pr <- oxcal_trap_upl_multiple[[i]]$posteriorRange
	pr <- pr[,c(1,2)]
	values <- as.vector(t(pr))
	## Put max min CIs, although this is wrong. Just for the recrod
	upl_df_mult[i,2] <- max(pr)
	upl_df_mult[i,3] <- min(pr)
	upl_df_mult[i,c(14:(13+length(values)))] <- values
}

ind <- 1
for (i in (length(oxcal_trap_upl_multiple)+1):nrow(upl_df_mult)){
	pr <- oxcal_unif_upl_multiple[[ind]]$posteriorRange
	pr <- pr[,c(1,2)]
	values <- as.vector(t(pr))
	## Put max min CIs, although this is wrong. Just for the recrod
	upl_df_mult[i,2] <- max(pr)
	upl_df_mult[i,3] <- min(pr)
	upl_df_mult[i,c(14:(13+length(values)))] <- values
	ind <- ind + 1
}

upl_df <- rbind(upl_df_1,upl_df_mult)

## All together
full_output <- rbind(hol_df,upl_df)
saveRDS(full_output,"../Results/full_output.rds")



