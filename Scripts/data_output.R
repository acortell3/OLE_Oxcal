
#######################################
######## This scripts homogeneises different data outputs into one single df
#######################################

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
				"simID_seed" = rep(NA,length(oxcal_trap_hol)))

oxcal_unif_hol_df <- data.frame("Estimate" = rep(NA,length(oxcal_unif_hol)),
				"upperCI" = sapply(oxcal_unif_hol, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_unif_hol, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_unif_hol, function(x) x$start_date),
				"r" = sapply(oxcal_unif_hol, function(x) x$r),
				"Sd" = sapply(oxcal_unif_hol, function(x) x$Sd),
				"ESS" = sapply(oxcal_unif_hol, function(x) x$ESS),
				"sample_size" = sapply(oxcal_unif_hol, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_unif_hol)))

hol_df <- rbind(CRIWM_hol,OLE_medians_hol,OLE_resamp_caldate_hol,OLE_resamp_norm_hol,OLE_resamp_unif_hol,oxcal_trap_hol_df,oxcal_unif_hol_df)

hol_df$method <- c(rep("CRIWM",nrow(CRIWM_hol)),rep("OLE_medians",nrow(OLE_medians_hol)),rep("OLE_resamp_caldate",nrow(OLE_resamp_caldate_hol)),rep("OLE_resamp_norm",nrow(OLE_resamp_norm_hol)),rep("OLE_resamp_unif",nrow(OLE_resamp_unif_hol)),rep("oxcal_trap",nrow(oxcal_trap_hol_df)),rep("oxcal_unif",nrow(oxcal_unif_hol_df)))

hol_df$Period <- rep("Holocene",nrow(hol_df))

## The error had turned some columns to strings. Now we put all numerics back to numeric
hol_df[] <- lapply(hol_df, function(x){
			   y <- suppressWarnings(as.numeric(x))
			   if (all(is.na(x) == is.na(y))) y else x})
## Discard NAs
#hol_df <- hol_df[complete.cases(hol_df[,c(2:8)]),]
## Remove if there are negative values
#hol_df <- hol_df[apply(hol_df[,c(2:8)] >= 0, 1, all),]


### PLEISTOCENE
## Load rds
CRIWM_upl <- readRDS(paste0(path,"CRIWM_upl.rds"))
OLE_medians_upl <- readRDS(paste0(path,"OLE_medians_upl.rds"))
OLE_resamp_caldate_upl <- readRDS(paste0(path,"OLE_resamp_caldate_upl.rds"))
OLE_resamp_norm_upl <- readRDS(paste0(path,"OLE_resamp_norm_upl.rds"))
OLE_resamp_unif_upl <- readRDS(paste0(path,"OLE_resamp_unif_upl.rds"))
oxcal_trap_upl <- readRDS(paste0(path,"oxcal_trap_upl.rds"))
oxcal_unif_upl <- readRDS(paste0(path,"oxcal_unif_upl.rds"))

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
				"simID_seed" = rep(NA,length(oxcal_trap_upl)))

oxcal_unif_upl_df <- data.frame("Estimate" = rep(NA,length(oxcal_unif_upl)),
				"upperCI" = sapply(oxcal_unif_upl, function(x) x$posteriorRange[1,1]),
				"lowerCI" = sapply(oxcal_unif_upl, function(x) x$posteriorRange[1,2]),
				"start_date" = sapply(oxcal_unif_upl, function(x) x$start_date),
				"r" = sapply(oxcal_unif_upl, function(x) x$r),
				"Sd" = sapply(oxcal_unif_upl, function(x) x$Sd),
				"ESS" = sapply(oxcal_unif_upl, function(x) x$ESS),
				"sample_size" = sapply(oxcal_unif_upl, function(x) x$ESS),
				"simID_seed" = rep(NA,length(oxcal_unif_upl)))

upl_df <- rbind(CRIWM_upl,OLE_medians_upl,OLE_resamp_caldate_upl,OLE_resamp_norm_upl,OLE_resamp_unif_upl,oxcal_trap_upl_df,oxcal_unif_upl_df)

upl_df$method <- c(rep("CRIWM",nrow(CRIWM_upl)),rep("OLE_medians",nrow(OLE_medians_upl)),rep("OLE_resamp_caldate",nrow(OLE_resamp_caldate_upl)),rep("OLE_resamp_norm",nrow(OLE_resamp_norm_upl)),rep("OLE_resamp_unif",nrow(OLE_resamp_unif_upl)),rep("oxcal_trap",nrow(oxcal_trap_upl_df)),rep("oxcal_unif",nrow(oxcal_unif_upl_df)))

upl_df$Period <- rep("Pleistocene",nrow(upl_df))

## All together
full_output <- rbind(hol_df,upl_df)
saveRDS(full_output,"../Results/full_output.rds")



