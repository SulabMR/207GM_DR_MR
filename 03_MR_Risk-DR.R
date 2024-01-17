# pos /share2/pub/zhenggw/zhenggw/TwoSampleMR/eyedisease/NC/207GM-4DR/Riskresult/
# finnGen datasets H7_RETINOPATHYDIAB (labelled as DR1), DM_RETINOPATHY_EXMORE (labelled as DR2) and DM_RETINOPATHY (labelled as DR3)

# env MR
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
setwd("/share2/pub/zhenggw/zhenggw/TwoSampleMR/eyedisease/NC/207GM-4DR/Riskresult/")
riskfactors_ids<-as.character(unlist(read.table("/share2/pub/zhenggw/zhenggw/TwoSampleMR/eyedisease/NC/207GM-4DR/supp_riskfactors/riskfactors_id.txt",header=F)))
out_gwas_list <- c("finn-b-DM_RETINOPATHY","finn-b-H7_RETINOPATHYDIAB","finn-b-DM_RETINOPATHY_EXMORE")

for(i in riskfactors_ids){
	riskfactors <- gsub("_","-",i)
	expo_dat <- extract_instruments(outcomes=i,p1=1e-5,r2 = 0.1,kb = 1000,access_token = NULL)
	for(j in out_gwas_list){
			DR <- gsub("_","-",j)
		    res = tryCatch({
    		outcome_dat<- extract_outcome_data(snps = expo_dat$SNP,outcomes = j)
			dat <- harmonise_data(exposure_dat = expo_dat,outcome_dat = outcome_dat)
			dat2 <- dat[dat$mr_keep,]
				if(length(dat2$SNP) == 1){
					res <- mr(dat, method_list = c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				} else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
					res <- mr(dat, method_list = c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				} else if (length(dat2$SNP) > 3){
					res <- mr(dat, method_list = c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				}
	        }, warning = function(w){
			outcome_dat<- extract_outcome_data(snps = expo_dat$SNP,outcomes = j)
			dat <- harmonise_data(exposure_dat = expo_dat,outcome_dat = outcome_dat)
				if(length(dat2$SNP) == 1){
					res <- mr(dat, method_list = c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				} else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
					res <- mr(dat, method_list = c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				} else if (length(dat2$SNP) > 3){
					res <- mr(dat, method_list = c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
				}   
	        }, error = function(e){print("error exist!")},finally = {print("completed!")}
	        )

	  filename <- paste(riskfactors,"_",DR,"_res.txt",sep='')
	  write.table(res, file = filename, sep = "\t",col.names = TRUE,row.names = FALSE)
	  print(filename)
	}
}