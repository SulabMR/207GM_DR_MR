# combined results
# MR PRESSO + Egger_intercept (Pval)+Q_statistics (Pval)+Heterogeneity_I2
# 可以在/share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps 读入处理好的GM数据

library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)

GMlist <- c("g__Collinsella","g__Burkholderiales_noname",
"f__Burkholderiales_noname","s__Burkholderiales_bacterium_1_1_47",
"s__Collinsella_aerofaciens","s__Streptococcus_salivarius","s__Bacteroides_faecis",
"s__Ruminococcus_torques")

DRlist <- c(
        "finn-b-H7_RETINOPATHYDIAB",
        "finn-b-DM_RETINOPATHY_EXMORE",
        "finn-b-DM_RETINOPATHY")

data2 <- c()
for(i in GMlist){
	file1 <- paste0(i,"_1kr0.1.txt")
	exposure <- read.table(file1)
	exposure$id.exposure <- i
	for(j in DRlist){
				outcome <- extract_outcome_data(snps = exposure$SNP,outcomes = j)
      	dat <- harmonise_data(exposure_dat = exposure,outcome_dat = outcome)
      	dat2 <- dat[dat$mr_keep,]
        if(length(dat2$SNP) == 1){
        	method_list1 <- c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
			res <- mr(dat, method_list = method_list1)
       		het <- mr_heterogeneity(dat,method_list = method_list1[1])
       		plt <- mr_pleiotropy_test(dat)
        } else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
        	method_list2 <- c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
         	res <- mr(dat, method_list = method_list2)
        	het <- mr_heterogeneity(dat,method_list = method_list2[1])
       		plt <- mr_pleiotropy_test(dat)
        } else if (length(dat2$SNP) > 3){
        	method_list3 <- c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
         	res <- mr(dat, method_list = method_list3)
        	het <- mr_heterogeneity(dat,method_list = method_list3[1])
       		plt <- mr_pleiotropy_test(dat)
        }

		if(dim(het)[2] == 8){
			het <- het[,5:8]
			aa <- merge(res[1,],het,by.x="method",by.y="method",all.x=TRUE)
		}else{
			het <- data.frame(res$method[1],NA,NA,NA)
			colnames(het) <- c("method","Q","Q_df","Q_pval")
			aa <- merge(res[1,],het,by.x="method",by.y="method",all.x=TRUE)
		}


		if(dim(plt)[2] == 7){
			plt <- plt[,c(5,7)]
			plt$method <- het$method
			colnames(plt) <- lapply(colnames(plt), FUN=function(x) paste0("plt_",x))
			aa <- merge(aa,plt,by.x="method",by.y="plt_method",all.x=TRUE)
		}else{
			plt <- data.frame(res$method[1],NA,NA)
			colnames(plt) <- c("het_method","plt_egger_intercept","plt_pval")
			aa <- merge(aa,plt,by.x="method",by.y="het_method",all.x=TRUE)
		}
			aa <- aa[,-c(4,5)]
	
			res_pre = tryCatch({res_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdOutcome = "se.outcome", SdExposure = "se.exposure",OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat2,NbDistribution = 1000,  SignifThreshold = 0.05)
							 data.frame("mr_presso_beta" = res_presso$`Main MR results`$`Causal Estimate`[1],
									    "mr_presso_sd" = res_presso$`Main MR results`$`Sd`[1],
									    "mr_presso_P-value" = res_presso$`Main MR results`$`P-value`[1],
									    "GlobalTest_Pval" = res_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
							}, warning = function(w) {
							    message("warning")
							}, error = function(e) {
							    data.frame("mr_presso_beta" = NA,"mr_presso_sd" = NA,"mr_presso_P-value" = NA,"GlobalTest_Pval" = NA)
							}, finally = {

							})
		aa <- cbind(aa,res_pre)
		data2 <- rbind(data2,aa)
	}
}
data2 <- data2 %>%
	mutate(Heterogeneity_I2=ifelse(data2$Q - data2$Q_df > 0 , paste0(100*(data2$Q - data2$Q_df)/data2$Q,"%") , paste0(0,"%") ))
data2$`Diabetic Retinopathy`[data2$id.outcome == "finn-b-H7_RETINOPATHYDIAB"] <- "DR_1"
data2$`Diabetic Retinopathy`[data2$id.outcome == "finn-b-DM_RETINOPATHY_EXMORE"] <- "DR_2"
data2$`Diabetic Retinopathy`[data2$id.outcome == "finn-b-DM_RETINOPATHY"] <- "DR_3"


files <- list.files(pattern = "*_res.txt$")
data <- c()
for(i in files){
	indexs <- strsplit(i, "_finn", fixed= T)
	output <- read.table(i,header=T)
	if(dim(output)[2] == 9){
		output$outcome <- NULL
		output$exposure <- NULL
		output$id.exposure <- indexs[[1]][1]
		data <- rbind(data,output)
	}else{
		indexs1 <- strsplit(i, "_finn", fixed= T)
		output <- data.frame(NA,NA,NA,NA,NA,NA,NA)
		colnames(output) <- c("id.exposure","id.outcome","method","nsnp","b","se","pval")
		output$id.exposure <- indexs1[[1]][1]

		indexs1 <- strsplit(i, "_finn-b-", fixed= T)
		indexs1 <- gsub("_res.txt","",indexs1[[1]][2])
		indexs1 <- gsub("-","_",indexs1)
		output$id.outcome <- paste0("finn-b-",indexs1)
		output$nsnp <- as.numeric(0)
		data <- rbind(data,output)
	}
}
DR <- c("finn-b-H7_RETINOPATHYDIAB","finn-b-DM_RETINOPATHY_EXMORE","finn-b-DM_RETINOPATHY")
data <- data[which(data$id.outcome %in% DR),]
data$Zscore <- data$b/data$se
method <- c("Inverse variance weighted","MR Egger","Maximum likelihood","Weighted median","Weighted mode")
dat <- data[which(!(data$method %in% method)),]

disease <- unique(dat$id.outcome)

extract <- function(dat,gwascode){
  dat2 <- dat[which(dat$id.outcome == gwascode),]
  dat2 <- dat2[order(dat2$pval),]
  dat2$padj <- p.adjust(dat2$pval,method = "fdr")
  dat2
}


dat_adj <- data.frame()
for(i in disease) 
{
  tmp <- extract(dat,i)
  dat_adj <-rbind(dat_adj,tmp)
}


dat_adj <- dat_adj %>%
  mutate(match_id = paste(id.exposure,id.outcome,method, sep='_'))
data2 <- data %>%
  mutate(match_id = paste(id.exposure,id.outcome,method, sep='_'))
data_f <- data2 %>%
  left_join(dat_adj %>% dplyr::select(match_id, padj), by = 'match_id')

method <- c("Inverse variance weighted","MR Egger","Maximum likelihood","Weighted median","Weighted mode")
for(m in method){
	dat <- data[which(data$method %in% m),]
	disease <- unique(dat$id.outcome)

	dat_adj <- data.frame()
	for(i in disease) 
	{
	  tmp <- extract(dat,i)
	  dat_adj <-rbind(dat_adj,tmp)
	}

	dat_adj <- dat_adj %>%
	  mutate(match_id = paste(id.exposure,id.outcome,method, sep='_'))
	data_f <- data_f %>%
	  left_join(dat_adj %>% dplyr::select(match_id, padj), by = 'match_id')

	data_f$padj.x <- ifelse(is.na(data_f$padj.x), data_f$padj.y, data_f$padj.x)
	data_f$padj.y <- NULL
	colnames(data_f)[10] <- c("padj")
}

data_f$match_id <- NULL
data_f$`Diabetic Retinopathy`[data_f$id.outcome == "finn-b-H7_RETINOPATHYDIAB"] <- "DR_1"
data_f$`Diabetic Retinopathy`[data_f$id.outcome == "finn-b-DM_RETINOPATHY_EXMORE"] <- "DR_2"
data_f$`Diabetic Retinopathy`[data_f$id.outcome == "finn-b-DM_RETINOPATHY"] <- "DR_3"