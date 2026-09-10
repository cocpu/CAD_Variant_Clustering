########
#part learn from Inferring multi-organ genetic connections using imaging and clinical data through Mendelian randomization
########
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(mr.divw)
library(mr.raps)
library(GRAPPLE)

library(openxlsx)
library(optparse)

option_list <- list(
  make_option(c("-j","--job"), type = "character", default = FALSE,
              help="")
)

opt_parser = OptionParser(
  usage = "usage: %prog [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL , 
  description = ""
)

opt = parse_args(opt_parser)
i <- opt$job

setwd("/mnt/data/chenmingyang/project/CAD_pathway_analysis/")

snp_cluster <- read.table("2_bNMF/snp_cluster.txt", header = T)
snp_cluster$idxall <- 1
cad_cleaned <- readRDS("0_dataprepare/CAD_noUKB.rds")
cad_cleaned$N <- 129014
colnames(cad_cleaned)[3] <- "POS"
cad_cleaned <- cad_cleaned[cad_cleaned$SNP %in% snp_cluster$SNP,]
cad_cleaned <- cad_cleaned[as.numeric(cad_cleaned$P) < 1e-5,]
cad_cleaned <- cad_cleaned %>%
  dplyr::rename(REF=A1) %>%
  dplyr::rename(ALT=A2) %>%
  dplyr::rename(BP=POS) %>%
  mutate(Risk_Allele = ifelse(BETA>=0, ALT, REF)) %>%
  dplyr::select(SNP,CHR,BP,REF,ALT,BETA,SE,P,Risk_Allele,N) %>%
  dplyr::rename(P_VALUE = P)
cad_cleaned$VAR_ID <- paste0(cad_cleaned$CHR,"_",cad_cleaned$BP,"_",cad_cleaned$REF,"_",cad_cleaned$ALT)
temp <- readRDS("0_dataprepare/CAD_small.rds")
cad_cleaned <- left_join(cad_cleaned, temp[,c("VAR_ID","FRQ")])

trait_path <- fread("/mnt/data/chenmingyang/project/CAD_pathway_analysis/3_MR/trait.txt")

trait_gwas_full <- fread(trait_path[trait_path$trait == i,]$full_path)
trait_gwas_full$N <- as.numeric(trait_gwas_full$N_PH)
trait_gwas_full[is.na(trait_gwas_full$N),]$N <- trait_path[trait_path$trait == i,]$samplesize
trait_gwas_full$SNP <- NULL
trait_gwas_full <- left_join(trait_gwas_full, cad_cleaned[,c("SNP","CHR","BP","REF","ALT","VAR_ID")])
trait_gwas <- trait_gwas_full[trait_gwas_full$VAR_ID %in% snp_cluster$VAR_ID,]

exposure <- format_data( as.data.frame(cad_cleaned), type="exposure", snp_col = "SNP", se_col = "SE",
                         beta_col = "BETA",pval_col = "P_VALUE",min_pval = 1e-200,
                         chr_col = "CHR",pos_col = "POS",effect_allele_col = "ALT", 
                         other_allele_col = "REF",eaf_col = "FRQ",samplesize_col ="N")
exposure$exposure <- "CAD"
outcome <- format_data( as.data.frame(trait_gwas), type="outcome", snp_col = "SNP", se_col = "SE",
                        beta_col = "BETA",pval_col = "P_VALUE",min_pval = 1e-200,
                        chr_col = "CHR",pos_col = "POS",effect_allele_col = "ALT", 
                        other_allele_col = "REF",eaf_col = "FRQ",samplesize_col ="N")
outcome$outcome <- i
mydata_all <- harmonise_data(exposure_dat = exposure, 
                             outcome_dat = outcome)
mydata_all$eaf.outcome <- mydata_all$eaf.exposure

mydata_all$vg.exposure <- (mydata_all$beta.exposure)^2 * mydata_all$eaf.exposure * (1-mydata_all$eaf.exposure)
mydata_all$rsq.exposure <- mydata_all$vg.exposure/(mydata_all$vg.exposure + 3.29)

if(trait_path[trait_path$trait == i,]$type == "continuous"){
  mydata_all$rsq.outcome <- (get_r_from_bsen(mydata_all$beta.outcome,mydata_all$se.outcome, mydata_all$samplesize.outcome))^2
}else{
  mydata_all$vg.outcome <- (mydata_all$beta.outcome)^2 * mydata_all$eaf.outcome * (1-mydata_all$eaf.outcome)
  mydata_all$rsq.outcome <- mydata_all$vg.outcome/(mydata_all$vg.outcome + 3.29)
}

mydata_all$effective_n.exposure <- mydata_all$samplesize.exposure
mydata_all$effective_n.outcome <- mydata_all$samplesize.outcome

mydata_all <- steiger_filtering(mydata_all)
mydata_all <- mydata_all[mydata_all$steiger_dir == T,]

result_all <- data.frame()
loo_all <- data.frame()
for(j in c("idx2","idx4","idx7","idx8","idxall")){
  if(j != "idxall"){
    snps_need <- snp_cluster[snp_cluster[[j]] > 0,]$SNP
  }else{
    snps_need <- snp_cluster$SNP
  }
  
  mydata <- mydata_all[mydata_all$SNP %in% snps_need,]
  
  mr_methods <- c('mr_ivw_fe', 'mr_ivw_mre',"mr_egger_regression", 
                  'mr_weighted_mode','mr_weighted_median')
  
  res_mr <- mr(mydata, method_list = mr_methods)
  res_mr$snps <- paste(mydata$SNP,collapse = ",")
  
  res_DIVW <- mr.divw(mydata$beta.exposure, mydata$beta.outcome, mydata$se.exposure, mydata$se.outcome)
  res_DIVW$snps <- paste(mydata$SNP,collapse = ",")
  
  
  if (tryCatch({res_raps <- mr.raps.simple.robust(mydata$beta.exposure, mydata$beta.outcome, mydata$se.exposure, mydata$se.outcome,diagnostics = F);
  res_raps$snps <- paste(mydata$SNP,collapse = ",");
  print(FALSE)},
  error = function(e){print(TRUE)})){
    res_raps <- NULL
  }
  
  data_grapple <- rename(mydata,selection_pvals = pval.exposure,
                         gamma_exp1=beta.exposure, se_exp1=se.exposure,
                         gamma_out=beta.outcome, se_out=se.outcome)
  
  pv_grapple = as.numeric(1e-5)
  er <- tryCatch({diagnosis_grapple <-findModes(data_grapple,p.thres=pv_grapple,map.marker = FALSE);print(FALSE)},
                 error = function(e){print(TRUE)})
  if (er){
    diagnosis_grapple = NULL
  }
  
  source("3_MR/modified_grapple_function.R")
  if (tryCatch({res_grapple <- grappleRobustEst(data_grapple,p.thres = pv_grapple);
  res_grapple$snps <- paste(mydata$SNP,collapse = ",");
  print(FALSE)},
  error = function(e){print(TRUE)})){
    res_grapple <- NULL
  }
  
  res_extra <- data.frame(
    id.exposure = rep(res_mr$id.exposure[1],3),
    id.outcome = rep(res_mr$id.outcome[1],3),
    outcome = rep(res_mr$outcome[1],3),
    exposure = rep(res_mr$exposure[1],3),
    method = c('DIVW','GRAPPLE','mr_raps'),
    nsnp = c(res_DIVW$n.IV,nrow(data_grapple),nrow(mydata)),
    b = c(res_DIVW$beta.hat,
          ifelse(is.null(res_grapple),NA,res_grapple$beta.hat),
          ifelse(is.null(res_raps),NA,res_raps$beta.hat)),
    se = c(res_DIVW$beta.se,
           ifelse(is.null(res_grapple),NA,sqrt(res_grapple$beta.var)),
           ifelse(is.null(res_raps),NA,res_raps$beta.se)),
    pval = c(1-pnorm(abs(res_DIVW$beta.hat/res_DIVW$beta.se)),
             ifelse(is.null(res_grapple),NA,res_grapple$beta.p.value),
             ifelse(is.null(res_raps),NA,res_raps$beta.p.value)),
    snps = c(res_DIVW$snps, 
             ifelse(is.null(res_grapple),NA,res_grapple$snps),
             ifelse(is.null(res_raps),NA,res_raps$snps))
  )
  
  res <- rbind(res_mr, res_extra)
  res <- select(res, -c('id.exposure', 'id.outcome'))
  
  res$cluster <- j
  
  heterogeneity_test <- mr_heterogeneity(mydata)
  pleiotropy_test <- mr_pleiotropy_test(mydata)
  I2_Egger <- (heterogeneity_test$Q[1] - heterogeneity_test$Q_df[1] + 1) / heterogeneity_test$Q[1]
  I2_IVW <- (heterogeneity_test$Q[2] - heterogeneity_test$Q_df[2] + 1) / heterogeneity_test$Q[2]
  
  make_conclusion <- function(p_val){ifelse(p_val < 0.05, "Warning","Pass")}
  test_result <- data.frame(
    exposure = rep(res$exposure[1],10),
    outcome = rep(res$outcome[1],10),
    name = c('heterogenity', 'heterogenity','heterogenity', 'heterogenity','heterogenity',
             'directional pleiotropy',
             'directional pleiotropy',
             'directional pleiotropy',
             'outlier',
             'wrong direction'),
    test_methods = c('IVW re v.s IVW fixed', 
                     "MR-Egger Cochran's Q","IVW Cochran's Q",
                     "MR-Egger Higgin's I2","IVW Higgin's I2",
                     "MR-Egger intercept",
                     "funnel plot",
                     "number of mode in grapple",
                     'outlier detected in grapple',
                     'mode in grapple'),
    test = c(paste(res$b[res$method == "Inverse variance weighted (fixed effects)"],
                   res$b[res$method == "Inverse variance weighted (multiplicative random effects)"],
                   sep = ','),
             heterogeneity_test$Q_pval[1],
             heterogeneity_test$Q_pval[2],
             I2_Egger,
             I2_IVW,
             pleiotropy_test$pval,
             "See funnel plot",
             ifelse(is.null(res_grapple),NA,paste(diagnosis_grapple$modes,sep = " ")),
             ifelse(is.null(res_grapple),NA,nrow(res_grapple$outliers)),
             ifelse(is.null(res_grapple),NA,paste(diagnosis_grapple$modes,sep = " "))),
    suggestion = c("if differs too much, possible directional peiotropy and heterogenity",
                   make_conclusion(heterogeneity_test$Q_pval[1]),
                   make_conclusion(heterogeneity_test$Q_pval[2]),
                   "0-100, the higher the more heterogenity",
                   "0-100, the higher the more heterogenity",
                   make_conclusion(pleiotropy_test$pval),
                   "if asymmetric, then warning ",
                   ifelse(is.null(res_grapple),NA,ifelse(length(diagnosis_grapple$modes)>1 , 'warning','Pass')),
                   ifelse(is.null(res_grapple),NA,ifelse(nrow(res_grapple$outliers)>=1, 'warning', 'Pass')),
                   ifelse(is.null(res_grapple),NA,ifelse(length(diagnosis_grapple$modes)>1, 'warning','Pass'))
    )
  )
  
  test_tmp <- as.data.frame(matrix(test_result$test, nrow = 1))
  names(test_tmp) <- test_result$test_methods
  res_comb <- cbind(res, test_tmp)
  
  loo_res <- mr_leaveoneout(mydata, parameters = default_parameters(), method = mr_ivw_mre)
  loo_res$cluster <- j
  
  result_all <- rbind(result_all, res_comb)
  loo_all <- rbind(loo_all, loo_res)
}

write.csv(result_all, paste0("/mnt/data/chenmingyang/project/CAD_pathway_analysis/3_MR/result/",i,"_mr_result.csv"), row.names = F)
write.csv(loo_all, paste0("/mnt/data/chenmingyang/project/CAD_pathway_analysis/3_MR/result/",i,"_mr_loo.csv"), row.names = F)
