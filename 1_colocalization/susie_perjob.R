library(coloc)
library(data.table)
library(dplyr)
library(tidyr)
library(optparse)

option_list <- list(
  make_option(c("-j","--job"), type = "numeric", default = FALSE,
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

setwd("/mnt/data/chenmingyang/project/CAD_pathway_analysis")

cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")
gwas_trait <- readxl::read_excel("0_dataprepare/gwas_trait_raw.xlsx")
region_coloc <- read.table("0_dataprepare/region_coloc.txt",header=T)
valid_snp <- read.table("1_choose_trait_by_coloc/valid_snp_coloc.txt",header=F)

job_id <- opt$job
if(job_id < floor(nrow(region_coloc)/1)+1){
  job_vars <- region_coloc[(1*(job_id-1)+1):(1*job_id),]
}else{
  job_vars <- region_coloc[(1*(job_id-1)+1):nrow(region_coloc),]
}

trait_file <- list.files("gwas_regional/")
trait_file <- paste0("gwas_regional/",trait_file)

for(i in 1:nrow(job_vars)){
  result_susie <- data.frame()
  chr <- job_vars[i,]$chr
  start <- job_vars[i,]$start
  end <- job_vars[i,]$end
  cad_tested <- cad_cleaned[as.numeric(cad_cleaned$CHR) == as.numeric(chr) & as.numeric(cad_cleaned$BP) >= as.numeric(start) & 
                              as.numeric(cad_cleaned$BP) <= as.numeric(end), ]
  
  tmp <- read.table(paste0("1_colocalization/coloc_LD_matrix/region_coloc_",job_vars[i,]$No,".txt"), header = T)
  tmp <- as.matrix(tmp)
  rownames(tmp) <- colnames(tmp)
  tmp <- tmp[, !colSums(is.na(tmp)) == nrow(tmp)]
  tmp <- tmp[!rowSums(is.na(tmp)) == ncol(tmp), ]
  
  cad_tested <- cad_tested[cad_tested$SNP %in% rownames(tmp), ]
  cad_tested <- cad_tested[cad_tested$SNP %in% valid_snp$V1,]
  #pre converge
  tmp3 <- tmp[rownames(tmp) %in% cad_tested$SNP, colnames(tmp) %in% cad_tested$SNP]
  
  D1 <- list(beta = cad_tested$BETA,
             varbeta = (cad_tested$SE)^2,
             snp = cad_tested$SNP,
             position = cad_tested$BP,
             type = "cc",
             N = 530185,
             LD = tmp3)
  S1 <- try(runsusie(D1,maxit=1000, repeat_until_convergence = F), silent = T)
  if(class(S1) == "try-error"){
    tmp_info <- data.frame(
      job_id = job_vars[i,]$No,
      S1 = "no_converge",
      trait = NA,
      S2 = NA
    )
    write.table(tmp_info,paste0("1_colocalization/susie_result/susie_region_",job_vars[i,]$No,".log"),quote=F,row.names=F,col.names=T)
    break
  }else if(is.null(summary(S1)$cs)){
    tmp_info <- data.frame(
      job_id = job_vars[i,]$No,
      S1 = "no_credible_sets",
      trait = NA,
      S2 = NA
    )
    write.table(tmp_info,paste0("1_colocalization/susie_result/susie_region_",job_vars[i,]$No,".log"),quote=F,row.names=F,col.names=T)
    break
  }
  
  log_info <- data.frame()
  
  for(j in trait_file){
    trait_gwas_raw <- fread(j) 
    trait_gwas<- trait_gwas_raw[trait_gwas_raw$VAR_ID %in% cad_tested$VAR_ID,]
    trait_gwas$FRQ <- NULL
    trait_gwas <- inner_join(trait_gwas,cad_tested[,c("SNP","VAR_ID","FRQ")], by=c("VAR_ID","SNP"))
    
    trait_gwas <- trait_gwas %>%
      separate(VAR_ID, into=c("CHR","BP","REF","ALT"),sep="_",remove = F)
    
    final_tested <- cad_tested[cad_tested$VAR_ID %in% trait_gwas$VAR_ID,]
    
    if(nrow(final_tested)==0){next}
    
    trait_name <- sub(".txt.gz","",j)
    trait_name <- sub("gwas_regional/","",trait_name)
    
    tmp2 <- tmp[rownames(tmp) %in% final_tested$SNP, colnames(tmp) %in% final_tested$SNP]
    
    D1 <- list(beta = final_tested$BETA,
               varbeta = (final_tested$SE)^2,
               snp = final_tested$SNP,
               position = final_tested$BP,
               type = "cc",
               N = 530185,
               LD = tmp2)
    S1 <- try(runsusie(D1,maxit=1000, repeat_until_convergence = F), silent = T)
    
    if(class(S1) == "try-error"){
      tmp_info <- data.frame(
        job_id = job_vars[i,]$No,
        S1 = "no_converge",
        trait = trait_name,
        S2 = NA
      )
      log_info <- rbind(log_info, tmp_info)
      next
    }else if(is.null(summary(S1)$cs)){
      tmp_info <- data.frame(
        job_id = job_vars[i,]$No,
        S1 = "no_credible_sets",
        trait = trait_name,
        S2 = NA
      )
      log_info <- rbind(log_info, tmp_info)
      next
    }
    
    D2 <- list(beta = trait_gwas$BETA,
               varbeta = (trait_gwas$SE)^2,
               snp = trait_gwas$SNP,
               position = trait_gwas$BP,
               type = "quant",
               N = as.integer(gwas_trait[gwas_trait$trait == trait_name, ]$samplesize),
               MAF = ifelse(trait_gwas$FRQ > 0.5, 1- trait_gwas$FRQ, trait_gwas$FRQ),
               LD = tmp2)
    
    S2 <- try(runsusie(D2,maxit=1000, repeat_until_convergence = F), silent = T)
    
    if(class(S2) == "try-error"){
      tmp_info <- data.frame(
        job_id = job_vars[i,]$No,
        S1 = "converged",
        trait = trait_name,
        S2 = "no_converge"
      )
      log_info <- rbind(log_info, tmp_info)
      next
    }else if(is.null(summary(S2)$cs)){
      tmp_info <- data.frame(
        job_id = job_vars[i,]$No,
        S1 = "converged",
        trait = trait_name,
        S2 = "no_credible_sets"
      )
      log_info <- rbind(log_info, tmp_info)
      next
    }else{
      tmp_info <- data.frame(
        job_id = job_vars[i,]$No,
        S1 = "converged",
        trait = trait_name,
        S2 = "converged"
      )
      log_info <- rbind(log_info, tmp_info)
    }
    
    res=coloc.susie(S1,S2)
    tmp2 <- res$summary
    tmp2 <- as.data.frame(tmp2)
    tmp2$trait <- trait_name
    idx1 <- match(tmp2$hit1, final_tested$SNP)
    idx2 <- match(tmp2$hit2, trait_gwas$SNP)
    tmp2$p1 <- final_tested[idx1,]$P_VALUE
    tmp2$p2 <- trait_gwas[idx2,]$P_VALUE
    tmp2$beta1 <- final_tested[idx1,]$BETA
    tmp2$beta2 <- trait_gwas[idx2,]$BETA
    tmp2$se1 <- final_tested[idx1,]$SE
    tmp2$se2 <- trait_gwas[idx2,]$SE
    tmp2$region <- job_vars[i,]$No
    result_susie <- rbind(result_susie, tmp2)
  }
  write.table(result_susie,paste0("1_colocalization/susie_result/susie_region_",job_vars[i,]$No,".txt"),quote=F,row.names=F,col.names=T)
  write.table(log_info,paste0("1_colocalization/susie_result/susie_region_",job_vars[i,]$No,".log"),quote=F,row.names=F,col.names=T)
}
