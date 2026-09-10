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

job_id <- opt$job
if(job_id < floor(nrow(region_coloc)/1)+1){
  job_vars <- region_coloc[(1*(job_id-1)+1):(1*job_id),]
}else{
  job_vars <- region_coloc[(1*(job_id-1)+1):nrow(region_coloc),]
}

trait_file <- list.files("gwas_regional/")
trait_file <- paste0("gwas_regional/",trait_file)

for(i in 1:nrow(job_vars)){
  result_coloc <- data.frame()
  chr <- job_vars[i,]$chr
  start <- job_vars[i,]$start
  end <- job_vars[i,]$end
  cad_tested <- cad_cleaned[as.numeric(cad_cleaned$CHR) == as.numeric(chr) & as.numeric(cad_cleaned$BP) >= as.numeric(start) & 
                              as.numeric(cad_cleaned$BP) <= as.numeric(end), ]
  
  for(j in trait_file){
    trait_gwas_raw <- fread(j) 
    trait_gwas<- trait_gwas_raw[trait_gwas_raw$VAR_ID %in% cad_tested$VAR_ID,]
    trait_gwas$FRQ <- NULL
    trait_gwas <- inner_join(trait_gwas,cad_tested[,c("SNP","VAR_ID","FRQ")], by=c("VAR_ID","SNP"))
    
    trait_gwas <- trait_gwas %>%
      separate(VAR_ID, into=c("CHR","BP","REF","ALT"),sep="_",remove = F)
    
    if(nrow(trait_gwas)==0){next;}
    
    trait_name <- sub(".txt.gz","",j)
    trait_name <- sub("gwas_regional/","",trait_name)
    
    final_tested <- cad_tested[cad_tested$VAR_ID %in% trait_gwas$VAR_ID,]
    
    D1 <- list(beta = final_tested$BETA,
               varbeta = (final_tested$SE)^2,
               snp = final_tested$SNP,
               position = final_tested$BP,
               type = "cc",
               N = 530185) #effective sample size of CAD GWAS
    
    D2 <- list(beta = trait_gwas$BETA,
               varbeta = (trait_gwas$SE)^2,
               snp = trait_gwas$SNP,
               position = trait_gwas$BP,
               type = "quant",
               N = as.integer(gwas_trait[gwas_trait$trait == trait_name, ]$samplesize),
               MAF = ifelse(as.numeric(trait_gwas$FRQ) > 0.5, 1- as.numeric(trait_gwas$FRQ), as.numeric(trait_gwas$FRQ)))
    
    coloc.res = coloc.abf(dataset1=D1,dataset2=D2)
    tmp2 <- coloc.res$summary
    tmp2$SNP.PP.H4 <- max(coloc.res$results$SNP.PP.H4)
    tmp2$SNP.H4 <- coloc.res$results$snp[coloc.res$results$SNP.PP.H4 == max(coloc.res$results$SNP.PP.H4)]
    tmp2 <- as.data.frame(tmp2)
    tmp2$trait <- trait_name
    idx <- match(tmp2$SNP.H4, final_tested$SNP)
    idx2 <- match(tmp2$SNP.H4, trait_gwas$SNP)
    tmp2$p1 <- final_tested[idx,]$P_VALUE
    tmp2$p2 <- trait_gwas[idx2,]$P_VALUE
    tmp2$beta1 <- final_tested[idx,]$BETA
    tmp2$beta2 <- trait_gwas[idx2,]$BETA
    tmp2$se1 <- final_tested[idx,]$SE
    tmp2$se2 <- trait_gwas[idx2,]$SE
    tmp2$region <- job_vars[i,]$No
    result_coloc <- rbind(result_coloc, tmp2)
  }
  write.table(result_coloc,paste0("1_colocalization/coloc_result/coloc_region_",job_vars[i,]$No,".txt"),quote=F,row.names=F,col.names=T)
}
