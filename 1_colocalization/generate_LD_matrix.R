library(data.table)
library(dplyr)
library(tidyr)
library(optparse)
library(ieugwasr)

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
#the region need LD matrix
region_coloc <- read.table("0_dataprepare/region_coloc.txt", header = T)
#a smaller version for fast reading
cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")

#The check of SNP allele should be done previously
valid_snp <- read.table("1_colocalization/valid_snp_coloc.txt")

i<-opt$job
chr <- region_coloc[i,]$chr
start <- region_coloc[i,]$start
end <- region_coloc[i,]$end
cad_tested <- cad_cleaned[as.numeric(cad_cleaned$CHR) == as.numeric(chr) & as.numeric(cad_cleaned$BP) >= as.numeric(start) & 
                            as.numeric(cad_cleaned$BP) <= as.numeric(end), ]
cad_tested <- cad_tested[cad_tested$SNP %in% valid_snp$V1,]
if(nrow(cad_tested) < 100){
  print(paste0(i," need extra"))
}
ld <- ld_matrix(cad_tested$SNP,
                plink_bin = "/mnt/data/chenmingyang/tools/plink1.9/plink",
                bfile = "/mnt/data/chenmingyang/resource/1KG/hg19_1KGP3_plink/EUR",
                with_alleles=F)
write.table(ld,paste0("1_colocalization/coloc_LD_matrix/region_coloc_",i,".txt"), quote = F,row.names = F,col.names = T)