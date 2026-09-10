library(dplyr)
library(optparse)
library(data.table)
library(MungeSumstats)

#VAR_ID	Effect_Allele_PH	P_VALUE	BETA	SE	N_PH FRQ

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default=FALSE,help="target file"),
  make_option(c("--trait"), type = "character", default=FALSE,help="trait"),
  make_option(c("--dic"), type = "character", default=FALSE,help="target dic"),
  make_option(c("--rsid"),type = "character",default=NULL, help="rs id"),
  make_option(c("--chr"),type = "character",default=NULL, help="chr"),
  make_option(c("--pos"),type = "character",default=NULL, help="position"),
  make_option(c("--ref"),type = "character",default=NULL, help="ref"),
  make_option(c("--alt"),type = "character",default=NULL, help="alt"),
  make_option(c("--beta"),type = "character",default=NULL, help="beta"),
  make_option(c("--p"),type = "character",default=NULL, help="p value"),
  make_option(c("--se"),type = "character",default=NULL, help="se"),
  make_option(c("--frq"),type = "character",default=NULL, help="frequency"),
  make_option(c("--n"),type = "character",default=NULL, help="n"),
  make_option(c("--grch"),type= "character", default=NULL, help="grch")
)

opt_parser = OptionParser(
  usage = "usage: %prog [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL , 
  description = "This Script is to clean the gwas result."
)

opt = parse_args(opt_parser)


gwas <- fread(opt$file)
trait_name <- opt$trait

gwas <- gwas %>% select(all_of(c(opt$rsid, opt$chr, opt$pos, opt$ref, opt$alt, opt$beta, opt$p, opt$se, opt$frq, opt$n)))
col_name <- c("SNP" = opt$rsid, "CHR" = opt$chr,"BP" = opt$pos,"other_allele" = opt$ref, "effect_allele" = opt$alt, "BETA" = opt$beta,
              "P" = opt$p, "SE" = opt$se, "FRQ" = opt$frq, "N" = opt$n)
colnames(gwas) <- names(col_name)

chain_file <- "/mnt/data/chenmingyang/resource/chain_file/hg38ToHg19.over.chain.gz"
file.copy(from = chain_file, to = file.path(tempdir(), basename(chain_file)), overwrite = TRUE)
if(is.null(opt$grch)){
  gwas_cleaned <- MungeSumstats::format_sumstats(path = gwas,convert_ref_genome="GRCh37", dbSNP=144,
                                                 strand_ambig_filter = T,
                                                 return_data = T, return_format = "data.table")
}else if(opt$grch == "GRCH37"){
  gwas_cleaned <- MungeSumstats::format_sumstats(path = gwas,convert_ref_genome="GRCh37", dbSNP=144,ref_genome="GRCh37",
                                                 strand_ambig_filter = T,
                                                 return_data = T, return_format = "data.table")
}else if(opt$grch == "GRCH38"){
  gwas_cleaned <- MungeSumstats::format_sumstats(path = gwas,convert_ref_genome="GRCh37", dbSNP=144,ref_genome="GRCh38",
                                                 strand_ambig_filter = T,
                                                 return_data = T, return_format = "data.table")
}

target_gwas <- gwas_cleaned %>% select(all_of(c("SNP","CHR","BP","A1","A2","P","BETA","SE")))
if("N" %in% colnames(gwas_cleaned)){
  target_gwas$N_PH <- gwas_cleaned$N
}else{
  target_gwas$N_PH <- NA
}

if("FRQ" %in% colnames(gwas_cleaned)){
  target_gwas$FRQ <- gwas_cleaned$FRQ
}else{
  target_gwas$FRQ <- NA
}

target_gwas <- filter(target_gwas, nchar(A1) == 1 & nchar(A2) == 1)
target_gwas$VAR_ID <- paste0(target_gwas$CHR,"_",target_gwas$BP,"_",target_gwas$A1,"_",target_gwas$A2)
target_gwas$Effect_Allele_PH <- ifelse(target_gwas$BETA> 0, target_gwas$A2, target_gwas$A1 )
target_gwas <- select(target_gwas,-c("CHR","BP","A1","A2"))
target_gwas <- dplyr::rename(target_gwas,P_VALUE = P)

write.table(target_gwas,paste0(opt$dic,trait_name,".txt"),quote=F,row.names = F)
R.utils::gzip(paste0(opt$dic,trait_name,".txt"))