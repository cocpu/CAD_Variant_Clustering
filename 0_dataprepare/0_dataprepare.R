library(MungeSumstats)
library(data.table)
library(tidyverse)

cad <- fread("raw_data/GCST90132314.h.tsv.gz")
cad <- dplyr::select(cad, all_of(names(cad)[colSums(is.na(cad)) < nrow(cad)]))

cad <- cad[,c("rsid","chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","effect_allele_frequency",
              "p_value","cases","n")]

chain_file <- "raw_data/hg38ToHg19.over.chain.gz"
file.copy(from = chain_file, to = file.path(tempdir(), basename(chain_file)), overwrite = TRUE)
cad_cleaned <- MungeSumstats::format_sumstats(path = cad,convert_ref_genome="GRCh37", dbSNP=144,
                                              strand_ambig_filter = T,
                                              return_data = T, return_format = "data.table")

cad_cleaned <- cad_cleaned[,c("SNP","CHR","BP","A1","A2","P","FRQ","BETA","SE","CASES","N")]
cad_cleaned <- filter(cad_cleaned, nchar(A1) == 1 & nchar(A2) == 1)
saveRDS(cad_cleaned, "0_dataprepare/CAD.rds")

cad_noUKB <- fread("raw_data/26343387-GCST003116-EFO_0000378.h.tsv.gz")
cad_noUKB <- cad_noUKB[,c("hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value")]
colnames(cad_noUKB) <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","p_value")
chain_file <- "raw_data/Genome/hg38ToHg19.over.chain.gz"
file.copy(from = chain_file, to = file.path(tempdir(), basename(chain_file)), overwrite = TRUE)
cad_noUKB_cleaned <- MungeSumstats::format_sumstats(path = cad_noUKB,convert_ref_genome="GRCh37", dbSNP=144,
                                                    strand_ambig_filter = T,
                                                    return_data = T, return_format = "data.table")
cad_noUKB_cleaned  <- filter(cad_noUKB_cleaned , nchar(A1) == 1 & nchar(A2) == 1)
saveRDS(cad_noUKB_cleaned, "0_dataprepare/CAD_noUKB.rds")
