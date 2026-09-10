########
#bNMF pipeline learn from Multi-ancestry polygenic mechanisms of type 2 diabetes
########

pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               rstudioapi, DT, kableExtra, GenomicRanges)
library(ggplot2)
library(ggtree)
library(aplot)
library(openxlsx)
library(patchwork)
library(ggdendro)
library(legendry)

cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")

#######prepare rsID map########
rsID_map <- cad_cleaned[,c("SNP","VAR_ID")]
colnames(rsID_map) <- c("rsID","VAR_ID")

rsID_map <- cad_cleaned[,c("SNP","VAR_ID")]
colnames(rsID_map) <- c("rsID","VAR_ID")

write.table(rsID_map, "2_bNMF/rsID_map.txt",quote = F,row.names = F,col.names = T)

rsID_map <- cad_cleaned[as.numeric(cad_cleaned$P_VALUE) < 1e-5,c("SNP","VAR_ID")]
colnames(rsID_map) <- c("rsID","VAR_ID")
write.table(rsID_map, "2_bNMF/rsID_map_1e5.txt",quote = F,row.names = F,col.names = T)

######data prepare######
source("bNMF_function/choose_variants.R")
source("bNMF_function/prep_bNMF.R")
source("bNMF_function/run_bNMF.R")
user_token = '7e42f0c43fb2'
#A1:REF,A2:ALT
cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")

data_dir <- "D:/R_project/CAD_pathway_analysis/2_bNMF"
rsID_map_file <-"D:/R_project/CAD_pathway_analysis/2_bNMF/rsID_map.txt"
rsID_map_1e5 <-"D:/R_project/CAD_pathway_analysis/2_bNMF/rsID_map_1e5.txt"
gwas_traits <- readxl::read_excel("D:/R_project/CAD_pathway_analysis/0_dataprepare/gwas_trait_raw.xlsx")
main_ss_filepath <- "D:/R_project/CAD_pathway_analysis/0_dataprepare/CAD_small.txt"
trait_ss_files <- setNames(gwas_traits$full_path, gwas_traits$trait)
trait_ss_size <- setNames(gwas_traits$samplesize, gwas_traits$trait)

##########prepare variants for clustering#######
result_all <- fread("1_colocalization/result_all.txt",header = T) #to extract variants showing colocalization
length(unique(result_all$trait2)) #125

PVCUTOFF = 5e-8
cad_var <- cad_cleaned %>%
  filter(as.numeric(P_VALUE) <= PVCUTOFF)

cad_var2 <- cad_cleaned[cad_cleaned$SNP %in% c(result_all$hit1),]

cad_var <- rbind(cad_var, cad_var2)

cad_var <- cad_var %>%
  filter(!duplicated(SNP))
print(paste("No. unique SNPs:",nrow(cad_var))) #14224

cad_var <- cad_var %>% dplyr::rename(PVALUE = P_VALUE)
gwas_variants <- cad_var$VAR_ID
var_nonmissingness <- count_traits_per_variant(gwas_variants,trait_ss_files)

df <- data.frame(
  VAR_ID = names(var_nonmissingness),
  missing_rate = as.vector(var_nonmissingness),
  stringsAsFactors = FALSE
)
cad_var <- left_join(cad_var,df)

#pruning consider for P value and missing rate
ld_prune(df_snps = cad_var[,c("VAR_ID","PVALUE","Risk_Allele","missing_rate")],
         pop = "EUR",
         output_dir = "2_bNMF/",
         r2 = 0.05,
         maf=0.01,
         my_token = user_token,
         chr = c(1:22))
print("Combining SNP.clip results...")
ld_files <- list.files(path = "2_bNMF/",
                       pattern = "^snpClip_results",
                       full.names = T)
df_clipped_res = data.frame("RS_Number"=as.character(),
                            "Position"=as.character(),
                            "Alleles"= as.character(),
                            "Details"=as.character())

rename_cols_clipped <- c(RS_Number="RS Number",
                         Position_grch37="Position")

for (ld_file in ld_files){
  df <- fread(ld_file, stringsAsFactors = F, data.table = F) %>%
    dplyr::rename(any_of(rename_cols_clipped))
  df_clipped_res <- rbind(df_clipped_res, df)
}

df_clipped_kept <- df_clipped_res %>%
  filter(Details=="Variant kept.")

#"SNP" "CHR" "BP" "REF" "ALT" "BETA" "SE" "PVALUE" "Risk_Allele" "VAR_ID" "ChrPos"
pruned_vars <- cad_var %>%
  mutate(ChrPos = paste0("chr", CHR, ":", BP)) %>%
  filter(SNP %in% df_clipped_kept$RS_Number)

#########
print("Searching for variants in trait GWAS...")
gwas_variants <- pruned_vars$VAR_ID

var_nonmissingness <- count_traits_per_variant(gwas_variants,trait_ss_files)

print("Identifying variants needing proxies...")
proxies_needed_df <- find_variants_needing_proxies(pruned_vars,
                                                   var_nonmissingness,
                                                   rsID_map_file,
                                                   missing_cutoff = 0.8)
print("Searching for proxies with TopLD API...")
proxy_search_results <- choose_proxies(need_proxies = proxies_needed_df,
                                       method="LDlink",
                                       LDlink_token = user_token,
                                       topLD_path = NA,
                                       rsID_map_file = rsID_map_1e5,
                                       trait_ss_files = trait_ss_files,
                                       pruned_variants = pruned_vars,
                                       population="EUR",
                                       r2_num = 0.8)

df_proxies <- proxy_search_results %>%
  dplyr::select(VAR_ID, proxy_VAR_ID) %>%
  dplyr::inner_join(pruned_vars[,c("VAR_ID")], by="VAR_ID") %>%
  mutate(Risk_Allele=NA, PVALUE=NA) 

#####Fetch summary statistics for SNPs in trait GWAS
print("Prepping input for fetch_summary_stats...")
df_orig_snps <- pruned_vars %>%
  filter(!VAR_ID %in% proxies_needed_df$VAR_ID)

df_input_snps <- df_orig_snps[,c("VAR_ID","PVALUE", "Risk_Allele")] %>%
  arrange(PVALUE) %>%
  filter(!duplicated(VAR_ID))

cat(sprintf("\n%i original SNPs...\n", nrow(df_orig_snps)))
cat(sprintf("\n%i proxy SNPs...\n", nrow(df_proxies)))
cat(sprintf("\n%i total unique SNPs!\n", nrow(df_input_snps)))

initial_zscore_matrices <- fetch_summary_stats(
  df_input_snps,
  cad_cleaned,
  trait_ss_files,
  pval_cutoff=0.05
)
#saveRDS(initial_zscore_matrices,"2_bNMF/initial_zscore_matrices.rds")
initial_zscore_matrices <- readRDS("2_bNMF/initial_zscore_matrices.rds")
system(sprintf("mv alignment_GWAS_summStats.csv %s", data_dir))

for( i in colnames(initial_zscore_matrices[["N_mat"]])){
  initial_zscore_matrices[["N_mat"]][,i] <- ifelse(is.na(initial_zscore_matrices[["N_mat"]][,i]), trait_ss_size[i],initial_zscore_matrices[["N_mat"]][,i])
}

print("Getting rsIDs for final snps and saving to results...")
z_mat <- initial_zscore_matrices$z_mat
N_mat <- initial_zscore_matrices$N_mat

df_var_ids <- df_input_snps %>%
  separate(VAR_ID, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
  mutate(ChrPos=paste(Chr,Pos,sep = ":")) %>%
  subset(ChrPos %in% rownames(z_mat))

df_rsIDs <- fread(rsID_map_file,header = F,col.names = c("rsID","VAR_ID")) %>%
  filter(VAR_ID %in% df_var_ids$VAR_ID)

print(sprintf("rsIDs found for %i of %i SNPs...", nrow(df_rsIDs), nrow(df_var_ids)))

write_delim(x=df_rsIDs,
            file = file.path(data_dir, "rsID_map_small.txt"),
            col_names = T)

#########Fill missing data in z-score and N matrices####
df_snps <- df_input_snps %>%
  inner_join(df_rsIDs, by="VAR_ID") %>%
  data.frame()

print("Searching for cover proxies for missing z-scores...")
initial_zscore_matrices_final <- fill_missing_zscores(initial_zscore_matrices,
                                                      df_snps,
                                                      trait_ss_files,
                                                      trait_ss_size,
                                                      main_ss_filepath,
                                                      rsID_map_file,
                                                      method_fill="proxy",
                                                      population="EUR",
                                                      LDlink_token = user_token)

prep_z_output <- prep_z_matrix(z_mat = initial_zscore_matrices_final$z_mat,
                               N_mat = initial_zscore_matrices_final$N_mat,
                               corr_cutoff = 0.85)
final_zscore_matrix <- prep_z_output$final_z_mat
#saveRDS(final_zscore_matrix, file="2_bNMF/final_zscore_matrix.rds")
final_zscore_matrix <- readRDS("2_bNMF/final_zscore_matrix.rds")

df_traits_filtered <- prep_z_output$df_traits
write_csv(x = df_traits_filtered,
          file = file.path(data_dir,"df_traits.csv"))

print(sprintf("Final matrix: %i SNPs x %i traits",
              nrow(final_zscore_matrix),
              ncol(final_zscore_matrix)/2))

###bNMF iteration#######
bnmf_reps <- run_bNMF(final_zscore_matrix,
                      n_reps=25)
summarize_bNMF(bnmf_reps, dir_save=data_dir)
k <- NULL
if (is.null(k)){
  html_filename <- "results_for_maxK.html"
} else {
  html_filename <- sprintf("results_for_K_%i.html", k)
}

rmarkdown::render(
  'bNMF_function/format_bNMF_results.Rmd',
  output_file = html_filename,
  params = list(main_dir = "2_bNMF/",
                k = k,
                loci_file="query",
                GTEx=F,
                my_traits=gwas_traits)
)

##########validation#########
cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")
gwas_traits <- readxl::read_excel("0_dataprepare/gwas_trait_raw.xlsx")
snp_cluster <- read.table("2_bNMF/snp_cluster.txt", header = T) #extract from bNMF result
snp_cluster$used <- rowSums(snp_cluster[, c("idx1", "idx2", "idx3", "idx4", "idx5","idx6","idx7","idx8","idx9")]) #index indicating whether snps were in the clusters

cad_snp_used <- cad_cleaned[cad_cleaned$SNP %in% snp_cluster[snp_cluster$used > 0,]$SNP,]
flip <- ifelse(cad_snp_used$BETA > 0,1,-1)
flip_df <- data.frame(SNP = cad_snp_used$SNP, flip = flip, check.names = FALSE)

assoc_result <- data.frame("trait"=character(),"BETA1"=numeric(),"P1"=numeric(),"Z1"=numeric(),"FRAC1"=numeric(),
                           "BETA2"=numeric(),"P2"=numeric(),"Z2"=numeric(),"FRAC2"=numeric(),
                           "BETA3"=numeric(),"P3"=numeric(),"Z3"=numeric(),"FRAC3"=numeric(),
                           "BETA4"=numeric(),"P4"=numeric(),"Z4"=numeric(),"FRAC4"=numeric(),
                           "BETA5"=numeric(),"P5"=numeric(),"Z5"=numeric(),"FRAC5"=numeric(),
                           "BETA6"=numeric(),"P6"=numeric(),"Z6"=numeric(),"FRAC6"=numeric(),
                           "BETA7"=numeric(),"P7"=numeric(),"Z7"=numeric(),"FRAC7"=numeric(),
                           "BETA8"=numeric(),"P8"=numeric(),"Z8"=numeric(),"FRAC8"=numeric(),
                           "BETA9"=numeric(),"P9"=numeric(),"Z9"=numeric(),"FRAC9"=numeric())

for(i in 1:nrow(gwas_traits)){
  trait_name <- gwas_traits[i,]$trait
  
  trait_gwas <- fread(gwas_traits[i,]$full_path)
  trait_gwas <- trait_gwas[trait_gwas$SNP %in% snp_cluster[snp_cluster$used > 0,]$SNP,]
  trait_gwas <- left_join(trait_gwas, snp_cluster)
  trait_gwas <- left_join(trait_gwas, flip_df, by = "SNP")
  trait_gwas$BETA_flip <- trait_gwas$BETA * trait_gwas$flip 
  trait_gwas <- trait_gwas  %>% mutate(weights = 1 / (SE^2))
  trait_gwas <- trait_gwas  %>% mutate(vi = (SE^2))
  
  temp_result <- c(trait_name)
  for(j in 1:9){
    cluster_col <- paste0("idx", j)
    model_data <- trait_gwas[trait_gwas[[cluster_col]] == 1, ]
    beta = sum(model_data$BETA_flip * model_data$weights) / sum(model_data$weights)
    se = 1 / sqrt(sum(model_data$weights))
    Z=beta/se
    p = 2 * pnorm(-abs(Z))
    n_pos <- sum(model_data[as.numeric(model_data$P_VALUE) < 0.05,]$BETA_flip > 0, na.rm = TRUE)
    n_neg <- sum(model_data[as.numeric(model_data$P_VALUE) < 0.05,]$BETA_flip < 0, na.rm = TRUE)
    if(beta > 0){
      fraction <- n_pos/length(model_data$BETA_flip)
    }else{
      fraction <- n_neg/length(model_data$BETA_flip)
    }
    temp_result <- c(temp_result, beta, p, Z,fraction)
  }
  temp_result <- as.data.frame(t(temp_result))
  colnames(temp_result) <- c("trait","BETA1","P1","Z1","FRAC1",
                             "BETA2","P2","Z2","FRAC2",
                             "BETA3","P3","Z3","FRAC3",
                             "BETA4","P4","Z4","FRAC4",
                             "BETA5","P5","Z5","FRAC5",
                             "BETA6","P6","Z6","FRAC6",
                             "BETA7","P7","Z7","FRAC7",
                             "BETA8","P8","Z8","FRAC8",
                             "BETA9","P9","Z9","FRAC9")
  assoc_result <- rbind(assoc_result, temp_result)
}

writexl::write_xlsx(assoc_result, path ="2_bNMF/result_summary/genetic_assoc.xlsx")

#######LOO##########
cad_cleaned <- readRDS("0_dataprepare/CAD_small.rds")
gwas_traits <- readxl::read_excel("0_dataprepare/gwas_trait_raw.xlsx")
snp_cluster <- read.table("2_bNMF/snp_cluster.txt", header = T)
snp_cluster$used <- rowSums(snp_cluster[, c("idx1", "idx2", "idx3", "idx4", "idx5","idx6","idx7","idx8","idx9")])

cad_snp_used <- cad_cleaned[cad_cleaned$SNP %in% snp_cluster[snp_cluster$used > 0,]$SNP,]
flip <- ifelse(cad_snp_used$BETA > 0,1,-1)
flip_df <- data.frame(SNP = cad_snp_used$SNP, flip = flip, check.names = FALSE)

assoc_result_sensi <- data.frame("trait"=character(),"BETA"=numeric(),"P"=numeric(),"Z"=numeric(),"cluster"=numeric(),"SNP"=character())

for(i in 1:nrow(gwas_traits)){
  trait_name <- gwas_traits[i,]$trait
  
  trait_gwas <- fread(gwas_traits[i,]$full_path)
  trait_gwas <- trait_gwas[trait_gwas$SNP %in% snp_cluster[snp_cluster$used > 0,]$SNP,]
  trait_gwas <- left_join(trait_gwas, snp_cluster)
  trait_gwas <- left_join(trait_gwas, flip_df, by = "SNP")
  trait_gwas$BETA_flip <- trait_gwas$BETA * trait_gwas$flip 
  trait_gwas <- trait_gwas  %>% mutate(weights = 1 / (SE^2))
  
  for(j in 1:9){
    cluster_col <- paste0("idx", j)
    model_data <- trait_gwas[trait_gwas[[cluster_col]] == 1, ]
    for(k in c(model_data$SNP)){
      temp_result <- c(trait_name)
      model_data_sensi <- model_data[model_data$SNP!=k,]
      beta = sum(model_data_sensi$BETA_flip * model_data_sensi$weights) / sum(model_data_sensi$weights)
      se = 1 / sqrt(sum(model_data_sensi$weights))
      Z=beta/se
      p = 2 * pnorm(-abs(Z))
      temp_result <- as.data.frame(t(c(temp_result, beta, p, Z,j,k)))
      colnames(temp_result) <- c("trait","BETA","P","z","cluster","SNP")
      assoc_result_sensi <- rbind(assoc_result_sensi,temp_result)
    }
  }
}

writexl::write_xlsx(assoc_result_sensi, path ="2_bNMF/result_summary/genetic_assoc_sensi.xlsx")
