##########
#example for protein mediation
#learn from Integrative analysis of the plasma proteome and polygenic risk of cardiometabolic diseases
#######
library(data.table)
library(tidyverse)
library(medflex)

protein <- fread("olink_data.txt.gz")
protein_id <- fread("protein_id.tsv")
protein_id$protein_name <- ""
for(i in 1:nrow(protein_id)){
  protein_id[i,]$protein_name <- strsplit(as.character(protein_id[protein_id$coding ==i,]$meaning), ";")[[1]][1]
}
results_all_mediation <- data.frame()
for(i in c("1","2","4","8","9","all")){
  protein_pass <- get(paste0("protein_pass_step3_",i)) #the protein used in mediation
  protein_id_pass <- protein_id[protein_id$protein_name %in% protein_pass,]$coding
  for(j in protein_id_pass){
    protein_target <- protein[protein$protein_id== j,]
    protein_target <- left_join(cad_survival_protein, protein_target[,c("eid","result")],by="eid")
    protein_name <- protein_id[protein_id$coding == j,]$protein_name
    vars_needed <- c(
      "eid", "group", "result", "age", "sex", "genotype_batch","sample_age", 
      "sample_season", "fasting_time", 
      paste0("PC", 1:10)
    )
    protein_target <- protein_target[complete.cases(protein_target[, vars_needed]), ]
    formula <- paste0("group~PRS_",i,"+result+age+sex+genotype_batch+sample_age+sample_season+fasting_time+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
    formula <- as.formula(formula)
    expData <- neImpute(formula, family = binomial("logit"), data=protein_target)
    
    formula <- paste0("group~PRS_",i,"0+PRS_",i,"1+age+sex+genotype_batch+sample_age+sample_season+fasting_time+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
    formula <- as.formula(formula)
    nMod <- neModel(formula, family = binomial("logit"), expData=expData, se="robust")
    lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c(paste0("PRS_",i,"0=0"), paste0("PRS_",i,"1=0"), paste0("PRS_",i,"0+PRS_",i,"1=0"))))))
    
    results <- data.frame(
      protein = protein_name,
      PRS = i,

      EST_DIRECT = lht$Estimate[1],
      SE_DIRECT = lht$`Std. Error`[1],
      P_DIRECT = lht$`Pr(>|z|)`[1],

      EST_INDIRECT = lht$Estimate[2],
      SE_INDIRECT = lht$`Std. Error`[2],
      P_INDIRECT = lht$`Pr(>|z|)`[2],

      EST_TOTAL = lht$Estimate[3],
      SE_TOTAL = lht$`Std. Error`[3],
      P_TOTAL = lht$`Pr(>|z|)`[3]
    )
    results_all_mediation <- rbind(results_all_mediation, results)
  }
}
saveRDS(results_all_mediation, "4_Association/association_protein/results_all_mediation.rds")
