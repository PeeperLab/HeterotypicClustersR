##################################################################
# TCR expansion - statistical enrichment analysis per patient
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds from 1.Main/Tcell_Analysis_Code.R
# # TCR_filtered_NA from 2.Plots/Tcell_Plots.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(tidyverse)
library(dplyr)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

S_pat_all_t_cells = readRDS("Results_Tcells/Tcells_Final.Rds")
md = S_pat_all_t_cells@meta.data
md$cell_id = rownames(md)

# add TCR info
md$chain_cdr3 = NA
md$CTnt = NA
md$CTstrict = NA

# fix barcodes and match cells to combine with TCR info
sample_names_table = read.table("./TCR_filtered_NA/TCR_samplenames_barcodes_filNA.txt",header=T)
for (sample_x in sample_names_table$V2){
  sample_tcr_table = read.table(paste0("./TCR_filtered_NA/TCR_combined_",sample_x,"_filNA.txt"),header=T)
  sample_tcr_table[,1] = gsub(paste0(sample_x,"_"),"",sample_tcr_table[,1])
  md[md$cell_id %in% sample_tcr_table[,1],]$chain_cdr3 = sample_tcr_table[,11]
  md[md$cell_id %in% sample_tcr_table[,1],]$CTnt = sample_tcr_table[,10]
  md[md$cell_id %in% sample_tcr_table[,1],]$CTstrict = sample_tcr_table[,12]
}

# only cells with a TCRs
md_tcr <- md[!is.na(md$chain_cdr3),]

set.seed(10)

md_tcr$chain_cdr3_pat = paste0(md_tcr$chain_cdr3,"_",md_tcr$patient_id)
md_freq <- md_tcr %>% 
  group_by(patient_id, sample_id, chain_cdr3_pat) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n), total=sum(n))
md_freq <- md_freq[order(md_freq$freq,decreasing = T),]

md_freq$top_clones <- "none"
for (y in c("P1","P2","P3","P4","P5")){
  for (x in unique(md_freq$sample_id)){
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][1,]$top_clones <- "top_1"
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][2:5,]$top_clones <- "top_2_5"
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][6:15,]$top_clones <- "top_6_15"
  }
}

##TCR Stats
md_all<-merge(md_tcr,md_freq,by=c("patient_id","sample_id","chain_cdr3_pat"),all.x = TRUE)
md_all$categories<-md_all$top_clones

fisher_test<-function(table_data){
  patients <- unique(table_data$patient_id)
  states <- unique(table_data$categories)
  interaction_groups <- c("DB_Tumor_Tcell", "DB_APC_Tcell")
  singlet_group <- "SG_Singlets"
  
   results <- data.frame(
    patient = character(),
    group = character(),
    comparison = character(),
    p_value = numeric()
  )
  
  for (patient in patients) {
    print(paste("Processing patient:", patient))
    patient_data <- table_data[table_data$patient_id == patient, ]
    
    for (state in states) {
      print(paste("  Processing state:", state))
      state_data <- patient_data[patient_data$categories == state, ]
      print(paste("  Rows in state_data for state", state, ":", nrow(state_data)))
      
      if (nrow(state_data) == 0) {
        next
      }
      
      for (group in interaction_groups) {
        subset_data <- patient_data[patient_data$sample_id %in% c(group, singlet_group), ]
        
        if (nrow(subset_data) == 0) {
          next
        }
        
        subset_data$Interaction <- ifelse(subset_data$sample_id == group, "Interaction", "Singlet")
        subset_data$State <- ifelse(subset_data$categories == state, "InState", "NotInState")
        
        if (length(unique(subset_data$Interaction)) < 2 || length(unique(subset_data$State)) < 2) {
          next
        }
        
        contingency <- table(subset_data$Interaction, subset_data$State)
        fisher_result <- fisher.test(contingency)
        
        results <- rbind(
          results,
          data.frame(
            patient = patient,
            group = state,
            comparison = paste("SG_Singlets vs", group),
            p_value = fisher_result$p.value
          )
        )
      }
    }
  }
return(results)
}

tcr_results<-fisher_test(md_all)
if (nrow(tcr_results) > 0) {
  tcr_results$adjusted_p_value <- p.adjust(tcr_results$p_value, method = "fdr")
}
write.csv(tcr_results,"Results_stats/Top_TCR_stats_per_patient.txt")


