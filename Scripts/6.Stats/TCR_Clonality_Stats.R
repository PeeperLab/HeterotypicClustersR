##################################################################
# TCR expansion - statistical enrichment analysis
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds from 1.Main/Tcell_Analysis_Code.R
# # TCR_filtered_NA from 2.Plots/Tcell_Plots.R

#setwd("/YOUR/PATH/")
library(readxl)
library(car)
library(nlme)
library(lmerTest)
library(mclogit)
library(data.table)
library(ggplot2)
library(tibble)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

# reformat the data into a counts table
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
# quick save, with cluster names and frequency
set.seed(10)

md_tcr$chain_cdr3_pat = paste0(md_tcr$chain_cdr3,"_",md_tcr$patient_id)
md_freq <- md_tcr %>% 
  group_by(patient_id, sample_id, chain_cdr3_pat) %>% 
  dplyr::summarize(n=n()) %>% 
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

md_freq_bar_sample <- md_freq %>% 
  group_by(patient_id, sample_id, top_clones) %>% 
  dplyr::summarize(
    n = sum(n)  # Sum the `n` column for each group
  )

data<-md_freq_bar_sample

metadata_to_sig_enrichment=function(dataset.complete){
  
  data.newtot <- NULL
  for (k in unique(dataset.complete$patient_id)){
    dataset <- dataset.complete[dataset.complete$patient_id==k,]
    for(l in unique(dataset$sample_id)){
      data.seur <- dataset[dataset$sample_id==l,]
      data.seur$tot2 <- sum(data.seur$n)
      data.newtot <- rbind(data.newtot,data.seur)
    }
  }
  dataset.complete <- data.newtot
  # first run for tumor
  dataset1 <- dataset.complete[which(dataset.complete$sample_id %in% c('SG_Singlets','DB_Tumor_Tcell')),]
  
  results_Tumor <- NULL
  for(k in unique(dataset1$top_clones)){
    print(k)
    data_cluster <- dataset1[dataset1$top_clones==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    data_cluster$sample_id = factor(data_cluster$sample_id, levels=c('SG_Singlets','DB_Tumor_Tcell'))
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_Tumor <- rbind(results_Tumor,res$coefficients[2,])
  }
  
  results_Tumor <- data.frame(results_Tumor)
  results_Tumor$seurat_cluster <- unique(dataset1$top_clones)
  results_Tumor$comparison <- rep(paste('SG','DB_Tumor',sep='-'),nrow(results_Tumor))
  
  # now for APC
  dataset2 <- dataset.complete[which(dataset.complete$sample_id %in% c('SG_Singlets','DB_APC_Tcell')),]
  
  results_APC <- NULL
  for(k in unique(dataset2$top_clones)){
    data_cluster <- dataset2[dataset2$top_clones==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    data_cluster$sample_id = factor(data_cluster$sample_id, levels=c('SG_Singlets','DB_APC_Tcell'))
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_APC <- rbind(results_APC,res$coefficients[2,])
  }
  
  results_APC <- data.frame(results_APC)
  results_APC$seurat_cluster <- unique(dataset2$top_clones)
  results_APC$comparison <- rep(paste('SG','DB_APC',sep='-'),nrow(results_APC))
  
  # combine and get the padj
  total_results=rbind(results_Tumor,results_APC)
  total_results$p.adj <- p.adjust(total_results$Pr...z.., method ='bonferroni')
  total_results$significant <- ifelse(total_results$p.adj <0.05, 'Yes', 'No')
  return(total_results)
}

results=metadata_to_sig_enrichment(data)

write.csv(results,"Results_stats/TCR_Clonality_Stats.csv")
