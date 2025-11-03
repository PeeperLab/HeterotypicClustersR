##################################################################
# APCs - statistical enrichment analysis
##################################################################

# required input data in the working directory:
# # apc_data_filtered_fvf_corr.rds from 1.Main/APC_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(readxl)
library(car)
library(nlme)
library(lmerTest)
library(mclogit)
library(data.table)
library(ggplot2)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

# Load the objects of interest
all_apcs=readRDS("Results_APCs/apc_data_filtered_fvf_corr.rds")

metadata=setDT(data.frame(all_apcs@meta.data))
metadata=metadata[sample_id!="DB_Tumor_Tcell",] # remove 
data=metadata[, .N, by=c("patient_id", "sample_id", "hi_res_clus", "major_classes")]
setnames(data, "N","n")

# load the function
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
  
  # now with the classic APC doublets
  dataset2 <- dataset.complete[which(dataset.complete$sample_id %in% c('SG_Singlets','DB_APC_Tcell')),]
  
  results_APC <- NULL
  for(k in unique(dataset2$hi_res_clus)){
    data_cluster <- dataset2[dataset2$hi_res_clus==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_APC <- rbind(results_APC,res$coefficients[2,])
  }
  
  results_APC <- data.frame(results_APC)
  results_APC$seurat_cluster <- unique(dataset2$hi_res_clus)
  results_APC$comparison <- rep(paste('SG','DB_APC',sep='-'),nrow(results_APC))
  
  # combine and get the padj
  results_APC$p.adj <- p.adjust(results_APC$Pr...z.., method ='bonferroni')
  results_APC$significant <- ifelse(results_APC$p.adj <0.05, 'Yes', 'No')
  return(results_APC)
}

mono_macro_results=metadata_to_sig_enrichment(data[major_classes=="mono_mac",])

# export
fwrite(mono_macro_results, "Results_stats/pvals_mono_mac_res.txt")

# # # AFTER USING ONLY PATs OF INTEREST

dc_results_filt=metadata_to_sig_enrichment(data[data$major_classes=="DCs" & patient_id %in% c("P2","P3","P5")])
bcell_results_filt=metadata_to_sig_enrichment(data[data$major_classes=="B_plasma_cells" & patient_id %in% c("P3","P5")])

# export
fwrite(dc_results_filt, "Results_stats/pvals_dc_res.txt")
fwrite(bcell_results_filt, "Results_stats/pvals_bcell_res.txt")
