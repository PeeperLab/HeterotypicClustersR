##################################################################
# T cells - statistical enrichment analysis
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds from 1.Main/Tcell_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(car)
library(nlme)
library(lmerTest)
library(mclogit)
library(data.table)
library(ggplot2)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

# reformat the data into a counts table
sc.combined<-readRDS("Results_Tcells/Tcells_Final.Rds")
metadata<-sc.combined@meta.data

#########
setDT(metadata)
dataset.complete=metadata[, .N, by=c("patient_id", "sample_id", "annotated_clusters_labels_final")]
setnames(dataset.complete, "N","n")

# load the function
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
  for(k in unique(dataset2$annotated_clusters_labels_final)){
    data_cluster <- dataset2[dataset2$annotated_clusters_labels_final==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_APC <- rbind(results_APC,res$coefficients[2,])
  }
  
  results_APC <- data.frame(results_APC)
  results_APC$seurat_cluster <- unique(dataset2$annotated_clusters_labels_final)
  results_APC$comparison <- rep(paste('SG','DB_APC',sep='-'),nrow(results_APC))
  
  # now with the classic Tum doublets
  dataset3 <- dataset.complete[which(dataset.complete$sample_id %in% c('SG_Singlets','DB_Tumor_Tcell')),]
  
  results_tum <- NULL
  for(k in unique(dataset3$annotated_clusters_labels_final)){
    data_cluster <- dataset3[dataset3$annotated_clusters_labels_final==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_tum <- rbind(results_tum,res$coefficients[2,])
  }
  
results_tum <- data.frame(results_tum)
results_tum$seurat_cluster <- unique(dataset3$annotated_clusters_labels_final)
results_tum$comparison <- rep(paste('SG','DB_Tum',sep='-'),nrow(results_APC))
  
results_APC=rbind(results_APC,results_tum)
  
# combine and get the padj
results_APC$p.adj <- p.adjust(results_APC$Pr...z.., method ='bonferroni')
results_APC$significant <- ifelse(results_APC$p.adj <0.05, 'Yes', 'No')

results_APC$stars <- ifelse(results_APC$p.adj < 0.00001, "*****",
                            ifelse(results_APC$p.adj < 0.0001, "****",
                                   ifelse(results_APC$p.adj < 0.001, "***",
                                          ifelse(results_APC$p.adj < 0.01, "**",
                                                 ifelse(results_APC$p.adj < 0.05, "*",
                                                        ifelse(results_APC$p.adj < 0.1, ".", ""))))))
fwrite(results_APC,"Results_stats/Tcells_Phenotype_clusters_statistics.csv")
  
  
   
