##################################################################
# Tumor cells - statistical enrichment analysis
##################################################################

# required input data in the working directory:
# # Tumor_Annotated.Rds from 1.Main/Tumor_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(readxl)
library(car)
library(nlme)
library(lmerTest)
library(mclogit)
library(data.table)
library(ggplot2)
library(dplyr)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

# reformat the data into a counts table
sc.combined<-readRDS("Tumor_Annotated.Rds")

sc.combined_subset<-subset(x = sc.combined, subset = annotated_clusters %in% c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific","Stress(hypoxia response)","Stress(p53 response)",
                                                                               "Neural Crest Like"))

metadata<-sc.combined_subset@meta.data
md_freq_bar_sample = metadata %>% 
  group_by(patient_id,sample_id,annotated_clusters) %>% 
  dplyr::summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
data<-md_freq_bar_sample
x<-data[which(data$sample_id %in% c('SG_Singlets','DB_Tumor_Tcell')),]
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
  dataset1 <- dataset.complete[which(dataset.complete$sample_id %in% c('SG_Singlets','DB_Tumor_Tcell')),]
  
  results_Tumor <- NULL
  for(k in unique(dataset1$annotated_clusters)){
    print(k)
    data_cluster <- dataset1[dataset1$annotated_clusters==k,]
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    model <- glmer(response ~ 1 + sample_id  + (1|patient_id) ,data= data_cluster,family='binomial')
    res <- summary(model)
    results_Tumor <- rbind(results_Tumor,res$coefficients[2,])
    
  }
  colnames(data)
  results_Tumor <- data.frame(results_Tumor)
  results_Tumor$seurat_cluster <- unique(dataset1$annotated_clusters)
  results_Tumor$comparison <- rep(paste('SG','DB_multi',sep='-'),nrow(results_Tumor))
  # combine and get the padj
  results_Tumor$p.adj <- p.adjust(results_Tumor$Pr...z.., method ='bonferroni')
  results_Tumor$significant <- ifelse(results_Tumor$p.adj <0.05, 'Yes', 'No')
  return(results_Tumor)
}

tumor_results=metadata_to_sig_enrichment(data)


write.csv(tumor_results,"Results_stats/Tumor_Phenotype_clusters_statistics.csv")
