library(ggplot2)
library(ggpubr)
library(data.table)
library(ggbeeswarm)
#setwd("/DATA/j.simon/JSN/revisions_round_3/")
set.seed(150799)

# Load AUC scores and reformat the metadata
# # created in script s1_phenotypes_to_barras.R
auc_scored=fread("Results_barras_CD39/new_auc_scored_barras.txt")
pre_auc_scored=auc_scored[timepoint=="T0",]

pre_auc_scored[, response:= ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1))) %in% c(2,3,7,8,9,13), "R","NR")]
pre_auc_scored[, response2:= ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1))) %in% c(2,3), "CR",
                                ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1))) %in% c(8,9,13,7), "PR",
                                       ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1))) %in% c(1,4,10), "SD",
                                              ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1))) %in% c(5,6,11,12), "PD","??"))))]

pre_auc_scored[,pat:=stringr::str_split_i(patient_id, "_",1)]
pre_auc_scored[, tot_pat:= .N, by="pat"]

# # DB_100
pre_auc_scored[, db_100_split:= ifelse(db_100< quantile(db_100, .33), "low",
                                           ifelse(db_100>=quantile(db_100, .67), "high","medium"))]
pre_auc_scored[, db_100_split:= factor(db_100_split, levels=c("high","medium","low"))]
pre_auc_scored[, freq_pat_db_100:= .N/tot_pat, by=c("pat", "db_100_split")]

# # DB_30
pre_auc_scored[, db_30_split:= ifelse(db_30< quantile(db_30, .33), "low",
                                       ifelse(db_30>=quantile(db_30, .67), "high","medium"))]
pre_auc_scored[, db_30_split:= factor(db_30_split, levels=c("high","medium","low"))]
pre_auc_scored[, freq_pat_db_30:= .N/tot_pat, by=c("pat", "db_30_split")]

# # Oliveira (Oliveira_tumor_specific)
pre_auc_scored[, oliveira_split:= ifelse(Oliveira_tumor_specific< quantile(Oliveira_tumor_specific, .33), "low",
                                      ifelse(Oliveira_tumor_specific>=quantile(Oliveira_tumor_specific, .67), "high","medium"))]
pre_auc_scored[, oliveira_split:= factor(oliveira_split, levels=c("high","medium","low"))]
pre_auc_scored[, freq_pat_oliveira:= .N/tot_pat, by=c("pat", "oliveira_split")]

# # Meng (Offringa_TR_9_samples)
pre_auc_scored[, meng_split:= ifelse(Offringa_TR_9_samples< quantile(Offringa_TR_9_samples, .33), "low",
                                      ifelse(Offringa_TR_9_samples>=quantile(Offringa_TR_9_samples, .67), "high","medium"))]
pre_auc_scored[, meng_split:= factor(meng_split, levels=c("high","medium","low"))]
pre_auc_scored[, freq_pat_meng:= .N/tot_pat, by=c("pat", "meng_split")]

# # Lowery (NeoTCR)
pre_auc_scored[, neotcr_split:= ifelse(NeoTCR_cd8< quantile(NeoTCR_cd8, .33), "low",
                                      ifelse(NeoTCR_cd8>=quantile(NeoTCR_cd8, .67), "high","medium"))]
pre_auc_scored[, neotcr_split:= factor(neotcr_split, levels=c("high","medium","low"))]
pre_auc_scored[, freq_pat_neotcr:= .N/tot_pat, by=c("pat", "neotcr_split")]

# # plot the results --> db_30 to main
signatures_of_interest=c("db_30","db_100","meng","oliveira", "neotcr")

plot_sig_pat_freqs=function(sig_name){
  SIGNATURE_SPLIT=paste0(sig_name,"_split")
  SIGNATURE_PAT_FREQUENCY=paste0("freq_pat_",sig_name)
  mavg_SIGNATURE_pat=pre_auc_scored[, lapply(.SD, mean),
                                   by=c("pat","response", "response2", SIGNATURE_SPLIT),
                                   .SDcols=SIGNATURE_PAT_FREQUENCY]
  
  ggplot(mavg_SIGNATURE_pat, aes(response, .data[[SIGNATURE_PAT_FREQUENCY]]))+
    geom_boxplot(aes(fill = response), outliers = F)+
    geom_beeswarm(aes(colour = response2), size=2.5, alpha=1,cex=6)+
    stat_compare_means(method="t", size=2.75, label.y.npc = 1, method.args = list(var.equal = FALSE))+
    scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                                "PR"="#3b756e","CR"="#6d8754"))+
    scale_fill_manual(values=c("NR"="grey60", "R"="#c1de87"))+
    theme_bw()+ggtitle(paste0("Per pat -",sig_name))+theme(panel.grid=element_blank())+
    facet_wrap(as.formula(paste0("~", SIGNATURE_SPLIT)))

  ggsave(paste0("./output_figs/rev3_M3E_",sig_name, ".pdf"), width=6,height=4)
  write.csv(mavg_SIGNATURE_pat, paste0("./output_figs/rev3_M3E_",sig_name, ".csv"))
  return(mavg_SIGNATURE_pat)
}

# Concatenate the "high" results for all other signatures (supps)
signature_list=lapply(signatures_of_interest,plot_sig_pat_freqs)
sig_1=signature_list[[1]]
merged_sig=sig_1
colnames(merged_sig)[4]="split"

for(i in 2:length(signature_list)){
  sig_i=signature_list[[i]]
  colnames(sig_i)[4]="split"
  merged_sig=merge(merged_sig,sig_i, by=c("pat","response","response2","split"))
}

merged_sig[, freq_pat_db_30:=NULL]
m_merged_sig=melt(merged_sig, id.vars=c("pat","response","response2","split"))
m_merged_sig[, variable:= gsub("freq_pat_","",variable)]

ggplot(m_merged_sig[split=="high"], aes(response, value))+
  geom_boxplot(aes(fill = response))+
  geom_beeswarm(aes(colour = response2), size=2.5, alpha=1,cex=7)+
  stat_compare_means(method="t", size=2.75, label.y.npc = 1, method.args = list(var.equal = FALSE))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))+
  scale_fill_manual(values=c("NR"="grey60", "R"="#c1de87"))+
  theme_bw()+ylab("Freq in top tertile")+theme(panel.grid=element_blank())+
  facet_wrap(~variable, nrow=1)
ggsave("./output_figs/rev3_S5E_aggregated_top_tertiles_supp.pdf", width=7,height=4)
write.csv(m_merged_sig[split=="high"], "./output_figs/rev3_S5E_aggregated_top_tertiles_supp.txt")

