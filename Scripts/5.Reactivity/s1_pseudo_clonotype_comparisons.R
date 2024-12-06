##################################################################
# Clonotypes per interaction status analysis
##################################################################

# required input data in the working directory:
# # TCR_filtered_five_pats.txt from Tcell_Plots.R
# # auc_numbers_merged.rds from s3_AUC_object.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

fig_path = "Results_clonotype_comparison"
if(!dir.exists(fig_path)) dir.create(fig_path)

tcr_ids=fread("Meta_data_TCRs_Tcells.txt")
tcr_ids$barcode<-paste0(tcr_ids$combined_id,"_",tcr_ids$cell_id)
tcr_ids<-tcr_ids[,c("barcode","combined_id","patient_id","chain_cdr3")]
colnames(tcr_ids)<-c("barcode", "sample" ,"patient_id", "TCR" )

tcr_ids[, cell_id:= gsub(paste0(sample,"_"), "", barcode), by="sample"]

auc57=readRDS("auc_numbers_merged.rds")
setDT(auc57)
auc57$cell_id=gsub("\\.","-",auc57$cell_id)

auc_dt=merge(auc57, tcr_ids, by=c("cell_id", "patient_id"))

signature_panel=c("neotcr_cd8","offringa_TR_9_samples",
                  "oliv_Virus.specific","oliv_Tumor.specific")
auc_mat=auc_dt[,colnames(auc_dt) %in% signature_panel]

levels=c("Tn", "Tn/Tm", "Early Tem", "Tem", "Tem-NK like", "Tc17 MAIT", "ISG+",
         "GZMK hi Tex", "TCF7+ stem-like Tex", "TOX hi Tex", "LAG3 hi Tex",
         "MKI67+ Tex/Tprol", "MKI67 hi Tex/Tprol", "MKI67+ Tem-NK like/Tprol")

auc_dt[, sample_id:=factor(sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell","DB_APC_Tcell"))]
auc_dt[, annotated_clusters_labels_final:=factor(annotated_clusters_labels_final, levels=levels)]

# tcr/interaction parameters - create groups per clonotype and sample_id
auc_dt[, id_TCR:= paste0(patient_id, "_", sample_id, "_", TCR)]
auc_dt[, n_id_TCR:= .N, by="id_TCR"]

tcr_scor=auc_dt[,lapply(.SD, mean), by=c("id_TCR", "patient_id", "sample_id","n_id_TCR"),
              .SDcols = c(signature_panel, "Tirosh_Mel_Exh","Yost_CD8.Exh",
                          "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY", "db_30","db_100",
                          "WP_TCELL_RECEPTOR_TCR_SIGNALING_PATHWAY",
                          "REACTOME_CTLA4_INHIBITORY_SIGNALING")]


# # without splitting by patient
mavg_sample=melt(tcr_scor, 
                 id.vars=c("id_TCR", "n_id_TCR", "patient_id","sample_id"),
                 measure.vars = c(signature_panel, "Tirosh_Mel_Exh","Yost_CD8.Exh",
                                  "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY", "db_30","db_100",
                                  "WP_TCELL_RECEPTOR_TCR_SIGNALING_PATHWAY",
                                  "REACTOME_CTLA4_INHIBITORY_SIGNALING"))

setDT(mavg_sample)
mavg_sample[, variable:= gsub("Tirosh_Mel_Exh", "Exh. Tirosh",variable)]
mavg_sample[, variable:= gsub("offringa_TR_9_samples", "React., Offringa",variable)]
mavg_sample[, variable:= gsub("neotcr_cd8", "React., Lowery",variable)]
mavg_sample[, variable:= gsub("oliv_Virus.specific", "Virus sp, Oliv.",variable)]
mavg_sample[, variable:= gsub("oliv_Tumor.specific", "Tum sp, Oliv.",variable)]
mavg_sample[, variable:= gsub("WP_TCELL_RECEPTOR_TCR_SIGNALING_PATHWAY", "TCR signaling",variable)]
mavg_sample[, variable:= gsub("REACTOME_CTLA4_INHIBITORY_SIGNALING", "CTLA4 signaling",variable)]
mavg_sample[, variable:= gsub("REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY", "CD28 signaling",variable)]
mavg_sample[, variable:= gsub("Yost_CD8.Exh", "Exh. Yost",variable)]

mavg_sample[,variable:= factor(variable, levels = c("Virus sp, Oliv.","Tum sp, Oliv.", "React., Offringa","React., Lowery",
                                             "Cytotoxicity","TCR signaling","db_30","db_100", "Exh. Yost",
                                             "Exh. Tirosh", "CD28 signaling", "CTLA4 signaling"))]

ggplot(mavg_sample[variable %in% c("Virus sp, Oliv.","Tum sp, Oliv.", "React., Offringa","React., Lowery"),], aes(sample_id, value, fill=sample_id))+
  geom_violin()+geom_boxplot(alpha=.4)+
  stat_compare_means(size=3,comparisons = list(c("DB_APC_Tcell", "DB_Tumor_Tcell"),
                                               c("DB_Tumor_Tcell","SG_Singlets"),
                                               c("DB_APC_Tcell","SG_Singlets")))+
  facet_wrap(~variable, nrow=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),
        axis.text.y = element_text(size=6),
        legend.position = "none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  ylab("Avg signature score per TCR group")+xlab("")+
  scale_fill_manual(values=c(DB_APC_Tcell="#009052",
                             DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))
ggsave("Results_clonotype_comparison/heatmap_panel_across_all_pseudo_clonotypes_all_pats_tcrs_2.pdf", width=22, height=5)

