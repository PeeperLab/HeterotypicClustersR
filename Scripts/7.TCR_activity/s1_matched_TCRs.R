##################################################################
# Expanded/ matched TCR comparison across interacting partners
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds from 1.Main/Tcell_Analysis_Code.R
# # Meta_data_TCRs_Tcells.txt from 2.Plots/Tcell_Plots.R
# # auc_numbers_merged.rds from 3.Signatures/s3_AUC_object.R

#setwd("/YOUR/PATH/")
library(data.table)
library(Seurat)
set.seed(150799)

five_pat=readRDS("Tcells_Final.Rds")

tcr_ids=fread("Meta_data_TCRs_Tcells.txt")
tcr_ids$barcode<-paste0(tcr_ids$combined_id,"_",tcr_ids$cell_id)
tcr_ids<-tcr_ids[,c("barcode","combined_id","patient_id","chain_cdr3")]
colnames(tcr_ids)<-c("barcode", "sample" ,"patient_id", "TCR" )


auc57=readRDS("auc_numbers_merged.rds")
setDT(auc57)
auc57$cell_id=gsub("\\.","-",auc57$cell_id)

# 1/ TCR selection
meta_tcells=data.frame(five_pat@meta.data)
meta_tcells$cell_id=rownames(meta_tcells)

tcr_ids[, cell_id:= gsub(paste0(sample,"_"), "", barcode), by="sample"]
tcr_ids=merge(meta_tcells, tcr_ids, by=c("cell_id","patient_id"))
setDT(tcr_ids)

tcr_ids[, patient_TCR:= paste0(patient_id, "_", TCR)]
tcr_ids[, sample_TCR:= paste0(sample, "_", TCR)]

tcr_ids[, patient_counts:= .N, by="patient_id"]
tcr_ids[, sample_counts:= .N, by="sample"]
tcr_ids[, tcr_patient_counts:= .N, by="patient_TCR"]
tcr_ids[, tcr_sample_counts:= .N, by="sample_TCR"]

tcr_ids[, patient_freqs:= tcr_patient_counts/patient_counts, by="patient_id"]
tcr_ids[, sample_freqs:= tcr_sample_counts/sample_counts, by="sample"]

tcr_ids[, count_TUM:=sum(sample_id=='DB_Tumor_Tcell'), by="patient_TCR"]
tcr_ids[, count_APC:=sum(sample_id=='DB_APC_Tcell'), by="patient_TCR"]

TCR_selection=tcr_ids[count_TUM>=10 & count_APC>=10,] #expanded only

# parameters per TCR depending on the setting
sample_clonotype=TCR_selection[, lapply(.SD, mean), by=c("patient_TCR", "patient_id"),
                               .SDcols = c("tcr_patient_counts")]

# we are NOT taking draws
# IF we do it per patient

pat_list=split(sample_clonotype, by="patient_id")

order_pat_and_get_top10=function(x){
  x[, random_number:= sample(x=1:nrow(x), size=nrow(x),replace = F)]
  setorder(x, -tcr_patient_counts, random_number)
  x=x[1:10,] # whoever is at 10th stays
  return(x)
}

top_pat=do.call(rbind, lapply(pat_list, order_pat_and_get_top10))

# 2/ matched comparisons
meta_sel=tcr_ids[patient_TCR %in% top_pat$patient_TCR]
meta_sel=meta_sel[sample_id!="SG_Singlets",]

meta_sel$sample_id=factor(meta_sel$sample_id,
                          levels=c("SG_Singlets","DB_Tumor_Tcell","DB_APC_Tcell"))


# test phenotypes
library(ggplot2)
library(ggpubr)
setDT(meta_sel)
tcr_auc=merge(auc57, meta_sel[,.(cell_id, patient_TCR)], by=c("cell_id"))
setDT(tcr_auc)
tcr_auc[, sample_id:= factor(tcr_auc$sample_id, levels=c("DB_Tumor_Tcell", "DB_APC_Tcell"))]

avg_tcr_auc=tcr_auc[, lapply(.SD, mean), .SDcols=c("Tirosh_Mel_Exh",
                                                   "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY",
                                                   "REACTOME_CTLA4_INHIBITORY_SIGNALING"),
                    by=c("patient_id", "sample_id", "patient_TCR")]

all_auc_tcr=melt(avg_tcr_auc, id.vars=c("patient_id","sample_id", "patient_TCR"))
setDT(all_auc_tcr)
setorder(all_auc_tcr, patient_TCR, sample_id)

all_auc_tcr[, variable:= gsub("Tirosh_Mel_Exh", "Exh., Tirosh",variable)]
all_auc_tcr[, variable:= gsub("REACTOME_CTLA4_INHIBITORY_SIGNALING", "CTLA4 signaling",variable)]
all_auc_tcr[, variable:= gsub("REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY", "CD28 signaling",variable)]
all_auc_tcr[, variable:= factor(variable, levels=c("Exh., Tirosh", "CD28 signaling","CTLA4 signaling"))]

ggplot(all_auc_tcr[variable %in% c("Exh., Tirosh",
                                   "CD28 signaling","CTLA4 signaling"),], aes(sample_id, value))+
  geom_boxplot(outlier.shape =NA)+geom_point(aes(col=patient_id), alpha=.7, size=2)+
  facet_wrap(~variable,nrow=1, scales="free_y")+
  stat_compare_means(method="wilcox.test", paired = TRUE,
                     method.args = list(alternative = "two.sided"), size=2.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),panel.grid = element_blank(),
        axis.text.y = element_text(size=6),
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  ylab("AUC score per pat_TCR")+
  geom_line(aes(group=patient_TCR, col=patient_id))+xlab("")+
  scale_color_manual(values=c("P2"="#3b528b", "P3"="#21918c", "P4"="#5ec962", "P5"="#fde725"))
ggsave("matched_tcrs_exhaustion_costimulation.pdf", width=5.5, height=4)



