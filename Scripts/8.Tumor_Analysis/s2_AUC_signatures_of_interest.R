##################################################################
# Interacting tumor cells phenotype highlight
##################################################################

# required input data in the working directory:
# # Tumor_Annotated.Rds from 1.Main/Tumor_Analysis_Code.R
# # c2.cp.v7.5.1.symbols.gmt from data/signatures_tcells
# # auc_numbers_merged.rds from 3.Signatures/s3_AUC_object.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(ggrepel)
library(scales)
library(data.table)
library(AUCell)
set.seed(150799)

fig_path = "Results_tumor_focused_analysis"
if(!dir.exists(fig_path)) dir.create(fig_path)

five_tum=readRDS("Tumor_Annotated.Rds")
DefaultAssay(five_tum)="RNA"

five_tum=subset(five_tum, annotated_clusters != "Low Gene Count Cells")

exprMatrix= five_tum@assays$RNA@counts
exprMatrix <- as(exprMatrix, "dgCMatrix")

gene.set=fgsea::gmtPathways("c2.cp.v7.5.1.symbols.gmt") 
gene.set=gene.set[names(gene.set) %in% c(
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "PID_HIF1_TFPATHWAY",
  "WP_AEROBIC_GLYCOLYSIS",#
  "WP_PHOTODYNAMIC_THERAPYINDUCED_HIF1_SURVIVAL_SIGNALING")]

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(gene.set, cells_rankings)

set.seed(333)
auc_numbers=t(data.frame(cells_AUC@assays@data$AUC))
auc_numbers=data.frame(auc_numbers)
auc_numbers$cell_id=rownames(auc_numbers)

meta=five_tum@meta.data
meta$cell_id=gsub("-",".",rownames(meta))
auc_numbers_merged=merge(auc_numbers, meta, by="cell_id")
saveRDS(auc_numbers_merged,"auc_numbers_merged_filt.rds")

setDT(auc_numbers_merged)
auc_numbers_merged[, sample_id:= factor(sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell"))]
auc_avg=auc_numbers_merged[, lapply(.SD, mean), by=c("patient_id", "sample_id"),
                                      .SDcols=names(gene.set)]

ggplot(auc_numbers_merged, aes(patient_id, fill=sample_id,y = PID_HIF1_TFPATHWAY))+
  geom_boxplot(alpha=1, width=1)+
  stat_compare_means(size=4, method = "wilcox.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  ylab("AUC score: HIF1 TF pathway")+
  theme_bw()+theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  scale_fill_manual(values=c(DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))+xlab("")
ggsave("Results_tumor_focused_analysis/HIF_score_per_pat.pdf", height=4.5, width=5)

ggplot(auc_numbers_merged, aes(patient_id, fill=sample_id,y = REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC))+
  geom_boxplot(alpha=1, width=1)+
  stat_compare_means(size=4, method = "wilcox.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  ylab("AUC score: Antigen presentation,\n folding assembly and\n peptide loading of class I MHC")+
  theme_bw()+theme( panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  scale_fill_manual(values=c(DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))+xlab("")
ggsave("Results_tumor_focused_analysis/MHC_score_per_pat.pdf", height=4.5, width=5)

ggplot(auc_numbers_merged, aes(patient_id, fill=sample_id,y = REACTOME_INTERFERON_ALPHA_BETA_SIGNALING))+
  geom_boxplot(alpha=1, width=1)+
  stat_compare_means(size=4, method = "wilcox.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  ylab("AUC score: IFN A/B pathway")+
  theme_bw()+theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  scale_fill_manual(values=c(DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))+xlab("")
ggsave("Results_tumor_focused_analysis/IFNAB_score_per_pat.pdf", height=4.5, width=5)

