##################################################################
# Unbiased interacting vs singlet tumor cells
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
set.seed(150799)

fig_path = "Results_tumor_unbiased_analysis"
if(!dir.exists(fig_path)) dir.create(fig_path)

five_tum=readRDS("Results_tumor/Tumor_Annotated.Rds")
DefaultAssay(five_tum)="RNA"

# DGE at cell level
# using MAST corrected for patient
five_tum=subset(five_tum, annotated_clusters != "Low Gene Count Cells") 
meta=five_tum@meta.data

Idents(five_tum)="sample_id"
res=FindMarkers(five_tum, ident.1="DB_Tumor_Tcell", "SG_Singlets", 
                logfc.threshold = 0, min.pct=0.1, 
                test.use = "MAST",
                latent.vars = "patient_id")
res=data.frame(res)
res$gene=rownames(res)

res=setDT(res)
res[, class:=ifelse(p_val_adj<0.05 & avg_log2FC<(-.05), "SG",
                    ifelse(p_val_adj<0.05 & avg_log2FC>(.05),"DB","NS"))]

res[, enrichment:=ifelse(p_val_adj<0.05 & avg_log2FC<(-0), "SG",
                    ifelse(p_val_adj<0.05 & avg_log2FC>(0),"DB","NS"))]

fwrite(res,"Results_tumor_unbiased_analysis/tum_DB_vs_tum_SG_filt.txt")

# Plotting the volcano
panel=c("HLA-A","HLA-B","HLA-C","HLA-E","TAP1","PSMB9",
        "ANGPTL4","VEGFA","ALDOA","SOD2","CXCL10",'PLIN2',
        "PLIN2","ENO2","CCL5","A2M",'IL32',
        "PKM","IFITM3","FN1",
        "NRG3", "APOD", "PDE3B", "DCT", "PMEL",
        "MLANA", "ZNF106","TYR","CD74")

ggplot(res, aes(avg_log2FC, -log10(p_val_adj)))+
  geom_point(aes(col=class))+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, col="black")+
  scale_color_manual("Significance", values=c("SG"="#7cd3f7", "NS"="grey", "DB"="#825d9e"))+
  scale_fill_manual("Significance", values=c("SG"="#7cd3f7", "NS"="grey", "DB"="#825d9e"))+
  geom_label_repel(data=res[p_val_adj<0.05 & abs(avg_log2FC)>.05 & gene %in% panel,],force = .3,label.padding = .2, max.overlaps = 12,
                   aes(label=gene, fill=class),col="white", size=5, segment.colour = "black")+
  xlim(c(-.3,.3))
ggsave("Results_tumor_unbiased_analysis/sc_tum_volcano_big_labels.pdf", width=8, height = 8)

write.csv(res,"Results_tumor_unbiased_analysis/sc_tum_volcano_big_labels.csv")
# Running GSEA analysis
gene.set=fgsea::gmtPathways("extdata/signatures_tcells/c2.cp.v7.5.1.symbols.gmt") 

setDT(res)
res_filt=res[,stats:=ifelse(avg_log2FC>0,-log10(p_val),log10(p_val))] # weight the pval
setorder(res_filt, -stats)
# create stats
stats=res_filt$stats
names(stats)=res_filt$gene
print(head(stats))
print(tail(stats))
#run results
set.seed(150799)
res_fgsea=fgsea::fgseaMultilevel(gene.set, stats, nPermSimple = 10000, minSize = 5, maxSize = 300)

setDT(res_fgsea)
res_fgsea=na.omit(res_fgsea)
res_fgsea[, sig:= ifelse(padj<0.05, "padj<0.05", "NS")]

setorder(res_fgsea, -NES)

fwrite(res_fgsea, "Results_tumor_unbiased_analysis/fgsea_results_filt.txt")

ff_res=res_fgsea[c(1:7, (nrow(res_fgsea)-6):nrow(res_fgsea))]

ggplot(ff_res, aes(NES, reorder(gsub("_"," ",pathway), NES), fill=NES))+
  geom_col(position="identity") +
  theme_bw()+
  theme(axis.text.y = element_text(size = 6), legend.key.width = unit(4,"mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(labels = label_wrap(width = 35))+
  ylab("")+scale_fill_gradient2("NES", high="#825d9e", 
                                low="#7cd3f7")+
  geom_point(data=ff_res[padj<0.05], shape=8)
ggsave("Results_tumor_unbiased_analysis/top.bottom_gsea.pdf", width=5, height=4.5)
df <- do.call(rbind,ff_res)
write.csv(df,"Results_tumor_unbiased_analysis/top.bottom_gsea.csv")
