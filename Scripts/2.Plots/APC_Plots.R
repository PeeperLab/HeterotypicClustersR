##################################################################
# APC phenotypical validations
##################################################################

# required input data in the working directory:
# # Objects output from APC_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(harmony)

fig_path = "Results_APC_plots"
if(!dir.exists(fig_path)) dir.create(fig_path) 

# 1st load the objects
apc_data_filtered=readRDS("Results_APCs/apc_data_filtered_fvf_corr.rds")
mm_filt=readRDS("Results_APCs/final_mm_fvf_corr.rds")
dc_filt=readRDS("Results_APCs/final_dc_fvf_corr.rds")
bcell_filt=readRDS("Results_APCs/final_bcell_fvf_corr.rds")
palette_apc=readRDS("Results_APCs/palette_apc_fvf_corr.rds") 

# 2nd: plot
# # first: by signatures
pdf("./Results_APC_plots/dotplot_UCell.signatures.pdf", width=5.5, height=5.5)
DotPlot(apc_data_filtered, features=c("Monocyte_UCell","Macrophage_UCell",
                                      "pDC_UCell","cDC1_UCell","cDC2_UCell","DC3_UCell",
                                      "PanBcell_UCell","Bcell_UCell","Plasma_cell_UCell"),
        group.by="hi_res_clus")+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=10),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=8))
dev.off()

# # second gene panel
pdf("./Results_APC_plots/dotplot_all_clusters.pdf", width=14.4, height=5.5)
DotPlot(apc_data_filtered, features=c(
  "FCN1","S100A8","S100A9", #mono cd14
  "PILRA", "FCGR3A", "HK3", #mono cd16
  "VCAN","TREM1","VEGFA","OLR1", # vcan
  "ISG15","CCL2","IL1RN",
  "DOCK4",# dock 4
  "C1QA", "C1QB", "C1QC", #C1Q
  "SLCO2B1","APOC1", #C1QC zhang
  "RPS23","RPS27","RPL39",# adding ribosomal genes
  "CXCL9", "CXCL10", #isg15/inflammatory
  "CD74","HLA-A","HLA-B","HLA-DRB5", #other inflamm_filtatory/AP
  "STAB1","TREM2","DAB2","FOLR2", # new lipid associated
  "SPP1", "MMP9","GPNMB","APOE", #spp1
  "LILRA4", "GZMB", "IL3RA",#pDCs
  "CLEC9A", "FLT3", "IDO1", #cDC1s
  "SIGLEC6", "PPP1R14A","AXL", # SIGLEC6 DCs
  "CD1C", "FCER1A", "HLA-DQA1", #cDC2s
  "CLEC10A","CD1E", #melanoma cDC2s
  "LAMP3", "CCR7", "FSCN1","IL4I1",# DC3s
  "CD19","FCER2","CD27","CD38",# bcell markers
  "SDC1","IGHA1","IGHG1","IGHG3", "IGKC", #immunoglobulins
  "MKI67","TOP2A" # proliferation markers
),
group.by="hi_res_clus")+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        legend.title = element_text(size=8))
dev.off()


pdf("./Results_APC_plots/major_classes_UMAP.pdf", width=6, height=5)
DimPlot(apc_data_filtered, group.by="major_classes")+
  scale_color_manual(values=c("B_plasma_cells"="palegreen3",
                              "DCs"="salmon","mono_mac"="#5D88E6"))+
  ggtitle("")
dev.off()

pdf("./Results_APC_plots/clusters_UMAP.pdf", width=6, height=5)
DimPlot(apc_data_filtered, group.by="hi_res_clus")+
  scale_color_manual(values=palette_apc)+
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7))+ggtitle("")
dev.off()


# # per compartment
# # # MACROPHAGES

pdf("./Results_APC_plots/mm_UMAP_per_sample_id_singlets.pdf", width=3.5, height=3.5)
mm_filt$sample_id=factor(mm_filt$sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell", "DB_APC_Tcell"))
DimPlot(subset(mm_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("SG_Singlets"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="grey",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "#7cd3f7"))
dev.off()
pdf("./Results_APC_plots/mm_UMAP_per_sample_DB_APCs.pdf", width=3.5, height=3.5)
mm_filt$sample_id=factor(mm_filt$sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell", "DB_APC_Tcell"))
DimPlot(subset(mm_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("DB_APC_Tcell"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="#009052",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "grey"))
dev.off()

pdf("./Results_APC_plots/mm_UMAP_annotated.pdf", width=6, height=4)
DimPlot(mm_filt, group.by = "hi_res_clus")+
  scale_color_manual(values=palette_apc)+
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7))
dev.off()

# # # DENDRITIC CELLS
pdf("./Results_APC_plots/DC_UMAP_per_sample_id_singlets.pdf", width=3.5, height=3.5)
DimPlot(subset(dc_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("SG_Singlets"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="grey",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "#7cd3f7"))
dev.off()

pdf("./Results_APC_plots/DC_UMAP_per_sample_id_DB_APCs.pdf", width=3.5, height=3.5)
DimPlot(subset(dc_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("DB_APC_Tcell"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="#009052",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "grey"))
dev.off()
pdf("./Results_APC_plots/DC_UMAP.pdf_annotated", width=6, height=5)
DimPlot(dc_filt, group.by = "hi_res_clus")+
  scale_color_manual(values=palette_apc)+
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7))
dev.off()

# # Bcell object

pdf("./Results_APC_plots/Bcells_UMAP_per_sample_id_singlets.pdf", width=3.5, height=3.5)
bcell_filt$sample_id=factor(bcell_filt$sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell", "DB_APC_Tcell"))
DimPlot(subset(bcell_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("SG_Singlets"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="grey",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "#7cd3f7"))
dev.off()

pdf("./Results_APC_plots/Bcells_UMAP_per_sample_id_DB_APCs.pdf", width=3.5, height=3.5)
DimPlot(subset(bcell_filt, subset = sample_id!="DB_Tumor_Tcell"),order=c("DB_APC_Tcell"), group.by = "sample_id")+
  theme(axis.text = element_text(size=6),
        legend.position="none",
        strip.text = element_text(size=8),
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_color_manual(values=c(DB_APC_Tcell="#009052",
                              DB_Tumor_Tcell="grey",
                              SG_Singlets = "grey"))
dev.off()

pdf("./Results_APC_plots/B_Bcells_UMAP_annotated.pdf", width=6, height=5)
DimPlot(bcell_filt, group.by = "hi_res_clus")+
  scale_color_manual(values=palette_apc)+
  theme(legend.text = element_text(size=9),
        legend.title = element_text(size=9))
dev.off()

# extra parameters
apc_data_filtered=subset(apc_data_filtered, subset= sample_id!="DB_Tumor_Tcell")
meta_apc=data.table(apc_data_filtered@meta.data)

pdf("./Results_APC_plots/counts_per_patient.pdf", width=4.5, height=4)
ggplot(meta_apc, aes(hi_res_clus, fill=patient_id))+geom_bar(position="stack")+
  facet_wrap(~major_classes, nrow=1, scales="free_x")+theme_bw()+
  scale_fill_viridis_d()+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                               strip.text = element_text(size=8),
                               strip.background = element_rect( linewidth = 0, fill="white"))+xlab("")
dev.off()
table_counts <- meta_apc %>%
  count(major_classes, hi_res_clus, patient_id)
write.csv(table_counts,"Results_APC_plots/counts_per_patient.csv")
pdf("./Results_APC_plots/counts_per_sample_id.pdf", width=4, height=4)
ggplot(meta_apc, aes(hi_res_clus, fill=sample_id))+geom_bar(position="stack")+
  facet_wrap(~major_classes, nrow=1, scales="free_x")+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.text = element_text(size=8),legend.position = "none",
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_fill_manual(values=c(DB_APC_Tcell="#009052",
                             DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))+xlab("")
dev.off()
table_counts2 <- meta_apc %>%
  count(major_classes, hi_res_clus, sample_id)

head(table_counts2)
# Wide format (matrix style, one column per sample_id)
table_wide2 <- table_counts2 %>%
  pivot_wider(names_from = sample_id, values_from = n, values_fill = 0)
write.csv(table_wide2,"Results_APC_plots/counts_per_sample_id.csv")
head(table_wide2)
pdf("./Results_APC_plots/freq_per_sample_id.pdf", width=4, height=4)
ggplot(meta_apc, aes(hi_res_clus, fill=sample_id))+geom_bar(position="fill")+
  facet_wrap(~major_classes, nrow=1, scales="free_x")+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.text = element_text(size=8),legend.position = "none",
        strip.background = element_rect( linewidth = 0, fill="white"))+
  scale_fill_manual(values=c(DB_APC_Tcell="#009052",
                             DB_Tumor_Tcell="#825d9e",
                             SG_Singlets = "#7cd3f7"))+xlab("")+ylab("Total frequency")
dev.off()

write.csv(meta_apc,"Results_APC_plots/APC_metadata.csv")
