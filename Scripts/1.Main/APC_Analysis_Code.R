##################################################################
#Subsetting only for APCs and APC characterization
##################################################################

# required input data in the working directory:
# # Combined_All_Cells.Rds from Allcells_Analysis_Code.R
# # MitoCarta2_human.txt from data folder

#setwd("/YOUR/PATH/")
library(Matrix)
library(Seurat)
library(ggrepel)
library(tidyverse)
library(harmony)
library(RSpectra)#0.16-2
library(DoubletFinder)#2.0.4
library(scGate)

set.seed(2000)
fig_path="Results_APCs"
if(!dir.exists(fig_path)) dir.create(fig_path) 

# object all patients from intial script
S_all = readRDS("Results_Main/Combined_All_Cells.Rds")
S_all<-SetIdent(S_all, value = S_all@meta.data$seurat_clusters)

# select APC compartment
pdf(file=paste0(fig_path,"/Features_UMAP_genes_APC.pdf"), width=14, height=10)
FeaturePlot(S_all,features = c("CD8A","SPI1","CD40","CD68","CD14","IGKC","S100A9","IRF7","CLEC9A","CD1C"))
dev.off()

S_all@meta.data$select_APCs = "none"
DimPlot(S_all,reduction="umap",label=TRUE)
S_all@meta.data[S_all@meta.data$seurat_clusters %in% c( 4,8,13,19,21,24,25,26,27,20),]$select_APCs = "yes"

pdf(file=paste0(fig_path,"/UMAP_clusters_all_patients_APC_selection.pdf"), width=7, height=7)
DimPlot(S_all, reduction = "umap",group.by = "select_APCs")+
  theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

apc_subset<-subset(x = S_all, subset = seurat_clusters %in% c( 4,8,13,19,21,24,25,26,27,20))

# # # # # # # # # # # # # # # # # #    PRE-FILTERING    # # # # # # # # # # # # # # # # # #

# load dimensionality processing of the object to remove doublets
apc_subset <- NormalizeData(apc_subset, normalization.method = "LogNormalize", scale.factor = 10000)
apc_subset <- FindVariableFeatures(apc_subset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(apc_subset)
apc_subset <- ScaleData(object = apc_subset, features =all.genes)
ElbowPlot(apc_subset, ndims = 20, reduction = "pca")
top.genes<-VariableFeatures(apc_subset)
apc_subset <- RunPCA(apc_subset, npcs = 20, verbose = FALSE,features=top.genes)
apc_subset <- RunHarmony(apc_subset, group.by.vars = "patient_id")
apc_subset <- RunUMAP(apc_subset, reduction = "harmony", dims = 1:20)
apc_subset <- FindNeighbors(apc_subset, reduction = "harmony", dims = 1:20)
apc_subset <- FindClusters(apc_subset)


# # some QCs
Idents(apc_subset) <- "combined_id" # includes hashtag info
apc_subset <- PercentageFeatureSet(apc_subset, "^RP[SL]", col.name = "percent_ribo")
apc_subset <- PercentageFeatureSet(apc_subset, pattern = "^MT-",col.name = "percent_mito")

FeaturePlot(apc_subset,features = "percent_ribo")+ theme_bw() +
  theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FeaturePlot(apc_subset,features = "percent_mito")+ theme_bw() +
  theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
FeaturePlot(apc_subset,features = "nFeature_RNA")+ theme_bw() +
  theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# # markers for seq doublets/ contaminations
pdf(file=paste0(fig_path,"/tumor_check_UMAP_genes_Tumor.pdf"), width=10, height=12)
FeaturePlot(apc_subset,features = c("KIT","PTPRC", "MLANA","MCAM",
                                    "SERPINE2","CD8A"), label=T)
dev.off()

# # Doublets check
# annotations <- apc_subset@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)
# nExp_poi <- round(0.075*nrow(apc_subset@meta.data))  ## Assuming 7.5% doublet formation rate 
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# apc_subset <- doubletFinder(apc_subset, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# 
# pdf(file=paste0(fig_path,"/APC_Doublets_Prefilter.pdf"), width=10, height=12)
# options(repr.plot.width=15, repr.plot.height=15)
# DimPlot(apc_subset, reduction = "umap",group.by="DF.classifications_0.25_0.09_1000",order = FALSE)
# dev.off()

# # Remove the problematic clusters
seurat_clusters_oi = unique(apc_subset@meta.data$seurat_clusters)
seurat_clusters_oi = as.character(seurat_clusters_oi[!(seurat_clusters_oi %in% c(7,17,8,16,18,19,21))])
apc_subset_sel <- subset(x = apc_subset, subset = seurat_clusters %in% seurat_clusters_oi)

# # Reprocess as validation
apc_subset_sel <- NormalizeData(apc_subset_sel, normalization.method = "LogNormalize", scale.factor = 10000)
apc_subset_sel <- FindVariableFeatures(apc_subset_sel, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(apc_subset_sel)
apc_subset_sel <- ScaleData(object = apc_subset_sel, features =all.genes)
ElbowPlot(apc_subset_sel, ndims = 20, reduction = "pca")
top.genes<-VariableFeatures(apc_subset_sel)
apc_subset_sel <- RunPCA(apc_subset_sel, npcs = 20, verbose = FALSE,features=top.genes)
apc_subset_sel <- RunHarmony(apc_subset_sel, group.by.vars = "patient_id")
apc_subset_sel <- RunUMAP(apc_subset_sel, reduction = "harmony", dims = 1:20)
apc_subset_sel <- FindNeighbors(apc_subset_sel, reduction = "harmony", dims = 1:20)
apc_subset_sel <- FindClusters(apc_subset_sel)

# # inspect the new object
pdf(file=paste0(fig_path,"/results.first.filter.pdf"), width=10, height=10)
DimPlot(apc_subset_sel, reduction = "umap",group.by = "seurat_clusters",label=T)+
  theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

pdf(file=paste0(fig_path,"/second_filter_16_expresses_cd8a.pdf"), width=5, height=5)
VlnPlot(apc_subset_sel,features = c("CD8A"))
dev.off()

# # cluster 16 expresses CD8A, remove and reprocess
apc_subset_sel <- subset(x = apc_subset_sel, subset = seurat_clusters !="16")
apc_subset_sel <- NormalizeData(apc_subset_sel, normalization.method = "LogNormalize", scale.factor = 10000)
apc_subset_sel <- FindVariableFeatures(apc_subset_sel, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(apc_subset_sel)
apc_subset_sel <- ScaleData(object = apc_subset_sel, features =all.genes)
top.genes<-VariableFeatures(apc_subset_sel)
apc_subset_sel <- RunPCA(apc_subset_sel, npcs = 20, verbose = FALSE,features=top.genes)
apc_subset_sel <- RunHarmony(apc_subset_sel, group.by.vars = "patient_id")
apc_subset_sel <- RunUMAP(apc_subset_sel, reduction = "harmony", dims = 1:20)
apc_subset_sel <- FindNeighbors(apc_subset_sel, reduction = "harmony", dims = 1:20)
apc_subset_sel <- FindClusters(apc_subset_sel)


pdf(file=paste0(fig_path,"/results_after_second_filter_out.pdf"), width=10, height=10)
DimPlot(apc_subset_sel, reduction = "umap",group.by = "seurat_clusters",label=T)+
  theme_bw() + theme(legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
apc_subset_sel$new_clusters<-as.character(apc_subset_sel$seurat_clusters)
dev.off()
# # this is the final object

# # # # # # # # # # # # # # # # # #    SIGNATURE SCORING    # # # # # # # # # # # # # # # # # #
scGate_models_DB <- get_scGateDB()
models.TME <- scGate_models_DB$human$TME_HiRes
apc_subset_sel <- scGate(apc_subset_sel, model=models.TME,reduction = "harmony")

write.table(apc_subset_sel@meta.data, paste0(fig_path,"/scGate_TME_HiRes_APCs.txt"),sep="\t",quote=F)
output_scgate = read.table( paste0(fig_path,"/scGate_TME_HiRes_APCs.txt"),sep="\t")

pdf(file=paste0(fig_path,"/scGate_multi_associations.pdf"), width=10, height=10)
DimPlot(apc_subset_sel, group.by = "scGate_multi") # cluster 8 contains all three populations
dev.off()

# # re-cluster 8 (prolif)
norm_and_clust_20=function(x,res,dims){
  set.seed(150799)
  x <- NormalizeData(x) 
  y <- VariableFeatures(x)
  x <- ScaleData(x, features = rownames(x)) %>% RunPCA(verbose = FALSE, features=y, npcs = 20)
  x <- RunHarmony(x, group.by.vars = "patient_id")
  print(ElbowPlot(x, ndims = 20, reduction = "pca"))
  x <- RunUMAP(x, reduction = "harmony", dims = 1:dims)
  x <- FindNeighbors(x, reduction = "harmony", dims = 1:dims)
  x <- FindClusters(x, resolution = res, reduction = "harmony")
  return(x)
}
cluster_mixed=subset(apc_subset_sel, subset = seurat_clusters %in% c("8"))
cluster_mixed=norm_and_clust_20(cluster_mixed, res=0.1, dims=20) # low res since the differences should be big
DimPlot(cluster_mixed, label=T)

pdf(file=paste0(fig_path,"/proliferating_cluster_subdiv.pdf"), width=6, height=4.5)
DotPlot(cluster_mixed, features=c("Monocyte_UCell","Macrophage_UCell","Monocyte_not_cDC2_UCell",
                                  "pDC_UCell","cDC1_UCell","cDC2_UCell","DC3_UCell",
                                  "PanBcell_UCell","Bcell_UCell","Plasma_cell_UCell"))+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow1")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
dev.off()

apc_subset_sel=AddMetaData(apc_subset_sel,
                           ifelse(apc_subset_sel$seurat_clusters%in% c("8") &
                                    colnames(apc_subset_sel) %in% colnames(cluster_mixed)[cluster_mixed$RNA_snn_res.0.1=="0"],
                                  "8_MM", ifelse(apc_subset_sel$seurat_clusters%in% c("8") &
                                                   colnames(apc_subset_sel) %in% colnames(cluster_mixed)[cluster_mixed$RNA_snn_res.0.1==c("2")],
                                                 "8_DC",
                                                 ifelse(apc_subset_sel$seurat_clusters%in% c("8") &
                                                          colnames(apc_subset_sel) %in% colnames(cluster_mixed)[cluster_mixed$RNA_snn_res.0.1=="1"],
                                                        "8_B",as.character(apc_subset_sel$seurat_clusters)))),
                           col.name="subclustered_8")

# # # # # # # # # # # # # # # # # #    SEPARATE CLASSES    # # # # # # # # # # # # # # # # # #
pdf(file=paste0(fig_path,"/final_class_calls.pdf"), width=6, height=12)
DotPlot(apc_subset_sel, features=c("PanBcell_UCell","Bcell_UCell","Plasma_cell_UCell",
                                   "Myeloid_UCell","Monocyte_not_cDC2_UCell",
                                   "Monocyte_UCell",'TAM_UCell',"Macrophage_UCell","MoMacDC_UCell",
                                   "cDC1_UCell","cDC2_UCell","DC3_UCell","pDC_UCell","cDC2b_UCell" ),
        group.by="subclustered_8")+
  scale_color_gradient2(high="red4", low="blue4", mid="lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
dev.off()

# # add the per class annotation based on signatures
apc_subset_sel=AddMetaData(apc_subset_sel, ifelse(apc_subset_sel$subclustered_8 %in% c("11","13","9","8_B"),"B_plasma_cells",
                                                  ifelse(apc_subset_sel$subclustered_8 %in% c("8_DC","2","15","14","12","10"),"DCs",
                                                         ifelse(apc_subset_sel$subclustered_8 %in% c("0","7","6","5","4","3","16","1",
                                                                                                     "8_MM"),"mono_mac","OOPS"))),
                           col.name = "major_classes")
DimPlot(apc_subset_sel, group.by="subclustered_8", label=T)
DimPlot(apc_subset_sel, group.by="major_classes")+
  scale_color_manual(values=c("B_plasma_cells"="palegreen3","DCs"="salmon","mono_mac"="#5D88E6"))

# # # # # # # # # # # # # # # # # #    ANALYSE INDIVIDUAL CLASSES    # # # # # # # # # # # # # # # # # #
# # higher dimensionality analysis + conflictive gene removal, for highest accuracy

# # genes to be removed from clustering analysis
gene_names_data = rownames(apc_subset_sel@assays$RNA)
mt_cands = grep("^MT-|^MTRN|^MTAT|^MTND|^MRP", gene_names_data, v=T, perl=T)
mito = data.table::fread("MitoCarta2_human.txt", header=T, sep="\t", stringsAsFactors=F)
mt_both = intersect(mt_cands, mito$Symbol)
mt_cands = setdiff(mt_cands, mt_both)
mitocarta = setdiff(mito$Symbol, mt_both)
mt_nms =  c(mt_cands, mitocarta, mt_both)
ig_nms = grep("^IGK|^IGL|^IGJ|^IGH|^IGBP|^IGSF", gene_names_data, v=T, perl=T)
ncrna_nms = c('MALAT1', 'XIST', 'NEAT1', 'hsa-mir-6723')
ncrp_nms = grep("^RP.[0-9]+", gene_names_data, v=T, perl=T)
tcr_genes = grep("^TRB[V|D|J|]|^TRA[V|D|J|]|^TRD[V|D|J|]|^TRG[V|D|J|]", gene_names_data, v=T, perl=T)

exclude_genes = c(mt_nms,ig_nms,ncrna_nms,ncrp_nms,
                  gene_names_data[grep("HIST", gene_names_data)],
                  gene_names_data[grep("HSP", gene_names_data)],tcr_genes)


norm_and_clust_2=function(x,res,dims){
  set.seed(150799)
  x <- NormalizeData(x) 
  y <- setdiff(VariableFeatures(x),exclude_genes) 
  x <- ScaleData(x, features = rownames(x)) %>% RunPCA(verbose = FALSE, features=y, npcs = 30)
  x <- RunHarmony(x, group.by.vars = "patient_id")
  print(ElbowPlot(x, ndims = 30, reduction = "pca"))
  x <- RunUMAP(x, reduction = "harmony", dims = 1:dims)
  x <- FindNeighbors(x, reduction = "harmony", dims = 1:dims)
  x <- FindClusters(x, resolution = res, reduction = "harmony")
  return(x)
}
# # bcells and plasma cells
apc_data_filtered=apc_subset_sel
bcell_filt=subset(apc_data_filtered, major_classes=="B_plasma_cells")
bcell_filt=norm_and_clust_2(bcell_filt, res=0.45, dims=20)
DimPlot(bcell_filt, label = T)

pdf(paste0(fig_path,"/bcell_breakdown.pdf"), height=5, width=5)
FeaturePlot(bcell_filt, c("CD27","FCER2"), label = T, order = T, pt.size = .5)
dev.off()

# cluster 1 contains CD27 and CD23 in separate populations
bcell_filt_1=norm_and_clust_2(subset(bcell_filt, RNA_snn_res.0.45=="1"), res=.2, dims=10)
DimPlot(bcell_filt_1, label = T)
FeaturePlot(bcell_filt_1, c("CD27","FCER2"))

pdf(paste0(fig_path,"/bcell_clus1_breakdown.pdf"), height=5, width=5)
VlnPlot(bcell_filt_1, c("CD27","FCER2"))
dev.off()

bcell_filt=AddMetaData(bcell_filt,
                       ifelse(bcell_filt$RNA_snn_res.0.45%in% c("1") &
                                colnames(bcell_filt) %in% colnames(bcell_filt_1)[bcell_filt_1$RNA_snn_res.0.2%in%c("1")],
                              "1_fcer2",ifelse(bcell_filt$RNA_snn_res.0.45%in% c("1") &
                                                 colnames(bcell_filt) %in% colnames(bcell_filt_1)[bcell_filt_1$RNA_snn_res.0.2%in%c("2")],
                                               "1_dn", ifelse(bcell_filt$RNA_snn_res.0.45%in% c("1") &
                                                                colnames(bcell_filt) %in% colnames(bcell_filt_1)[bcell_filt_1$RNA_snn_res.0.2%in%c("0")],
                                                              "1_cd27",as.character(bcell_filt$RNA_snn_res.0.45)))),
                       col.name="sub_1")

DotPlot(bcell_filt, features=c("CD19","FCER2","CD27","CD38","MKI67",
                               "SDC1","IGHA1","IGHG1","IGHG3", "IGKC"), group.by = "sub_1")+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow1")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=8),
        axis.text.y = element_text(size=8),
        legend.title = element_text(size=8))

DimPlot(bcell_filt, label = T, group.by = "sub_1")

bcell_filt=AddMetaData(bcell_filt, ifelse(bcell_filt$sub_1 %in% c("1_fcer2"), "CD27-CD23+ B cells",
                                          ifelse(bcell_filt$sub_1 %in% c("2"), "plasma cells",
                                                 ifelse(bcell_filt$sub_1 %in% c("0","1_dn"), "CD27-CD23- B cells",
                                                        ifelse(bcell_filt$sub_1 %in% c("1_cd27","4"), "CD27+CD23- B cells",
                                                               ifelse(bcell_filt$sub_1 == "3", "prol B cells","OOPS"))))),
                       col.name = "hi_res_clus")

bcell_filt$hi_res_clus=factor(bcell_filt$hi_res_clus, levels=c(
  "plasma cells", "prol B cells", "CD27+CD23- B cells",
  "CD27-CD23+ B cells", "CD27-CD23- B cells"))

DimPlot(bcell_filt, group.by = "hi_res_clus")
VlnPlot(bcell_filt, features = "FCER2", group.by = "hi_res_clus")+geom_boxplot()
VlnPlot(bcell_filt, features = "CD27", group.by = "hi_res_clus")+geom_boxplot()
VlnPlot(bcell_filt, features = "MKI67", group.by = "hi_res_clus")+geom_boxplot()

# # dendritic cells
dc_filt=subset(apc_data_filtered, major_classes=="DCs")
dc_filt=norm_and_clust_2(dc_filt, res=0.45, dims=20)

pdf(paste0(fig_path,"/DC_clus.pdf"), height=5, width=5)
DimPlot(dc_filt, label=T)
dev.off()

pdf(paste0(fig_path,"/DC_clus_dotplot.pdf"), height=8, width=11)
DotPlot(dc_filt, features=c(
  "LILRA4", "GZMB", "IL3RA","pDC_UCell",#pDCs
  "TCF4","SIGLEC6", "PPP1R14A","AXL", #AS
  "CD1C", "FCER1A", "HLA-DQA1", #cDC2s
  "CLEC10A","CD1E","cDC2_UCell",#melanoma cDC2s,
  "HSPA6","HSPH1","RPS27A", "RPLP1",
  "CLEC9A", "FLT3", "IDO1","MKI67","cDC1_UCell", #cDC1s
  "LAMP3", "CCR7", "FSCN1","IL4I1","DC3_UCell"# DC3s
))+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=8),
        axis.text.y = element_text(size=8),
        legend.title = element_text(size=8))
dev.off()

# # # add the final clusters
dc_filt=AddMetaData(dc_filt, ifelse(dc_filt$RNA_snn_res.0.45 %in% c("0","1","2","5"), "cDC2",
                                    ifelse(dc_filt$RNA_snn_res.0.45 %in% c("8"), "SIGLEC6+ DC",  
                                           ifelse(dc_filt$RNA_snn_res.0.45 =="3", "cDC1",
                                                  ifelse(dc_filt$RNA_snn_res.0.45 %in% c("4","9"), "mreg DC",
                                                         ifelse(dc_filt$RNA_snn_res.0.45 == "6", "pDC",
                                                                ifelse(dc_filt$RNA_snn_res.0.45 == "7", "prol cDC1","forgot_this_one")))))),
                    col.name = "hi_res_clus")

dc_filt$hi_res_clus=factor(dc_filt$hi_res_clus, levels=c("mreg DC", "cDC2","SIGLEC6+ DC", "prol cDC1", "cDC1", "pDC"))
DimPlot(dc_filt, group.by = "hi_res_clus")


# # monocytes and macrophages
mm_filt=subset(apc_data_filtered, major_classes=="mono_mac")
mm_filt=norm_and_clust_2(mm_filt, res = 0.6, dims=20)

pdf(paste0(fig_path,"/mm_clus.pdf"), height=6, width=6)
DimPlot(mm_filt, label = T)
dev.off()

pdf(paste0(fig_path,"/mm_clus_dotplot.pdf"), height=8, width=11)
DotPlot(mm_filt, group.by="RNA_snn_res.0.6", features=c(
  "FCN1","S100A8","S100A9", #mono cd14
  "PILRA", "FCGR3A", "HK3", #mono cd16
  "VCAN","TREM1","OLR1","VEGFA", # vcan
  "ISG15","CCL2","IL1RN",
  "DOCK4",# dock 4
  "C1QA", "C1QB", "C1QC", #C1Q
  "SLCO2B1","APOC1", #C1QC zhang
  "MKI67", #prol
  "RPS23","RPS27","RPL39",# adding ribosomal genes
  "CXCL9", "CXCL10", #isg15/inflamm_filtatory 
  "CD74","HLA-A","HLA-B","HLA-DRB5", #other inflamm_filtatory/AP
  "STAB1","TREM2","DAB2","FOLR2", # new lipid associated
  "SPP1", "MMP9","GPNMB","APOE" #spp1
))+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=8),
        axis.text.y = element_text(size=8),
        legend.title = element_text(size=8))
dev.off()

mm_filt$hi_res_clus=ifelse(mm_filt$RNA_snn_res.0.6=="6","FCN1hi monocyte-like cells",
                           ifelse(mm_filt$RNA_snn_res.0.6=="8","CD16hi monocyte-like cells",
                                  ifelse(mm_filt$RNA_snn_res.0.6=='4', "DOCK4hi macrophages",
                                         ifelse(mm_filt$RNA_snn_res.0.6=="0","ISG15hi mono/macrophages",
                                                ifelse(mm_filt$RNA_snn_res.0.6=="2","VCANhi mono/macrophages",
                                                       ifelse(mm_filt$RNA_snn_res.0.6=="5","C1Qhi inflammatory macrophages",
                                                              ifelse(mm_filt$RNA_snn_res.0.6=="7","C1Qhi MKI67+ prol macrophages",
                                                                     ifelse(mm_filt$RNA_snn_res.0.6=="9","SPP1hi macrophages",
                                                                            ifelse(mm_filt$RNA_snn_res.0.6 %in% c("3"),"C1Qhi Ribohi macrophages",
                                                                                   ifelse(mm_filt$RNA_snn_res.0.6 %in% c("1"),"C1Qhi LA macrophages","OOPS"
                                                                                   ))))))))))

mm_filt$hi_res_clus=factor(mm_filt$hi_res_clus, levels=c(
  "SPP1hi macrophages","C1Qhi LA macrophages","C1Qhi inflammatory macrophages",
  "C1Qhi Ribohi macrophages", "C1Qhi MKI67+ prol macrophages", 
  "DOCK4hi macrophages", "ISG15hi mono/macrophages",
  "VCANhi mono/macrophages","CD16hi monocyte-like cells",
  "FCN1hi monocyte-like cells"
))
DimPlot(mm_filt, group.by = "hi_res_clus")

# 3rd: create the overall object annotation
m_sum <- data.frame(
  og_clus = mm_filt@meta.data$`RNA_snn_res.0.6`,
  hi_res_clus = mm_filt@meta.data$hi_res_clus,
  row.names = rownames(mm_filt@meta.data))

dc_sum <- data.frame(
  og_clus = dc_filt@meta.data$`RNA_snn_res.0.45`,
  hi_res_clus = dc_filt@meta.data$hi_res_clus,
  row.names = rownames(dc_filt@meta.data))

bc_sum <- data.frame(
  og_clus = bcell_filt@meta.data$sub_1,
  hi_res_clus = bcell_filt@meta.data$hi_res_clus,
  row.names = rownames(bcell_filt@meta.data))

meta_2=do.call(rbind, list(m_sum,bc_sum,dc_sum))
meta_2$og_clus=paste0(meta_2$hi_res_clus,"_",meta_2$og_clus)

# # add the metadata
apc_data_filtered=AddMetaData(apc_data_filtered, meta_2)
DimPlot(apc_data_filtered, group.by = "hi_res_clus", label = T)


levels=c(
  "plasma cells", "prol B cells", "CD27+CD23- B cells",
  "CD27-CD23+ B cells","CD27-CD23- B cells",
  "mreg DC", "cDC2", "SIGLEC6+ DC","prol cDC1","cDC1", "pDC",
  "SPP1hi macrophages", "C1Qhi LA macrophages","C1Qhi inflammatory macrophages",
  "C1Qhi Ribohi macrophages","C1Qhi MKI67+ prol macrophages",
  "DOCK4hi macrophages",
  "ISG15hi mono/macrophages","VCANhi mono/macrophages",
  "CD16hi monocyte-like cells", "FCN1hi monocyte-like cells"
)

apc_data_filtered$hi_res_clus=factor(apc_data_filtered$hi_res_clus, levels=levels)

palette_apc=c(
  "plasma cells"="#032507",
  "prol B cells"="#14743E",
  "CD27+CD23- B cells"="#33A46B",
  "CD27-CD23+ B cells"="palegreen3",
  "CD27-CD23- B cells"="palegreen1",
  
  "mreg DC"="firebrick4",
  "cDC2"="salmon1",
  "SIGLEC6+ DC"="salmon4",
  "prol cDC1"="lightpink",
  "cDC1"="#D3457B", 
  "pDC"="#7900C6",
  
  "SPP1hi macrophages"= "#FF66FF",
  "C1Qhi MKI67+ prol macrophages"= "#910099",
  "C1Qhi Ribohi macrophages"= "#A284F1",
  "C1Qhi inflammatory macrophages"= "#0000BF",
  "C1Qhi LA macrophages"= "#739FFF",
  "DOCK4hi macrophages"= "steelblue4",
  "ISG15hi mono/macrophages"= "#E4D3FF",
  "VCANhi mono/macrophages"= "#B9F2FF",
  "CD16hi monocyte-like cells"= "#00FFFF",
  "FCN1hi monocyte-like cells"= "#008B8B"
)

DimPlot(apc_data_filtered, group.by = "hi_res_clus", label=T)+scale_color_manual(values=palette_apc)+theme(legend.position = "none")
DimPlot(apc_data_filtered, group.by = "og_clus", label=T)+theme(legend.position = "none")

##### export the objects
saveRDS(apc_data_filtered, paste0(fig_path,"/apc_data_filtered_fvf_corr.rds"))
saveRDS(mm_filt, paste0(fig_path,"/final_mm_fvf_corr.rds"))
saveRDS(dc_filt, paste0(fig_path,"/final_dc_fvf_corr.rds"))
saveRDS(bcell_filt,paste0(fig_path, "/final_bcell_fvf_corr.rds"))
saveRDS(palette_apc, paste0(fig_path,"/palette_apc_fvf_corr.rds"))

#apc_data_filtered_old<-readRDS("/Users/a.george/Documents/CD39_revision_code_check/all_processed_objects/apc_data_filtered_fvf_corr.rds")

#table(apc_data_filtered_old$major_classes)
#table(apc_data_filtered$major_classes)

