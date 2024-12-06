##################################################################
#Subsetting only for T cells and T cell characterization
##################################################################

# required input data in the working directory:
# # Combined_All_Cells.Rds from Allcells_Analysis_Code.R
# # MitoCarta2_human.txt from data folder

#R version 4.3.3
#setwd("/YOUR/PATH/")
library(Seurat)#Version 4.4.0 
library(ggplot2)
library(dplyr)
library(harmony)# Version 1.2.1
library(RSpectra)
library(Nebulosa)
set.seed(2000)

fig_path = "Results_Tcells"
if(!dir.exists(fig_path)) dir.create(fig_path) 


pbmc.combined<-readRDS("Combined_All_Cells.Rds")
#Further analysis of Tcells by same harmony integration by filtering out genes in the "MitoCarta2_human.txt" file related to mitochondria and other like TCR and ribo
tcell.combined<-subset(x = pbmc.combined, subset = seurat_clusters %in% c( 0,1,3,6,12,14,15 ))


tcell.combined <- NormalizeData(tcell.combined, normalization.method = "LogNormalize", scale.factor = 10000)
tcell.combined <- FindVariableFeatures(tcell.combined, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(tcell.combined)
tcell.combined <- ScaleData(object = tcell.combined, 
                            do.par = TRUE, num.cores = 8, features =all.genes)

# filter variable genes
gene_names_data = rownames(tcell.combined@assays$RNA)
## filter feature genes based on following sets:
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

## save all exluding genes
exclude_genes = c(mt_nms,ig_nms,ncrna_nms,ncrp_nms,
                  gene_names_data[grep("HIST", gene_names_data)],
                  gene_names_data[grep("HSP", gene_names_data)],tcr_genes)


top.genes<-VariableFeatures(tcell.combined)
ex_genes<-top.genes[!(top.genes %in% exclude_genes)]

tcell.combined <- RunPCA(tcell.combined, npcs = 30, verbose = FALSE,features=ex_genes)

ElbowPlot(tcell.combined, ndims = 20, reduction = "pca")

tcell.combined <- RunHarmony(tcell.combined, group.by.vars = "patient_id")

tcell.combined <- RunUMAP(tcell.combined, reduction = "harmony", dims = 1:15)

tcell.combined <- FindNeighbors(tcell.combined, reduction = "harmony", dims = 1:15)
tcell.combined <- FindClusters(tcell.combined, resolution = 2)
pdf("Results_Tcells/Tcells_without_clustering.pdf",width=13,height=12)

DimPlot(tcell.combined, reduction = "umap",order = TRUE, label = TRUE,pt.size=1.2,label.size = 12)
dev.off()
pdf("Results_Tcells/Tcells_filter_Markers.pdf")
options(repr.plot.width=15, repr.plot.height=15)
FeaturePlot(tcell.combined,features=c("FOXP3","ITGAX","CD4"),pt.size =0.1,order=T)
dev.off()

#saveRDS(tcell.combined,"Tcells_prefilter.Rds")

#Not including Treg cells, CD4 cells and  APCs from the further analysis marked by FOXP3, CD4, and ITGAX
#remove 3 seurat clusters in total two from bottom and one from top 
tcell.combined_subset<-subset(x = tcell.combined, subset = seurat_clusters %in% c( 0,1,2,3,4,5,6,7,8,9,10,12,13,15,16,17,18,19,20,
                                                                                   21,22,23,24,26))

tcell.combined_subset <- NormalizeData(tcell.combined_subset, normalization.method = "LogNormalize", scale.factor = 10000)
tcell.combined_subset <- FindVariableFeatures(tcell.combined_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(tcell.combined_subset)
tcell.combined_subset <- ScaleData(object = tcell.combined_subset, 
                                   do.par = TRUE, num.cores = 8, features =all.genes)

top.genes<-VariableFeatures(tcell.combined_subset)
ex_genes<-top.genes[!(top.genes %in% exclude_genes)]

tcell.combined_subset <- RunPCA(tcell.combined_subset, npcs = 30, verbose = FALSE,features=ex_genes)

ElbowPlot(tcell.combined_subset, ndims = 20, reduction = "pca")

tcell.combined_subset <- RunHarmony(tcell.combined_subset, group.by.vars = "patient_id")

tcell.combined_subset <- RunUMAP(tcell.combined_subset, reduction = "harmony", dims = 1:15)

tcell.combined_subset <- FindNeighbors(tcell.combined_subset, reduction = "harmony", dims = 1:15)
tcell.combined_subset <- FindClusters(tcell.combined_subset, resolution = 1)

DimPlot(tcell.combined_subset,reduction = "umap",group.by = "seurat_clusters",label=T)

#Cluster 9 and 2 needs to be reclustered to higher resolution, to recover MAIT cells and Early Tem and Tem cells, so repeat the below code for 9 and 2, then identify the cells from subsetted object
tcell.combined_subset$new_clusters<-as.character(tcell.combined_subset$seurat_clusters)

tcell.combined_cluster<-subset(x = tcell.combined_subset, subset = seurat_clusters %in% c( 9))

tcell.combined_cluster <- NormalizeData(tcell.combined_cluster, normalization.method = "LogNormalize", scale.factor = 10000)
tcell.combined_cluster <- FindVariableFeatures(tcell.combined_cluster, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(tcell.combined_cluster)
tcell.combined_cluster <- ScaleData(object = tcell.combined_cluster, 
                                    do.par = TRUE, num.cores = 8, features =all.genes)

top.genes<-VariableFeatures(tcell.combined_cluster)
ex_genes<-top.genes[!(top.genes %in% exclude_genes)]

tcell.combined_cluster <- RunPCA(tcell.combined_cluster, npcs = 30, verbose = FALSE,features=ex_genes)

ElbowPlot(tcell.combined_cluster, ndims = 20, reduction = "pca")

tcell.combined_cluster <- RunHarmony(tcell.combined_cluster, group.by.vars = "patient_id")

tcell.combined_cluster <- RunUMAP(tcell.combined_cluster, reduction = "harmony", dims = 1:10)

tcell.combined_cluster <- FindNeighbors(tcell.combined_cluster, reduction = "harmony", dims = 1:10)
tcell.combined_cluster <- FindClusters(tcell.combined_cluster, resolution = 0.09)

tcell.combined_subset_0<-subset(x = tcell.combined_cluster, subset = seurat_clusters %in% c( 0))

tcell.combined_subset_1<-subset(x = tcell.combined_cluster, subset = seurat_clusters %in% c( 1))

tcell.combined_subset@meta.data$new_clusters[WhichCells(tcell.combined_subset) %in% rownames(tcell.combined_subset_0@meta.data)]<-"9_0"

tcell.combined_subset@meta.data$new_clusters[WhichCells(tcell.combined_subset) %in% rownames(tcell.combined_subset_1@meta.data)]<-"9_1"

############## Subset 2 and part of 9 

tcell.combined_cluster<-subset(x = tcell.combined_subset, subset = new_clusters %in% c( "9_0","2"))

tcell.combined_cluster <- NormalizeData(tcell.combined_cluster, normalization.method = "LogNormalize", scale.factor = 10000)
tcell.combined_cluster <- FindVariableFeatures(tcell.combined_cluster, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(tcell.combined_cluster)
tcell.combined_cluster <- ScaleData(object = tcell.combined_cluster, 
                                    do.par = TRUE, num.cores = 8, features =all.genes)

all.genes <- rownames(tcell.combined_cluster)
tcell.combined_cluster <- ScaleData(object = tcell.combined_cluster, 
                                    do.par = TRUE, num.cores = 8, features =all.genes)

top.genes<-VariableFeatures(tcell.combined_cluster)
ex_genes<-top.genes[!(top.genes %in% exclude_genes)]

tcell.combined_cluster <- RunPCA(tcell.combined_cluster, npcs = 30, verbose = FALSE,features=ex_genes)

ElbowPlot(tcell.combined_cluster, ndims = 20, reduction = "pca")


tcell.combined_cluster <- RunHarmony(tcell.combined_cluster, group.by.vars = "patient_id")

tcell.combined_cluster <- RunUMAP(tcell.combined_cluster, reduction = "harmony", dims = 1:10)

tcell.combined_cluster <- FindNeighbors(tcell.combined_cluster, reduction = "harmony", dims = 1:10)
tcell.combined_cluster <- FindClusters(tcell.combined_cluster, resolution = 0.09)
DimPlot(tcell.combined_cluster, reduction = "umap",order = TRUE, label = TRUE) 

tcell.combined_subset_0<-subset(x = tcell.combined_cluster, subset = seurat_clusters %in% c( 0))

tcell.combined_subset_1<-subset(x = tcell.combined_cluster, subset = seurat_clusters %in% c( 1))

tcell.combined_subset@meta.data$new_clusters[WhichCells(tcell.combined_subset) %in% rownames(tcell.combined_subset_0@meta.data)]<-"9_2_0"

tcell.combined_subset@meta.data$new_clusters[WhichCells(tcell.combined_subset) %in% rownames(tcell.combined_subset_1@meta.data)]<-"9_2_1"

pdf("Results_Tcells/Tcells_new_clusters.pdf",width=15.5,height=13)
DimPlot(tcell.combined_subset, reduction = "umap",group.by="new_clusters",order = TRUE, label = TRUE,pt.size=1.2,label.size = 12) 
dev.off()
##########################


tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(0)]<-"LAG3 hi Tex"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(1)]<-"GZMK hi Tex"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(3)]<-"TOX hi Tex"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(4)]<-"Tn"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(5)]<-"Tn/Tm"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(6)]<-"TCF7+ stem-like Tex"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(7)]<-"MKI67 hi Tex/Tprol"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(8)]<-"MKI67+ Tex/Tprol"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(10)]<-"MKI67 hi Tex/Tprol"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(11)]<-"ISG+"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(12)]<-"Tem-NK like"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(13)]<-"MKI67 hi Tex/Tprol"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c(14)]<-"MKI67+ Tem-NK like/Tprol"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c("9_1")]<-"Tc17 MAIT"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c("9_2_0")]<-"Early Tem"
tcell.combined_subset@meta.data$annotated_clusters_labels_final[tcell.combined_subset@meta.data$new_clusters==c("9_2_1")]<-"Tem"



color_combo<-c("Tn" = "#85965f",
               "Tn/Tm" =	"#3b8e87"	,
               "Early Tem"=	"#EE9572"	,
               "Tem"=	"#b35f42"	,
               "Tem-NK like"= "#CD3700"	,	
               "TOX hi Tex"="#8e477f",
               "TCF7+ stem-like Tex"="#d379af",	
               "LAG3 hi Tex"=	"#c7abdd"	 ,
               "GZMK hi Tex"=	"pink",
               "MKI67+ Tex/Tprol" ="#a3dbdf"	,
               "MKI67 hi Tex/Tprol"="#6495ed"	,
               "MKI67+ Tem-NK like/Tprol"= "#104E8B",	
               "ISG+"=  "#dac565",
               "Tc17 MAIT" =	"#d8d5b2" )	

pdf("Results_Tcells/Tcells_Annotated_UAMP.pdf",width=16.5,height=13)
DimPlot(tcell.combined_subset,reduction="umap",group.by="annotated_clusters_labels_final",pt.size=3,cols=color_combo,order=c("ISG+","Tc17 MAIT" ,"Tn","Tn/Tm","Tem","Early Tem","Tem-NK like","LAG3 hi Tex","GZMK hi Tex","TOX hi Tex","TCF7+ stem-like Tex","MKI67+ Tem-NK like/Tprol","MKI67+ Tex/Tprol",
                                                                                                                             "MKI67 hi Tex/Tprol"))+ggtitle(" ")+guides(color = guide_legend(override.aes = list(size=6), ncol=1) )+theme(axis.text=element_text(size=25),text=element_text(size=25))+scale_color_manual("CD8+ T cell states:",breaks = c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol","MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol"),
                                                                                                                                                                                                                                                                                                                         values = c("#85965f",	"#3b8e87"	,"#EE9572"	,"#b35f42"	, "#CD3700"	,		"#d8d5b2", "#dac565","pink","#d379af",	"#8e477f","#c7abdd"	 ,"#a3dbdf"	,"#6495ed"	,"#104E8B"), 
                                                                                                                                                                                                                                                                                                                         labels = c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol",
                                                                                                                                                                                                                                                                                                                                    "MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol" ))
dev.off()

#object
saveRDS(tcell.combined_subset,"Results_Tcells/Tcells_Final.Rds")
