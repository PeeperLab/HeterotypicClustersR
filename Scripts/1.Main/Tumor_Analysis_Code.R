##################################################################
#Subsetting only for Tumor cells and Tumor cells characterization
##################################################################

# required input data in the working directory:
# # Combined_All_Cells.Rds from Allcells_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(Seurat) 
library(ggplot2)
library(dplyr)
library(data.table)
library(RSpectra)

set.seed(2000)

fig_path = "Results_tumor"
if(!dir.exists(fig_path)) dir.create(fig_path) 

#Read All cells integrated  object
mal.combined<-readRDS("Combined_All_Cells.Rds")

#Filter out the Tumor cells and T:Tumor cells and Singlets
mal.combined_db_subset<-subset(x = mal.combined, subset = seurat_clusters %in% c( 2,5,7,10,11,9,16,17,18,20,22,23))
mal.combined_db_subset<-subset(x = mal.combined_db_subset, subset = sample_id %in% c( "SG_Singlets","DB_Tumor_Tcell"))

DimPlot(mal.combined_db_subset,reduction = "umap")
###Seurat Integration

sc.list <- SplitObject(mal.combined_db_subset,split.by="patient_id")

sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})

sc.anchors <- FindIntegrationAnchors(object.list = sc.list, dims = 1:20)

pbmc.combined <- IntegrateData(anchorset = sc.anchors, dims = 1:20)


DefaultAssay(pbmc.combined) <- "integrated"
sc.combined <- ScaleData(pbmc.combined, verbose = FALSE)

sc.combined <- RunPCA(sc.combined, npcs = 50, verbose = FALSE)

ElbowPlot(sc.combined, ndims = 20, reduction = "pca")
DefaultAssay(sc.combined)<-"integrated"
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15, k.param = 30)
sc.combined <- FindClusters(sc.combined, resolution = 0.8)

DimPlot(sc.combined, reduction = "umap", group.by="seurat_clusters",label = TRUE,order=TRUE)

FeaturePlot(sc.combined,features="PTPRC",order=T)
FeaturePlot(sc.combined,features="CD8A",order=T)
FeaturePlot(sc.combined,features="CD4",order=T)
FeaturePlot(sc.combined,features="ITGAX",order=T)

#saveRDS(sc.combined,"Tumor_prefilter.Rds")

##Subset clusters based on QC and reintegrate on the previous object 
seurat_clusters = unique(sc.combined@meta.data$seurat_clusters)
seurat_clusters_include = as.character(seurat_clusters[!(seurat_clusters %in% c(6,10,13,14))])

sc.combined_subset<-subset(x = sc.combined, subset = seurat_clusters %in% seurat_clusters_include)

DimPlot(sc.combined_subset, reduction = "umap", group.by="seurat_clusters",label = TRUE,order=TRUE)

cell_ids1 <- rownames(sc.combined_subset@meta.data)
mal.combined_db_subset@meta.data <- mal.combined_db_subset@meta.data %>% mutate(mal.combined_db_subset@meta.data$non_tumor <- ifelse((WhichCells(mal.combined_db_subset) %in% cell_ids1),"Present",NA))
mal.combined_db_subset<-AddMetaData(object = mal.combined_db_subset, metadata = mal.combined_db_subset@meta.data['... <- NULL'], col.name = 'tumor')
mal.combined_db_subset_trim<-subset(x = mal.combined_db_subset, subset = tumor %in% c( "Present"))
mal.combined_db_subset_trim


######################################################################################################

options(future.globals.maxSize = 8000 * 1024^2)
mal.combined_db_subset_trim <- SCTransform(mal.combined_db_subset_trim,verbose = FALSE, vars.to.regress = c("percent.mito", 'nFeature_RNA')) 

sc.list <- SplitObject(mal.combined_db_subset_trim,split.by="patient_id")
sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
  x
})

mito_genes <- grep("^MT-", rownames(mal.combined_db_subset_trim), value = TRUE)
top_features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 4000)
genes.to.keep <- setdiff(top_features, mito_genes)



sc.list <- PrepSCTIntegration( sc.list, anchor.features = genes.to.keep)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = genes.to.keep,dims = 1:20,normalization.method = "SCT",reduction="cca")
pbmc.combined <- IntegrateData(anchorset = sc.anchors,dims = 1:20, new.assay.name = "CCA",normalization.method = "SCT")
DefaultAssay(pbmc.combined) <- "CCA"
pbmc.combined <- RunPCA(pbmc.combined, npcs = 50,features=genes.to.keep,  verbose = FALSE)
ElbowPlot(sc.combined, ndims = 20, reduction = "pca")
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca",dims = 1:15)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:15, k.param = 30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 1)

DefaultAssay(pbmc.combined)<-"RNA"
DimPlot(pbmc.combined, reduction = "umap", group.by="seurat_clusters",label = TRUE,order=TRUE)
FeaturePlot(pbmc.combined,features="NGFR",order=T)

pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("16","1")]<-"Immune Response"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("11","14")]<-"Mitotic"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("4","5")]<-"Melanocytic"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("7","15")]<-"Stress(p53 response)" 
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("10")]<-"Low Gene Count Cells"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("13","6","8","9","0")]<-"Transitory Melanocytic"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("2")]<-"Patient Specific"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("3")]<-"Stress(hypoxia response)"
pbmc.combined@meta.data$annotated_clusters[pbmc.combined@meta.data$seurat_clusters %in% c("12")]<-"Neural Crest Like"


color_combo2<-c("Immune Response"="#99CC66",
                "Mitotic"="#006633",
                "Melanocytic"="#990000",
                "Stress(hypoxia response)"="#dac565",
                "Transitory Melanocytic"="#FF3300",
                "Stress(p53 response)"="#0099CC",
                "Low Gene Count Cells"="#CCCCCC",
                "Neural Crest Like"="#003366",
                "Patient Specific"="#FF9933")




pdf("Results_tumor/Tumor_Annotated_UMAP.pdf",width=13,height=10)
DimPlot(pbmc.combined, reduction = "umap",group.by = "annotated_clusters", label = F,order=c( "Immune Response", "Stress(p53 response)" ,
                                                                                              "Mitotic", "Low Gene Count Cells", "Melanocytic",
                                                                                              "Transitory Melanocytic" , "Neural Crest Like",
                                                                                              "Stress(hypoxia response)","Patient Specific" ),
        pt.size = 1.8,label.size = 25,cols=color_combo2)+ggtitle(" ")+
  guides(color = guide_legend(override.aes = list(size=6), ncol=1) )+
  theme(text=element_text(size=25))+
  scale_color_manual("Tumor cell phenotypes:",breaks = c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic",
                                                         "Patient Specific", "Stress(hypoxia response)", "Stress(p53 response)",
                                                         "Neural Crest Like","Low Gene Count Cells"),
                     values = c("#99CC66","#006633","#990000","#FF3300","#FF9933","#dac565","#0099CC","#003366","#CCCCCC"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               labels = c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific", "Stress(hypoxia response)", "Stress(p53 response)","Neural Crest Like","Low Gene Count Cells"))
dev.off()

saveRDS(pbmc.combined,"Results_tumor/Tumor_Annotated.Rds")

