
###################################################################
#scRNA seq Analysis of clusters data from 5 patients
##################################################################

# required input data in the working directory:
# # subfolders in data/sequence_data/
# # a results folder will be created in this directory

#R version 4.3.3


set.seed(2000)
#setwd("/YOUR/PATH/")
library(Seurat)# Version 4.4.0 
library(ggplot2)
library(dplyr)
library(harmony)# Version 1.2.1
library(RSpectra)
library(Nebulosa)
library(Matrix)# Version 1.6-5

fig_path = "Results_Main"
if(!dir.exists(fig_path)) dir.create(fig_path)
#Patient 1 from the samples named as Patient 1 in code,
#Patient 2 from the samples named as Patient 2 in code,
#Patient 8 from the samples named as Patient 3 in code,
#Patient 9 from the samples named as Patient 4 in code,
#Patient 15 from the samples named as Patient 5 in code

##Loading cellranger data for 5 patients from folders

#Patient 1 folders S1_TTD,S2_TDD,S3_STCD
#Patient 2 folders S1_DB_2, S3_SG_2
#Patient 3 folders DB_3, SG_3
#Patient 4 folders DB1_4, DB2_4, SG_4
#Patient 5 folders S5_TT, S5_AT, S5_SG

#Patient 1 
ids<-c("S1_TTD","S2_TDD","S3_STCD")
for (i in ids){
  
  pbmc.data <- Read10X(file.path(i,"/"))
  
  pbmc <- CreateSeuratObject(counts = pbmc.data,  project = "scRNA", min.cells = 10, min.features = 200,
                             names.field = 2, names.delim = "\\-")
  assign(paste0("pbmc_", i), pbmc)
}

#Patient 2 contains the hashtagged clusters samples and hashtaged are present in "Antibody Capture" segment of the data 
pbmc.data <- Read10X(file.path("S1_DB_2/sample_filtered_feature_bc_matrix/"))

pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`,  project = "scRNA", min.cells = 10, min.features = 200,
                           names.field = 2, names.delim = "\\-")

pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.data[["Antibody Capture"]][, colnames(x = pbmc)])

pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")

DefaultAssay(pbmc) <- 'ADT'

#Histogram of Hashtag expression
h1<-VlnPlot(pbmc,features="HTO1")
h1<-h1$data
hist(h1$HTO1,breaks=1000)

h2<-VlnPlot(pbmc,features="HTO2")
h2<-h2$data
hist(h2$HTO2,breaks=1000)
#HTO1 and HTO2 cutoff chosen based on local minimum
pbmc.ht01<-subset(x = pbmc, subset = HTO1 > 0.5 )
pbmc.ht02<-subset(x = pbmc, subset = HTO2 > 1 )

DefaultAssay(pbmc.ht01) <- 'RNA'
DefaultAssay(pbmc.ht02) <- 'RNA'
cell_ids1 <- rownames(pbmc.ht01@meta.data)
cell_ids2 <- rownames(pbmc.ht02@meta.data)

# Find overlapping cell IDs
overlapping_cell_ids <- intersect(cell_ids1, cell_ids2)

length(overlapping_cell_ids)

pbmc.ht01a <- pbmc.ht01[, !colnames(pbmc.ht01) %in% overlapping_cell_ids]
pbmc.ht02a <- pbmc.ht02[, !colnames(pbmc.ht02) %in% overlapping_cell_ids]

pbmc.ht01a<-AddMetaData(pbmc.ht01a, "DB_Tumor_Tcell", col.name = "sample_id")
pbmc.ht02a<-AddMetaData(pbmc.ht02a, "DB_APC_Tcell", col.name = "sample_id")

clone.S1_DB_2<- read.csv("S1_DB_2/filtered_contig_annotations.csv") 

clone.S1_DB_2_TT<-clone.S1_DB_2[clone.S1_DB_2$barcode %in% colnames(pbmc.ht01a),]
clone.S1_DB_2_AT<-clone.S1_DB_2[clone.S1_DB_2$barcode %in% colnames(pbmc.ht02a),]

write.csv(clone.S1_DB_2_TT,"S1_DB_2/clone.S1_DB_2_TT.csv")
write.csv(clone.S1_DB_2_AT,"S1_DB_2/clone.S1_DB_2_AT.csv")

# Create Seurat object for singlets patient 2

pbmc.data <- Read10X(data.dir = "S3_SG_2/filtered_feature_bc_matrix")

pbmc.s3_2 <- CreateSeuratObject(counts = pbmc.data, project = "scRNA", assay = "RNA", 
                                min.cells = 10, min.features = 200)

pbmc.s3_2<-AddMetaData(pbmc.s3_2, "SG_Singlets", col.name = "sample_id")

#Patient 3 contains the hashtagged clusters samples and hashtaged are present in "Antibody Capture" segment of the data 
pbmc.data3 <- Read10X(file.path("DB_3/sample_filtered_feature_bc_matrix/"))

pbmc3 <- CreateSeuratObject(counts = pbmc.data3$`Gene Expression`,  project = "scRNA", min.cells = 10, min.features = 200,
                            names.field = 2, names.delim = "\\-")
pbmc3[["ADT"]] <- CreateAssayObject(counts = pbmc.data3[["Antibody Capture"]][, colnames(x = pbmc3)])

pbmc3 <- NormalizeData(pbmc3, assay = "ADT", normalization.method = "CLR")

DefaultAssay(pbmc3) <- 'ADT'

#Histogram of Hashtag expression
h1<-VlnPlot(pbmc3,features="HTO1")
h1<-h1$data
hist(h1$HTO1,breaks=1000)

h2<-VlnPlot(pbmc3,features="HTO2")
h2<-h2$data
hist(h2$HTO2,breaks=1000)


pbmc3.ht01<-subset(x = pbmc3, subset = HTO1 > 0.4 )
pbmc3.ht02<-subset(x = pbmc3, subset = HTO2 > 1 )

DefaultAssay(pbmc3.ht01) <- 'RNA'
DefaultAssay(pbmc3.ht02) <- 'RNA'
cell_ids1 <- rownames(pbmc3.ht01@meta.data)
cell_ids2 <- rownames(pbmc3.ht02@meta.data)


# Find overlapping cell IDs
overlapping_cell_ids <- intersect(cell_ids1, cell_ids2)


length(overlapping_cell_ids)

pbmc3.ht01a <- pbmc3.ht01[, !colnames(pbmc3.ht01) %in% overlapping_cell_ids]
pbmc3.ht02a <- pbmc3.ht02[, !colnames(pbmc3.ht02) %in% overlapping_cell_ids]

pbmc3.ht01a<-AddMetaData(pbmc3.ht01a, "DB_Tumor_Tcell", col.name = "sample_id")

pbmc3.ht02a<-AddMetaData(pbmc3.ht02a, "DB_APC_Tcell", col.name = "sample_id")


clone.S1_DB_2<- read.csv("DB_3/filtered_contig_annotations.csv") 

clone.S1_DB_2_TT<-clone.S1_DB_2[clone.S1_DB_2$barcode %in% colnames(pbmc3.ht01a),]
clone.S1_DB_2_AT<-clone.S1_DB_2[clone.S1_DB_2$barcode %in% colnames(pbmc3.ht02a),]

write.csv(clone.S1_DB_2_TT,"DB_3/clone.S1_DB_3_TT.csv")
write.csv(clone.S1_DB_2_AT,"DB_3/clone.S1_DB_3_AT.csv")

# Create Seurat object for singlets patient 3
pbmc3.data <- Read10X(data.dir = "SG_3/filtered_feature_bc_matrix")

pbmc3.s3_2 <- CreateSeuratObject(counts = pbmc3.data, project = "scRNA", assay = "RNA", 
                                 min.cells = 10, min.features = 200)

pbmc3.s3_2<-AddMetaData(pbmc3.s3_2, "SG_Singlets", col.name = "sample_id")

#Patient 4
ids4<-c("DB1_4","DB2_4","SG_4")

for (i in ids4){
  
  pbmc.data <- Read10X(file.path(i,"/filtered_feature_bc_matrix"))
  
  pbmc <- CreateSeuratObject(counts = pbmc.data,  project = "scRNA", min.cells = 10, min.features = 200,
                             names.field = 2, names.delim = "\\-")
  assign(paste0("pbmc_", i), pbmc)
}

#Patient 5
ids5<-c("S5_TT","S5_AT","S5_SG")

for (i in ids5){
  
  pbmc.data <- Read10X(file.path(i,"/filtered_feature_bc_matrix"))
  
  pbmc <- CreateSeuratObject(counts = pbmc.data,  project = "scRNA", min.cells = 10, min.features = 200,
                             names.field = 2, names.delim = "\\-")
  assign(paste0("pbmc_", i), pbmc)
}

##Estimate the Mitochondrial percentages

for(i in ids){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^MT-",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.mito, col.name = "percent.mito"))}

mt.genes <- rownames(pbmc.ht01a)[grep("^MT-",rownames(pbmc.ht01a))]
C<-GetAssayData(object = pbmc.ht01a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.ht01a <- AddMetaData(pbmc.ht01a, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc.ht02a)[grep("^MT-",rownames(pbmc.ht02a))]
C<-GetAssayData(object = pbmc.ht02a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.ht02a <- AddMetaData(pbmc.ht02a, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc.s3_2)[grep("^MT-",rownames(pbmc.s3_2))]
C<-GetAssayData(object = pbmc.s3_2, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.s3_2 <- AddMetaData(pbmc.s3_2, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc3.ht01a)[grep("^MT-",rownames(pbmc3.ht01a))]
C<-GetAssayData(object = pbmc3.ht01a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.ht01a <- AddMetaData(pbmc3.ht01a, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc3.ht02a)[grep("^MT-",rownames(pbmc3.ht02a))]
C<-GetAssayData(object = pbmc3.ht02a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.ht02a <- AddMetaData(pbmc3.ht02a, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc3.s3_2)[grep("^MT-",rownames(pbmc3.s3_2))]
C<-GetAssayData(object = pbmc3.s3_2, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.s3_2 <- AddMetaData(pbmc3.s3_2, percent.mito, col.name = "percent.mito")


for(i in ids4){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^MT-",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.mito, col.name = "percent.mito"))}

for(i in ids5){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^MT-",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.mito, col.name = "percent.mito"))}

#Apply the general cutoff of 200< nFeature_RNA <8000 and mitochondrial percentage of below 15%
for(i in ids){
  
  
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA >= 200)))
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA <= 8000)))
  assign(paste0("pbmc_",i), subset(x = get(paste0("pbmc_",i)), 
                                   subset= (percent.mito < 15)))}



pbmc.ht01a<-subset(x = pbmc.ht01a, 
                   subset= (nFeature_RNA >= 200))
pbmc.ht01a<-subset(x = pbmc.ht01a, 
                   subset= (nFeature_RNA <= 8000))
pbmc.ht01a<-subset(x = pbmc.ht01a, 
                   subset= (percent.mito < 15))

pbmc.ht02a<-subset(x = pbmc.ht02a, 
                   subset= (nFeature_RNA >= 200))
pbmc.ht02a<-subset(x = pbmc.ht02a, 
                   subset= (nFeature_RNA <= 8000))
pbmc.ht02a<-subset(x = pbmc.ht02a, 
                   subset= (percent.mito < 15))

pbmc.s3_2<-subset(x = pbmc.s3_2, 
                  subset= (nFeature_RNA >= 200))
pbmc.s3_2<-subset(x = pbmc.s3_2, 
                  subset= (nFeature_RNA <= 8000))
pbmc.s3_2<-subset(x = pbmc.s3_2, 
                  subset= (percent.mito < 15))


pbmc3.ht01a<-subset(x = pbmc3.ht01a, 
                    subset= (nFeature_RNA >= 200))
pbmc3.ht01a<-subset(x = pbmc3.ht01a, 
                    subset= (nFeature_RNA <= 8000))
pbmc3.ht01a<-subset(x = pbmc3.ht01a, 
                    subset= (percent.mito < 15))

pbmc3.ht02a<-subset(x = pbmc3.ht02a, 
                    subset= (nFeature_RNA >= 200))
pbmc3.ht02a<-subset(x = pbmc3.ht02a, 
                    subset= (nFeature_RNA <= 8000))
pbmc3.ht02a<-subset(x = pbmc3.ht02a, 
                    subset= (percent.mito < 15))

pbmc3.s3_2<-subset(x = pbmc3.s3_2, 
                   subset= (nFeature_RNA >= 200))
pbmc3.s3_2<-subset(x = pbmc3.s3_2, 
                   subset= (nFeature_RNA <= 8000))
pbmc3.s3_2<-subset(x = pbmc3.s3_2, 
                   subset= (percent.mito < 15))

for(i in ids4){
  
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA >= 200)))
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA <= 8000)))
  assign(paste0("pbmc_",i), subset(x = get(paste0("pbmc_",i)), 
                                   subset= (percent.mito < 15)))}

for(i in ids5){
  
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA >= 200)))
  assign(paste0("pbmc_",i),subset(x = get(paste0("pbmc_",i)), 
                                  subset= (nFeature_RNA <= 8000)))
  assign(paste0("pbmc_",i), subset(x = get(paste0("pbmc_",i)), 
                                   subset= (percent.mito < 15)))}

for(i in ids){
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), i, col.name = "sample_id"))
}

for(i in ids){
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), i, col.name = "sample_id"))
}

#Add Sample specific metadata informations

ids<-c("S1_TTD","S2_TDD","S3_STCD")

pbmc_S1_TTD<-AddMetaData(pbmc_S1_TTD, "DB_Tumor_Tcell", col.name = "sample_id")
pbmc_S2_TDD<-AddMetaData(pbmc_S2_TDD, "DB_APC_Tcell", col.name = "sample_id")
pbmc_S3_STCD<-AddMetaData(pbmc_S3_STCD, "SG_Singlets", col.name = "sample_id")


ids4<-c("DB1_4","DB2_4","SG_4")

pbmc_DB1_4<-AddMetaData(pbmc_DB1_4, "DB_Tumor_Tcell", col.name = "sample_id")
pbmc_DB2_4<-AddMetaData(pbmc_DB2_4, "DB_APC_Tcell", col.name = "sample_id")
pbmc_SG_4<-AddMetaData(pbmc_SG_4, "SG_Singlets", col.name = "sample_id")

pbmc_S5_TT<-AddMetaData(pbmc_S5_TT, "DB_Tumor_Tcell", col.name = "sample_id")
pbmc_S5_AT<-AddMetaData(pbmc_S5_AT, "DB_APC_Tcell", col.name = "sample_id")
pbmc_S5_SG<-AddMetaData(pbmc_S5_SG, "SG_Singlets", col.name = "sample_id")

for(i in ids){
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), "P1", col.name = "patient_id"))
}

pbmc.ht01a<-AddMetaData(pbmc.ht01a, "P2", col.name = "patient_id")
pbmc.ht02a<-AddMetaData(pbmc.ht02a, "P2", col.name = "patient_id")
pbmc.s3_2<-AddMetaData(pbmc.s3_2, "P2", col.name = "patient_id")

pbmc3.ht01a<-AddMetaData(pbmc3.ht01a, "P3", col.name = "patient_id")
pbmc3.ht02a<-AddMetaData(pbmc3.ht02a, "P3", col.name = "patient_id")
pbmc3.s3_2<-AddMetaData(pbmc3.s3_2, "P3", col.name = "patient_id")

for(i in ids4){
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), "P4", col.name = "patient_id"))
}

for(i in ids5){
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), "P5", col.name = "patient_id"))
}

pbmc_S1_TTD<-AddMetaData(pbmc_S1_TTD, "P1_DB_Tumor_Tcell", col.name = "combined_id")
pbmc_S2_TDD<-AddMetaData(pbmc_S2_TDD, "P1_DB_APC_Tcell", col.name = "combined_id")
pbmc_S3_STCD<-AddMetaData(pbmc_S3_STCD, "P1_SG_Singlets", col.name = "combined_id")

pbmc.ht01a<-AddMetaData(pbmc.ht01a, "P2_DB_Tumor_Tcell", col.name = "combined_id")
pbmc.ht02a<-AddMetaData(pbmc.ht02a, "P2_DB_APC_Tcell", col.name = "combined_id")
pbmc.s3_2<-AddMetaData(pbmc.s3_2, "P2_SG_Singlets", col.name = "combined_id")

pbmc3.ht01a<-AddMetaData(pbmc3.ht01a, "P3_DB_Tumor_Tcell", col.name = "combined_id")
pbmc3.ht02a<-AddMetaData(pbmc3.ht02a, "P3_DB_APC_Tcell", col.name = "combined_id")
pbmc3.s3_2<-AddMetaData(pbmc3.s3_2, "P3_SG_Singlets", col.name = "combined_id")

pbmc_DB1_4<-AddMetaData(pbmc_DB1_4, "P4_DB_Tumor_Tcell", col.name = "combined_id")
pbmc_DB2_4<-AddMetaData(pbmc_DB2_4, "P4_DB_APC_Tcell", col.name = "combined_id")
pbmc_SG_4<-AddMetaData(pbmc_SG_4, "P4_SG_Singlets", col.name = "combined_id")

pbmc_S5_TT<-AddMetaData(pbmc_S5_TT, "P5_DB_Tumor_Tcell", col.name = "combined_id")
pbmc_S5_AT<-AddMetaData(pbmc_S5_AT, "P5_DB_APC_Tcell", col.name = "combined_id")
pbmc_S5_SG<-AddMetaData(pbmc_S5_SG, "P5_SG_Singlets", col.name = "combined_id")

pbmc_S1_TTD@meta.data$barcode<-rownames(pbmc_S1_TTD@meta.data)
pbmc_S2_TDD@meta.data$barcode<-rownames(pbmc_S2_TDD@meta.data)
pbmc_S3_STCD@meta.data$barcode<-rownames(pbmc_S3_STCD@meta.data)

pbmc.ht01a@meta.data$barcode<-rownames(pbmc.ht01a@meta.data)
pbmc.ht02a@meta.data$barcode<-rownames(pbmc.ht02a@meta.data)
pbmc.s3_2@meta.data$barcode<-rownames(pbmc.s3_2@meta.data)

pbmc3.ht01a@meta.data$barcode<-rownames(pbmc3.ht01a@meta.data)
pbmc3.ht02a@meta.data$barcode<-rownames(pbmc3.ht02a@meta.data)
pbmc3.s3_2@meta.data$barcode<-rownames(pbmc3.s3_2@meta.data)

pbmc_DB1_4@meta.data$barcode<-rownames(pbmc_DB1_4@meta.data)
pbmc_DB2_4@meta.data$barcode<-rownames(pbmc_DB2_4@meta.data)
pbmc_SG_4@meta.data$barcode<-rownames(pbmc_SG_4@meta.data)

pbmc_S5_TT@meta.data$barcode<-rownames(pbmc_S5_TT@meta.data)
pbmc_S5_AT@meta.data$barcode<-rownames(pbmc_S5_AT@meta.data)
pbmc_S5_SG@meta.data$barcode<-rownames(pbmc_S5_SG@meta.data)

#Append the cell id to remove the duplicate cellid names

pbmc_S1_TTD <- RenameCells(object = pbmc_S1_TTD, add.cell.id = "DB1_1_")
pbmc_S2_TDD <- RenameCells(object = pbmc_S2_TDD, add.cell.id = "DB2_1_")
pbmc_S3_STCD <- RenameCells(object = pbmc_S3_STCD, add.cell.id = "SG_1_")

pbmc.ht01a <- RenameCells(object = pbmc.ht01a, add.cell.id = "DB1_2_")
pbmc.ht02a <- RenameCells(object = pbmc.ht02a, add.cell.id = "DB2_2_")
pbmc.s3_2 <- RenameCells(object = pbmc.s3_2, add.cell.id = "SG_2_")

pbmc3.ht01a <- RenameCells(object = pbmc3.ht01a, add.cell.id = "DB1_3_")
pbmc3.ht02a <- RenameCells(object = pbmc3.ht02a, add.cell.id = "DB2_3_")
pbmc3.s3_2 <- RenameCells(object = pbmc3.s3_2, add.cell.id = "SG_3_")

pbmc_DB1_4 <- RenameCells(object = pbmc_DB1_4, add.cell.id = "DB1_4_")
pbmc_DB2_4 <- RenameCells(object = pbmc_DB2_4, add.cell.id = "DB2_4_")
pbmc_SG_4 <- RenameCells(object = pbmc_SG_4, add.cell.id = "SG_4_")

pbmc_S5_TT <- RenameCells(object = pbmc_S5_TT, add.cell.id = "DB1_5_")
pbmc_S5_AT <- RenameCells(object = pbmc_S5_AT, add.cell.id = "DB2_5_")
pbmc_S5_SG <- RenameCells(object = pbmc_S5_SG, add.cell.id = "SG_5_")



pbmc.combined<-merge(pbmc_S1_TTD ,
                     y=c(pbmc_S2_TDD,pbmc_S3_STCD,pbmc.ht01a,pbmc.ht02a,pbmc.s3_2,pbmc3.ht01a,pbmc3.ht02a,pbmc3.s3_2,pbmc_DB1_4,pbmc_DB2_4,pbmc_SG_4,pbmc_S5_TT,pbmc_S5_AT,pbmc_S5_SG),project = "scRNA_merge")


pbmc.combined

#Using harmony to perfom integration on all the cells
pbmc.combined <- NormalizeData(pbmc.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
pbmc.combined <- RunHarmony(pbmc.combined, group.by.vars = "patient_id")
ElbowPlot(pbmc.combined, ndims = 50, reduction = "pca")
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "harmony", dims = 1:20)


pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "harmony", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.8)


pdf("Results_Main/Allcells_seurat_clusters.pdf")
options(repr.plot.width=17, repr.plot.height=15)
DimPlot(pbmc.combined, reduction = "umap",order = TRUE, label = TRUE,label.size = 7,pt.size = 1.2)
dev.off()


pbmc.combined<-SetIdent(pbmc.combined, value = pbmc.combined@meta.data$seurat_clusters)

#Now based on the expression of CD8A, ITGAX and PMEL seprate out possible Tcells, APCs and Tumor cells
pbmc.combined@meta.data$source_id[pbmc.combined@meta.data$seurat_clusters %in% c(2,5,7,10,11,9,16,17,18,20,22,23,28)]<-"Tumor cells"
pbmc.combined@meta.data$source_id[pbmc.combined@meta.data$seurat_clusters %in% c(0,1,3,6,12,14,15)]<-"T cells"  
pbmc.combined@meta.data$source_id[pbmc.combined@meta.data$seurat_clusters %in% c( 4,8,13,19,21,24,25,26,27)]<-"APCs"


color_combo_sample<-c("Tumor cells" = "#990033",
                      "T cells" =	"#509abd",
                      "APCs"="#dac565" )

p <- DimPlot(
  pbmc.combined,
  reduction = "umap",
  group.by = "source_id",
  pt.size = 2.5,
  cols = color_combo_sample, # Ensure color_combo is defined properly
  order = c(
    "T cells","APCs","Tumor cells"
  )
) + labs(title = NULL)+
  guides(
    color = guide_legend(override.aes = list(size = 6), ncol = 1)
  ) + 
  theme(
    text = element_text(size = 20)
  ) + 
  scale_color_manual(" ",
                     values = c(
                       "#509abd", "#990033", "#dac565" 
                     ),
                     breaks = c(
                       "T cells", "Tumor cells" ,"APCs"
                     ),
                     labels = c(
                       "T cells", "Tumor cells" ,"APCs"
                     )
  )
pdf("Results_Main/Allcells_clusters.pdf", width = 13, height = 12)
p
dev.off()



pdf(paste0("Results_Main/","Allcells_CD19_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "CD19",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("CD19")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_CD8A_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "CD8A",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("CD8A")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_ITGAX_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "ITGAX",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("ITGAX")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_PMEL_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "PMEL",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("PMEL")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_MCAM_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "MCAM",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("MCAM")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_CD19_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "CD19",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("CD19")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()
pdf(paste0("Results_Main/","Allcells_CD3D_Markers_Nebulosa.pdf"),width=12,height=10)
plot_density(pbmc.combined, "CD3D",size=2,adjust = 20,method = "wkde",slot="counts") +ggtitle("CD3D")+theme(plot.title = element_text(hjust = 0.5),title=element_text(size=35),axis.text=element_text(size=12),axis.title=element_text(size=14))
dev.off()


saveRDS(pbmc.combined,"Results_Main/Combined_All_Cells.Rds")

