


###################################################################

#scRNA seq Analysis of clusters data from 5 patients

##################################################################


#R version 4.3.3

set.seed(2000)

library(Seurat)# Version 4.4.0 
library(ggplot2)
library(dplyr)
library(harmony)# Version 1.2.1
library(RSpectra)
library(Nebulosa)
library(Matrix)# Version 1.6-5

fig_path = "Results_Main_CD39"
if(!dir.exists(fig_path)) dir.create(fig_path)

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

##CD39#################################################################################################################################################################################################
#Patient 2 contains the hashtagged clusters samples and hashtaged are present in "Antibody Capture" segment of the data 
pbmc.data <- Read10X(file.path("CD39/sample_filtered_feature_bc_matrix/"))

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
pbmc.ht02<-subset(x = pbmc, subset = HTO2 > 1.2 )

DefaultAssay(pbmc.ht01) <- 'RNA'
DefaultAssay(pbmc.ht02) <- 'RNA'
cell_ids1 <- rownames(pbmc.ht01@meta.data)
cell_ids2 <- rownames(pbmc.ht02@meta.data)

# Find overlapping cell IDs
overlapping_cell_ids <- intersect(cell_ids1, cell_ids2)


length(overlapping_cell_ids)

pbmc_cd39.ht01a <- pbmc.ht01[, !colnames(pbmc.ht01) %in% overlapping_cell_ids]
pbmc_cd39.ht02a <- pbmc.ht02[, !colnames(pbmc.ht02) %in% overlapping_cell_ids]

pbmc_cd39.ht01a<-AddMetaData(pbmc_cd39.ht01a, "P3_CD39", col.name = "sample_id")
pbmc_cd39.ht02a<-AddMetaData(pbmc_cd39.ht02a, "P5_CD39", col.name = "sample_id")

pbmc_cd39.ht01a<-AddMetaData(pbmc_cd39.ht01a, "P3", col.name = "patient_id")
pbmc_cd39.ht02a<-AddMetaData(pbmc_cd39.ht02a, "P5", col.name = "patient_id")

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

mt.genes <- rownames(pbmc_cd39.ht01a)[grep("^MT-",rownames(pbmc_cd39.ht01a))]
C<-GetAssayData(object = pbmc_cd39.ht01a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc_cd39.ht01a <- AddMetaData(pbmc_cd39.ht01a, percent.mito, col.name = "percent.mito")

mt.genes <- rownames(pbmc_cd39.ht02a)[grep("^MT-",rownames(pbmc_cd39.ht02a))]
C<-GetAssayData(object = pbmc_cd39.ht02a, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc_cd39.ht02a <- AddMetaData(pbmc_cd39.ht02a, percent.mito, col.name = "percent.mito")



################Ribo

for(i in ids){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^RP[SL]",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.ribo, col.name = "percent.ribo"))}

mt.genes <- rownames(pbmc.ht01a)[grep("^RP[SL]",rownames(pbmc.ht01a))]
C<-GetAssayData(object = pbmc.ht01a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.ht01a <- AddMetaData(pbmc.ht01a, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc.ht02a)[grep("^RP[SL]",rownames(pbmc.ht02a))]
C<-GetAssayData(object = pbmc.ht02a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.ht02a <- AddMetaData(pbmc.ht02a, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc.s3_2)[grep("^RP[SL]",rownames(pbmc.s3_2))]
C<-GetAssayData(object = pbmc.s3_2, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc.s3_2 <- AddMetaData(pbmc.s3_2, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc3.ht01a)[grep("^RP[SL]",rownames(pbmc3.ht01a))]
C<-GetAssayData(object = pbmc3.ht01a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.ht01a <- AddMetaData(pbmc3.ht01a, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc3.ht02a)[grep("^RP[SL]",rownames(pbmc3.ht02a))]
C<-GetAssayData(object = pbmc3.ht02a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.ht02a <- AddMetaData(pbmc3.ht02a, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc3.s3_2)[grep("^RP[SL]",rownames(pbmc3.s3_2))]
C<-GetAssayData(object = pbmc3.s3_2, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc3.s3_2 <- AddMetaData(pbmc3.s3_2, percent.ribo, col.name = "percent.ribo")


for(i in ids4){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^RP[SL]",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.ribo, col.name = "percent.ribo"))}

for(i in ids5){
  mt.genes <- rownames(get(paste0("pbmc_",i)))[grep("^RP[SL]",rownames(get(paste0("pbmc_",i))))]
  C<-GetAssayData(object = get(paste0("pbmc_",i)), slot = "counts")
  percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  assign(paste0("pbmc_",i), AddMetaData(get(paste0("pbmc_",i)), percent.ribo, col.name = "percent.ribo"))}

mt.genes <- rownames(pbmc_cd39.ht01a)[grep("^RP[SL]",rownames(pbmc_cd39.ht01a))]
C<-GetAssayData(object = pbmc_cd39.ht01a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc_cd39.ht01a <- AddMetaData(pbmc_cd39.ht01a, percent.ribo, col.name = "percent.ribo")

mt.genes <- rownames(pbmc_cd39.ht02a)[grep("^RP[SL]",rownames(pbmc_cd39.ht02a))]
C<-GetAssayData(object = pbmc_cd39.ht02a, slot = "counts")

percent.ribo <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
pbmc_cd39.ht02a <- AddMetaData(pbmc_cd39.ht02a, percent.ribo, col.name = "percent.ribo")
################

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

pbmc_cd39.ht01a<-subset(x = pbmc_cd39.ht01a, 
                   subset= (nFeature_RNA >= 200))
pbmc_cd39.ht01a<-subset(x = pbmc_cd39.ht01a, 
                   subset= (nFeature_RNA <= 8000))
pbmc_cd39.ht01a<-subset(x = pbmc_cd39.ht01a, 
                   subset= (percent.mito < 15))

pbmc_cd39.ht02a<-subset(x = pbmc_cd39.ht02a, 
                   subset= (nFeature_RNA >= 200))
pbmc_cd39.ht02a<-subset(x = pbmc_cd39.ht02a, 
                   subset= (nFeature_RNA <= 8000))
pbmc_cd39.ht02a<-subset(x = pbmc_cd39.ht02a, 
                   subset= (percent.mito < 15))

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

pbmc_cd39.ht01a<-AddMetaData(pbmc_cd39.ht01a, "P3_CD39", col.name = "combined_id")
pbmc_cd39.ht02a<-AddMetaData(pbmc_cd39.ht02a, "P5_CD39", col.name = "combined_id")


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

pbmc_cd39.ht01a@meta.data$barcode<-rownames(pbmc_cd39.ht01a@meta.data)
pbmc_cd39.ht02a@meta.data$barcode<-rownames(pbmc_cd39.ht02a@meta.data)

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

pbmc_cd39.ht01a <- RenameCells(object = pbmc_cd39.ht01a, add.cell.id = "CD39_3_")
pbmc_cd39.ht02a <- RenameCells(object = pbmc_cd39.ht02a, add.cell.id = "CD39_5_")


pbmc.combined<-merge(pbmc_S1_TTD ,
                     y=c(pbmc_S2_TDD,pbmc_S3_STCD,pbmc.ht01a,pbmc.ht02a,pbmc.s3_2,pbmc3.ht01a,pbmc3.ht02a,pbmc3.s3_2,pbmc_DB1_4,pbmc_DB2_4,pbmc_SG_4,pbmc_S5_TT,pbmc_S5_AT,pbmc_S5_SG,pbmc_cd39.ht01a,pbmc_cd39.ht02a),project = "scRNA_merge")


pbmc.combined
VlnPlot(pbmc.combined,features="ENTPD1",group.by="combined_id")
#Using harmony to perfom integration on all the cells
pbmc.combined <- NormalizeData(pbmc.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
pbmc.combined <- RunHarmony(pbmc.combined, group.by.vars = "patient_id")
ElbowPlot(pbmc.combined, ndims = 50, reduction = "pca")
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "harmony", dims = 1:20)


pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "harmony", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.8)


pdf("Results_Main_CD39/Allcells_seurat_clusters.pdf")
options(repr.plot.width=17, repr.plot.height=15)
DimPlot(pbmc.combined, reduction = "umap",order = TRUE, label = TRUE,label.size = 7,pt.size = 1.2)
dev.off()

DimPlot(pbmc.combined, reduction = "umap",group.by="sample_id",order = TRUE, label = TRUE,label.size = 7,pt.size = 1.2)

pbmc.combined<-SetIdent(pbmc.combined, value = pbmc.combined@meta.data$seurat_clusters)

##################################################################
#Subsetting only for T cells and T cell characterization
##################################################################
#R version 4.3.3
library(Seurat)#Version 4.4.0 
library(ggplot2)
library(dplyr)
#library(data.table)
library(harmony)# Version 1.2.1
library(RSpectra)
library(Nebulosa)
set.seed(2000)

fig_path = "Results_Tcells_CD39"
if(!dir.exists(fig_path)) dir.create(fig_path) 


#saveRDS(pbmc.combined,"Combined_All_Cells_CD39.Rds")
#Further analysis of Tcells by same harmony integration by filtering out genes in the "MitoCarta2_human.txt" file related to mitochondria and other like TCR and ribo
tcell.combined<-subset(x = pbmc.combined, subset = seurat_clusters %in% c( 12,10,19,4,0,1,3))


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

pdf("Results_Tcells_CD39/Tcells_withoutclustering.pdf",width=13,height=12)

DimPlot(tcell.combined, reduction = "umap",order = TRUE, label = TRUE,pt.size=1.2,label.size = 12)
dev.off()
pdf("Results_Tcells_CD39/Tcells_filter_Markers.pdf")
options(repr.plot.width=15, repr.plot.height=15)
FeaturePlot(tcell.combined,features=c("FOXP3","ITGAX","CD4"),pt.size =0.1,order=T)
dev.off()



#saveRDS(tcell.combined,"Tcells_prefilter_CD39.Rds")
#tcell.combined<-readRDS("Tcells_prefilter_CD39.Rds")
#Not including Treg cells, CD4 cells and  APCs from the further analysis marked by FOXP3, CD4, and ITGAX
#remove 3 seurat clusters in total two from bottom and one from top 
tcell.combined_subset<-subset(x = tcell.combined, subset = seurat_clusters %in% c( 0,1,2,3,4,5,6,7,8,9,10,
                                                                                   11,12,13,16,17,18,20,
                                                                                   21,22,23,26))



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

saveRDS(tcell.combined_subset,"Results_Tcells_CD39/Tcells_CD39_harmony_integrated_per_patient.Rds")
#tcell.combined_subset2<-tcell.combined_subset
#Cluster 9 and 2 needs to be reclustered to higher resolution, to recover MAIT cells and Early Tem and Tem cells, so repeat the below code for 9 and 2, then identify the cells from subsetted object
#####################
#remotes::install_version("harmony", version = "1.2.1")
#saveRDS(tcell.combined_subset,"Tcells_CD39_harmony_integrated_P3P5_singlets.Rds")
#tcell.combined_subset<-readRDS("Tcells_CD39_harmony_integrated_P3P5_singlets.Rds")

# 
# 
# 
# color_combo<-c("Tn" = "#85965f",
#                "Tn/Tm" =	"#3b8e87"	,
#                "Early Tem"=	"#EE9572"	,
#                "Tem"=	"#b35f42"	,
#                "Tem-NK like"= "#CD3700"	,	
#                "TOX hi Tex"="#8e477f",
#                "TCF7+ stem-like Tex"="#d379af",	
#                "LAG3 hi Tex"=	"#c7abdd"	 ,
#                "GZMK hi Tex"=	"pink",
#                "MKI67+ Tex/Tprol" ="#a3dbdf"	,
#                "MKI67 hi Tex/Tprol"="#6495ed"	,
#                "MKI67+ Tem-NK like/Tprol"= "#104E8B",	
#                "ISG+"=  "#dac565",
#                "Tc17 MAIT" =	"#d8d5b2" )	
# 
# pdf("Results_Tcells/Tcells_Annotations_Main_B.pdf",width=16.5,height=13)
# DimPlot(tcell.combined_subset,reduction="umap",group.by="annotated_clusters_labels_final",pt.size=3,cols=color_combo,order=c("ISG+","Tc17 MAIT" ,"Tn","Tn/Tm","Tem","Early Tem","Tem-NK like","LAG3 hi Tex","GZMK hi Tex","TOX hi Tex","TCF7+ stem-like Tex","MKI67+ Tem-NK like/Tprol","MKI67+ Tex/Tprol",
#                                                                                                                              "MKI67 hi Tex/Tprol"))+ggtitle(" ")+guides(color = guide_legend(override.aes = list(size=6), ncol=1) )+theme(axis.text=element_text(size=25),text=element_text(size=25))+scale_color_manual("CD8+ T cell states:",breaks = c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol","MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol"),
#                                                                                                                                                                                                                                                                                                                          values = c("#85965f",	"#3b8e87"	,"#EE9572"	,"#b35f42"	, "#CD3700"	,		"#d8d5b2", "#dac565","pink","#d379af",	"#8e477f","#c7abdd"	 ,"#a3dbdf"	,"#6495ed"	,"#104E8B"), 
#                                                                                                                                                                                                                                                                                                                          labels = c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol",
#                                                                                                                                                                                                                                                                                                                                     "MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol" ))
# dev.off()
# 
# 
# 
# #tcell.combined_subset<-readRDS("Tcells_CD39_seurat_integrated_filtered.Rds")
# pdf("CD39_Clusters_Tcells_annotated.pdf",width=15,height = 12)
# DimPlot(tcell.combined_subset, reduction = "umap",group.by="annotated_clusters_labels_final",order = TRUE, label = TRUE,pt.size=1.2,label.size = 8) 
# dev.off()
# DefaultAssay(tcell.combined_subset)<-"RNA"
# FeaturePlot(tcell.combined_subset,features="ITGAX",order=T)
# DimPlot(tcell.combined_subset, reduction = "umap",group.by="annotated_clusters_labels_final",order = TRUE, label = TRUE,pt.size=1.2,label.size = 12) 
# 
# 
# gene_panel<-c("SELL","IL7R","CCR7",#(Naïve/memory makers)
#               "LEF1", "TCF7",  "ZNF683", 
#               "GZMK", "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", "NKG7", "XCL1", "XCL2", #(Effector cytokines)
#               "KLRC3","KIR2DL1",#Tem NK Like
#               "KLRB1", "ZBTB16", # (MAIT markers))# 
#               "IFIT1","STAT1",  "ISG15",# (Interferon stimulated genes)
#               "TOX","PDCD1", "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4", "TNFRSF9", "CXCL13",#(Exhaustion or inhibitory molecules)
#               "MKI67", "TOP2A"# (transcription factors),
#               
# )
# 
# gene_panel_detailed<- c("SELL","IL7R","CCR7","LEF1", "TCF7", # Naïve/memory makers
#                 "ZNF683","GZMK", # pre? Exh
#                 "IFIT1","IFIT3", "IFIT2","STAT1","ISG15",# ISG
#                 "HSPA1A","HSPA6", # HSP
#                 "KLRC3","KIR2DL1",#Tem NK Like
#                 "KLRB1", "ZBTB16", # MAIT
#                 "CXCR6", "GNLY", # eff 1
#                 "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", # eff 2
#                 "AREG","ANKRD28", # Doublet-rich exh
#                 "XCL1", "XCL2", "TNFRSF9", "CXCL13", # Active
#                 "CCL4","IFNG","PDCD1", # exh 1
#                 "TOX","VAV3","HIPK2", # exh 2
#                 "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4", # Exh classic
#                 "MKI67", "TOP2A","CDC20","FOXP3","HIST1H3B") # prol subclass
#         
# tcell.combined$annotated_clusters_labels_final = factor(tcell.combined$annotated_clusters_labels_final, levels=c("MKI67+ Tem-NK like/Tprol", "MKI67 hi Tex/Tprol", "MKI67+ Tex/Tprol", "LAG3 hi Tex", "TOX hi Tex", 
#                                                                                                                  "TCF7+ stem-like Tex", "GZMK hi Tex", "ISG+", "Tc17 MAIT", "Tem-NK like", "Tem", 
#                                                                                                                  "Early Tem", "Tn/Tm", "Tn"))
# 
# h<-DotPlot(tcell.combined_subset, features = gene_panel_detailed, group.by = "annotated_clusters_labels_final", 
#            dot.scale = 8, col.min = -2, col.max = 2) + ylab(" ")+xlab(" ")+
#   scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
#   theme(axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5, hjust = 1),
#         panel.spacing = unit(0.001, "lines"), # Reduce space further
#         axis.text.y = element_text(size = 16)) + 
#   theme(text = element_text(size = 18), legend.title = element_text(size = 15))
# 
# pdf(paste0("Dotplot_Tcell_Markers_Supps_CD39_clusters_annotated.pdf"),width=15,height=6)
# h
# dev.off()
# 
# data_total<-tcell.combined_subset@meta.data
# 
# 
# ggplot(data_total) + 
#   geom_bar(aes(x = seurat_clusters, fill = sample_id), position =  "fill") +  # Adjusting y aesthetic
#   theme_bw() + ylab("Frequency")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=16),
#         legend.key = element_rect(fill = NA),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())  
# 
# 
# p<-ggplot(data_total[data_total$patient_id %in% c("P5"),], aes(x =sample_id , fill = annotated_clusters_labels_final)) + 
#   geom_bar(position = "fill") + 
#   theme_bw() + 
#   ylab("Frequency") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         text = element_text(size = 16),
#         legend.key = element_rect(fill = NA),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# p
# pdf(paste0("boxplot_Tcell_Markers_Supps_CD39_clusters2.pdf"),width=12,height=6)
# p
# dev.off()
# 
# Idents(tcell.combined_subset)="annotated_clusters_labels_final"
# res=FindAllMarkers(tcell.combined_subset,
#                 logfc.threshold = 0, min.pct=0.1, 
#                 test.use = "MAST",
#                 latent.vars = "patient_id")
# res=data.frame(res)
# 
# res$gene<-rownames(res)
# 
# write.csv(res,"CD39_Cluster_DEGs2.csv")
# 
# table(tcell.combined_subset$annotated_clusters_labels_final,tcell.combined_subset$patient_id)
# 
# ####Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis########
# ####Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis########
# ####Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis############Signature Analysis########
# 
# library(readxl)
# tcr_sig = read_excel("tcell_signatures_of_interest.xlsx")
# tcr_sig
# library(tidyr)
# library(dplyr)
# all_signatures <- tcr_sig %>%
#   pivot_longer(cols = everything(), 
#                names_to = "signature", 
#                values_to = "gene") %>%
#   filter(!is.na(gene))  %>% arrange(signature)
# # Optional: drop NAs if any
# 
# print(long_df)
# head(all_signatures)
# write.csv(all_signatures,"all_signatures.csv")
# geneSets <- list(neoTCRCD4=tcr_sig$NeoTCR4,neoTCRCD8=tcr_sig$NeoTCR8)
# 
# 
# library(AUCell)
# 
# 
# cell_types <-  unique(all_signatures$signature)
# cell_type_list <- c()
# for (cell_x in cell_types){
#   to_sigs_sel <- all_signatures[all_signatures$signature == cell_x,]
#   cell_type_list <- c(cell_type_list,list(to_sigs_sel$gene))
# }
# names(cell_type_list) <- cell_types
# 
# 
# exprMatrix <- as.matrix(tcell.combined_subset@assays$RNA@counts)
# 
# 
# cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
# cells_AUC <- AUCell_calcAUC(cell_type_list,cell_rank)
# 
# output_data <- cells_AUC@assays@data@listData$AUC
# output_data <- as.data.frame(t(output_data))
# 
# md <- tcell.combined_subset@meta.data
# check_data_plot <- merge(output_data,md,by=0)
# check_data_plot
# #tcell.combined_subset<-AddMetaData(tcell.combined_subset,check_data_plot)
# colnames(check_data_plot_trimmed)
# check_data_plot_trimmed<-check_data_plot[,c(2:48,54)]
# check_data_plot_trimmed
# check_data_plot_trimmed<-check_data_plot[,c("Rambow Immune", "Pozniak Antigen_presentation",
#                                             "Pozniak Mitotic","Rambow mitosis",
#                                             "Wouters Melanocytic cell state",  "Tsoi Melanocytic","Rambow MITFtargets",  "Rambow pigmentation", 
#                                             "Tsoi Transitory-Melanocytic","Tsoi Transitory",
#                                             "Pozniak Stress (p53 response)",
#                                             "Rambow Neuro", "Tsoi Neural crest-like","Pozniak Neural_Crest_like",
#                                             "Pozniak Stress (hypoxia response)",
#                                             "annotated_clusters")]
# 
# 
# aggregated_AUCell <- check_data_plot_trimmed %>%
#   group_by(sample_id) %>%
#   summarise(across(everything(), mean, na.rm = TRUE))
# 
# AUCell_matrix <- as.matrix(aggregated_AUCell[, -1])
# rownames(AUCell_matrix) <- aggregated_AUCell$sample_id
# 
# 
# order<-c("Immune Response" ,"Mitotic" ,"Patient Specific" ,"Melanocytic" ,"Transitory Melanocytic",
#          "Stress(p53 response)" ,"Neural Crest Like" ,"Stress(hypoxia response)")
# AUCell_matrix_order_row<-AUCell_matrix[order,]
# 
# transposed_matrix <- t(AUCell_matrix)
# 
# Sig_order<-c("Rambow Immune", "Pozniak Antigen_presentation",
#              "Pozniak Mitotic","Rambow mitosis",
#              "Wouters Melanocytic cell state",  "Tsoi Melanocytic","Rambow MITFtargets",  "Rambow pigmentation", 
#              "Tsoi Transitory-Melanocytic","Tsoi Transitory",
#              "Pozniak Stress (p53 response)",
#              "Rambow Neuro", "Tsoi Neural crest-like","Pozniak Neural_Crest_like",
#              "Pozniak Stress (hypoxia response)"
# )
# 
# transposed_matrix_ordered<-transposed_matrix[Sig_order,]
# 
# scaled_matrix <- t(scale(t(transposed_matrix), center = TRUE, scale = TRUE))
# library(ComplexHeatmap)
# library(colorRamp2)
# ht7 <- Heatmap(
#   scaled_matrix, 
#   name = "Z Score", 
#   cluster_rows = FALSE, 
#   cluster_columns = FALSE,
#   column_names_rot = 45,  
#   column_names_gp = gpar(fontsize = 45),  
#   row_names_gp = gpar(fontsize = 45),     
#   heatmap_width = unit(32, "cm"),         
#   heatmap_height = unit(60, "cm"),        
#   col = colorRamp2(breaks=c(-2,0,1.5,2), colors = c("darkblue","white", "red1","red4")),  
#   heatmap_legend_param = list(
#     title_gp = gpar(fontsize = 25),        
#     labels_gp = gpar(fontsize = 20),       
#     legend_height = unit(4, "cm"),        
#     legend_width = unit(2, "cm")          
#   )
# )
# pdf("HeatMap_Tcell_Signatures_Supps_All.pdf",width=50,height=40)
# draw(ht7, heatmap_legend_side = "left")  
# dev.off()



