##################################################################
# Infer CNV
##################################################################

# required input data in the working directory:
# # genes_location_sorted.txt from /data/tumor_data
# # objects generated in /1.Main/ (Tcells_Final.Rds/apc_data_filtered_fvf_corr.rds/Tumor_Annotated.Rds)

#setwd("/YOUR/PATH/")
library(Seurat)
library(infercnv)
library(dplyr)

tcells<-readRDS("Tcells_Final.Rds")
tumorcells<-readRDS("Tumor_Annotated.Rds")
apc<-readRDS("apc_data_filtered_fvf_corr.rds")

tcells$cells<-"T cells"
apc$cells<-"APC cells"
tumorcells$cells<-tumorcells$patient_id


pbmc.combined<-merge(tcells , y=c(tumorcells,apc),project = "scRNA_merge")


counts_matrix <- as.matrix(GetAssayData(pbmc.combined, assay = "RNA", slot = "counts"))


annotations <- data.frame(
  Cell = colnames(pbmc.combined),
  Cluster = pbmc.combined$cells
)

rm(pbmc.combined)
rm(tcells)
rm(tumorcells)
rm(apc)
write.table(annotations, "annotations.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

rm(annotations)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = "annotations.txt",
  gene_order_file = "genes_location_sorted.txt",  # Path to the gene ordering file
  ref_group_names = c("T cells","APC cells")
)

options(scipen = 100)
rm(counts_matrix)


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir="output_directory",
                              cluster_by_groups=TRUE,
                              denoise=TRUE,output_format="pdf",
                              HMM=FALSE,
                              num_threads=30,#choose threads as per the server capability
                              analysis_mode="subclusters")
