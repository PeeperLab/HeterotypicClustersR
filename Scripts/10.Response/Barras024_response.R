##################################################################
# Barras et al., reactivity association to response
##################################################################

# required input data in the working directory:
# # external dataset
# # signature output from /3.Signatures/s1.create_signature.R

#setwd("/YOUR/PATH/")
set.seed(150799)
library(Seurat)
library(ggpubr)
library(stringr)
library(data.table)
library(RSpectra)
library(ggplot2)
library(dplyr)
library(AUCell)

# external dataset
counts_files=list.files("./GSE221553_RAW/")[grepl("counts",list.files("./GSE221553_RAW/"))]

# aggregate counts by simply merging
setwd("./GSE221553_RAW/")
objs_all_pats=lapply(counts_files, function(x){
  print(100*which(counts_files==x)/length(counts_files))
  counts=as(as.matrix(read.delim(x)), "dgCMatrix")
  obj=CreateSeuratObject(counts, assay = "RNA")
  return(obj)
})

seu_obj=scCustomize::Merge_Seurat_List(
  objs_all_pats,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject"
)
rm(objs_all_pats); gc()
saveRDS(seu_obj, "../merged_raw_counts.rds")
print("step1- complete")

# normalise and process
# seu_obj=readRDS("../merged_raw_counts.rds")
seu_obj=NormalizeData(seu_obj, scale.factor = 10000)
seu_obj[["percent.mt"]] <- PercentageFeatureSet(object = seu_obj, pattern = "^MT-")
seu_obj[["percent.ribo"]] <- PercentageFeatureSet(object = seu_obj,
                                                  features = rownames(seu_obj)[grepl("^RPL", rownames(seu_obj))|
                                                                                 grepl("^RPS", rownames(seu_obj))])
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu_obj=CellCycleScoring(seu_obj, s.features = s.genes, g2m.features = g2m.genes)
print("step2- checkpoint")

library(future)
# plan(strategy = "multicore", workers = 3)
# options(future.globals.maxSize= 50*2048 * 1024^2)
seu_obj=ScaleData(seu_obj, vars.to.regress = c("percent.mt","percent.ribo","S.Score", "G2M.Score","nCount_RNA"))
saveRDS(seu_obj, '../norm_scaled_seurat_object_new.rds')
print("step2- complete")

seu_obj<-readRDS("norm_scaled_seurat_object_new.rds")
mt<- seu_obj@meta.data
unique(mt$patient_id)
# associate the cells to their patients of origin
seu_obj$patient_id= str_split_i(colnames(seu_obj), "\\.",2)
seu_obj$cd45_pos= ifelse(grepl("45",colnames(seu_obj)), "CD45","other")
seu_obj$timepoint= ifelse(grepl("T30",colnames(seu_obj)), "T30","T0")



# process towards main types
set.seed(150799)
seu_obj= FindVariableFeatures(seu_obj)
seu_obj <- RunPCA(seu_obj, npcs = 75)
seu_obj <- RunTSNE(seu_obj, method="Rtsne")
seu_obj <- RunUMAP(seu_obj, dims = 1:75, min.dist = 0.75)

seu_obj <- FindNeighbors(seu_obj, dims = 1:75,nn.eps=0.5)

seu_obj <- FindClusters(seu_obj, resolution = .5)


pdf("class_division_dotplot_res05.pdf", height=8, width=6)
DotPlot(seu_obj, group.by="RNA_snn_res.0.5",
        features=c(
          "MLANA", "PRAME", "SOX10", "S100B", #melanoma
          "PTPRC", #immune cells
          "CD3E", "CD8A", "CD8B", "CD4", # T cells
          "CD79A", "MS4A1", # B cells
          "DCN", "FAP", # CAFs
          "PECAM1", "VWF", # endothelial cells
          "IGKC", "IGHG1", "IGHA1", #plasma cells
          "CD68", "HLA-DRA", "LYZ", "CD86" # myeloid
        ))+
  scale_color_gradient2(high="red4", low="blue4", mid = "lightyellow2")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=8),
        axis.text.y = element_text(size=8),
        legend.title = element_text(size=8))
dev.off()

# get only the group of interest and re-integrate

seu_obj$first_split=ifelse(seu_obj$RNA_snn_res.0.5 %in% 
                             c('9',"6","7","22","19","15","13","11","0",
                               "24", "2","16","1"),
                           "cd8_nk_cd4", "other")


seu_obj=subset(seu_obj, first_split=="cd8_nk_cd4")
saveRDS(seu_obj, "after_cd8a_cd4_nk_filtering.rds")


################### here finishes the class selection ################### 

############# here starts the cd4 and gamma delta filtering  #############

seu_obj <- FindClusters(seu_obj, resolution = 10)

# first classify on per cell basis
count_matrix=t(data.frame(seu_obj@assays$RNA$data[rownames(seu_obj@assays$RNA$data) %in% c("CD8A", "CD4"),]))
count_matrix=data.frame(count_matrix)
count_matrix$cell_id=rownames(count_matrix)
setDT(count_matrix)
count_matrix[, dn_dp:= ifelse(CD8A+CD4==0, "dn",
                              ifelse(CD8A>0 & CD4>0, "dp", 
                                     ifelse(CD8A>0 & CD4==0, "cd8", 
                                            ifelse(CD8A==0 & CD4>0, "cd4", "HUH?"))))]
table(count_matrix$dn_dp)

# second get per cluster res10 classif
meta=data.frame(seu_obj@meta.data)
meta$cell_id=rownames(meta)
setDT(meta)

meta_counts=merge(meta[, cell_id:= paste0("X",cell_id)],count_matrix, by="cell_id")
meta_counts[, perc75_cd8:= quantile(CD8A, probs=0.75), by="RNA_snn_res.10"]
meta_counts[, perc75_cd4:= quantile(CD4, probs=0.75), by="RNA_snn_res.10"]
meta_counts[, perc_call:= ifelse(perc75_cd8>perc75_cd4, "cd8_q",
                                 ifelse(perc75_cd8<perc75_cd4, "cd4_q",
                                        ifelse( perc75_cd8+perc75_cd4==0, "dn", "HUH?")))]


# combine both parameters
meta_counts[, cd8_call:= ifelse((perc_call=="cd8_q" & dn_dp %in% c("cd8", "dn"))|dn_dp=="cd8","cd8",
                                ifelse(dn_dp=="dp", "dp",
                                       ifelse((perc_call=="cd4_q" & dn_dp %in% c("cd4", "dn"))|dn_dp=="cd4","cd4",
                                              ifelse(perc_call=="dn" & dn_dp=='dn', 'dn', "other"))))]

colnames(seu_obj)=paste0("X", colnames(seu_obj))
meta_counts=data.frame(meta_counts)
rownames(meta_counts)=meta_counts$cell_id
seu_obj=AddMetaData(seu_obj, meta_counts)

identical(as.character(seu_obj$cell_id), colnames(seu_obj))

VlnPlot(seu_obj, "nCount_RNA", group.by = "cd8_call")+geom_boxplot()+scale_y_log10() # their DP theory could make sense
VlnPlot(seu_obj, "rna_CD8A", group.by = "cd8_call")+geom_boxplot()
VlnPlot(seu_obj, "rna_CD4", group.by = "cd8_call")+geom_boxplot()

# try to identify gamma deltas
seu_obj[["gamma_delta"]] <- PercentageFeatureSet(object = seu_obj, pattern = "^TR[GD]")
hist(seu_obj$gamma_delta)

count_matrix_GD=t(data.frame(seu_obj@assays$RNA$data[grepl("^TR[GD]", rownames(seu_obj@assays$RNA$data)),]))
count_matrix_GD=data.frame(count_matrix_GD)
count_matrix_GD$TRGD_agg=rowMeans(count_matrix_GD)
count_matrix_GD$cell_id=rownames(count_matrix_GD)

hist(count_matrix_GD$TRGD_agg)
table(count_matrix_GD$TRGD_agg>0.5)

seu_obj$gamma_delta=ifelse(seu_obj$cell_id %in% count_matrix_GD[count_matrix_GD$TRGD_agg>0.5,]$cell_id, "GD","other")
table(seu_obj$gamma_delta)

# # try to identify NKs based on markers
# 
seu_obj$final_call_with_nks=ifelse(seu_obj$cd8_call=="cd8"  & seu_obj$gamma_delta=="other", "CD8pos", "other")
table(seu_obj$final_call)

VlnPlot(seu_obj, "GNLY", group.by = "nk")+geom_boxplot()
saveRDS(seu_obj, "new_cd8_annotation.rds")

##### actual filtering - with the papers indications
seu_obj=subset(seu_obj, final_call_with_nks=="CD8pos")
#####

seu_obj_list=SplitObject(seu_obj, split.by="patient_id")
all_genes=rownames(seu_obj)
rm(seu_obj)

# process and integrate
library(dplyr)
seu_obj_list=seu_obj_list[names(seu_obj_list)!="patient5_T30"]
seu_obj_list=lapply(seu_obj_list, function(x){
  x<- NormalizeData(x) %>% FindVariableFeatures() 
})

sif=SelectIntegrationFeatures(seu_obj_list, nfeatures = 800)
tcr_genes = grep("^TRB[V|D|J|]|^TRA[V|D|J|]|^TRD[V|D|J|]|^TRG[V|D|J|]", all_genes, v=T, perl=T)
sif_noTCR=setdiff(sif, tcr_genes)

seu_obj_list=lapply(seu_obj_list, function(x){
  x<- ScaleData(x, vars.to.regress = c("percent.mt","percent.ribo","S.Score", "G2M.Score","nCount_RNA"),sif_noTCR) %>%
    RunPCA(features=sif) 
})

anchors=FindIntegrationAnchors(seu_obj_list, reduction = "rpca", dims = 1:30)

final_obj=IntegrateData(anchors)

# reduce dimensionality and recluster

final_obj=ScaleData(final_obj, assay = "integrated", vars.to.regress = c("percent.mt","percent.ribo","S.Score", "G2M.Score","nCount_RNA")) 
final_obj <- RunPCA(final_obj, npcs = 20, assay = "integrated")
final_obj <- RunTSNE(final_obj, method="Rtsne", assay = "integrated")
final_obj <- RunUMAP(final_obj, dims = 1:20, min.dist = 0.75, assay = "integrated")

final_obj <- FindNeighbors(final_obj, dims = 1:20,n.eps=0.5)
final_obj <- FindClusters(final_obj, resolution = .5, assay = "integrated")
saveRDS(final_obj, "processed_cd8_nk_final.rds")

pdf("quality_control_dot.plot_with_nk.pdf", width=12, height=5)
DotPlot(final_obj, assay = "RNA",  c("CD3E", "rna_CD8A", "rna_CD4",
                                     "TPT1","FOS","IL7R",
                                     "GZMK", "GZMA", "CCL5", 
                                     "CXCL13", "CRTAM", "FABP5",
                                     "LAG3", "HAVCR2", "CCL3",
                                     "HSPA1A","HSPA6", "DNAJB1",
                                     "FOXP3", "CTLA4", "IL32", 
                                     "FGFBP2", "CX3CR1", "GZMH",
                                     "ISG15","MX1", "IFI6",
                                     "MALAT1", "NEAT1","CD164",
                                     "ZNF683", "TRDV1","LEF1",
                                     "TYROBP", "FCER1G", "KLRB1",
                                     "GNLY", "KLRD1","KLRC1",
                                     "TRDV2", "TRGV1"), group.by = "integrated_snn_res.0.5")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_color_gradient2(high="red4", low="darkblue", mid="lightyellow2")
dev.off()


# next step: remove nk/ score auc/ produce figures
DimPlot(final_obj, group.by = "integrated_snn_res.0.5", label = T)

pdf("NK_selection.pdf", width=10, height=10)
FeaturePlot(final_obj,c("KLRD1","GNLY", "KLR","TYROBP"), label=T)
dev.off()

pdf("quality_control_feature.plots.pdf", width=16, height=14)
FeaturePlot(final_obj,c("IL7R", "GZMA","FABP5","CRTAM","HAVCR2", "DNAJB1","CTLA4","GZMH","ISG15","NEAT1","TCF7",'TOX'))
dev.off()


#final_obj<-readRDS("final_obj.Rds")

meta=final_obj@meta.data
db_30=scan("top30_genes_DB_signature.txt",what="character")
db_100=scan("top100_genes_DB_signature.txt",what="character")

filt_obj=subset(final_obj,integrated_snn_res.0.5!="9" )
exprMatrix= filt_obj@assays$RNA@counts
exprMatrix <- as(exprMatrix, "dgCMatrix")

set.seed(333)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(list("db_30"=db_30, "db_100"=db_100), cells_rankings)

auc_numbers=t(data.frame(cells_AUC@assays@data$AUC))
auc_numbers=data.frame(auc_numbers)
auc_numbers$cell_id=rownames(auc_numbers)
auc_numbers_merged=merge(auc_numbers, meta, by="cell_id")

setDT(auc_numbers_merged)
auc_numbers_merged=auc_numbers_merged[timepoint=="T0",]
auc_numbers_merged[, response:= ifelse(patient_id %in% c("patient2_T0_CD45", "patient3_T0_CD45", "patient7_T0_CD45",
                                                         "patient8_T0_CD45", "patient9_T0_CD45", "patient13_0K_T0_CD45","patient13_1Z_T0_CD45"), "R",
                                       ifelse(patient_id %in% c("patient1_T0_CD45", "patient4_T0_CD45", "patient10_B_T0_CD45","patient10_A_T0_CD45",
                                                                "patient5_T0_CD45", "patient6_T0_CD45", "patient11_T0_CD45", "patient12_T0_CD45"),"NR",
                                              "other"))]


auc_numbers_merged[, response2:= ifelse(patient_id %in% c("patient2_T0_CD45", "patient3_T0_CD45"),'CR',
                                        ifelse(patient_id %in% c( "patient7_T0_CD45","patient8_T0_CD45", "patient9_T0_CD45", "patient13_0K_T0_CD45","patient13_1Z_T0_CD45"), "PR",
                                               ifelse(patient_id %in% c("patient1_T0_CD45", "patient4_T0_CD45", "patient10_B_T0_CD45","patient10_A_T0_CD45"), "SD",
                                                      ifelse(patient_id %in% c("patient5_T0_CD45", "patient6_T0_CD45", "patient11_T0_CD45", "patient12_T0_CD45"),"PD",
                                                             "other"))))]

setDT(auc_numbers_merged)

auc_numbers_merged[ ,pat_unique:= ifelse(grepl("patient10", patient_id) & timepoint=="T0","patient10_T0_CD45",patient_id)]
auc_numbers_merged[ ,pat_unique:= ifelse(grepl("patient13", pat_unique) & timepoint=="T0","patient13_T0_CD45",pat_unique)]
fwrite(auc_numbers_merged, "auc_scored_cd8s.txt")
auc_numbers_merged[, response2:= factor(response2, levels=c("CR","PR","SD","PD"))]
auc_numbers_merged[, response:= factor(response, levels=c("R","NR"))]

ggplot(auc_numbers_merged[timepoint=="T0",],
       aes(response, db_30))+geom_violin(draw_quantiles = c(.25,.5,.75), aes(fill=response))+geom_boxplot(alpha=.3)+
  stat_compare_means(comparisons = list(c("R","NR")),method = 't.test')+
  scale_fill_manual(values=c("NR"="grey10",
                             "R"="coral1"))+theme_bw()
ggsave("all_cells_response_T0.pdf", width = 4, height=4)


ggplot(auc_numbers_merged[timepoint=="T0",],
       aes(response, db_30))+geom_violin(draw_quantiles = c(.25,.5,.75), aes(fill=response2))+
  stat_compare_means(comparisons = list(c("R","NR")),method = 't.test')+
  scale_fill_manual(values=c("SD"="grey60","PD"="grey10",
                             "PR"="coral4","CR"="coral1"))+theme_bw()
ggsave("all_cells_response2_T0.pdf", width = 4.3, height=4)

mavg_2_10s=auc_numbers_merged[, lapply(.SD, mean), by=c("patient_id","response", "response2"), .SDcols=c("db_30","db_100")]
mavg=auc_numbers_merged[, lapply(.SD, mean), by=c("pat_unique","response", "response2"), .SDcols=c("db_30","db_100")]




# # divide reactive cells into 3 groups
hist(auc_numbers_merged$db_30)
hist(auc_numbers_merged$db_100)

auc_numbers_merged[, db_100_split:= ifelse(db_100< quantile(db_100, .33), "low",
                                           ifelse(db_100>=quantile(db_100, .67), "high",
                                                  "medium"))]
auc_numbers_merged[, db_100_split:= factor(db_100_split, levels=c("high","medium","low"))]

auc_numbers_merged[, db_30_split:= ifelse(db_30< quantile(db_30, .33), "low",
                                          ifelse(db_30>=quantile(db_30, .67), "high",
                                                 "medium"))]
auc_numbers_merged[, db_30_split:= factor(db_30_split, levels=c("high","medium","low"))]

auc_numbers_merged[, tot_pat_unique:= .N, by="pat_unique"]
auc_numbers_merged[, tot_site:= .N, by="patient_id"]

auc_numbers_merged[, freq_cat_pat_100:= .N/tot_pat_unique, by=c("pat_unique", "db_100_split")]
auc_numbers_merged[, freq_cat_pat_30:= .N/tot_pat_unique, by=c("pat_unique", "db_30_split")]
auc_numbers_merged[, freq_cat_site_100:= .N/tot_site, by=c("patient_id", "db_100_split")]
auc_numbers_merged[, freq_cat_site_30:= .N/tot_site, by=c("patient_id", "db_30_split")]



mavg_db100_pat=auc_numbers_merged[, lapply(.SD, mean),
                                  by=c("pat_unique","response", "response2", "db_100_split"),
                                  .SDcols=c("db_100","freq_cat_pat_100")]

mavg_db100_site=auc_numbers_merged[, lapply(.SD, mean),
                                   by=c("patient_id","response", "response2", "db_100_split"),
                                   .SDcols=c("db_100","freq_cat_site_100")]

mavg_db30_pat=auc_numbers_merged[, lapply(.SD, mean),
                                 by=c("pat_unique","response", "response2", "db_30_split"),
                                 .SDcols=c("db_30","freq_cat_pat_30")]

mavg_db30_site=auc_numbers_merged[, lapply(.SD, mean),
                                  by=c("patient_id","response", "response2", "db_30_split"),
                                  .SDcols=c("db_30","freq_cat_site_30")]

ggplot(mavg_db30_pat, aes(response, freq_cat_pat_30))+
  geom_boxplot(aes(fill = response), alpha=.4)+
  geom_point(aes(colour = response2), size=3, alpha=1,
             position = position_jitter(width=.3))+
  stat_compare_means(method="t.test")+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="coral4","CR"="coral1"))+
  scale_fill_manual(values=c("NR"="grey60", "R"="coral1"))+
  theme_bw()+ggtitle("Per pat")+
  facet_wrap(~db_30_split)
ggsave("reactive_freq_per.pat_db30_t.test.pdf", width=6, height=4)


ggplot(mavg_db100_pat, aes(response, freq_cat_pat_100))+
  geom_boxplot(aes(fill = response), alpha=.4)+
  geom_point(aes(colour = response2), size=3, alpha=1,
             position = position_jitter(width=.3))+
  stat_compare_means(method="t.test")+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="coral4","CR"="coral1"))+
  scale_fill_manual(values=c("NR"="grey60", "R"="coral1"))+
  theme_bw()+ggtitle("Per pat")+
  facet_wrap(~db_100_split)
ggsave("reactive_freq_per.pat_db100_t.test.pdf", width=6, height=4)

