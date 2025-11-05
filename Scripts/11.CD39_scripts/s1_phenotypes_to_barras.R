library(data.table)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(AUCell)
library(harmony)

# # # # 
# preceeding AUSTIN CODE TO GET ALL SAMPLES TOGETHER
# # # #
fig_path = "Results_barras_CD39"
if(!dir.exists(fig_path)) dir.create(fig_path)

# object integration and clustering
#five_pat=readRDS("Results_Tcells_CD39/Tcells_CD39_harmony_integrated_per_patient.Rds")
five_pat=readRDS("/Users/a.george/Documents/Sofia_combined/Tcells_CD39_harmony_integrated_per_sample.Rds")


table(five_pat$patient_id)
five_pat$sample_id_notpat=ifelse(five_pat$sample_id %in% c("P3_CD39", "P5_CD39"),
                                 "CD39", "previous")

set.seed(1507)

five_pat <- NormalizeData(five_pat, normalization.method = "LogNormalize", scale.factor = 10000)
five_pat <- FindVariableFeatures(five_pat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(five_pat)
five_pat <- ScaleData(object = five_pat, 
                      do.par = TRUE, num.cores = 8, features =all.genes)

# filter variable genes
gene_names_data = rownames(five_pat@assays$RNA)
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
top.genes<-VariableFeatures(five_pat)
ex_genes<-top.genes[!(top.genes %in% exclude_genes)]
five_pat <- RunPCA(five_pat, npcs = 30, verbose = FALSE,features=ex_genes)

five_pat <- RunHarmony(five_pat, group.by.vars = c("patient_id"))
ElbowPlot(five_pat, ndims = 50, reduction = "pca")
five_pat <- RunUMAP(five_pat, reduction = "harmony", dims = 1:25, min.dist = 0.2)
five_pat <- FindNeighbors(five_pat, reduction = "harmony", dims = 1:25)
five_pat <- FindClusters(five_pat, resolution = 1)
# high dimensionality and resolution, because in the previous analysis we had to subset and recluster = we already knew there were more subpopulations than initially shown.
# clear differences in the former TCF7+ stem group subdivided populations deserve to be annotated and reported.
# identification of expected TC17 and TEM populations
# established a min.dist that shows this clusters as their own populations.

DimPlot(five_pat, label=T)
FeaturePlot(five_pat, order=T, "TCF7")

Seurat::DotPlot(object=five_pat,features = c("SELL","IL7R","CCR7",#(Naïve/memory makers)
                                             "LEF1", "TCF7",  "ZNF683", 
                                             "GZMK", "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", "NKG7", "XCL1", "XCL2", #(Effector cytokines)
                                             "KLRC3","KIR2DL1",#Tem NK Like
                                             "HSPA6", "HSPA1A",
                                             "KLRB1", "ZBTB16", # (MAIT markers))# 
                                             "IFIT1","STAT1",  "ISG15",# (Interferon stimulated genes)
                                             "TOX","PDCD1", "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4", "TNFRSF9", "CXCL13",#(Exhaustion or inhibitory molecules)
                                             "MKI67", "TOP2A"),# (transcription factors),
                assay = "RNA",dot.scale = 6, col.min = -2, col.max = 2)+
  ylab(" ")+xlab(" ")+
  scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.001, "lines"), # Reduce space further
        axis.text.y = element_text(size = 10)) + 
  theme(text = element_text(size = 12), legend.title = element_text(size = 10))

meta=five_pat@meta.data
meta=data.frame(meta)
setDT(meta)
meta[, s_id:= ifelse(grepl("CD39", sample_id), "CD39", ifelse(grepl("DB", sample_id), "DB", "SG"))]
ggplot(meta, aes(s_id, fill=seurat_clusters==10))+geom_bar(position = "fill")+facet_wrap(~patient_id, scale="free")

color_vector <- c(
  "Tn" = "#85965f",
  "Tn/Tm" = "#3b8e87",
  "Early Tem" = "#EE9572",
  "Tem" = "#b35f42",
  "Tem-NK like" = "#CD3700",
  "Tc17 MAIT" = "#d8d5b2",
  "ISG+" = "#dac565",
  "GZMK hi Tex" = "pink",
  "TCF7+ stem-like Tex" = "#d379af",
  "TOX hi Tex" = "#8e477f",
  "LAG3 hi Tex" = "#c7abdd",
  "MKI67+ Tex/Tprol" = "#a3dbdf",
  "MKI67 hi Tex/Tprol" = "#6495ed",
  "MKI67+ Tem-NK like/Tprol" = "#104E8B"
)

# Add the final annotation column

five_pat$final_clusters=factor(ifelse(five_pat$seurat_clusters==2, "LAG3 hi Tex", 
                                      ifelse(five_pat$seurat_clusters==3, "TOX hi Tex", 
                                             ifelse(five_pat$seurat_clusters==7, "TCF7+ stem-like Tex",
                                                    ifelse(five_pat$seurat_clusters %in% c(6), "MKI67+ Tex/Tprol",  
                                                           ifelse(five_pat$seurat_clusters %in% c(8,11,12), "MKI67 hi Tex/Tprol", 
                                                                  ifelse(five_pat$seurat_clusters==9, "ISG+", 
                                                                         ifelse(five_pat$seurat_clusters==14, "HSP", 
                                                                                ifelse(five_pat$seurat_clusters==5, "Tn/Tm",
                                                                                       ifelse(five_pat$seurat_clusters==0, "GZMK hi Tex", 
                                                                                              ifelse(five_pat$seurat_clusters==13, "Tem", 
                                                                                                     ifelse(five_pat$seurat_clusters==15, "Tc17 MAIT",
                                                                                                            ifelse(five_pat$seurat_clusters==1, "Early Tem",
                                                                                                                   ifelse(five_pat$seurat_clusters==4, "Tn",
                                                                                                                          ifelse(five_pat$seurat_clusters %in% c(10), "CD137 hi Tex",
                                                                                                                                 as.character(five_pat$seurat_clusters))))))))))))))), # change
                               levels=c("MKI67 hi Tex/Tprol","MKI67+ Tex/Tprol","LAG3 hi Tex","TOX hi Tex",
                                        "CD137 hi Tex","TCF7+ stem-like Tex","GZMK hi Tex","ISG+","HSP",
                                        "Tc17 MAIT","Tem","Early Tem",
                                        "Tn/Tm", "Tn"))

# plot nice figures
color_vector_2 <- c(
  "Tn" = "#85965f",
  "Tn/Tm" = "#3b8e87",
  "Early Tem" = "#EE9572",
  "Tem" = "#b35f42",
  "Tem-NK like" = "#CD3700",
  "Tc17 MAIT" = "#d8d5b2",
  "ISG+" = "#dac565",
  "GZMK hi Tex" = "pink",
  "TCF7+ stem-like Tex" = "#d379af",
  "TOX hi Tex" = "#8e477f",
  "LAG3 hi Tex" = "#c7abdd",
  "MKI67+ Tex/Tprol" = "#a3dbdf",
  "MKI67 hi Tex/Tprol" = "#6495ed",
  "CD137 hi Tex" = "#907193",
  "HSP" = "gold4"
)

DimPlot(five_pat, group.by = "final_clusters", label = T, repel = T)+scale_color_manual(values = color_vector_2)

# plot figures to characterise the object

# # UMAP
pdf("Results_barras_CD39/new_UMAP_including_CD39_samples.pdf", width=6, height=4)
DimPlot(five_pat, group.by = "final_clusters", pt.size = .01)+scale_color_manual(values = color_vector_2)
dev.off()

# # Dotplot
pdf("Results_barras_CD39/new_dotplot_including_CD39_samples.pdf", width=9, height=3.5)
Seurat::DotPlot(object=five_pat,features = c("SELL","IL7R","CCR7",#(Naïve/memory makers)
                                             "LEF1", "TCF7",  "ZNF683", 
                                             "GZMK", "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", "NKG7",  #(Effector cytokines)
                                             #"KLRC3","KIR2DL1",#Tem NK Like
                                             "KLRB1", "ZBTB16", # (MAIT markers))# 
                                             "HSPA6","HSPA1A",
                                             "IFIT1","STAT1",  "ISG15",# (Interferon stimulated genes)
                                             "TOX",
                                             "TNFRSF9", "CXCL13","XCL1", "XCL2",#(Exhaustion or inhibitory molecules)
                                             "PDCD1", "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4",
                                             "MKI67", "TOP2A"),# (transcription factors),
                assay = "RNA",dot.scale = 6, col.min = -2, col.max = 2, group.by = "final_clusters")+
  ylab(" ")+xlab(" ")+
  scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.001, "lines"), # Reduce space further
        axis.text.y = element_text(size = 10)) + 
  theme(text = element_text(size = 12), legend.title = element_text(size = 10))
dev.off()

# # Sample distribution
five_pat$simple_sample=gsub("P3_|P5_","",five_pat$sample_id)
five_pat$simple_sample_2=factor(ifelse(grepl("DB",five_pat$simple_sample),"DB", five_pat$simple_sample),
                                levels=c("SG_Singlets", "DB", 'CD39'))

# # # umaps
pdf("Results_barras_CD39/new_per_sample_umaps.pdf", width=7, height=8)
DimPlot(five_pat, group.by = "simple_sample_2", pt.size = 0.001)+
  scale_color_manual(values=c("SG_Singlets"="steelblue1",
                              "DB"="red",
                              "CD39"="gold2"))+
  facet_grid(five_pat$patient_id~five_pat$simple_sample_2)
dev.off()

# # # barplots (pat avgs)
meta=data.frame(five_pat@meta.data)
setDT(meta)

meta[, clus_n:= .N, by=c("simple_sample_2", "patient_id","final_clusters")]
meta[, pat_n:= .N, by=c("simple_sample_2", "patient_id")]
meta[, Averaged_Frequency:= clus_n/pat_n]
meta_avg=meta[, lapply(.SD, mean), by=c("simple_sample_2", "patient_id","final_clusters"),
              .SDcols = "Averaged_Frequency" ]

ggplot(meta_avg, aes(patient_id, Averaged_Frequency, fill=final_clusters))+
  geom_col(position="stack")+facet_grid(~simple_sample_2, space="free", scale="free_x")+
  scale_fill_manual(values=color_vector_2)+theme_bw()+xlab("")+ylab("Frequency")

meta_p3p5=meta_avg[patient_id %in% c("P3", "P5")]
ggplot(meta_p3p5, aes(patient_id, Averaged_Frequency, fill=final_clusters))+
  geom_col(position="stack")+facet_grid(~simple_sample_2, space="free", scale="free_x")+
  scale_fill_manual(values=color_vector_2)+theme_bw()+xlab("")+ylab("Frequency")

meta_avg_p3p5=meta_p3p5[, lapply(.SD, sum), by=c("simple_sample_2","final_clusters"),
                        .SDcols = "Averaged_Frequency"]
meta_avg_p3p5[, Averaged_Frequency:=  Averaged_Frequency/2]

pdf("Results_barras_CD39/new_avg_pat_barplot.pdf", width=3.5, height=4.5)
ggplot(meta_avg_p3p5, aes(simple_sample_2, Averaged_Frequency, fill=final_clusters))+
  geom_col(position="stack")+
  scale_fill_manual(values=color_vector_2)+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid = element_blank())+xlab("")
dev.off()

# Add statistics to the enrichment with Mixed model
library(car)
library(nlme)
library(lmerTest)
meta_p3_p5=meta[patient_id %in% c("P3", "P5"),]

metadata_to_sig_enrichment=function(dataset.complete,level_1,level_2){
  data.newtot <- NULL
  for (k in unique(dataset.complete$patient_id)){
    dataset <- dataset.complete[dataset.complete$patient_id==k,]
    for(l in unique(dataset$simple_sample_2)){
      data.seur <- dataset[dataset$simple_sample_2==l,]
      data.seur$tot2 <- sum(data.seur$n)
      data.newtot <- rbind(data.newtot,data.seur)
    }
  }
  dataset.complete <- data.newtot
  
  # now with the classic APC doublets
  dataset2 <- dataset.complete[which(dataset.complete$simple_sample_2 %in% c(level_1,level_2)),]
  
  results_DT <- NULL
  for(k in unique(dataset2$final_clusters)){
    data_cluster <- dataset2[dataset2$final_clusters==k,]
    if (length(unique(data_cluster$simple_sample_2)) < 2) next # added functionality
    response <- cbind(data_cluster$n, data_cluster$tot2 -data_cluster$n)
    print(response)
    model <- glmer(response ~ 1 + simple_sample_2  + (1|patient_id), data=data_cluster, family='binomial')
    res <- summary(model)
    results_DT <- rbind(results_DT, res$coefficients[2,])
  }
  
  results_DT <- data.frame(results_DT)
  results_DT$seurat_cluster <- unique(dataset2$final_clusters)
  results_DT$comparison <- rep(paste(level_1,level_2,sep='-'),nrow(results_DT))
  
  # combine and get the padj
  results_DT$p.adj <- p.adjust(results_DT$Pr...z.., method ='bonferroni')
  results_DT$significant <- ifelse(results_DT$p.adj <0.05, 'Yes', 'No')
  return(results_DT)
}

data=meta_p3_p5[, .N, by=c("patient_id", "simple_sample_2", "final_clusters")]

blank_0s = CJ(patient_id=unique(meta_p3_p5$patient_id),
              simple_sample_2=unique(meta_p3_p5$simple_sample_2),
              final_clusters=unique(meta_p3_p5$final_clusters))
data=merge(blank_0s,data, all.x=T, by=c("patient_id", "simple_sample_2", "final_clusters"))
data[is.na(data)]=0

data[, total:= sum(N), by=c("patient_id", "simple_sample_2")]
setnames(data, "N","n")
data$simple_sample_2=factor(data$simple_sample_2, levels = c("DB","SG_Singlets", "CD39"))
db_cd39=metadata_to_sig_enrichment(data, "DB", "CD39")
fwrite(db_cd39, "Results_barras_CD39/new_DB_cd39_enrichment.txt")
db_sg=metadata_to_sig_enrichment(data, "DB", "SG_Singlets")
fwrite(db_sg, "Results_barras_CD39/new_DB_sg_enrichment.txt")

ggplot(db_cd39, aes(Estimate, -log10(p.adj)))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=seurat_cluster))

db_comps=rbind(db_cd39,db_sg)
setDT(db_comps)
db_comps=dcast(db_comps, seurat_cluster~comparison, value.var = "Estimate")
ggplot(db_comps, aes(-`DB-CD39`,`DB-SG_Singlets`))+geom_point()+
  ggrepel::geom_text_repel(aes(label=seurat_cluster))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

db_comps=rbind(db_cd39,db_sg)
setDT(db_comps)
db_comps[, group:=ifelse(comparison=="DB-CD39" & Estimate<0, "DB-CD39-hi",
                         ifelse(comparison=="DB-CD39" & Estimate>0, "DB-CD39-low",
                                ifelse(comparison=="DB-SG_Singlets" & Estimate>0, "DB-SG-hi",
                                       ifelse(comparison=="DB-SG_Singlets" & Estimate<0, "DB-SG-low","?"))))]
ggvenn::ggvenn(split(db_comps[p.adj<0.05]$seurat_cluster,
                     f=db_comps$group[db_comps$p.adj<0.05]),
               set_name_size = 3, show_percentage = F,
               fill_color = c("firebrick4","gold1","firebrick4","steelblue1"))+ 
  theme(plot.margin = margin(10, 40, 10, 40)) + 
  coord_cartesian(clip = "off")
ggsave("Results_barras_CD39/new_enrichment_venn_diagram.pdf", height=4.5, width=6.5)
write.csv(db_comps, "Results_barras_CD39/S11E.csv")

# Compare frequency among most expanded TCRs
#tcr_data=fread("/DATA/j.simon/DB_revisions_repository/required_input/TCR_filtered_five_pats.txt")
#tcr_data=fread("Results_Tcells_Plots/Meta_data_TCRs_Tcells.txt")
tcr_data=fread("TCR_filtered_five_pats.txt")

p3_cd39_tcrs=fread("CD39/clone.CD39_P3.csv")
p5_cd39_tcrs=fread("CD39/clone.CD39_P5.csv")

tcr_data[, pat_TCR:= paste0(sample, "_", TCR )]
tcr_data[, patient_id:= NULL]
tcr_data[, barcode:= paste0(sample, "__", stringr::str_split_i(barcode, "__",2))]

meta[, barcode:= paste0(combined_id, "__", barcode)]

d_p3_cd39_tcrs=dcast(p3_cd39_tcrs, barcode~chain, value.var = "cdr3", fun.aggregate = function(x){paste0(x, collapse = ";")})
d_p3_cd39_tcrs[, sample:= "P3_CD39"]
d_p3_cd39_tcrs[, pat_TCR:= paste(sample,TRA,TRB, sep="_")]
d_p3_cd39_tcrs[, barcode:= paste0(sample, "__", barcode)]

d_p5_cd39_tcrs=dcast(p5_cd39_tcrs, barcode~chain, value.var = "cdr3", fun.aggregate = function(x){paste0(x, collapse = ";")})
d_p5_cd39_tcrs[, sample:= "P5_CD39"]
d_p5_cd39_tcrs[, pat_TCR:= paste(sample,TRA,TRB, sep="_")]
d_p5_cd39_tcrs[, barcode:= paste0(sample, "__", barcode)]

tcr_data_new=rbindlist(list(tcr_data[,.(barcode,pat_TCR)],
                            d_p3_cd39_tcrs[,.(barcode,pat_TCR)],
                            d_p5_cd39_tcrs[,.(barcode,pat_TCR)]))

tcr_data_new[, pat_TCR:= gsub("_APC|_Tumor","", pat_TCR)]

meta_tcr=merge(meta, tcr_data_new, by="barcode")
meta_tcr[, count_complex_TCR:= .N, by=c("pat_TCR")]
meta_tcr[, freq_tcf7_tex:= sum(grepl("stem", final_clusters)), by="pat_TCR"]
meta_tcr[, freq_cd137_tex:= sum(final_clusters %in% c("CD137 hi Tex")), by="pat_TCR"]
meta_tcr[, freq_tn_tntm_tnEtem_tex:= sum(final_clusters %in% c("Early Tem", "Tn","Tn/Tm")), by="pat_TCR"]
meta_tcr[, freq_lag3_tex:= sum(grepl("LAG3", final_clusters, ignore.case = T)), by="pat_TCR"]

meta_tcr[, freq_tcf7_tex:= freq_tcf7_tex/count_complex_TCR]
meta_tcr[, freq_tn_tntm_tnEtem_tex:= freq_tn_tntm_tnEtem_tex/count_complex_TCR]
meta_tcr[, freq_lag3_tex:= freq_lag3_tex/count_complex_TCR]
meta_tcr[, freq_cd137_tex:= freq_cd137_tex/count_complex_TCR]

exp_meta_tcr=meta_tcr[count_complex_TCR>=10,]

avg_exp_meta_tcr=exp_meta_tcr[, lapply(.SD, mean), by=c('patient_id',"pat_TCR",
                                                        "simple_sample_2",
                                                        "count_complex_TCR"),
                              .SDcols=c("freq_tcf7_tex","freq_lag3_tex",
                                        "freq_tn_tntm_tnEtem_tex","freq_cd137_tex")]

setorder(avg_exp_meta_tcr, -count_complex_TCR)
avg_exp_meta_tcr=avg_exp_meta_tcr[, lapply(.SD, head, n=15), by=c('patient_id', "simple_sample_2")] # select top 15 per group

mexp=melt(avg_exp_meta_tcr, id.vars=c('patient_id',
                                      "simple_sample_2"),
          measure.vars = c("freq_tcf7_tex","freq_lag3_tex",
                           "freq_tn_tntm_tnEtem_tex","freq_cd137_tex"))

mexp[, variable:= ifelse(variable == "freq_tcf7_tex", "TCF7+ stem-like Tex",
                         ifelse(variable == "freq_cd137_tex", "CD137 hi Tex",
                         ifelse(variable == "freq_lag3_tex", "LAG3 hi Tex",
                                ifelse(variable == "freq_tn_tntm_tnEtem_tex", "Naive/Early Tem","??"))))]
ggplot(mexp[patient_id %in% c("P3","P5")], aes(simple_sample_2, value))+
  geom_boxplot(aes(fill=simple_sample_2), outlier.shape = NA)+facet_wrap(~variable, nrow=1, scale="free_y")+
  stat_compare_means(comparisons=list(c("DB","SG_Singlets"),
                                      c("DB","CD39")), method='t', size=3)+
  scale_fill_manual(values=c("SG_Singlets"="steelblue1",
                             "DB"="red",
                             "CD39"="gold2"))+
  theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
                   legend.position="none", panel.grid = element_blank())+
  xlab("")+ylab("Frequency in top15 CTs per group\n(N>=10)")+
  ggbeeswarm::geom_beeswarm(size=1, cex = 2, alpha=.5)
ggsave("Results_barras_CD39/new_per_pheno_expanded_CTs.pdf", width=7.2, height=4.3)
# # # we extend on this observation in the second script.

# save the scRNA object
saveRDS(five_pat,"Results_barras_CD39/new_reprocessed_five_pat_cd39.rds")


# LOAD EXTERNAL DATASET WITH RESONSE INFORMATION
barras_obj=readRDS("extdata/final_cd8_anno.rds") 
# not the one used in previous submission
# the only difference is that the clusters were annotated
# can use the previous object, and change NK in subtypes for clusters_int_0.5!=9
barras_obj=subset(barras_obj, subset = subtypes !="NK")

# reprocessing the object to mimic the conditions in our dataset (reference)
DefaultAssay(barras_obj)='RNA'
barras_obj <- NormalizeData(barras_obj, normalization.method = "LogNormalize", scale.factor = 10000)
barras_obj <- FindVariableFeatures(barras_obj, selection.method = "vst", nfeatures = 2000)
barras_obj <- ScaleData(object = barras_obj,
                        do.par = TRUE, num.cores = 8, features =all.genes)

# create new harmony object to fix lack of feature loadings
# https://github.com/immunogenomics/harmony/pull/144
five_pat[['harmony_2']] <- CreateDimReducObject(embeddings = five_pat[['harmony']]@cell.embeddings,
                                                key = "harmony2_",
                                                loadings = five_pat[['pca']]@feature.loadings,
                                                assay = "RNA")

tcell_anchors <- FindTransferAnchors(reference = five_pat, query = barras_obj, dims = 1:25,
                                     reference.reduction = "harmony_2", query.assay ="RNA",
                                     normalization.method = "LogNormalize",
                                     k.filter = 30) # enough anchors found to be more strict 
saveRDS(tcell_anchors, "Results_barras_CD39/new_tcell_anchors.rds")

# # mapquery function
five_pat_ref <- RunUMAP(five_pat, dims = 1:25, reduction = "harmony_2", return.model = TRUE, min.dist = 0.2)
# a1=DimPlot(five_pat_ref)+theme(legend.position = "none")
# a2=DimPlot(five_pat)+theme(legend.position = "none")
# plot(a1 + a2)
queried_red <- MapQuery(anchorset = tcell_anchors, reference = five_pat_ref, query = barras_obj,
                        refdata = list(celltype = "final_clusters"), reference.reduction = "harmony_2", reduction.model = "umap")

p1 <- DimPlot(five_pat, reduction = "umap", group.by = "final_clusters", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")+
  scale_color_manual(values=color_vector_2)

p2 <- DimPlot(queried_red, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Barras et al., transferred labels")+
  scale_color_manual(values=color_vector_2)

pdf("Results_barras_CD39/new_label_transfer_umaps.pdf", height = 4, width = 7.5)
plot(p1 + p2 )
dev.off()

p3 <- FeaturePlot(five_pat, reduction = "umap","TCF7",order=T) + NoLegend() + ggtitle("Reference annotations")
p4 <- FeaturePlot(subset(queried_red, queried_red$predicted.celltype.score>.4), reduction = "ref.umap","TCF7",order = T) +  NoLegend() + ggtitle("Barras et al., transferred labels")

pdf("Results_barras_CD39/new_tcf7_label_transfer_umaps.pdf", height = 4, width = 7.5)
p3 + p4 + patchwork::plot_layout(nrow = 1, ncol=2)
dev.off()

pdf("Results_barras_CD39/new_score_umaps.pdf", height = 4, width = 4.5)
FeaturePlot(queried_red, reduction = "ref.umap", feature="predicted.celltype.score", pt.size = .1) +
  scale_color_gradient2(high="palegreen3", mid="yellow",low="red4", midpoint = 0.5)
dev.off()

# check: similar patterns?
saveRDS(queried_red,"Results_barras_CD39/new_projected_coukos_data.rds")
queried_red$predicted.celltype=factor(queried_red$predicted.celltype,
                                      levels=c("MKI67 hi Tex/Tprol","MKI67+ Tex/Tprol","LAG3 hi Tex","TOX hi Tex",
                                               "CD137 hi Tex","TCF7+ stem-like Tex","GZMK hi Tex","ISG+",'Mem_like',"HSP",
                                               "Tc17 MAIT","Tem","Early Tem",
                                               "Tn/Tm", "Tn"))

pdf("Results_barras_CD39/new_coukos_after_transfer_dotplot.pdf", width=9, height=3.7)
Seurat::DotPlot(object=subset(queried_red, queried_red$predicted.celltype.score>.4),
                features = c("SELL","IL7R","CCR7",#(Naïve/memory makers)
                                                "LEF1", "TCF7",  "ZNF683", 
                                                "GZMK", "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", "NKG7",  #(Effector cytokines)
                                                #"KLRC3","KIR2DL1",#Tem NK Like
                                                "KLRB1", "ZBTB16", # (MAIT markers))# 
                                                "HSPA6","HSPA1A",
                                                "IFIT1","STAT1",  "ISG15",# (Interferon stimulated genes)
                                                "TOX",
                                                "TNFRSF9", "CXCL13","XCL1", "XCL2",#(Exhaustion or inhibitory molecules)
                                                "PDCD1", "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4",
                                                "MKI67", "TOP2A"),# (transcription factors),
                assay = "RNA",dot.scale = 6, col.min = -2, col.max = 2, group.by = "predicted.celltype")+
  ylab(" ")+xlab(" ")+
  scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.001, "lines"), # Reduce space further
        axis.text.y = element_text(size = 10)) + 
  theme(text = element_text(size = 12), legend.title = element_text(size = 10))
dev.off()

# associate the clusters to annotations
coketa=data.frame(queried_red@meta.data)
ggplot(coketa[coketa$timepoint=="T0",], aes(predicted.celltype, fill=subtypes))+geom_bar(position="fill")+coord_flip()
ggplot(coketa[,], aes(subtypes, fill=predicted.celltype))+geom_bar(position="fill")+coord_flip()+
  scale_fill_manual(values=color_vector_2)

# # some redundant checks on phenotype prediction (this chunk is only for internal consumption)
predictions=data.frame(t(data.frame(queried_red@assays$prediction.score.celltype@data)))
identical(rownames(predictions), rownames(coketa))
co_pred=cbind(coketa, predictions)
setDT(co_pred)
co_pred=co_pred[, lapply(.SD, mean), .SDcols=colnames(predictions), by="subtypes"]
library(ComplexHeatmap)
dt_co_pred=data.frame(co_pred[,2:ncol(co_pred), with=F])
rownames(dt_co_pred)=co_pred$subtypes
pdf("Results_barras_CD39/new_coukos_transfer_pred_scores_heatmap.pdf", height = 4, width=4)
draw(Heatmap(as.matrix(dt_co_pred),col=circlize::colorRamp2(breaks = c(0,.25,.5,1), colors = c("white","coral1","firebrick4","black")), name="Avg\nprediction\nscore"))
dev.off()
draw(Heatmap(scale(dt_co_pred),col=circlize::colorRamp2(breaks = c(-4,0,1,2.5,4), colors = c("darkblue","white","salmon1","red2","firebrick4"))))

# calculate AUC score for neoTCR
library(AUCell)
# # get reads in a sparse matrix format 
exprMatrix= queried_red@assays$RNA@counts
exprMatrix <- as(exprMatrix, "dgCMatrix")
set.seed(333)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)

all_sigs=readxl::read_xlsx("extdata/tcell_signatures_of_interest.xlsx")
all_sigs=lapply(all_sigs, function(x){return(x)})
geneSets=all_sigs
geneSets=lapply(geneSets, function(x){return(x[which(!is.na(x))])})
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc=data.frame(t(data.frame(cells_AUC@assays@data$AUC)))
auc$cell_id=rownames(auc)

setDT(coketa)
coketa=coketa[subtypes!="NK"]
coketa_auc=merge(coketa, auc, by="cell_id")
fwrite(coketa_auc, "Results_barras_CD39/new_auc_scored_barras.txt")
coketa_auc=coketa_auc[predicted.celltype.score>0.4,] # NEW: removing low confidence predictions. 0.5 would remove 180/300 TCF7+.

setDT(coketa_auc)
coketa_auc[, response:= ifelse(as.numeric(gsub("patient","",
                                               stringr::str_split_i(patient_id,"_",1)))%in%
                                 c(2,3,7,8,9,13), "R","NR")]
coketa_auc[, response2:= ifelse(as.numeric(gsub("patient","",
                                                stringr::str_split_i(patient_id,"_",1)))%in%
                                  c(2,3), "CR",
                                ifelse(as.numeric(gsub("patient","",
                                                       stringr::str_split_i(patient_id,"_",1)))%in%
                                         c(8,9,13,7), "PR",
                                       ifelse(as.numeric(gsub("patient","",
                                                              stringr::str_split_i(patient_id,"_",1)))%in%
                                                c(1,4,10), "SD",
                                              ifelse(as.numeric(gsub("patient","",
                                                                     stringr::str_split_i(patient_id,"_",1)))%in%
                                                       c(5,6,11,12), "PD","??"))))]

coketa_auc[,pat:=stringr::str_split_i(patient_id, "_",1)]

# check association with response after NeoTCR correction
coketa_auc[, neotcr_cd8_BIN:= cut(NeoTCR_cd8,
                                  quantile(NeoTCR_cd8, probs = seq(0,1,1/3)),
                                  include.lowest=T,
                                  labels=c("NeoTCR_low","NeoTCR_mid" ,"NeoTCR_high")), by="timepoint"]
coketa_auc[, response:= factor(response, levels=c("R","NR"))]
coketa_auc[, clus_stem_freq_neo:= sum(predicted.celltype=='TCF7+ stem-like Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, tot_pat_neo:= .N, by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_stem_freq_neo:= 100*clus_stem_freq_neo/tot_pat_neo]

coketa_auc[, clus_lag3_freq_neo:= sum(predicted.celltype=='LAG3 hi Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_lag3_freq_neo:= 100*clus_lag3_freq_neo/tot_pat_neo]

coketa_auc[, clus_ne_freq_neo:= sum(predicted.celltype=='CD137 hi Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_ne_freq_neo:= 100*clus_ne_freq_neo/tot_pat_neo]

coketa_auc[, clus_gk_freq_neo:= sum(predicted.celltype=='GZMK hi Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_gk_freq_neo:= 100*clus_gk_freq_neo/tot_pat_neo]

avg_coketa_auc=coketa_auc[, lapply(.SD, mean),
                          .SDcols = c('clus_stem_freq_neo',"clus_lag3_freq_neo","clus_ne_freq_neo","clus_gk_freq_neo"),
                          by=c("pat","neotcr_cd8_BIN","timepoint","response","response2")]

avg_coketa_auc[,neotcr_cd8_BIN:=factor(neotcr_cd8_BIN, levels=c("NeoTCR_high","NeoTCR_mid" ,"NeoTCR_low"))]

pdf("Results_barras_CD39/S11I_t.test_transf_tcf7_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_stem_freq_neo/100))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width = .3,seed = 150799))+
  theme_bw()+ylab("Frequency TCF7+ stem-like Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()
write.csv(avg_coketa_auc[timepoint=="T0"],"Results_barras_CD39/output_s1/S11I_t.test.csv")

# # # # # None of the following figures are included in the manuscript: (remove from final script) # # # # # 
pdf("Results_barras_CD39/new_transf_lag3_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_lag3_freq_neo))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency LAG3 high Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()


pdf("Results_barras_CD39/new_transf_cd137_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_ne_freq_neo))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency CD137 high Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

pdf("Results_barras_CD39/new_transf_gzmk_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_gk_freq_neo))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency GZMK high Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

pdf("Results_barras_CD39/new_transf_stem_among_alltex_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, 100*(clus_stem_freq_neo/(clus_gk_freq_neo+clus_stem_freq_neo+clus_ne_freq_neo+clus_lag3_freq_neo))))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("%TCF7+ within Tex states")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

# not biasing by neoTCR
coketa_auc[, clus_stem_freq:= sum(predicted.celltype=='TCF7+ stem-like Tex'), by=c("pat","timepoint")]
coketa_auc[, tot_pat:= .N, by=c("pat","timepoint")]
coketa_auc[, clus_stem_freq:= 100*clus_stem_freq/tot_pat]

coketa_auc[, clus_lag3_freq:= sum(predicted.celltype=='LAG3 hi Tex'), by=c("pat","timepoint")]
coketa_auc[, clus_lag3_freq:= 100*clus_lag3_freq/tot_pat]

coketa_auc[, clus_gk_freq:= sum(predicted.celltype=='GZMK hi Tex'), by=c("pat","timepoint")]
coketa_auc[, clus_gk_freq:= 100*clus_gk_freq/tot_pat]

coketa_auc[, clus_ne_freq:= sum(predicted.celltype=='CD137 hi Tex'), by=c("pat","timepoint")]
coketa_auc[, clus_ne_freq:= 100*clus_ne_freq/tot_pat]

avg_coketa_pat=coketa_auc[, lapply(.SD, mean), .SDcols = c('clus_stem_freq',"clus_lag3_freq","clus_ne_freq","clus_gk_freq"),
                          by=c("pat","timepoint","response","response2")]

pdf("Results_barras_CD39/new_total_freq_tcf7_in_response.pdf", width=3.2, height=3.5)
ggplot(avg_coketa_pat[timepoint=="T0"],aes(response, clus_stem_freq))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  stat_compare_means(method="t")+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency TCF7+ stem-like Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

pdf("Results_barras_CD39/new_total_freq_GZMK_in_response.pdf", width=3.2, height=3.5)
ggplot(avg_coketa_pat[timepoint=="T0"],aes(response, clus_gk_freq))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  stat_compare_means(method="t")+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency GZMK high Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

pdf("Results_barras_CD39/new_total_freq_cd137_in_response.pdf", width=3.2, height=3.5)
ggplot(avg_coketa_pat[timepoint=="T0"],aes(response, clus_ne_freq))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  stat_compare_means(method="t")+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency CD137 high Tex")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

pdf("Results_barras_CD39/new_total_freq_all.ex_in_response.pdf", width=3.2, height=3.5)
ggplot(avg_coketa_pat[timepoint=="T0"],aes(response, clus_ne_freq+clus_gk_freq+clus_lag3_freq+clus_stem_freq))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  stat_compare_means(method="t")+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency all Tex states")+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()


pdf("Results_barras_CD39/new_total_freq_lag3_in_response.pdf", width=3.2, height=3.5)
ggplot(avg_coketa_pat[timepoint=="T0"],aes(response, clus_lag3_freq))+
  geom_boxplot(aes(fill=response), outlier.shape = NA)+
  stat_compare_means(method="t")+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  ylab("Frequency LAG3 hi Tex")+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()
