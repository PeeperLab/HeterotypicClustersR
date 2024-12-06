##################################################################
# Clonotypes per interaction status analysis
##################################################################

# required input data in the working directory:
# # TCR_filtered_five_pats.txt from Tcell_Plots.R
# # auc_numbers_merged.rds from s3_AUC_object.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

fig_path = "Results_clonotype_heatmap"
if(!dir.exists(fig_path)) dir.create(fig_path)

# Load generated objects
tcr_ids=fread("Meta_data_TCRs_Tcells.txt")
tcr_ids$barcode<-paste0(tcr_ids$combined_id,"_",tcr_ids$cell_id)
tcr_ids<-tcr_ids[,c("barcode","combined_id","patient_id","chain_cdr3")]
colnames(tcr_ids)<-c("barcode", "sample" ,"patient_id", "TCR" )

tcr_ids[, cell_id:= gsub(paste0(sample,"_"), "", barcode), by="sample"]


# GET STARTED
auc=readRDS("auc_numbers_merged.rds")
setDT(auc)
auc$cell_id=gsub("\\.","-",auc$cell_id)

auc_dt=merge(auc, tcr_ids, by=c("cell_id", "patient_id"))

signature_panel=c("neotcr_cd8","offringa_TR_9_samples",
                  "oliv_Virus.specific","oliv_Tumor.specific")
auc_mat=auc_dt[,colnames(auc_dt) %in% signature_panel, with=F]

levels=c("Tn", "Tn/Tm", "Early Tem", "Tem", "Tem-NK like", "Tc17 MAIT", "ISG+",
         "GZMK hi Tex", "TCF7+ stem-like Tex", "TOX hi Tex", "LAG3 hi Tex",
         "MKI67+ Tex/Tprol", "MKI67 hi Tex/Tprol", "MKI67+ Tem-NK like/Tprol")

auc_dt[, sample_id:=factor(sample_id, levels=c("SG_Singlets", "DB_Tumor_Tcell","DB_APC_Tcell"))]
auc_dt[, annotated_clusters_labels_final:=factor(annotated_clusters_labels_final, levels=levels)]

# tcr/interaction parameters - create groups per clonotype and sample_id
auc_dt[, id_TCR:= paste0(patient_id, "_", sample_id, "_", TCR)]
auc_dt[, n_id_TCR:= .N, by="id_TCR"]
tcr_scor=auc_dt[,lapply(.SD, mean), by=c("id_TCR", "patient_id", "sample_id","n_id_TCR"),
                .SDcols = signature_panel]


# # if using only tcr groups larger than 10 Tcells
tcr_scor_filt10=tcr_scor[n_id_TCR>=10,]
ha32_filt10=HeatmapAnnotation(df=tcr_scor_filt10[,.("Interactive state"=sample_id)],
                              col = list(
                                "Interactive state"=c("DB_APC_Tcell"="#009052",
                                                      "DB_Tumor_Tcell"="#825d9e",
                                                      "SG_Singlets" = "#7cd3f7")))

heat_filt10= tcr_scor_filt10[,names(tcr_scor_filt10) %in% signature_panel, with=F]

pdf("Results_clonotype_heatmap/reactivity_heatmap.pdf", width=10, height=3.6)
draw(Heatmap(t(scale(heat_filt10)), top_annotation = ha32_filt10,
             name="Row Z-scores", row_names_gp = gpar(cex=.6)), merge_legends = TRUE)
dev.off()

# heatmap to compile all per tcr comparisons by patient
m_scor=melt(tcr_scor, id.vars = c("id_TCR", "patient_id", "sample_id"), measure.vars = names(tcr_scor)[5:ncol(tcr_scor)])
setDT(m_scor)
m_scor_list=split(m_scor, by=c("patient_id", "variable"))

m_scor_agg=rbindlist(lapply(m_scor_list, function(x){
  
  # log2FCs
  x[, apc_mean:= mean(value[sample_id=="DB_APC_Tcell"])]
  x[, tum_mean:= mean(value[sample_id=="DB_Tumor_Tcell"])]
  x[, sg_mean:= mean(value[sample_id=="SG_Singlets"])]
  
  x[, log2fc_tum_sg:=log2(tum_mean/sg_mean)]
  x[, log2fc_apc_sg:=log2(apc_mean/sg_mean)]
  
  # pvals
  x[, pval_tum:= wilcox.test(value[sample_id=="DB_Tumor_Tcell"],value[sample_id=="SG_Singlets"])$p.value]
  x[, pval_apc:= wilcox.test(value[sample_id=="DB_APC_Tcell"],value[sample_id=="SG_Singlets"])$p.value]
  
  # aggregate
  x=x[, lapply(.SD, mean),by=c("patient_id", "variable"),
      .SDcols=c("log2fc_tum_sg", 'log2fc_apc_sg', 
                "pval_tum", "pval_apc")]
  return(x)
}))

# format lfc
m_scor_agg_lfc_tum=dcast(m_scor_agg, patient_id~ variable, value.var = "log2fc_tum_sg")
m_scor_agg_lfc_tum[, patient_id:= paste0(patient_id,"_","tum")]
m_scor_agg_lfc_apc=dcast(m_scor_agg, patient_id~ variable, value.var = "log2fc_apc_sg")
m_scor_agg_lfc_apc[, patient_id:= paste0(patient_id,"_","apc")]

m_scor_agg_lfc=rbind(m_scor_agg_lfc_tum,m_scor_agg_lfc_apc)

library(stringr)
m_scor_agg_lfc$comp=factor(str_split_i(m_scor_agg_lfc$patient_id,"_",2), levels=c("tum","apc"))
m_scor_agg_lfc$patient_id=str_split_i(m_scor_agg_lfc$patient_id,"_",1)

setorder(m_scor_agg_lfc, patient_id, comp)
mat_scor_agg_lfc=as.matrix(m_scor_agg_lfc[,2:(ncol(m_scor_agg_lfc)-1), with=F])
rownames(mat_scor_agg_lfc)=paste0(m_scor_agg_lfc$patient_id,"_",m_scor_agg_lfc$comp)

m_scor_agg_lfc_tum=dcast(m_scor_agg, patient_id~ variable, value.var = "log2fc_tum_sg")
m_scor_agg_lfc_tum[, patient_id:= paste0(patient_id,"_","tum")]
m_scor_agg_lfc_apc=dcast(m_scor_agg, patient_id~ variable, value.var = "log2fc_apc_sg")
m_scor_agg_lfc_apc[, patient_id:= paste0(patient_id,"_","apc")]

m_scor_agg_lfc=rbind(m_scor_agg_lfc_tum,m_scor_agg_lfc_apc)

library(stringr)
m_scor_agg_lfc$comp=factor(str_split_i(m_scor_agg_lfc$patient_id,"_",2), levels=c("tum","apc"))
m_scor_agg_lfc$patient_id=str_split_i(m_scor_agg_lfc$patient_id,"_",1)

setorder(m_scor_agg_lfc, patient_id, comp)
mat_scor_agg_lfc=as.matrix(m_scor_agg_lfc[,2:(ncol(m_scor_agg_lfc)-1), with=F])
rownames(mat_scor_agg_lfc)=paste0(m_scor_agg_lfc$patient_id,"_",m_scor_agg_lfc$comp)

# format significance
m_scor_agg_pval_tum=dcast(m_scor_agg, patient_id~ variable, value.var = "pval_tum")
m_scor_agg_pval_tum[, comp:= "tum"]
m_scor_agg_pval_apc=dcast(m_scor_agg, patient_id~ variable, value.var = "pval_apc")
m_scor_agg_pval_apc[, comp:= "apc"]

m_scor_agg_pval=rbind(m_scor_agg_pval_tum,m_scor_agg_pval_apc)
m_scor_agg_pval$comp
m_scor_agg_pval$comp=factor(m_scor_agg_pval$comp, levels=c("tum","apc"))

setorder(m_scor_agg_pval, patient_id, comp)
mat_scor_agg_pval=as.matrix(m_scor_agg_pval[,2:(ncol(m_scor_agg_pval)-1), with=F])
rownames(mat_scor_agg_pval)=paste0(m_scor_agg_pval$patient_id,"_",m_scor_agg_pval$comp)
identical(rownames(mat_scor_agg_pval),rownames(mat_scor_agg_lfc))

anno=rowAnnotation(df=data.frame(m_scor_agg_lfc[,.(patient_id, comp)]),
                   col=list(comp=c("apc"="#009052", "tum"="#825d9e"),
                            patient_id=c("P1"="#440154", "P2"="#3b528b", "P3"="#21918c", "P4"="#5ec962", "P5"="#fde725")))

mat_scor_agg_lfc=mat_scor_agg_lfc[,c(3,4,1,2)]
mat_scor_agg_lfc
library(circlize)
pdf("Results_clonotype_heatmap/per_pat_heatmap_stats_without_singlets.pdf", height=4, width=6)
Heatmap(mat_scor_agg_lfc, cluster_rows = F, cluster_columns=F,name = "Log2FC",
        left_annotation = anno,col=colorRamp2(breaks=c(-1.5,0,1,1.5), colors = c("darkblue","white", "red1","red4")),  
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(mat_scor_agg_pval[i, j] < 0.0001) {
            grid.text("****", x, y)
          } else if(mat_scor_agg_pval[i, j] < 0.001) {
            grid.text("***", x, y)
          } else if(mat_scor_agg_pval[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if(mat_scor_agg_pval[i, j] < 0.05) {
            grid.text("*", x, y)
          } else {
            grid.text("", x, y)}},
        row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list( # Explicit ticks for the color legend
          labels = c("-1.5", "-1", "0", "1", "1.5")
        ))
dev.off()


