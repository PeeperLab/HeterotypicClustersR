library(data.table)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(AUCell)
library(harmony)
library(readxl)

fig_path = "output_s2_clean"
if(!dir.exists(fig_path)) dir.create(fig_path)

five_pat=readRDS("Results_barras_CD39/new_reprocessed_five_pat_cd39.rds")
tcr_data=fread("TCR_filtered_five_pats.txt")
p3_cd39_tcrs=fread("CD39/clone.CD39_P3.csv")
p5_cd39_tcrs=fread("CD39/clone.CD39_P5.csv")
coketa_auc=fread("Results_barras_CD39/new_auc_scored_barras.txt")
all_sigs=readxl::read_xlsx("extdata/tcell_signatures_of_interest.xlsx")

meta=data.frame(five_pat@meta.data)
meta$cell_id=gsub("\\-",".",rownames(meta))

# 1/ Classify cells as CD39+/-
cd39_mat=t(data.frame(five_pat@assays$RNA@data[rownames(five_pat@assays$RNA@data) %in% c("ENTPD1", "CD69"),]))
cd39_mat=data.frame(cd39_mat)
cd39_mat$cell_id=rownames(cd39_mat)
auc_entpd1=merge(cd39_mat,meta,by="cell_id")

city_neighborhood<-FindNeighbors(five_pat,return.neighbor=TRUE)
society=rbindlist(lapply(colnames(five_pat), function(x){
  cul_de_sac<-TopNeighbors(city_neighborhood[["RNA.nn"]],cell=x,n=10)
  cul_de_sac_dt=data.table(cell=rep(x, n=10), neighbors=cul_de_sac)
  return(cul_de_sac_dt)}))

setDT(cd39_mat)
setnames(cd39_mat, 'cell_id', "neighbors")
society_cd39=merge(society[, neighbors:= gsub("-", "\\.", neighbors)], cd39_mat[,.(neighbors,ENTPD1)], by="neighbors")
society_cd39[, num_pos_neigh:= sum(ENTPD1>0.35), by="cell"] # bimodal threshold
society_cd39=society_cd39[, lapply(.SD,mean), .SDcols = c('ENTPD1',"num_pos_neigh"), by="cell"]
society_cd39[, cell_id:= gsub("-", "\\.", cell)]
setnames(society_cd39, "ENTPD1", "mean_ENTPD1")

auc_entpd1_n=merge(auc_entpd1,society_cd39,by="cell_id")
setDT(auc_entpd1_n)
auc_entpd1_n[, CD39_selection:= ifelse(ENTPD1<=.35 & mean_ENTPD1<=.38, "CD39neg", "CD39pos")] # bimodal thresholds

gene.plot=ggplot(auc_entpd1_n, aes(ENTPD1))+geom_density(fill="grey40")+
  geom_vline(xintercept = .35)+theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm") )+
  coord_flip()+xlab("")

gene.avg.plot=ggplot(auc_entpd1_n, aes(mean_ENTPD1))+geom_density(fill="grey40")+
  geom_vline(xintercept = .38)+theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm") )+
  xlab("")

scatter_plot=ggplot(auc_entpd1_n, aes(mean_ENTPD1,ENTPD1))+
  geom_point(aes(col=CD39_selection), alpha=.2, size=1)+
  geom_vline(xintercept = .38)+
  geom_hline(yintercept = .35)+theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"))+
  scale_color_manual(values=c("CD39neg"="pink2",
                              "CD39pos"="grey40"))

empty <- ggplot() + theme_void()

combined_plot <- cowplot::plot_grid(
  gene.avg.plot,empty, 
  scatter_plot,gene.plot,  
  ncol = 2, 
  nrow = 2,
  align = "hv",
  rel_widths = c(3, 1),
  rel_heights = c(1, 3))

pdf("output_s2_clean/new_cd39_neg_population_selection.pdf", height=5, width=5)
combined_plot
dev.off()

# load others
color_vector_2 <- c(
  "Tn" = "#85965f",
  "Tn/Tm" = "#3b8e87",
  "Early Tem" = "#EE9572",
  "Tem" = "#b35f42",
  "Tem-NK like" = "#CD3700",
  "Tc17 MAIT" = "#d8d5b2",
  "Mem_like"="#85965f",
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

auc_entpd1_n[, sample_id:= gsub("P3_|P5_","",sample_id)]
auc_entpd1_n[, s_id:= factor(gsub("APC_|Tumor_","",sample_id),
                             levels=c("SG_Singlets", "DB_Tcell", "CD39"))]

auc_entpd1_n[,final_clusters:=factor(final_clusters,levels=c("MKI67 hi Tex/Tprol","MKI67+ Tex/Tprol",
                                                             "TCF7+ stem-like Tex", "LAG3 hi Tex","TOX hi Tex",
                                                             "CD137 hi Tex","GZMK hi Tex","ISG+","HSP",
                                                             "Tc17 MAIT","Tem","Early Tem",
                                                             "Tn/Tm", "Tn"))] # change order of TCf7 and LAG3

# 2- Include TCR data

tcr_data[, pat_TCR:= paste0(sample, "_", TCR )]
tcr_data[, patient_id:= NULL]
tcr_data[, barcode:= paste0(sample, "__", stringr::str_split_i(barcode, "__",2))]

auc_entpd1_n[, barcode:= paste0(combined_id, "__", barcode)]

d_p3_cd39_tcrs=dcast(p3_cd39_tcrs, barcode~chain, value.var = "cdr3", fun.aggregate = function(x){paste0(x, collapse = ";")})
d_p3_cd39_tcrs[, sample:= "P3_CD39"]
d_p3_cd39_tcrs[, pat_TCR:= paste(sample,TRA,TRB, sep="_")]
d_p3_cd39_tcrs[, barcode:= paste0(sample, "__", barcode)]

d_p5_cd39_tcrs=dcast(p5_cd39_tcrs, barcode~chain, value.var = "cdr3", fun.aggregate = function(x){paste0(x, collapse = ";")})
d_p5_cd39_tcrs[, sample:= "P5_CD39"]
d_p5_cd39_tcrs[, pat_TCR:= paste(sample,TRA,TRB, sep="_")]
d_p5_cd39_tcrs[, barcode:= paste0(sample, "__", barcode)]

d_p3_cd39_tcrs[, TCR:=paste0(TRA,"_",TRB)]
d_p5_cd39_tcrs[, TCR:=paste0(TRA,"_",TRB)]

tcr_data_new=rbindlist(list(tcr_data[,.(barcode,pat_TCR,TCR)],
                            d_p3_cd39_tcrs[,.(barcode,pat_TCR,TCR)],
                            d_p5_cd39_tcrs[,.(barcode,pat_TCR,TCR)]))

tcr_data_new[, pat_TCR:= gsub("_APC|_Tumor","", pat_TCR)]

auc_tcr=merge(auc_entpd1_n, tcr_data_new, by="barcode")
auc_tcr[, simple_sample:= ifelse(grepl("DB", sample_id), "DB_Tcell", sample_id)]
auc_tcr[, count_complex_TCR:= .N, by=c("pat_TCR")]
auc_tcr[, sum_CD39pos:= sum(CD39_selection=="CD39pos"), by=c("pat_TCR")]
auc_tcr[, sum_CD39neg:= sum(CD39_selection=="CD39neg"), by=c("pat_TCR")]


mean_error=mean(c(nrow(auc_tcr[s_id=="CD39"& patient_id=="P3"& CD39_selection=="CD39neg"])/nrow(auc_tcr[s_id=="CD39"& patient_id=="P3"]),
                  nrow(auc_tcr[s_id=="CD39"& patient_id=="P5"& CD39_selection=="CD39neg"])/nrow(auc_tcr[s_id=="CD39"& patient_id=="P5"])))

auc_tcr_f=auc_tcr[,.(patient_id, simple_sample, sum_CD39neg, sum_CD39pos,
                     count_complex_TCR, pat_TCR, TCR)]
auc_tcr_f=auc_tcr_f[!duplicated(auc_tcr_f),]

setorder(auc_tcr_f, -count_complex_TCR)
exp_test=auc_tcr_f[count_complex_TCR>=10,]

# 3- Add signatures and parametrise expansion levels
exprMatrix= five_pat@assays$RNA@counts
exprMatrix <- as(exprMatrix, "dgCMatrix")
set.seed(333)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)

all_sigs=lapply(all_sigs, function(x){return(x)})
geneSets=all_sigs
geneSets=lapply(geneSets, function(x){return(x[which(!is.na(x))])})
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc=data.frame(t(data.frame(cells_AUC@assays@data$AUC)))
auc$cell_id=rownames(auc)
auc_tcr_new=merge(auc_tcr, auc, by="cell_id", all.x=T)

auc_tcr_new[, expansion:= factor(ifelse(count_complex_TCR>=10, "N>=10",
                                        ifelse(count_complex_TCR>=3, "3<=N<10",
                                               "N<3")), levels=c("N<3", "3<=N<10", "N>=10"))]
auc_tcr_new[, expansion_bin:= ifelse(count_complex_TCR>=10, "Expanded", "NonExp")]

# 4- Are there CD39- DBs? Is this population in reactive (expanded) clonotypes?
by_expansion=auc_tcr_new[,.N, by=c("expansion", "s_id", "patient_id","CD39_selection")]

by_expansion[,tot:=sum(N), by=c("expansion", "s_id", "patient_id")]
by_expansion[, perc_neg:= 100*N/tot]
blank_0s = CJ(expansion=unique(by_expansion$expansion),
              patient_id=unique(by_expansion$patient_id),
              CD39_selection=unique(by_expansion$CD39_selection),
              s_id=unique(by_expansion$s_id))
by_expansion=merge(blank_0s,by_expansion, all.x=T, by=c("expansion","patient_id","s_id","CD39_selection"))
by_expansion[is.na(by_expansion)]=0

by_expansion[, tot_tot:= sum(N), by=c("patient_id", "s_id")]
by_expansion[, tot_freq:= N/tot_tot]
by_expansion[, exp_cd39:= factor(paste(expansion,CD39_selection,sep = "_"),
                                 levels = c("N>=10_CD39pos","3<=N<10_CD39pos","N<3_CD39pos",
                                            "N>=10_CD39neg","3<=N<10_CD39neg","N<3_CD39neg"))]

ggplot(by_expansion[s_id=="DB_Tcell"], aes(patient_id, y=tot_freq, fill=exp_cd39))+
  geom_col(position = "stack")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid= element_blank())+
  scale_fill_manual(values=c("N<3_CD39neg"="pink1",
                             "N<3_CD39pos"="grey95",
                             "3<=N<10_CD39neg"="pink3",
                             "3<=N<10_CD39pos"="grey85",
                             "N>=10_CD39neg"="pink4",
                             "N>=10_CD39pos"="grey75"))+
  xlab("")+geom_hline(yintercept=mean_error, linetype = 'dashed')+ylab("Freq")
ggsave("output_s2_clean/S11B_clus_cell_distribution_CD39pos_neg_expansion.pdf", width=4, height=5)
write.csv(by_expansion[s_id=="DB_Tcell"], "output_s2_clean/S11b.csv")

export_totals=by_expansion[s_id=="DB_Tcell"][, .(s_id,patient_id, CD39_selection, expansion,
                                                 group_size=N,pat_total=tot_tot, group_pat_freq=tot_freq)]
setorder(export_totals, -patient_id, -expansion, -CD39_selection)
export_totals[, CD39_total:= sum(group_size), by="CD39_selection"]
fwrite(export_totals,"output_s2_clean/S11B_ordered.csv")

ggplot(by_expansion[s_id=="DB_Tcell"& CD39_selection=="CD39neg"], aes(patient_id, y=tot_freq, fill=exp_cd39))+
  geom_col(position = "fill")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid= element_blank())+
  scale_fill_manual(values=c("N<3_CD39neg"="pink1",
                             "N<3_CD39pos"="grey95",
                             "3<=N<10_CD39neg"="pink3",
                             "3<=N<10_CD39pos"="grey85",
                             "N>=10_CD39neg"="pink4",
                             "N>=10_CD39pos"="grey75"))+
  xlab("")+ylab("Freq")
ggsave("output_s2_clean/FigS11b_wihtin_CD39_neg.pdf", width=5, height=5)

selection=by_expansion[s_id=="DB_Tcell"& CD39_selection=="CD39neg"]
selection[, tot_CD39neg:= sum(N), by="patient_id"]
selection[, freq_within_CD39neg:= N/tot_CD39neg]
selection[, .(s_id,patient_id, CD39_selection, expansion,
              group_size=N,tot_CD39neg,freq_within_CD39neg)]
fwrite(selection,"output_s2_clean/S11B_within_neg.csv")


# 5- How much are the TCF7+ and LAG3+ dependent on this CD39- cells?

# # Across clusters + singlets
f_auc_tcr_new=auc_tcr_new[s_id!="CD39"]
pre_sum2=f_auc_tcr_new[, .N, by=c("patient_id", "CD39_selection", "final_clusters")]

blank_0s = CJ(final_clusters=unique(pre_sum2$final_clusters),
              patient_id=unique(pre_sum2$patient_id),
              CD39_selection=unique(pre_sum2$CD39_selection))
pre_sum2=merge(blank_0s,pre_sum2, all.x=T,
               by=c("final_clusters","patient_id","CD39_selection"))
pre_sum2[is.na(pre_sum2)]=0

pre_sum2[, tot:= sum(N), by=c("patient_id", "final_clusters")]
pre_sum2[, perc:= 100*N/tot]

ggplot(pre_sum2[final_clusters %in% c("TCF7+ stem-like Tex","LAG3 hi Tex") &
                  CD39_selection=="CD39neg"], aes(final_clusters, perc))+
  geom_boxplot(aes(fill=final_clusters))+
  geom_point(aes(colour = patient_id))+
  geom_line(aes(group=patient_id,colour = patient_id))+
  stat_compare_means(method="t", paired=T)+
  ylab("% CD39- cells in cell states")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid= element_blank())+
  scale_fill_manual(values=color_vector_2)
ggsave("output_s2_clean/S11F_left_perc_CD39neg_LAG3_TCF7_clus_sg.pdf", width=5, height=5)
fwrite(pre_sum2[final_clusters %in% c("TCF7+ stem-like Tex","LAG3 hi Tex") &
                  CD39_selection=="CD39neg"], "output_s2_clean/S11F_left.csv")

test_44=dcast(pre_sum2[final_clusters %in% c("TCF7+ stem-like Tex","LAG3 hi Tex") &
                         CD39_selection=="CD39neg"], patient_id~final_clusters)
test_44[, diff:= `TCF7+ stem-like Tex`-`LAG3 hi Tex`]
shapiro.test(test_44$diff)


# 6- Krishna signatures in the clusters

# # Across clusters + Singlets
pheno_auc_tcr_new2=auc_tcr_new[s_id!='CD39', lapply(.SD, mean),
                               .SDcols = c("Krishna_CD39neg.CD69neg"),
                               by=c("final_clusters", "patient_id")]

ggplot(pheno_auc_tcr_new2[final_clusters%in%c("TCF7+ stem-like Tex","LAG3 hi Tex")],
       aes(final_clusters, Krishna_CD39neg.CD69neg))+
  geom_boxplot(aes(fill=final_clusters))+geom_point(aes(col=patient_id))+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid= element_blank())+
  stat_compare_means(method="t", paired=T, size=3)+
  geom_line(aes(group=patient_id,col=patient_id))+
  ylab("CD39-CD69- signature score per patient")+
  scale_fill_manual(values=color_vector_2)
ggsave("output_s2_clean/S11F_right_krishna_sig_LAG3_TCF7_clus_sg.pdf", width=5, height=5)
write.csv(pheno_auc_tcr_new2[final_clusters%in%c("TCF7+ stem-like Tex","LAG3 hi Tex")], 
          "output_s2_clean/S11F_right.csv")

test22=dcast(pheno_auc_tcr_new2[final_clusters%in%c("TCF7+ stem-like Tex","LAG3 hi Tex")], patient_id~final_clusters, value.var = "Krishna_CD39neg.CD69neg")
test22[, diff:= `TCF7+ stem-like Tex`- `LAG3 hi Tex`]
shapiro.test(test22$diff)

# 7- Compare the phenotype presence in matched clonotypes
exp_test_clus=auc_tcr[count_complex_TCR>=10,]
exp_test_clus[, sum_tcf7_tcr:= sum(final_clusters=="TCF7+ stem-like Tex"), by=c("pat_TCR")]
exp_test_clus[, sum_lag3_tcr:= sum(final_clusters=="LAG3 hi Tex"), by=c("pat_TCR")]

exp_test_clus[, freq_tcf7_tcr:= sum_tcf7_tcr/count_complex_TCR]
exp_test_clus[, freq_lag3_tcr:= sum_lag3_tcr/count_complex_TCR]

matched=unique(intersect(intersect(exp_test_clus$TCR[exp_test_clus$s_id=="DB_Tcell"],
                                   exp_test_clus$TCR[exp_test_clus$s_id=="SG_Singlets"]),
                         exp_test_clus$TCR[exp_test_clus$s_id=="CD39"]))

matched_exp_test_clus=exp_test_clus[TCR %in% matched,lapply(.SD, mean),
                                    by=c("count_complex_TCR", "pat_TCR", "TCR","s_id","patient_id"),
                                    .SDcols=c("freq_tcf7_tcr","freq_lag3_tcr")]
m_matched_exp_test_clus=melt(matched_exp_test_clus, id.vars=c("count_complex_TCR", "pat_TCR", "TCR","s_id","patient_id"))

ggplot(m_matched_exp_test_clus[patient_id %in% c("P3","P5")], aes(s_id, value))+
  facet_wrap(~variable, scales="free", nrow=1)+
  geom_boxplot(aes(fill=s_id))+geom_point()+
  geom_line(aes(group=TCR), alpha=.5)+
  stat_compare_means(method="wilcox", paired=T,size=3,
                     comparisons=list(c("DB_Tcell","SG_Singlets"),
                                      c("DB_Tcell","CD39")))+
  scale_fill_manual(values=c("SG_Singlets"="lightblue",
                             "DB_Tcell"="salmon",
                             "CD39"="lightyellow2"))+
  theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
                   panel.grid = element_blank())+
  xlab("")+ylab("Frequency in matched expanded clonotypes (N>=10)")
ggsave("output_s2_clean/S11G_phenotypes_in_matched_CTs.pdf", height=5, width=5)
write.csv(m_matched_exp_test_clus[patient_id %in% c("P3","P5")],"output_s2_clean/S11G.csv" )
table(m_matched_exp_test_clus[patient_id %in% c("P3","P5")]$patient_id,
      m_matched_exp_test_clus[patient_id %in% c("P3","P5")]$s_id,
      m_matched_exp_test_clus[patient_id %in% c("P3","P5")]$variable)

lapply(unique(as.character(m_matched_exp_test_clus[patient_id %in% c("P3","P5")]$variable)),
       function(x){
         print(x)
         shapiro.test(m_matched_exp_test_clus[patient_id %in% c("P3","P5")][s_id %in% c("DB_Tcell")][variable==x]$value - 
                        m_matched_exp_test_clus[patient_id %in% c("P3","P5")][s_id %in% c("CD39")][variable==x]$value)$p.value})

lapply(unique(as.character(m_matched_exp_test_clus[patient_id %in% c("P3","P5")]$variable)),
       function(x){
         print(x)
         shapiro.test(m_matched_exp_test_clus[patient_id %in% c("P3","P5")][s_id %in% c("DB_Tcell")][variable==x]$value - 
                        m_matched_exp_test_clus[patient_id %in% c("P3","P5")][s_id %in% c("SG_Singlets")][variable==x]$value)$p.value})

# 8- response association
coketa_auc=coketa_auc[predicted.celltype.score>0.4,] 

setDT(coketa_auc)
coketa_auc[, response:= ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1)))%in%c(2,3,7,8,9,13), "R","NR")]
coketa_auc[, response2:= ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1)))%in%c(2,3), "CR",
                                ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1)))%in%c(8,9,13,7), "PR",
                                       ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1)))%in%c(1,4,10), "SD",
                                              ifelse(as.numeric(gsub("patient","",stringr::str_split_i(patient_id,"_",1)))%in%c(5,6,11,12), "PD","??"))))]

coketa_auc[,pat:=stringr::str_split_i(patient_id, "_",1)]
coketa_auc[, neotcr_cd8_BIN:= cut(NeoTCR_cd8,
                                  quantile(NeoTCR_cd8, probs = seq(0,1,1/3)),
                                  include.lowest=T,
                                  labels=c("NeoTCR_low","NeoTCR_mid" ,"NeoTCR_high")), by="timepoint"]

coketa_auc[, response:= factor(response, levels=c("R","NR"))]
coketa_auc[, tot_pat_neo:= .N, by=c("pat","neotcr_cd8_BIN","timepoint")]

coketa_auc[, clus_stem_freq_neo:= sum(predicted.celltype=='TCF7+ stem-like Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_stem_freq_neo:= 100*clus_stem_freq_neo/tot_pat_neo]

coketa_auc[, clus_lag3_freq_neo:= sum(predicted.celltype=='LAG3 hi Tex'), by=c("pat","neotcr_cd8_BIN","timepoint")]
coketa_auc[, clus_lag3_freq_neo:= 100*clus_lag3_freq_neo/tot_pat_neo]


avg_coketa_auc=coketa_auc[, lapply(.SD, mean),
                          .SDcols = c('clus_stem_freq_neo',"clus_lag3_freq_neo"),
                          by=c("pat","neotcr_cd8_BIN","timepoint","response","response2")]

avg_coketa_auc[,neotcr_cd8_BIN:=factor(neotcr_cd8_BIN, levels=c("NeoTCR_high","NeoTCR_mid" ,"NeoTCR_low"))]

pdf("output_s2_clean/S11I_transf_tcf7_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_stem_freq_neo/100))+
  geom_boxplot(aes(fill=response), alpha=.5, outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width = .3,seed = 150799))+
  theme_bw()+ylab("Frequency TCF7+ stem-like Tex")+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()
write.csv(avg_coketa_auc[timepoint=="T0"], "output_s2_clean/S11I.csv")

pdf("output_s2_clean/out_FigS11i_transf_lag3_response_tertiles.pdf", width=5.5, height=3.5)
ggplot(avg_coketa_auc[timepoint=="T0"],aes(response, clus_lag3_freq_neo/100))+
  geom_boxplot(aes(fill=response), alpha=.5, outlier.shape = NA)+
  facet_wrap(~neotcr_cd8_BIN)+stat_compare_means(method="t", size=3)+
  geom_point(aes(colour = response2), size=3, alpha=1,position=position_jitter(width=.2))+
  theme_bw()+ylab("Frequency LAG3 high Tex")+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(values=c("R"="#c1de87", "NR"="grey40"))+
  scale_color_manual(values=c("SD"="grey60","PD"="grey10",
                              "PR"="#3b756e","CR"="#6d8754"))
dev.off()

