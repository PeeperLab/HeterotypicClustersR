##################################################################
# General analysis of T cell populations
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds from Tcell_Analysis_Code.R
# # subfolders in data/sequence_data/

#R 4.3.3 
#setwd("/YOUR/PATH/")
library(Seurat)#Version 4.4.0 
library(ggplot2)
library(dplyr)
library(scRepertoire)

fig_path = "Results_Tcells_Plots"
if(!dir.exists(fig_path)) dir.create(fig_path) 

tcell.combined<-readRDS("Results_Tcells/Tcells_Final.Rds")

##########################################
#Dotplot with defined t cell states
##########################################

gene_panel<-c("SELL","IL7R","CCR7",#(Naïve/memory makers)
              "LEF1", "TCF7",  "ZNF683", 
              "GZMK", "GZMH", "GZMM", "GZMB", "GZMA", "PRF1", "NKG7", "XCL1", "XCL2", #(Effector cytokines)
              "KLRC3","KIR2DL1",#Tem NK Like
              "KLRB1", "ZBTB16", # (MAIT markers))# 
              "IFIT1","STAT1",  "ISG15",# (Interferon stimulated genes)
              "TOX","PDCD1", "LAG3", "ENTPD1", "HAVCR2", "TIGIT", "CTLA4", "TNFRSF9", "CXCL13",#(Exhaustion or inhibitory molecules)
              "MKI67", "TOP2A"# (transcription factors),
              
)

tcell.combined$annotated_clusters_labels_final = factor(tcell.combined$annotated_clusters_labels_final, levels=c("MKI67+ Tem-NK like/Tprol", "MKI67 hi Tex/Tprol", "MKI67+ Tex/Tprol", "LAG3 hi Tex", "TOX hi Tex", 
                                                                                                                               "TCF7+ stem-like Tex", "GZMK hi Tex", "ISG+", "Tc17 MAIT", "Tem-NK like", "Tem", 
                                                                                                                               "Early Tem", "Tn/Tm", "Tn"))

h<-DotPlot(tcell.combined, features = gene_panel, group.by = "annotated_clusters_labels_final", 
           dot.scale = 8, col.min = -2, col.max = 2) + ylab(" ")+xlab(" ")+
  scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
  theme(axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.001, "lines"), # Reduce space further
        axis.text.y = element_text(size = 16)) + 
  theme(text = element_text(size = 18), legend.title = element_text(size = 15))

pdf(paste0("Results_Tcells_Plots/","Dotplot_Tcell_Markers.pdf"),width=12,height=6)
h
dev.off()

##########################################
#Dotplot with defined t cell exhausted states
##########################################

features  <-c( "GZMA", "GZMH", "GZMK","IL32","TCF7",  "XCL1", "XCL2", "CXCL13", "TOX", "FYN", "HAVCR2", "ENTPD1","LAG3", "PDCD1","GZMB", "S100A4", "MKI67", "TOP2A")

tcell.combined_cluster<-subset(x = tcell.combined, subset = annotated_clusters_labels_final %in% c("TOX hi Tex","TCF7+ stem-like Tex","LAG3 hi Tex","GZMK hi Tex","MKI67+ Tex/Tprol", "MKI67 hi Tex/Tprol" ))

tcell.combined_cluster$annotated_clusters_labels_final = factor(tcell.combined_cluster$annotated_clusters_labels_final, levels=c("MKI67 hi Tex/Tprol","MKI67+ Tex/Tprol","LAG3 hi Tex","TOX hi Tex","TCF7+ stem-like Tex", "GZMK hi Tex" ))

h<-DotPlot(tcell.combined_cluster, features = features, group.by = "annotated_clusters_labels_final", 
           dot.scale = 8, col.min = -2, col.max = 2) + ylab(" ")+xlab(" ")+
  scale_color_gradient2(high = "red4", low = "blue4", mid = "lightyellow2") + 
  theme(axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.001, "lines"), # Reduce space further
        axis.text.y = element_text(size = 16)) + 
  theme(text = element_text(size = 18), legend.title = element_text(size = 15))

pdf(paste0("Results_Tcells_Plots/","Dotplot_Tex_Markers.pdf"),width=8.5,height=4)
h
dev.off()



#######################################################
#plots of t cell states 
#######################################################


##Phenotypes distribution between Clusters and Singlets##############################################################################



umap_coo = as.data.frame(tcell.combined@reductions$umap@cell.embeddings)
umap_coo_meta = merge(umap_coo, tcell.combined@meta.data,by=0)


color_combo<-c("Tn" = "#85965f",
               "Tn/Tm" =	"#3b8e87"	,
               "Early Tem"=	"#EE9572"	,
               "Tem"=	"#b35f42"	,##8B4513
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

#Aggregate the t cell states per patient and samples

data_total<-data.frame()
for(pt in unique(umap_coo_meta$patient_id)){
  for(sd in unique(umap_coo_meta$sample_id)){
    data <- umap_coo_meta %>%
      filter(patient_id == pt, sample_id == sd) %>%
      group_by(patient_id, sample_id, annotated_clusters_labels_final) %>%
      dplyr::summarise(n = n()) %>%
      group_by(patient_id, sample_id) %>%
      mutate(freq = n / sum(n))
    data_total<-rbind(data_total,data)
    
  }
  
}
data_total$annotated_clusters_labels_final = factor(data_total$annotated_clusters_labels_final, levels=c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol", "MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol" ))
data_total$sample_id =factor(data_total$sample_id,levels=c("SG_Singlets","DB_Tumor_Tcell","DB_APC_Tcell"))

data_total <- data_total %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Single T cells",
    sample_id == "DB_Tumor_Tcell" ~ "T cells from\n tumor clusters",
    sample_id == "DB_APC_Tcell" ~ "T cells from\n APC clusters"
  ))

data_total$clusters = factor(data_total$clusters, levels=c("Single T cells","T cells from\n tumor clusters","T cells from\n APC clusters"))
data_total <- data_total %>%
  mutate(patients = case_when(
    patient_id == "P1" ~ "Patient 1",
    patient_id == "P2" ~ "Patient 2",
    patient_id == "P3" ~ "Patient 3",
    patient_id == "P4" ~ "Patient 4",
    patient_id == "P5" ~ "Patient 5"
  ))
data_total$patients = factor(data_total$patients, levels=c("Patient 1","Patient 2","Patient 3","Patient 4","Patient 5"))

pdf(paste0("Results_Tcells_Plots/","Tcells_states_clusters_distribution_per_patient.pdf"),height=8,width=15)
ggplot(data_total) + 
  geom_bar(aes(x = clusters, y = freq, fill = annotated_clusters_labels_final), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + ylab("Frequency")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=16),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(name ="CD8+ T cell states:",values = color_combo) + facet_grid(.~patients)+theme(text=element_text(size=16))
dev.off()
#write.csv(data_total,"Results_Tcells_Plots/Tcells_states_clusters_distribution_per_patient.csv")
##All Patient together phenotype distribution

mean_freq = data_total %>% 
  group_by(clusters,annotated_clusters_labels_final) %>% 
  summarize(Total_freq=(sum(freq)/5))

pdf(paste0("Results_Tcells_Plots/","Tcell_cellstates_distribution_across_clusters_and_singlets.pdf"),height=7,width=7)
ggplot(mean_freq) + 
  geom_bar(aes(x = clusters, y = Total_freq, fill = annotated_clusters_labels_final), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + ylab("Frequency")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=20),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(name ="CD8+ T cell states:",values = color_combo) 
dev.off()
#write.csv(mean_freq,"Results_Tcells_Plots/Tcell_cellstates_distribution_across_clusters_and_singlets.csv")

########### TCR Clonality distribution Among clusters ##############################################################################
#Here Initially we 
get_barcodes_md = function(md,patient,sample_types){
  save_barcodes = c()
  for (s_x in sample_types){
    md_sel = md[md$patient_id == patient & md$sample_id == s_x,]
    md_sel_barcode = rownames(md_sel[1,])
    md_sel_barcode = strsplit(md_sel_barcode,"__")[[1]][1]
    save_barcodes = c(save_barcodes,md_sel_barcode)
  }
  return(save_barcodes)
}


# patient 1
p1_S1_TTD_VDJ <- read.csv("S1_TTD_VDJ/filtered_contig_annotations.csv") #Tcell
p1_S2_TDD_VDJ <- read.csv("S2_TDD_VDJ/filtered_contig_annotations.csv") #APC
p1_S3_STCD_VDJ <- read.csv("S3_STCD_VDJ/filtered_contig_annotations.csv") #singlets

# patient 2
p2_S1_DB_2 <- read.csv("S1_DB_2/clone.S1_DB_2_TT.csv") #Tcell
p2_S1_DB_2 <- p2_S1_DB_2[,colnames(p2_S1_DB_2) != c("X")]
p2_S1_DB_2_2 <- read.csv("S1_DB_2/clone.S1_DB_2_AT.csv") #APC
p2_S1_DB_2_2 <- p2_S1_DB_2_2[,colnames(p2_S1_DB_2_2) != c("X")]
p2_S3_SG_2 <- read.csv("S3_SG_2/filtered_contig_annotations.csv") #singlets
p2_S3_SG_2 <- p2_S3_SG_2[,colnames(p2_S3_SG_2) != c("X")]

# patient 3
p3_DB_3 <- read.csv("DB_3/clone.S1_DB_3_TT.csv") #Tcell
p3_DB_3 <- p3_DB_3[,colnames(p3_DB_3) != c("X")]
p3_DB_3_2 <- read.csv("DB_3/clone.S1_DB_3_AT.csv") #APC
p3_DB_3_2 <- p3_DB_3_2[,colnames(p3_DB_3_2) != c("X")]
p3_SG_3 <- read.csv("SG_3/filtered_contig_annotations.csv") #singlets
p3_SG_3 <- p3_SG_3[,colnames(p3_SG_3) != c("X")]

# patient 4
p4_DB1_4 <- read.csv("DB1_4/filtered_contig_annotations.csv") #Tcell
p4_DB2_4 <- read.csv("DB2_4/filtered_contig_annotations.csv") #APC
p4_SG_4 <- read.csv("SG_4/filtered_contig_annotations.csv") #singlets

p5_DB1_5 <- read.csv("S5_TT_VDJ/filtered_contig_annotations.csv") #Tcell
p5_DB2_5 <- read.csv("S5_AT_VDJ/filtered_contig_annotations.csv") #APC
p5_SG_5 <- read.csv("S5_SG_VDJ/filtered_contig_annotations.csv") #singlets


#seurat_combined<-tcell.combined_subset
sample_types = c( "DB_Tumor_Tcell","DB_APC_Tcell","SG_Singlets" )
md <- tcell.combined@meta.data


# fix barcodes to match the scRNA data and scTCR seq data
barcodes_p1 = get_barcodes_md(md, "P1", sample_types)
p1_S1_TTD_VDJ$barcode <- paste0(barcodes_p1[1],"__",p1_S1_TTD_VDJ$barcode)
p1_S2_TDD_VDJ$barcode <- paste0(barcodes_p1[2],"__",p1_S2_TDD_VDJ$barcode)
p1_S3_STCD_VDJ$barcode <- paste0(barcodes_p1[3],"__",p1_S3_STCD_VDJ$barcode)

barcodes_p2 = get_barcodes_md(md, "P2", sample_types)
p2_S1_DB_2$barcode <- paste0(barcodes_p2[1],"__",p2_S1_DB_2$barcode)
p2_S1_DB_2_2$barcode <- paste0(barcodes_p2[2],"__",p2_S1_DB_2_2$barcode)
p2_S3_SG_2$barcode <- paste0(barcodes_p2[3],"__",p2_S3_SG_2$barcode)

barcodes_p3 = get_barcodes_md(md, "P3", sample_types)
p3_DB_3$barcode <- paste0(barcodes_p3[1],"__",p3_DB_3$barcode)
p3_DB_3_2$barcode <- paste0(barcodes_p3[2],"__",p3_DB_3_2$barcode)
p3_SG_3$barcode <- paste0(barcodes_p3[3],"__",p3_SG_3$barcode)

barcodes_p4 = get_barcodes_md(md, "P4", sample_types)
p4_DB1_4$barcode <- paste0(barcodes_p4[1],"__",p4_DB1_4$barcode)
p4_DB2_4$barcode <- paste0(barcodes_p4[2],"__",p4_DB2_4$barcode)
p4_SG_4$barcode <- paste0(barcodes_p4[3],"__",p4_SG_4$barcode)

barcodes_p5 = get_barcodes_md(md, "P5", sample_types)
p5_DB1_5$barcode <- paste0(barcodes_p5[1],"__",p5_DB1_5$barcode)
p5_DB2_5$barcode <- paste0(barcodes_p5[2],"__",p5_DB2_5$barcode)
p5_SG_5$barcode <- paste0(barcodes_p5[3],"__",p5_SG_5$barcode)


S_all_contig = rbind(p1_S1_TTD_VDJ,p1_S2_TDD_VDJ,p1_S3_STCD_VDJ,
                     p2_S1_DB_2,p2_S1_DB_2_2,p2_S3_SG_2,
                     p3_DB_3,p3_DB_3_2,p3_SG_3,
                     p4_DB1_4,p4_DB2_4,p4_SG_4,
                     p5_DB1_5,p5_DB2_5,p5_SG_5)

tcell.combined@meta.data$barcode_name = do.call(rbind, strsplit(rownames(tcell.combined@meta.data),"__"))[,1]
tcell.combined@meta.data$barcode_name = paste0(tcell.combined@meta.data$barcode_name,"__")




contig.list <- createHTOContigList(S_all_contig, tcell.combined,group.by = "barcode_name")

combined <- combineTCR(contig.list, 
                       samples = unique(tcell.combined@meta.data$combined_id),
                       filterMulti=F,removeNA = TRUE)

# continue with this one and save
save_path = "TCR_filtered_NA"
if(!dir.exists(save_path)) dir.create(save_path) 
for (x in unique(tcell.combined@meta.data$combined_id)){
  write.table(combined[x],paste0(save_path,"/TCR_combined_",x,"_filNA.txt"),sep="\t",quote=F,row.names =F)
}

bar_names = cbind(unique(tcell.combined@meta.data$barcode_name),
                  unique(tcell.combined@meta.data$combined_id))
write.table(bar_names,paste0(save_path,"/TCR_samplenames_barcodes_filNA.txt"),sep="\t",quote=F,row.names =F)

md = tcell.combined@meta.data
md$cell_id = rownames(md)

# add TCR info
md$chain_cdr3 = NA
md$CTnt = NA
md$CTstrict = NA

# fix barcodes and match cells to combine with TCR info
sample_names_table = read.table("./TCR_filtered_NA/TCR_samplenames_barcodes_filNA.txt",header=T)
for (sample_x in sample_names_table$V2){
  sample_tcr_table = read.table(paste0("./TCR_filtered_NA/TCR_combined_",sample_x,"_filNA.txt"),header=T)
  sample_tcr_table[,1] = gsub(paste0(sample_x,"_"),"",sample_tcr_table[,1])
  md[md$cell_id %in% sample_tcr_table[,1],]$chain_cdr3 = sample_tcr_table[,11]
  md[md$cell_id %in% sample_tcr_table[,1],]$CTnt = sample_tcr_table[,10]
  md[md$cell_id %in% sample_tcr_table[,1],]$CTstrict = sample_tcr_table[,12]
}

# only cells with a TCRs
md_tcr <- md[!is.na(md$chain_cdr3),]
# quick save, with cluster names and frequency
write.table(md_tcr,"Results_Tcells_Plots/Meta_data_TCRs_Tcells.txt",quote=F,row.names =F,sep="\t")


md_tcr_freq <- md_tcr %>% 
  group_by(combined_id, chain_cdr3) %>% 
  dplyr::summarize(n=n()) %>% 
  mutate(freq=n/sum(n), total=sum(n))
write.table(md_tcr_freq,"Results_Tcells_Plots/Meta_data_TCRs_Tcells_frequency.txt",quote=F,row.names =F,sep="\t")

set.seed(10)
library(tibble)

## TCR absolute cut off, top most frequent ##

md_tcr$chain_cdr3_pat = paste0(md_tcr$chain_cdr3,"_",md_tcr$patient_id)
md_freq <- md_tcr %>% 
  group_by(patient_id, sample_id, chain_cdr3_pat) %>% 
  dplyr::summarize(n=n()) %>% 
  mutate(freq=n/sum(n), total=sum(n))
md_freq <- md_freq[order(md_freq$freq,decreasing = T),]

md_freq$top_clones <- "none"
for (y in c("P1","P2","P3","P4","P5")){
  for (x in unique(md_freq$sample_id)){
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][1,]$top_clones <- "top_1"
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][2:5,]$top_clones <- "top_2_5"
    md_freq[md_freq$sample_id == x & md_freq$patient_id == y,][6:15,]$top_clones <- "top_6_15"
  }
}

md_freq_sum = md_freq %>%
  group_by(patient_id, sample_id, top_clones) %>% 
  mutate(freq_sum=sum(freq), n_sum = sum(n)) %>%
  select(patient_id,sample_id, top_clones, freq_sum, n_sum)
md_freq_sum = md_freq_sum[!duplicated(md_freq_sum),]
md_freq_sum$sample_id = factor(md_freq_sum$sample_id, levels=c("SG_Singlets","DB_Tumor_Tcell","DB_APC_Tcell"))


md_freq_sum <- md_freq_sum %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Single T cells",
    sample_id == "DB_Tumor_Tcell" ~ "T cells from\n tumor clusters",
    sample_id == "DB_APC_Tcell" ~ "T cells from\n APC clusters"
  ))
md_freq_sum$clusters = factor(md_freq_sum$clusters, levels=c("Single T cells","T cells from\n tumor clusters","T cells from\n APC clusters"))
md_freq_sum <- md_freq_sum %>%
  mutate(patients = case_when(
    patient_id == "P1" ~ "Patient 1",
    patient_id == "P2" ~ "Patient 2",
    patient_id == "P3" ~ "Patient 3",
    patient_id == "P4" ~ "Patient 4",
    patient_id == "P5" ~ "Patient 5"
  ))
md_freq_sum$patients = factor(md_freq_sum$patients, levels=c("Patient 1","Patient 2","Patient 3","Patient 4","Patient 5"))

md_freq_sum <- md_freq_sum %>%
  mutate(combined_tcr = paste(top_clones, sample_id, sep = "_"))


custom_colors <- c("top_1_SG_Singlets" = "#1F78B4", "top_1_DB_Tumor_Tcell" = "#6A3D9A", "top_1_DB_APC_Tcell" = "#009052",
                   "top_2_5_SG_Singlets" = "#6EB8E1", "top_2_5_DB_Tumor_Tcell" = "#9A78B8", "top_2_5_DB_APC_Tcell" = "#00b469",
                   "top_6_15_SG_Singlets" = "#A6CEE3", "top_6_15_DB_Tumor_Tcell" = "#CAB2D6", "top_6_15_DB_APC_Tcell" = "#00d691",
                   "none_SG_Singlets"="grey","none_DB_Tumor_Tcell"="grey","none_DB_APC_Tcell"="grey")


pdf(file=paste0("Results_Tcells_Plots/","Top_TCRs_across_clusters_per_patient_Supps_tcells1_H.pdf"), width=12.8, height=8)
ggplot(md_freq_sum,aes(x=clusters,y=freq_sum,fill=combined_tcr)) + 
  geom_bar(stat="identity") +
  theme_bw()  + ylab("Frequency")+ xlab(" ")+
  scale_fill_manual(values=custom_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=16),
        legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+
  facet_grid(.~patients) 
dev.off()
write.csv(md_freq_sum,"Results_Tcells_Plots/Top_TCRs_across_clusters_per_patient_Supps_tcells1_H.csv")


md_freq_all = md_freq_sum %>%
  group_by(sample_id, top_clones) %>% 
  mutate(freq_all=mean(freq_sum))%>%
  select(sample_id, top_clones, freq_all)

md_freq_all = md_freq_all[!duplicated(md_freq_all),]
md_freq_all$sample_id = factor(md_freq_all$sample_id, levels=c("SG_Singlets","DB_Tumor_Tcell","DB_APC_Tcell"))



md_freq_all <- md_freq_all %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Single T cells",
    sample_id == "DB_Tumor_Tcell" ~ "T cells from\n tumor clusters",
    sample_id == "DB_APC_Tcell" ~ "T cells from\n APC clusters"
  ))
md_freq_all$clusters = factor(md_freq_all$clusters, levels=c("Single T cells","T cells from\n tumor clusters","T cells from\n APC clusters"))

md_freq_all <- md_freq_all %>%
  mutate(combined_tcr = paste(top_clones, sample_id, sep = "_"))
custom_colors <- c("top_1_SG_Singlets" = "#1F78B4", "top_1_DB_Tumor_Tcell" = "#6A3D9A", "top_1_DB_APC_Tcell" = "#009052",
                   "top_2_5_SG_Singlets" = "#6EB8E1", "top_2_5_DB_Tumor_Tcell" = "#9A78B8", "top_2_5_DB_APC_Tcell" = "#00b469",
                   "top_6_15_SG_Singlets" = "#A6CEE3", "top_6_15_DB_Tumor_Tcell" = "#CAB2D6", "top_6_15_DB_APC_Tcell" = "#00d691",
                   "none_SG_Singlets"="grey","none_DB_Tumor_Tcell"="grey","none_DB_APC_Tcell"="grey")

pdf(paste0("Results_Tcells_Plots/","top_TCRs_across_clusters.pdf"), width=2.6, height=5)
ggplot(md_freq_all,aes(x=clusters,y=freq_all,fill=combined_tcr)) + 
  geom_bar(stat="identity") +
  theme_bw() +  xlab(" ")+ ylab("Frequency")+
  scale_fill_manual(values=custom_colors) + theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text=element_text(size=12),
        legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

#write.csv(md_freq_all,"Results_Tcells_Plots/top_TCRs_across_clusters.csv")


###Patients Distribution  in Phenotypes ###########################################################################
md_freq_bar = umap_coo_meta %>% 
  group_by(patient_id,annotated_clusters_labels_final) %>% 
  dplyr::summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
md_freq_bar$annotated_clusters_labels_final = factor(md_freq_bar$annotated_clusters_labels_final, levels=rev(c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol",
                                                                                                               "MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol" )))
md_freq_bar <- md_freq_bar %>%
  mutate(patients = case_when(
    patient_id == "P1" ~ "Patient 1",
    patient_id == "P2" ~ "Patient 2",
    patient_id == "P3" ~ "Patient 3",
    patient_id == "P4" ~ "Patient 4",
    patient_id == "P5" ~ "Patient 5"
  ))
md_freq_bar$patients = factor(md_freq_bar$patients, levels=c("Patient 1","Patient 2","Patient 3","Patient 4","Patient 5"))

color_combo<-c("Patient 1" = "#440154",
               "Patient 2" =	"#3b528b"	,
               "Patient 3"=	"#21918c"	,
               "Patient 4"=	"#5ec962"	,
               "Patient 5"="#fde725")

pdf(paste0("Results_Tcells_Plots/","Tcells_patient_distribution.pdf"),height=5,width=6.5)
ggplot(md_freq_bar) + 
  geom_bar(aes(x = n, y = annotated_clusters_labels_final, fill = patients), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + xlab("Counts")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size=12),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Patients:",values = color_combo)+ylab(" ")+xlab("Counts")
dev.off()
write.csv(md_freq_bar,"Results_Tcells_Plots/Tcells_patient_distribution.csv")

###Clusters Distribution  in CD8+ T cell Phenotypes ###########################################################################

md_freq_bar_sample = umap_coo_meta %>% 
  group_by(sample_id,annotated_clusters_labels_final) %>% 
  dplyr::summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
md_freq_bar_sample$annotated_clusters_labels_final = factor(md_freq_bar_sample$annotated_clusters_labels_final, levels=rev(c("Tn","Tn/Tm","Early Tem","Tem","Tem-NK like","Tc17 MAIT","ISG+","GZMK hi Tex","TCF7+ stem-like Tex","TOX hi Tex","LAG3 hi Tex","MKI67+ Tex/Tprol",
                                                                                                                             "MKI67 hi Tex/Tprol","MKI67+ Tem-NK like/Tprol" )))
md_freq_bar_sample <- md_freq_bar_sample %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Single T cells",
    sample_id == "DB_Tumor_Tcell" ~ "T cells from tumor clusters",
    sample_id == "DB_APC_Tcell" ~ "T cells from APC clusters"
  ))
md_freq_bar_sample$clusters = factor(md_freq_bar_sample$clusters, levels=c("Single T cells","T cells from tumor clusters","T cells from APC clusters"))


color_combo_sample<-c("Single T cells" = "#509abd",
                      "T cells from tumor clusters" =	"#825d9e"	,
                      "T cells from APC clusters"=	"#009052"	)
pdf(paste0("Results_Tcells_Plots/","Tcells_sample_distributions.pdf"),height=5,width=7)
ggplot(md_freq_bar_sample) + 
  geom_bar(aes(x = n, y = annotated_clusters_labels_final, fill = clusters), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + xlab("Counts")+ylab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text = element_text(size=12),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Clusters:",values = color_combo_sample) 
dev.off()
write.csv(md_freq_bar_sample,"Results_Tcells_Plots/Tcells_sample_distributions.csv")
tcell.combined@meta.data$sample_final[tcell.combined@meta.data$sample_id==c("SG_Singlets")]<-"Single T cells"
tcell.combined@meta.data$sample_final[tcell.combined@meta.data$sample_id==c("DB_Tumor_Tcell")]<-"T cells from tumor clusters"
tcell.combined@meta.data$sample_final[tcell.combined@meta.data$sample_id==c("DB_APC_Tcell")]<-"T cells from APC clusters"



color_combo_sample<-c("Single T cells" = "#509abd",
                      "T cells from tumor clusters" =	"grey"	,
                      "T cells from APC clusters"=	"grey"	)

pdf(paste0("Results_Tcells_Plots/","Umap_Tcells_singlets.pdf"), width = 16, height = 13)
DimPlot(tcell.combined, reduction = "umap", group.by = "sample_final",order=c("Single T cells","T cells from tumor clusters","T cells from APC clusters"), pt.size = 4, cols = color_combo_sample) +
  guides(color = guide_legend(override.aes = list(size = 5), reverse = TRUE)) + 
  theme(
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 20)
  )+ggtitle(" ")
dev.off()

color_combo_sample<-c("Single T cells" = "grey",
                      "T cells from tumor clusters" =	"#825d9e"	,
                      "T cells from APC clusters"=	"grey"	)

pdf(paste0("Results_Tcells_Plots/","Umap_Tcells_tumor_clusters.pdf"), width = 16, height = 13)
DimPlot(tcell.combined, reduction = "umap", group.by = "sample_final",order=c("T cells from tumor clusters","Single T cells","T cells from APC clusters"), pt.size = 4, cols = color_combo_sample) +
  guides(color = guide_legend(override.aes = list(size = 5), reverse = TRUE)) + 
  theme(
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 20)
  )+ggtitle(" ")
dev.off()

color_combo_sample<-c("Single T cells" = "grey",
                      "T cells from tumor clusters" =	"grey"	,
                      "T cells from APC clusters"=	"#009052"	)

pdf(paste0("Results_Tcells_Plots/","Umap_Tcells_APC_clusters.pdf"), width = 16, height = 13)
DimPlot(tcell.combined, reduction = "umap", group.by = "sample_final",order=c("T cells from APC clusters","Single T cells","T cells from tumor clusters"), pt.size = 4, cols = color_combo_sample) +
  guides(color = guide_legend(override.aes = list(size = 5), reverse = TRUE)) + 
  theme(
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 20)
  )+ggtitle(" ")
dev.off()

