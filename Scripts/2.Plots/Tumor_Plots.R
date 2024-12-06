##################################################################
# Tumor phenotypes description and enrichment analysis
##################################################################

# required input data in the working directory:
# # Tumor_Annotated.Rds output from Tumor_Analysis_Code.R
# # GSEA folder from data/tumor_data
# # Tumor_Signatures.csv from data/tumor_data

#setwd("/YOUR/PATH/")
library(Seurat) 
library(ggplot2)
library(dplyr)
library(data.table)
library(AUCell)
library(ComplexHeatmap)
library(circlize)

fig_path = "Results_Tumor_Plots"
if(!dir.exists(fig_path)) dir.create(fig_path) 

pbmc.combined<-readRDS("Tumor_Annotated.Rds")


#Remove Low Gene Count Cells group from further analysis
pbmc.combined_subset<-subset(x = pbmc.combined, subset = annotated_clusters %in% c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific","Stress(hypoxia response)","Stress(p53 response)",
                                                                                   "Neural Crest Like"))

color_combo2<-c("Immune Response"="#99CC66",
                "Mitotic"="#006633",
                "Melanocytic"="#990000",
                "Stress(hypoxia response)"="#dac565",
                "Transitory Melanocytic"="#FF3300",
                "Stress(p53 response)"="#0099CC",
                "Neural Crest Like"="#003366",
                "Patient Specific"="#FF9933")


umap_coo_meta<-pbmc.combined_subset@meta.data
data_total<-data.frame()
for(pt in unique(umap_coo_meta$patient_id)){
  for(sd in unique(umap_coo_meta$sample_id)){
    data <- umap_coo_meta %>%
      filter(patient_id == pt, sample_id == sd) %>%
      group_by(patient_id, sample_id, annotated_clusters) %>%
      dplyr::summarise(n = n()) %>%
      group_by(patient_id, sample_id) %>%
      mutate(freq = n / sum(n))
    data_total<-rbind(data_total,data)
    
  }
  
}

data_total <- data_total %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Single tumor\n cells",
    sample_id == "DB_Tumor_Tcell" ~ "Tumor cells \nfrom clusters"
  ))
data_total$clusters = factor(data_total$clusters, levels=c("Single tumor\n cells","Tumor cells \nfrom clusters"))
data_total$annotated_clusters = factor(data_total$annotated_clusters, levels=c("Immune Response",  "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific","Stress(hypoxia response)","Stress(p53 response)",
                                                                               "Neural Crest Like"))

data_total <- data_total %>%
  mutate(patients = case_when(
    patient_id == "P1" ~ "Patient 1",
    patient_id == "P2" ~ "Patient 2",
    patient_id == "P3" ~ "Patient 3",
    patient_id == "P4" ~ "Patient 4",
    patient_id == "P5" ~ "Patient 5"
  ))
data_total$patients = factor(data_total$patients, levels=c("Patient 1","Patient 2","Patient 3","Patient 4","Patient 5"))

pdf("Results_Tumor_Plots/Tumor_Annotated_clusters_distribution_per_patient.pdf",width=13,height=8)
ggplot(data_total) + 
  geom_bar(aes(x = clusters, y = freq, fill = annotated_clusters), 
           stat = "identity", 
           color = NA) +  # Removes borders around the bars
  theme_bw() + ylab("Frequency")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=20),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Tumor phenotypes:",values = color_combo2) +facet_grid(.~patients)
dev.off()


data_combined = data_total %>% 
  group_by(clusters,annotated_clusters) %>% 
  summarize(Total_freq=(sum(freq)/5))

pdf("Results_Tumor_Plots/Tumor_Annotated_clusters_distribution.pdf",height=10,width=9.5)
ggplot(data_combined) + 
  geom_bar(aes(x = clusters, y = Total_freq, fill = annotated_clusters), 
           stat = "identity", 
           color = NA) +  # Removes borders around the bars
  theme_bw() + ylab("Frequency")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=35),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Tumor phenotypes:", values = color_combo2)
dev.off()
########################################################
###TUMOR SIGNATURE ANALYSIS########
########################################################
Tumor_Signatures <- read.csv("Tumor_Signatures.csv", row.names=1)


cell_types <-  unique(Tumor_Signatures$Signature)
cell_type_list <- c()
for (cell_x in cell_types){
  to_sigs_sel <- Tumor_Signatures[Tumor_Signatures$Signature == cell_x,]
  cell_type_list <- c(cell_type_list,list(to_sigs_sel$Gene))
}
names(cell_type_list) <- cell_types


exprMatrix <- as.matrix(pbmc.combined_subset@assays$RNA@counts)


cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(cell_type_list,cell_rank)

output_data <- cells_AUC@assays@data@listData$AUC
output_data <- as.data.frame(t(output_data))

md <- pbmc.combined_subset@meta.data
check_data_plot <- merge(output_data,md,by=0)

check_data_plot_trimmed<-check_data_plot[,c("Rambow Immune", "Pozniak Antigen_presentation",
                                            "Pozniak Mitotic","Rambow mitosis",
                                            "Wouters Melanocytic cell state",  "Tsoi Melanocytic","Rambow MITFtargets",  "Rambow pigmentation", 
                                            "Tsoi Transitory-Melanocytic","Tsoi Transitory",
                                            "Pozniak Stress (p53 response)",
                                            "Rambow Neuro", "Tsoi Neural crest-like","Pozniak Neural_Crest_like",
                                            "Pozniak Stress (hypoxia response)",
                                            "annotated_clusters")]


aggregated_AUCell <- check_data_plot_trimmed %>%
  group_by(annotated_clusters) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

AUCell_matrix <- as.matrix(aggregated_AUCell[, -1])
rownames(AUCell_matrix) <- aggregated_AUCell$annotated_clusters


order<-c("Immune Response" ,"Mitotic" ,"Patient Specific" ,"Melanocytic" ,"Transitory Melanocytic",
         "Stress(p53 response)" ,"Neural Crest Like" ,"Stress(hypoxia response)")
AUCell_matrix_order_row<-AUCell_matrix[order,]

transposed_matrix <- t(AUCell_matrix_order_row)

Sig_order<-c("Rambow Immune", "Pozniak Antigen_presentation",
             "Pozniak Mitotic","Rambow mitosis",
             "Wouters Melanocytic cell state",  "Tsoi Melanocytic","Rambow MITFtargets",  "Rambow pigmentation", 
             "Tsoi Transitory-Melanocytic","Tsoi Transitory",
             "Pozniak Stress (p53 response)",
             "Rambow Neuro", "Tsoi Neural crest-like","Pozniak Neural_Crest_like",
             "Pozniak Stress (hypoxia response)"
)

transposed_matrix_ordered<-transposed_matrix[Sig_order,]

scaled_matrix <- t(scale(t(transposed_matrix_ordered), center = TRUE, scale = TRUE))

ht7 <- Heatmap(
  scaled_matrix, 
  name = "Z Score", 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_names_rot = 45,  
  column_names_gp = gpar(fontsize = 45),  
  row_names_gp = gpar(fontsize = 45),     
  heatmap_width = unit(32, "cm"),         
  heatmap_height = unit(60, "cm"),        
  col = colorRamp2(breaks=c(-2,0,1.5,2), colors = c("darkblue","white", "red1","red4")),  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 25),        
    labels_gp = gpar(fontsize = 20),       
    legend_height = unit(4, "cm"),        
    legend_width = unit(2, "cm")          
  )
)
pdf("Results_Tumor_Plots/HeatMap_Tumor_Signatures.pdf",width=50,height=40)
draw(ht7, heatmap_legend_side = "left")  
dev.off()
######################################################################
###Enriched GSEA pathway Heatmap##########
######################################################################
input_dir <- "./GSEA" 

# List all TSV files in the directory
tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
tsv_files

read_tsv_with_filename <- function(file) {
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  data$filename <- basename(file)
  return(data)
}

# Load all TSV files into a list of data frames
tsv_data_list <- lapply(tsv_files, read_tsv_with_filename)
tsv_data_list<-as.data.frame(tsv_data_list)
# Optionally combine all data frames into one
combined_data <- do.call(rbind, tsv_data_list)
combined_data<-as.data.frame(combined_data)

filtered_combined_data <- combined_data[grepl("^M", combined_data[["V1"]]), ]
filtered_combined_data$name<-rownames(filtered_combined_data)
filtered_combined_data<-filtered_combined_data[,c("V5","V17")]
colnames(filtered_combined_data)<-c("signature","genes")

#Trim some long gene set names
filtered_combined_data$signature <- ifelse(filtered_combined_data$signature == "Genes up-regulated by activation of WNT signaling through accumulation of beta catenin CTNNB1 [GeneID=1499].","Genes up-regulated by activation of WNT signaling", filtered_combined_data$signature)
filtered_combined_data$signature <- ifelse(filtered_combined_data$signature == "Antigen Presentation: Folding, assembly and peptide loading of class I MHC","Antigen Presentation: MHC class I", filtered_combined_data$signature)
filtered_combined_data$signature <- ifelse(filtered_combined_data$signature == "Genes encoding proteins involved in glycolysis and gluconeogenesis.","Genes encoding proteins involved in glycolysis", filtered_combined_data$signature)


split_genes <- function(signature, genes) {
  gene_list <- unlist(strsplit(genes, split = ","))
  data.frame(SignatureName = rep(signature, length(gene_list)), Gene = gene_list, stringsAsFactors = FALSE)
}

# Apply the function to each row and combine the results
expanded_data <- do.call(rbind, apply(filtered_combined_data, 1, function(row) split_genes(row["signature"], row["genes"])))

expanded_data <- expanded_data[!(expanded_data$Gene == "" | is.na(expanded_data$Gene) ), ]

expanded_data <- unique(expanded_data)

cell_types <-  unique(expanded_data$SignatureName)
cell_type_list <- c()
for (cell_x in cell_types){
  to_sigs_sel <- expanded_data[expanded_data$SignatureName == cell_x,]
  cell_type_list <- c(cell_type_list,list(to_sigs_sel$Gene))
}

names(cell_type_list) <- cell_types

exprMatrix <- as.matrix(pbmc.combined_subset@assays$RNA@counts)

cell_rank <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats = F)
cells_AUC <- AUCell_calcAUC(cell_type_list,cell_rank)

output_data <- cells_AUC@assays@data@listData$AUC
output_data <- as.data.frame(t(output_data))

md <- pbmc.combined_subset@meta.data
check_data_plot2 <- merge(output_data,md,by=0)

check_data_plot_trimmed<-check_data_plot2[,c("annotated_clusters","Interferon alpha/beta signaling","PD-1 signaling","Antigen Presentation: MHC class I", "MHC class II antigen presentation",
                                             "Cell cycle","AURKA Activation by TPX2","PKR-mediated signaling" ,"Genes up-regulated through activation of mTORC1 complex.","Genes up-regulated by activation of the PI3K/AKT/mTOR pathway.",
                                             "Melanin biosynthesis",
                                             "Translation","Genes encoding proteins involved in oxidative phosphorylation.",
                                             "Apoptosis" ,"VEGF signaling pathway","Genes up-regulated by activation of WNT signaling" ,"IGF-1 Signaling Pathway",
                                             "ECM-receptor interaction",
                                             "Genes up-regulated in response to low oxygen levels (hypoxia).","Genes encoding proteins involved in glycolysis"  
)
]

aggregated_AUCell <- check_data_plot_trimmed %>%
  group_by(annotated_clusters) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

AUCell_matrix <- as.matrix(aggregated_AUCell[, -1])
rownames(AUCell_matrix) <- aggregated_AUCell$annotated_clusters


order<-c("Immune Response" ,"Mitotic" ,"Patient Specific" ,"Melanocytic" ,"Transitory Melanocytic",
         "Stress(p53 response)" ,"Neural Crest Like" ,"Stress(hypoxia response)")
AUCell_matrix_order_row<-AUCell_matrix[order,]

transposed_matrix <- t(AUCell_matrix_order_row)

Pathway_order<-c("Interferon alpha/beta signaling","PD-1 signaling","Antigen Presentation: MHC class I", "MHC class II antigen presentation",
                 "Cell cycle","AURKA Activation by TPX2","PKR-mediated signaling" ,"Genes up-regulated through activation of mTORC1 complex.","Genes up-regulated by activation of the PI3K/AKT/mTOR pathway.",
                 "Melanin biosynthesis",
                 "Translation","Genes encoding proteins involved in oxidative phosphorylation.",
                 "Apoptosis" ,"VEGF signaling pathway","Genes up-regulated by activation of WNT signaling" ,"IGF-1 Signaling Pathway",
                 "ECM-receptor interaction",
                 "Genes up-regulated in response to low oxygen levels (hypoxia).","Genes encoding proteins involved in glycolysis"  
)

transposed_matrix_ordered<-transposed_matrix[Pathway_order,]

scaled_matrix <- t(scale(t(transposed_matrix_ordered), center = TRUE, scale = TRUE))

ht7 <- Heatmap(
  scaled_matrix, 
  name = "Z Score", 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_names_rot = 45,  
  column_names_gp = gpar(fontsize = 45),  
  row_names_gp = gpar(fontsize = 45),     
  heatmap_width = unit(32, "cm"),         
  heatmap_height = unit(60, "cm"),        
  col = colorRamp2(breaks=c(-2,0,1.5,2), colors = c("darkblue","white", "red1","red4")),  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 25),        
    labels_gp = gpar(fontsize = 20),       
    legend_height = unit(4, "cm"),         
    legend_width = unit(2, "cm")           
  )
)
pdf("Results_Tumor_Plots/HeatMap_Tumor_Pathways.pdf",width=50,height=40)
draw(ht7, heatmap_legend_side = "left")  
dev.off()

###Tumor Clusters/Singlets distribution in UMAP
pbmc.combined@meta.data$source_id[pbmc.combined@meta.data$sample_id == 'SG_Singlets'] <- "Tumor cells from singlets"
pbmc.combined@meta.data$source_id[pbmc.combined@meta.data$sample_id == 'DB_Tumor_Tcell'] <- "Tumor cells from tumor clusters"

color_combo_sample<-c("Tumor cells from singlets" = "#509abd",
                      "Tumor cells from tumor clusters" =	"grey"		)

pdf("Results_Tumor_Plots/Tumor_UMAP_clusters_SG.pdf",width=16,height=13)
DimPlot(pbmc.combined,reduction="umap",group.by="source_id",pt.size=2,order=c("Tumor cells from singlets","Tumor cells from tumor clusters"),cols=color_combo_sample)+labs(title=NULL)+guides(color = guide_legend(override.aes = list(size=6), ncol=1) )+theme(text=element_text(size=16))
dev.off()

color_combo_sample<-c("Tumor cells from singlets" = "grey",
                      "Tumor cells from tumor clusters" =	"#825d9e"		)

pdf("Results_Tumor_Plots/Tumor_UMAP_clusters_DB.pdf",width=16,height=13)
DimPlot(pbmc.combined,reduction="umap",group.by="source_id",pt.size=2,order=c("Tumor cells from tumor clusters","Tumor cells from singlets"),cols=color_combo_sample)+labs(title=NULL)+guides(color = guide_legend(override.aes = list(size=6), ncol=1) )+theme(text=element_text(size=16))
dev.off()


umap_coo_meta<-pbmc.combined_subset@meta.data


md_freq_bar_sample = umap_coo_meta %>% 
  group_by(annotated_clusters,sample_id) %>% 
  dplyr::summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
md_freq_bar_sample$annotated_clusters = factor(md_freq_bar_sample$annotated_clusters, levels=rev(c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific","Stress(hypoxia response)","Stress(p53 response)",
                                                                                                   "Neural Crest Like")))

md_freq_bar_sample <- md_freq_bar_sample %>%
  mutate(clusters = case_when(
    sample_id == "SG_Singlets" ~ "Tumor cells from singlets",
    sample_id == "DB_Tumor_Tcell" ~ "Tumor cells from tumor clusters"
  ))
md_freq_bar_sample$clusters = factor(md_freq_bar_sample$clusters, levels=c("Tumor cells from singlets","Tumor cells from tumor clusters"))


color_combo_sample<-c("Tumor cells from singlets" = "#509abd",
                      "Tumor cells from tumor clusters" =	"#825d9e"		)

pdf("Results_Tumor_Plots/Tumor_Sample_Distributions_Counts.pdf",height=6,width=12)
ggplot(md_freq_bar_sample) + 
  geom_bar(aes(x = n, y = annotated_clusters, fill = clusters), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + xlab("Counts")+ylab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=20),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Source:",values = color_combo_sample) 
dev.off()


#Patient Distribution

md_freq_bar = umap_coo_meta %>% 
  group_by(patient_id,annotated_clusters) %>% 
  dplyr::summarize(n=n()) %>% mutate(freq=n/sum(n), total=sum(n)) 
md_freq_bar$annotated_clusters = factor(md_freq_bar$annotated_clusters, levels=rev(c("Immune Response", "Mitotic","Melanocytic","Transitory Melanocytic","Patient Specific","Stress(hypoxia response)","Stress(p53 response)",
                                                                                     "Neural Crest Like")) )
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
pdf("Results_Tumor_Plots/Tumor_Patient_Distributions.pdf",height=6,width=9)
ggplot(md_freq_bar) + 
  geom_bar(aes(x = annotated_clusters, y = n, fill = patients), stat = "identity") +  # Adjusting y aesthetic
  theme_bw() + ylab("Counts")+xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text=element_text(size=20), 
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  scale_fill_manual(name ="Patients:",values = color_combo) +coord_flip()
dev.off()


