##################################################################
# LR interactions across the dataset (T cells as receivers)
##################################################################

# required input data in the working directory:
# # output from s0_extract_curated_pairs.R
# # objects generated in /1.Main/ (Tcells_Final.Rds/apc_data_filtered_fvf_corr.rds/Tumor_Annotated.Rds/Combined_All_Cells.Rds)
# # complete_comparison.txt from s1_create_signature.R

# setwd("/YOUR/PATH/)
library(Seurat)
library(data.table)
library(nichenetr)
library(dplyr)
library(tidyverse)
library(data.table)

fig_path = "Results_general_distribution"
if(!dir.exists(fig_path)) dir.create(fig_path)

# get the T cell compartment (Receivers) 
tcells=readRDS("Results_Tcells/Tcells_Final.Rds")
meta_t=tcells@meta.data
meta_t=meta_t[,c(5,5)] #only care about sample_id
colnames(meta_t)[2]="mid_divisions"
meta_t$mid_divisions=rep("cd8", nrow(meta_t))
meta_t$class=rep("cd8", nrow(meta_t))
meta_t$phenotypes=tcells$annotated_clusters_labels_final

# get the APCs (senders)
all_apcs=readRDS("Results_APCs/apc_data_filtered_fvf_corr.rds")

all_apcs=subset(all_apcs, subset= sample_id != "DB_Tumor_Tcell")
apc_groups=unique(all_apcs$major_classes)
meta_a=all_apcs@meta.data
meta_a=meta_a[,c(5,5,which(colnames(meta_a)=="major_classes"))] # get sample id and compartment
colnames(meta_a)[2]="mid_divisions"
meta_a$mid_divisions=rep("apc", nrow(meta_a))
colnames(meta_a)[3]="class"
meta_a$phenotypes=all_apcs$hi_res_clus

# get the Tumor (senders)
all_tum=readRDS("Results_tumor/Tumor_Annotated.Rds")
all_tum=subset(all_tum, sample_id != "DB_APC_Tcell")
all_tum=subset(all_tum, annotated_clusters != "Low Gene Count Cells")
DefaultAssay(all_tum) <- "RNA"
meta_tum=all_tum@meta.data
meta_tum=meta_tum[,c(5,5)] # get sample id and compartment
colnames(meta_tum)[2]="mid_divisions"
meta_tum$mid_divisions=rep("tum", nrow(meta_tum))
meta_tum$class=rep("tum", nrow(meta_tum))
meta_tum$phenotypes=all_tum$annotated_clusters

# get the combined object
all_cells=readRDS("Results_Main/Combined_All_Cells.Rds")
meta_tot=rbind(meta_t, meta_a)
meta_tot=rbind(meta_tot, meta_tum)
meta_tot$final_divisions=ifelse(meta_tot$mid_divisions=="cd8", paste0(meta_tot$sample_id,"_",meta_tot$mid_divisions),
                                paste0(stringr::str_split_i(meta_tot$sample_id,"_",1),"_",meta_tot$mid_divisions) )
all_cells=AddMetaData(all_cells, meta_tot)

DimPlot(all_cells, group.by = "mid_divisions")
DimPlot(all_cells, group.by = "final_divisions")

filt_cells=subset(all_cells, cells = colnames(all_cells)[!is.na(all_cells$final_divisions)])

filt_cells$interaction=ifelse(filt_cells$sample_id!="SG_Singlets", "db","sg")
Idents(filt_cells)="class"
sender_celltypes <- setdiff(as.character(meta_tot$class), as.character(meta_t$class))
filt_cells$phenotypes=ifelse(filt_cells$phenotypes %in% as.character(meta_t$phenotypes), "cd8", as.character(filt_cells$phenotypes))
#saveRDS(filt_cells, "s1_output/filt_cells_scRNA.rds")

# nichenet databases
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# custom, step by step nichenet
# # get restricted possibilities
curated_pairs= fread("Results_Nichenet/s0_potential_pairs.txt")

# # expressed ligands & receptors
Idents(filt_cells)="class"
expressed_genes_receiver <- get_expressed_genes("cd8", subset(filt_cells, interaction=="db"), pct = 0.1) 

# # # break the APCs into classes to avoid losing ligands specific to subtypes
expressed_genes_sender_tum <- get_expressed_genes("tum", subset(filt_cells, interaction=="db"), pct = 0.1) 
expressed_genes_sender_apc_bc <- get_expressed_genes("B_plasma_cells", subset(filt_cells, interaction=="db"), pct = 0.1) 
expressed_genes_sender_apc_dc <- get_expressed_genes("DCs", subset(filt_cells, interaction=="db"), pct = 0.1) 
expressed_genes_sender_apc_mm <- get_expressed_genes("mono_mac", subset(filt_cells, interaction=="db"), pct = 0.1) 
expressed_genes_sender=union(union(expressed_genes_sender_tum,expressed_genes_sender_apc_bc),
                             union(expressed_genes_sender_apc_dc,expressed_genes_sender_apc_mm))

potential_ligands=intersect(expressed_genes_sender, lr_network$from)
potential_receptors=intersect(intersect(expressed_genes_receiver, lr_network$to), 
                              curated_pairs$to) # only Rs curated pairs

potential_ligands <- intersect(potential_ligands,
                               curated_pairs[curated_pairs$to %in% potential_receptors,]$from %>% unique())

potential_ligands= intersect(potential_ligands, curated_pairs$from) # only Ls curated pairs

# # # background genes
background_expressed_genes <- expressed_genes_receiver %>%
  .[. %in% rownames(ligand_target_matrix)]

# predict activities
Idents(filt_cells)="interaction"

# # # get the DGE comparison from the interaction_signature.R
tcell_db_sg=fread("Results_create_signature/complete_comparison.txt") #should be generated previously
setorder(tcell_db_sg, -avg_log2FC)
goi=tcell_db_sg[p_val_adj<0.05 & avg_log2FC>0.1 & pct.1>0.05,]$gene # low pct threshold to include subtype specific
ligand_activities <-  predict_ligand_activities(geneset = goi,
                                                background_expressed_genes = background_expressed_genes,
                                                ligand_target_matrix = ligand_target_matrix,
                                                potential_ligands = potential_ligands, single = T)

ligand_activities %>% arrange(-aupr_corrected) 

best_upstream_ligands <- ligand_activities %>%
  top_n(100, aupr_corrected) %>% # broad representation of all posible pathways
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

# manually assign the ligands
db_obj=subset(filt_cells, interaction=="db")
avg_Exp=AverageExpression(db_obj, group.by = "class",assays = "RNA",slot = "data", features = best_upstream_ligands)
avg_Exp=data.frame(avg_Exp$RNA) %>% rownames_to_column(var="gene")
setDT(avg_Exp)
mavg=melt(avg_Exp, id.vars = "gene")
mavg[, mean:= mean(value), by="gene"]
mavg[, sd:= sd(value), by="gene"]
mavg[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]

dexp=dcast(mavg[!is.na(enriched),], gene~variable, value.var = "enriched")
dexp[, tot_senders:= rowSums(dexp[,.(B_plasma_cells, DCs, mono_mac, tum)])]
dexp[, tot_apc:= rowSums(dexp[,.(B_plasma_cells, DCs, mono_mac)])]

# # don't split the APCs
avg_Exp_bin=AverageExpression(subset(filt_cells, interaction=="db"), group.by = "mid_divisions",
                              assays = "RNA",slot = "data", features = best_upstream_ligands)
avg_Exp_bin=data.frame(avg_Exp_bin$RNA) %>% rownames_to_column(var="gene")
setDT(avg_Exp_bin)
mavg_alt=melt(avg_Exp_bin, id.vars = "gene")
mavg_alt[, mean:= mean(value), by="gene"]
mavg_alt[, sd:= sd(value), by="gene"]
mavg_alt[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]
dexp_alt=dcast(mavg_alt[!is.na(enriched),], gene~variable, value.var = "enriched")
dexp_alt[, tot_senders:= rowSums(dexp_alt[,.(apc, tum)])]
dexp_alt[, ligand_type:= ifelse(tot_senders!=1, "General",
                                ifelse(apc>0 & tum==0, "apc",
                                       ifelse(apc==0 & tum>0, "tum","QUE")))]

ligand_type_df_opt3=data.table(dexp_alt[, .(ligand_type,"ligand"=gene)]) 
ligand_colors <- c("apc" = "#009052","General" = "grey70","tum" = "#825d9e") 
cell_order <- c("tum", "General","apc") 

# manually assign the Receptors to a general T cell type
#tcells=readRDS("/DATA/j.simon/DB_revisions_repository/required_input/Tcells_Final.Rds")
tcells$annotated_clusters_final_red=ifelse(tcells$annotated_clusters_labels_final %in% c("TOX hi Tex", 'GZMK hi Tex',
                                                                                         "LAG3 hi Tex", "TCF7+ stem-like Tex"),"Exhausted_CD8",
                                           ifelse(tcells$annotated_clusters_labels_final %in% c("MKI67 hi Tex/Tprol", "MKI67+ Tex/Tprol", 'MKI67+ Tem-NK like/Tprol'),"Proliferating_CD8",
                                                  ifelse(tcells$annotated_clusters_labels_final %in% c("Tn/Tm", "Tn"),"Naive_Memory_CD8",
                                                         ifelse(tcells$annotated_clusters_labels_final %in% c("Early Tem", "Tem", "Tem-NK like"),"Tem_CD8",
                                                                ifelse(tcells$annotated_clusters_labels_final =="ISG+","ISGpos_CD8",
                                                                       ifelse(tcells$annotated_clusters_labels_final =="Tc17 MAIT","Tc17_MAIT_CD8","NS"))))))


# # assign receptors
t_avg_Exp=AverageExpression(tcells, group.by = "annotated_clusters_final_red",
                            assays = "RNA",slot = "data", features = potential_receptors)
t_avg_Exp=data.frame(t_avg_Exp$RNA) %>% rownames_to_column(var="receptor")
setDT(t_avg_Exp)
t_mavg=melt(t_avg_Exp, id.vars = "receptor")
t_mavg[, mean:= mean(value), by="receptor"]
t_mavg[, sd:= sd(value), by="receptor"]
t_mavg[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]
t_dexp=dcast(t_mavg[!is.na(enriched),], receptor~variable, value.var = "enriched")

t_dexp[, target_type := ifelse(rowSums(.SD == 1) == 1,  
                               names(.SD)[max.col(.SD == 1, ties.method = "first")],
                               "Unspecific"), .SDcols = !c("receptor")]

target_colors_pastel <- c(
  "Exhausted_CD8" = "orchid2",         
  "Naive_Memory_CD8" = "palegreen2",   
  "ISGpos_CD8" = "lightgoldenrod3",    
  "Unspecific" = "gainsboro",              
  "Proliferating_CD8" = "#183676", 
  "Tc17_MAIT_CD8" = "#ffe489",     
  "Tem_CD8" = "lightcoral"             
)

# create final object
lr_network_top_df <- get_weighted_ligand_receptor_links(best_upstream_ligands,
                                                        potential_receptors,
                                                        lr_network,
                                                        weighted_networks$lr_sig) %>%
  dplyr::rename(ligand=from, receptor=to)

setDT(ligand_activities)
lr_network_top_df=merge(lr_network_top_df, ligand_activities[, .("ligand"=test_ligand, aupr_corrected)], by="ligand")
lr_network_top_df=merge(lr_network_top_df, t_dexp[,.(receptor, target_type)], by="receptor")
lr_network_top_df=merge(lr_network_top_df, dexp_alt[,.(ligand=gene, ligand_type)], by="ligand")
lr_network_top_df$lr_pair=paste0(lr_network_top_df$ligand, "_", lr_network_top_df$receptor)
lr_network_top_df=lr_network_top_df[lr_network_top_df$lr_pair %in% curated_pairs$lr_pair,]

# # change confidence to activity
setDT(lr_network_top_df)
lr_network_top_df[, weight:= NULL]
setnames(lr_network_top_df, "aupr_corrected", "weight")
setnames(lr_network_top_df, "receptor", "target")

# support table to choose the best 30 across every step
supp_table=lr_network_top_df[, lapply(.SD, function(x){head(x, n=1)}), by="ligand"]
setorder(supp_table, -weight)
final50=supp_table$ligand[1:50]
lr_network_top_df=lr_network_top_df[ligand %in% final50,]
# # #

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
  ligand_colors =  ligand_colors,
  target_colors = target_colors_pastel,
  celltype_order = cell_order) 

# # load modified circos plot function with minimum transparency adjustment
source("Scripts/4.Nichenet/min_thres_circos.R")
library(circlize)
pdf("Results_general_distribution/confidence_pairs_all.groups.pdf", width=7.7, height=7.7)
min_thres_circos(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))
dev.off()

write.csv(lr_network_top_df,"Results_general_distribution/lr_network_top_df_all.groups.csv")
