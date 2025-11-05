##################################################################
# LR interactions across the dataset in APCs - visualization
##################################################################

# required input data in the working directory:
# # output from s0_extract_curated_pairs.R
# # objects generated in /1.Main/ (Tcells_Final.Rds/ final_bcell_fvf_corr.rds/final_dc_fvf_corr.rds/final_mm_fvf_corr.rds)
# # /Results_APC_nn from s2.1_APC_classes_nichenet_backbone.R

#setwd("/YOUR/PATH/")
library(nichenetr)
library(data.table)
library(ggplot2)
library(circlize)
library(Seurat)
library(dplyr)
library(tidyverse)

# # create the T cell annotation
tcells=readRDS("Results_Tcells/Tcells_Final.Rds")
tcells$annotated_clusters_final_red=ifelse(tcells$annotated_clusters_labels_final %in% c("TOX hi Tex", 'GZMK hi Tex',
                                                                                         "LAG3 hi Tex", "TCF7+ stem-like Tex"),"Exhausted_CD8",
                                           ifelse(tcells$annotated_clusters_labels_final %in% c("MKI67 hi Tex/Tprol", "MKI67+ Tex/Tprol", 'MKI67+ Tem-NK like/Tprol'),"Proliferating_CD8",
                                                  ifelse(tcells$annotated_clusters_labels_final %in% c("Tn/Tm", "Tn"),"Naive_Memory_CD8",
                                                         ifelse(tcells$annotated_clusters_labels_final %in% c("Early Tem", "Tem", "Tem-NK like"),"Tem_CD8",
                                                                ifelse(tcells$annotated_clusters_labels_final =="ISG+","ISGpos_CD8",
                                                                       ifelse(tcells$annotated_clusters_labels_final =="Tc17 MAIT","Tc17_MAIT_CD8","NS"))))))

target_colors_pastel <- c(
  "Exhausted_CD8" = "orchid2",         
  "Naive_Memory_CD8" = "palegreen2",   
  "ISGpos_CD8" = "lightgoldenrod3",    
  "Unspecific" = "gainsboro",              
  "Proliferating_CD8" = "#183676", 
  "Tc17_MAIT_CD8" = "#ffe489",     
  "Tem_CD8" = "lightcoral"             
)


curated_pairs= fread("Results_Nichenet/s0_potential_pairs.txt")
qc_ligands=curated_pairs$from
qc_receptors=curated_pairs$to

# # assign receptors
t_avg_Exp=AverageExpression(tcells, group.by = "annotated_clusters_final_red",
                            assays = "RNA",slot = "data", features = qc_receptors)
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

######################################## macrophages ######################################## 

mm=readRDS("Results_APCs/final_mm_fvf_corr.rds")
activity_ligands_mm=fread("Results_APC_nn/backbone_ligand_activitymono_mac.txt") 
predicted_pairs_mm=fread("Results_APC_nn/backbone_mono_mac_nn_results.txt")

# # assign ligands
avg_Exp=AverageExpression(mm, group.by = "hi_res_clus",assays = "RNA",slot = "data",
                          features = activity_ligands_mm$test_ligand)
avg_Exp=data.frame(avg_Exp$RNA) %>% rownames_to_column(var="ligand")
setDT(avg_Exp)
mavg=melt(avg_Exp, id.vars = "ligand")
mavg[, mean:= mean(value), by="ligand"]
mavg[, sd:= sd(value), by="ligand"]
mavg[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]

dexp=dcast(mavg[!is.na(enriched),], ligand~variable, value.var = "enriched")
dexp[, ligand_type := ifelse(rowSums(.SD == 1) == 1,  
                             names(.SD)[max.col(.SD == 1, ties.method = "first")],
                             "shared"), .SDcols = !c("ligand")]

dexp[, ligand_type := ifelse(ligand_type=="shared" & C1Qhi.inflammatory.macrophages>0 & C1Qhi.LA.macrophages>0, "Intersect_C1QhiLA_C1Q.infl",
                             ifelse(ligand_type=="shared" & C1Qhi.LA.macrophages>0, "C1Qhi.LA.macro_unspecific",
                                    ifelse(ligand_type=="shared" & C1Qhi.inflammatory.macrophages>0, "C1Qhi.infl.macro_unspecific",ligand_type)))]

dexp[, ligand_type := ifelse(grepl("LA", ligand_type)|grepl("inf", ligand_type, ignore.case = T),
                                                            ligand_type, "other")]
table(dexp$ligand_type)

lr_table_30=merge(activity_ligands_mm[,.("ligand"=test_ligand, "weight"=aupr_corrected)],predicted_pairs_mm )

# support table to choose the best 30 across every step
supp_table=lr_table_30[, lapply(.SD, function(x){head(x, n=1)}), by="ligand"]
setorder(supp_table, -weight)
final30=supp_table$ligand[1:30]
lr_table_30=lr_table_30[ligand %in% final30,]
# # #

lr_table_30=merge(lr_table_30, t_dexp[,.(receptor, target_type)], by="receptor")
lr_table_30=merge(lr_table_30,dexp[,.(ligand, ligand_type)],by="ligand")
setnames(lr_table_30, "receptor", "target")

table(lr_table_30$ligand_type)

values_mm= c( "C1Qhi.LA.macro_unspecific"="#7cCFFF",
              "C1Qhi.LA.macrophages"= "#739FFF",
              "Intersect_C1QhiLA_C1Q.infl"="red",
              "C1Qhi.inflammatory.macrophages"= "#0000BF",
              "C1Qhi.infl.macro_unspecific"="#0069ff",
              "other"="grey")
# # plot
vis_circos_receptor_obj <- prepare_circos_visualization(
  lr_table_30,
  ligand_colors = values_mm,
  target_colors = target_colors_pastel,
  celltype_order = names(values_mm)) 

source("Scripts/4.Nichenet/min_thres_circos.R")
write.csv(lr_table_30,"Results_APC_nn/lr_table_30_monomacro.csv")

pdf("Results_APC_nn/circos_plot_all_groups_nichenet_monomacro.pdf", width=7, height=7)
min_thres_circos(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))
dev.off()
write.csv(vis_circos_receptor_obj$links_circle,"Results_APC_nn/circos_plot_all_groups_nichenet_monomacro.csv")
######################################## DCs ######################################## 
dcs=readRDS("Results_APCs/final_dc_fvf_corr.rds")
activity_ligands_dcs=fread("Results_APC_nn/backbone_ligand_activityDCs.txt") # change to current folder
predicted_pairs_dcs=fread("Results_APC_nn/backbone_DCs_nn_results.txt")

# # assign ligands
avg_Exp=AverageExpression(subset(dcs, patient_id %in% c("P2", "P3", "P5")), group.by = "hi_res_clus",
                          assays = "RNA",slot = "data", features = activity_ligands_dcs$test_ligand)
avg_Exp=data.frame(avg_Exp$RNA) %>% rownames_to_column(var="ligand")
setDT(avg_Exp)
mavg=melt(avg_Exp, id.vars = "ligand")
mavg[, mean:= mean(value), by="ligand"]
mavg[, sd:= sd(value), by="ligand"]
mavg[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]

dexp=dcast(mavg[!is.na(enriched),], ligand~variable, value.var = "enriched")
dexp[, ligand_type := ifelse(rowSums(.SD == 1) == 1,  
                             names(.SD)[max.col(.SD == 1, ties.method = "first")],
                             "shared"), .SDcols = !c("ligand")]

dexp[, ligand_type := ifelse(ligand_type=="shared" & pDC>0 & mreg.DC>0, "Intersect_mreg_pDC",
                             ifelse(ligand_type=="shared" & mreg.DC>0, "mreg.DC_unspecific",
                                    ifelse(ligand_type=="shared"  & pDC>0, "pDC_unspecific",ligand_type)))]

dexp[, ligand_type := ifelse(grepl("pdc", ligand_type, ignore.case = T)|grepl("mreg", ligand_type, ignore.case = T),
                             ligand_type, "other")]

table(dexp$ligand_type)

lr_table_30=merge(activity_ligands_dcs[,.("ligand"=test_ligand, "weight"=aupr_corrected)],predicted_pairs_dcs )

# support table to choose the best 30 across every step
supp_table=lr_table_30[, lapply(.SD, function(x){head(x, n=1)}), by="ligand"]
setorder(supp_table, -weight)
final30=supp_table$ligand[1:30]
lr_table_30=lr_table_30[ligand %in% final30,]
# # #

lr_table_30=merge(lr_table_30, t_dexp[,.(receptor, target_type)], by="receptor")
lr_table_30=merge(lr_table_30, dexp[,.(ligand, ligand_type)], by="ligand")
setnames(lr_table_30, "receptor", "target")

ligand_colors_dc=c( 
  "mreg.DC_unspecific"="#E6A6A6","mreg.DC"="firebrick4",
  "Intersect_mreg_pDC"="red","pDC"="#7900C6","pDC_unspecific"="#DAB8F0",
                    "other"="grey")

# # plot
vis_circos_receptor_obj <- prepare_circos_visualization(
  lr_table_30,
  ligand_colors =   ligand_colors_dc,
  target_colors = target_colors_pastel,
  celltype_order = names(ligand_colors_dc)) 
write.csv(lr_table_30,"lr_table_30_DC.csv")
pdf("Results_APC_nn/circos_plot_all_groups_nichenet_DC.pdf", width=7, height=7)
min_thres_circos(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))
dev.off()
write.csv(vis_circos_receptor_obj$links_circle,"Results_APC_nn/circos_plot_all_groups_nichenet_DC.csv")
######################################## B cells ######################################## 
bcells=readRDS("Results_APCs/final_bcell_fvf_corr.rds")
activity_ligands_bcells=fread("Results_APC_nn/backbone_ligand_activityB_cells.txt") # change to current folder
predicted_pairs_bcells=fread("Results_APC_nn/backbone_B_cells_nn_results.txt")

activity_ligands_bcells=activity_ligands_bcells[test_ligand %in% predicted_pairs_bcells$ligand,]

# # assign ligands
avg_Exp=AverageExpression(subset(bcells, patient_id %in% c("P3", "P5")), group.by = "hi_res_clus",
                          assays = "RNA",slot = "data", features = activity_ligands_bcells$test_ligand)
avg_Exp=data.frame(avg_Exp$RNA) %>% rownames_to_column(var="ligand")
setDT(avg_Exp)
mavg=melt(avg_Exp, id.vars = "ligand")
mavg[, mean:= mean(value), by="ligand"]
mavg[, sd:= sd(value), by="ligand"]
mavg[, enriched:= ifelse((value- c(mean+sd))>0, 1,0)]

dexp=dcast(mavg[!is.na(enriched),], ligand~variable, value.var = "enriched")
dexp[, ligand_type := ifelse(rowSums(.SD == 1) == 1,  
                             names(.SD)[max.col(.SD == 1, ties.method = "first")],
                             "shared"), .SDcols = !c("ligand")]

dexp[, ligand_type := ifelse(ligand_type=="shared" & plasma.cells>0 , "plasma.cells_unspecific",ligand_type)]

dexp[, ligand_type := ifelse(grepl("plasma", ligand_type),ligand_type, "other")]

table(dexp$ligand_type)

lr_table_30=merge(activity_ligands_bcells[,.("ligand"=test_ligand, "weight"=aupr_corrected)],predicted_pairs_bcells )

# support table to choose the best 30 across every step
supp_table=lr_table_30[, lapply(.SD, function(x){head(x, n=1)}), by="ligand"]
setorder(supp_table, -weight)
final30=supp_table$ligand[1:30]
lr_table_30=lr_table_30[ligand %in% final30,]
# # #

lr_table_30=merge(lr_table_30, t_dexp[,.(receptor, target_type)], by="receptor")
lr_table_30=merge(lr_table_30, dexp[,.(ligand, ligand_type)], by="ligand")
setnames(lr_table_30, "receptor", "target")

ligand_colors_bcell=c( 
  "plasma.cells_unspecific"="forestgreen","plasma.cells"="#032507",
  "other"="grey")

# # plot
vis_circos_receptor_obj <- prepare_circos_visualization(
  lr_table_30,
  ligand_colors =   ligand_colors_bcell,
  target_colors = target_colors_pastel,
  celltype_order = names(ligand_colors_bcell)) 
write.csv(lr_table_30,"Results_APC_nn/lr_table_30_bcell_fvf_corr.csv")

pdf("Results_APC_nn/circos_plot_all_groups_nichenet_bcell_fvf_corr.pdf", width=7, height=7)
min_thres_circos(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))
dev.off()

