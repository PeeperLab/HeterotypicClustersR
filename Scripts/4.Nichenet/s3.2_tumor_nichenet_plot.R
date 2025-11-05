##################################################################
# LR interactions across the dataset in Tumor cells - visualization
##################################################################

# required input data in the working directory:
# # output from s0_extract_curated_pairs.R
# # objects generated in /1.Main/ (Tcells_Final.Rds/ Tumor_Annotated.Rds)
# # /Results_tumor_nn from s3.1_tumor_backbone.R

#setwd("/YOUR/PATH/")
library(nichenetr)
library(data.table)
library(ggplot2)
library(circlize)
library(Seurat)

fig_path = "Results_Tumor_nn"
if(!dir.exists(fig_path)) dir.create(fig_path)

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

######################################## tumor ######################################## 
all_tum=readRDS("Results_tumor/Tumor_Annotated.Rds")
all_tum=subset(all_tum, annotated_clusters != "Low Gene Count Cells") 
activity_ligands_tum=fread("./Results_Tumor_nn/backbone_ligand_activity_tum.txt") # change to current folder
predicted_pairs_tum=fread("./Results_Tumor_nn/backbone_tum_nn_results.txt")

# # assign ligands
avg_Exp=AverageExpression(all_tum, group.by = "annotated_clusters",assays = "RNA",slot = "data",
                          features = activity_ligands_tum$test_ligand)
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

dexp[, ligand_type := ifelse(ligand_type=="shared" & Immune.Response>0 & Stress.hypoxia.response.>0, "intersect_hyp_IR",
                             ifelse(ligand_type=="shared" & Stress.hypoxia.response.>0, "hypoxia_unspecific",
                                    ifelse(ligand_type=="shared" & Immune.Response>0, "IR_unspecific",ligand_type)))]

dexp[, ligand_type := ifelse(ligand_type %in% c("intersect_hyp_IR","hypoxia_unspecific","IR_unspecific",
                                                 "Immune.Response","Stress.hypoxia.response."),
                             ligand_type, "other")]
table(dexp$ligand_type)

lr_table_30=merge(activity_ligands_tum[,.("ligand"=test_ligand, "weight"=aupr_corrected)],predicted_pairs_tum )

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
write.csv(lr_table_30,"Results_Tumor_nn/lr_table_30_tum.csv")

values_tum= c( "IR_unspecific"= "#689a79","Immune.Response"= "#99CC66",
              "intersect_hyp_IR"="limegreen",
              "Stress.hypoxia.response."="#dac565",
              "hypoxia_unspecific"="#b09000",
              "other"="grey")
# # plot
vis_circos_receptor_obj <- prepare_circos_visualization(
  lr_table_30,
  ligand_colors = values_tum,
  target_colors = target_colors_pastel,
  celltype_order = names(values_tum)) 

source("Scripts/4.Nichenet/min_thres_circos.R")

pdf("./Results_Tumor_nn/circos_plot_all_groups_nichenet_tum.pdf", width=7, height=7)
min_thres_circos(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))
dev.off()

write.csv(vis_circos_receptor_obj$links_circle,"Results_Tumor_nn/circos_plot_all_groups_nichenet_tum.csv")
