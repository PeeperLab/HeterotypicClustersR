##################################################################
# LR interactions across the dataset in APCs - prediction
##################################################################

# required input data in the working directory:
# # output from s0_extract_curated_pairs.R
# # objects generated in /1.Main/ (Tcells_Final.Rds/apc_data_filtered_fvf_corr.rds/Combined_All_Cells.Rds/
# #                                 final_bcell_fvf_corr.rds/final_dc_fvf_corr.rds/final_mm_fvf_corr.rds)

# figures will be plotted in s2.2_APC_classes_nichenet_plots.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(nichenetr)
library(tidyverse)

fig_path = "Results_APC_nn"
if(!dir.exists(fig_path)) dir.create(fig_path)

# load the nichenet databases (https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper.md)
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds")) 
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

tcells=readRDS("Tcells_Final.Rds")
all_cells=readRDS("Combined_All_Cells.Rds")
curated_pairs= fread("s0_potential_pairs.txt")
qc_ligands=curated_pairs$from
qc_receptors=curated_pairs$to

nichenet_analysis_per_class=function(SELECTED_PATS,OBJECT,DEFINED_CLUSTERS,CLASS){
  # get the T cell compartment (Receivers) 
  tcells=subset(tcells, subset= patient_id %in% SELECTED_PATS)
  meta_t=tcells@meta.data
  meta_t=meta_t[,c(5,5)] #only care about sample_id
  colnames(meta_t)[2]="mid_divisions"
  meta_t$mid_divisions=rep("cd8", nrow(meta_t))
  meta_t$class=rep("cd8", nrow(meta_t))
  
  # get the APCs (senders)
  all_apcs=readRDS(OBJECT)
  all_apcs=subset(all_apcs, subset= sample_id != "DB_Tumor_Tcell")
  all_apcs=subset(all_apcs, subset= patient_id %in% SELECTED_PATS)
  apc_groups=unique(all_apcs$hi_res_clus)
  meta_a=all_apcs@meta.data
  meta_a=meta_a[,c(5,which(colnames(meta_a)=="hi_res_clus"))] # get sample id and compartment
  meta_a$class=rep("apc", nrow(meta_a))
  colnames(meta_a)[2]="mid_divisions"
  
  # get the combined object
  meta_tot=rbind(meta_t, meta_a)
  meta_tot$final_divisions=ifelse(meta_tot$class=="cd8", paste0(meta_tot$sample_id,"_",meta_tot$mid_divisions),
                                  paste0(stringr::str_split_i(meta_tot$sample_id,"_",1),"_",meta_tot$mid_divisions) )
  all_cells=AddMetaData(all_cells, meta_tot)
  
  filt_cells=subset(all_cells, cells = colnames(all_cells)[!is.na(all_cells$final_divisions) & all_cells$final_divisions!="DB_Tumor_Tcell_cd8"])
  rm(all_cells)
  rm(tcells)
  rm(all_apcs)
  
  # or analysis for background as in singlets vs doublets 
  # the goal is to predict ligands/ receptors and then see which cell types express them
  
  # pre-eliminary
  Idents(filt_cells)="final_divisions"
  
  # 1. Define a set of potential ligands for sender-focused approach
  receiver = "DB_APC_Tcell_cd8"
  expressed_genes_receiver <- get_expressed_genes(receiver, filt_cells, pct = .1) 
  # all other cells for background
  # i think this is not part of the final outcome
  
  all_receptors <- unique(lr_network$to)  
  expressed_receptors <- intersect(intersect(all_receptors, expressed_genes_receiver), qc_receptors)
  
  potential_ligands <- curated_pairs[curated_pairs$to %in% expressed_receptors]$from %>% unique()
  
  sender_celltypes = setdiff(filt_cells$final_divisions %>% unique,c("DB_APC_Tcell_cd8","SG_Singlets_cd8")) # APC major classes??
  sender_celltypes=sender_celltypes[!grepl("SG_", sender_celltypes)]
  print(sender_celltypes)
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, filt_cells, 0.35) #UPDATE METHODS - its by cluster
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  
  potential_ligands_focused <- intersect(intersect(potential_ligands, expressed_genes_sender), lr_network$from )
  
  Idents(filt_cells)="mid_divisions"
  
  # 2. Define the gene set of interest
  condition_oi_tcell <-  "DB_APC_Tcell_cd8"
  condition_reference_tcell <- "SG_Singlets_cd8"
  
  Idents(filt_cells)="final_divisions"
  seurat_obj_receiver <- subset(filt_cells, idents = c(condition_reference_tcell,condition_oi_tcell))
  
  DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,test.use = "MAST", latent.vars = "patient_id",
                                    ident.1 = condition_oi_tcell, ident.2 = condition_reference_tcell,
                                    group.by = "final_divisions",
                                    min.pct = 0.05, logfc.threshold = 0.1, only.pos = T) %>% rownames_to_column("gene")
  
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.1) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  # 3. Define the background genes 
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # 4. Perform NicheNet ligand activity analysis
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands_focused,
                                                 single=T)
  
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
  fwrite(ligand_activities, paste0("Results_APC_nn/backbone_ligand_activity",CLASS,".txt"))
  
  best_upstream_ligands <- ligand_activities %>% top_n(100, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) 
  # some are lost in the last curation, so we overcalculate and filter downstream
  
  lr_network <- lr_network %>% distinct(from, to)
  lr_network_filtered <-  lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors)
  
  # check the pairs in our confidence function
  setDT(lr_network_filtered)
  
  names(lr_network_filtered)[1:2]=c('ligand','receptor')
  lr_network_filtered[, lr_pair:= paste0(ligand,"_",receptor)]
  results=lr_network_filtered[lr_pair %in% curated_pairs$lr_pair,] #double check
  fwrite(results, paste0("Results_APC_nn/backbone_",CLASS,"_nn_results.txt"))
}


nichenet_analysis_per_class(SELECTED_PATS = c("P3","P5"),
                            OBJECT = "final_bcell_fvf_corr.rds",
                            CLASS = "B_cells")

nichenet_analysis_per_class(SELECTED_PATS = c("P2","P3","P5"),
                            OBJECT = "final_dc_fvf_corr.rds",
                            CLASS = "DCs")

nichenet_analysis_per_class(SELECTED_PATS = c("P1","P2","P3","P4","P5"),
                            OBJECT = "final_mm_fvf_corr.rds",
                            CLASS = "mono_mac")


