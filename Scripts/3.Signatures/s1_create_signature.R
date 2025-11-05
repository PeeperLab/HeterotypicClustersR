##################################################################
# Creating the clusters signatures
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds output from Tcell_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)

fig_path = "Results_create_signature"
if(!dir.exists(fig_path)) dir.create(fig_path) 

five_pat=readRDS("Results_Tcells/Tcells_Final.Rds")

Idents(five_pat)="sample_id"

# ALL DB vs SG
five_pat$interaction=stringr::str_split_i(five_pat$sample_id, "_", 1)
Idents(five_pat)="interaction"

res_DB_SG_long=FindMarkers(five_pat, ident.1 = "DB", ident.2="SG", logfc.threshold = 0,
                           latent.vars = "patient_id", test.use = "MAST", only.pos = F)
res_DB_SG_long=data.frame(res_DB_SG_long)
res_DB_SG_long$gene=rownames(res_DB_SG_long)
setDT(res_DB_SG_long)

fwrite(res_DB_SG_long, "Results_create_signature/complete_comparison.txt")

res_DB_SG_df=data.frame(res_DB_SG_long)
setDT(res_DB_SG_df)
setorder(res_DB_SG_df, -avg_log2FC)
filt_res_DB_SG=res_DB_SG_df[pct.1>.3 & -log10(p_val_adj)>150,] # high threshold for significance based on the morphology of the volcano plot

top50=filt_res_DB_SG$gene[1:50]
top30=filt_res_DB_SG$gene[1:30]
top100=filt_res_DB_SG$gene[1:100]

int=setdiff(top100, top30)

res_DB_SG_df[, rank:= ifelse(gene %in% top30, "top30",
                             ifelse(gene %in% int, "top100", "outside"))]

write(top100, "Results_create_signature/top100_genes_DB_signature.txt")
write(top30, "Results_create_signature/top30_genes_DB_signature.txt")
