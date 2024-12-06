##################################################################
# Extract high confidence pairs
##################################################################

# required input data in the working directory:
# # database_curation folder from /data/

library(data.table)
#setwd("/YOUR/PATH")

# # single cell signal R
load("database_curation/data_LRdb.rda")
setDT(LRdb)
LRdb[, lr_pair:= paste0(ligand,"_",receptor)]

# # celltalk_db
celltalk_db=fread("database_curation/human_lr_pair.txt")

# # cell chat - PPI
load("database_curation/PPI.human.rda")
PPI.human=as.matrix(PPI.human) 
PPI.human=data.frame(PPI.human)
PPI.human$prot1=rownames(PPI.human)
setDT(PPI.human)

PPI.human=melt(PPI.human, id.vars="prot1") 
setDT(PPI.human)
PPI.human=PPI.human[value!=0,]
PPI.human[, lr_pair:= paste0(prot1, "_",variable)]
PPI.human[, ppi.hu:=1]

# # cell chat - interactions
load("database_curation/CellChatDB.human.rda")
cc_int=CellChatDB.human$interaction
setDT(cc_int)

# # nichenet
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
nn=weighted_networks$lr_sig
setDT(nn)
nn[, lr_pair:= paste0(from, "_",to)]

# compile info
merged=nn
merged[, PPI.human:= ifelse(lr_pair %in% PPI.human$lr_pair, 1,0)] # experimental
merged[, cc_int:= ifelse(lr_pair %in% cc_int$interaction_name, 1,0)] #annotated, curated - hi trust
merged[, scsR:= ifelse(lr_pair %in% LRdb$lr_pair, 1,0)]#annotated, curated - mid trust
merged[, celltalk_db:= ifelse(lr_pair %in% celltalk_db$lr_pair, 1,0)]#annotated, curated - hi trust
merged[, curated_counts_2:= rowSums(merged[,.(scsR, celltalk_db)])]

merged=merged[weight>0.75,]
merged=merged[weight>0.9 | (weight>0.8 & curated_counts_2>1)|cc_int>0 |  PPI.human>1,]
curated_receptors=union(union(cc_int$receptor, LRdb$receptor), celltalk_db$receptor_gene_symbol)
final_pairs=merged[to %in% curated_receptors ,]

# add receptor localization info from uniprot
uniprot_membrane=fread("database_curation/uniprotkb_cell_membrane_OR_surface_AND_2024_11_07.tsv")
setnames(uniprot_membrane,"Subcellular location [CC]","location")
uniprot_membrane=uniprot_membrane[ grepl("surface", location)| grepl("membrane", location),]
uniprot_membrane_genes=uniprot_membrane$`Gene Names`
uniprot_membrane_genes=unique(unlist(stringr::str_split(uniprot_membrane_genes, " ")))

final_pairs=final_pairs[to %in% uniprot_membrane_genes,]
fwrite(final_pairs[,.(from, to, lr_pair)], "s0_potential_pairs.txt")




