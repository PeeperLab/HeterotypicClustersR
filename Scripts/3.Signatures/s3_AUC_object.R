##################################################################
# Score T cells for genesets of interest
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds output from Tcell_Analysis_Code.R
# # cluster signatures from s1_create_signature
# # files from data/signatures_tcells

library(Seurat)
library(data.table)
library(AUCell)

fig_path = "Results_AUC_object"
if(!dir.exists(fig_path)) dir.create(fig_path)

# load the t-cell object
set.seed(150799)
#setwd("/YOUR/PATH/")
five_pat=readRDS("Tcells_Final.Rds")
meta=five_pat@meta.data
meta$cell_id=rownames(meta)

#get reads in a sparse matrix format 
exprMatrix= five_pat@assays$RNA@counts
exprMatrix <- as(exprMatrix, "dgCMatrix")

# load all signatures
neotcr_cd8=scan("rosenberg_neoTCR_cd8.txt",what="character")
offringa_TR_9_samples=scan("offringa_TR_9_samples.txt",what="character")
Tirosh_Mel_Exh=scan("Tirosh_Mel_Exh.txt",what="character")
Yost_CD8.Exh=scan("Yost_CD8.Exh",what="character")
db_30=scan("top30_genes_DB_signature.txt",what="character")
db_100=scan("top100_genes_DB_signature.txt",what="character")

# now for the csv files
specificity_oliv=read.csv2("tumor_viral_specificity_oliveira.csv")
setDT(specificity_oliv)
specificity_oliv=split(specificity_oliv, f=factor(specificity_oliv$Signature))
specificity_oliv=lapply(specificity_oliv, function(x){
  return(as.character(x$gene))
})
names(specificity_oliv)=paste0("oliv_", names(specificity_oliv))

all_pathways=fgsea::gmtPathways("c2.cp.v7.5.1.symbols.gmt")

geneSets=c(list("neotcr_cd8"=neotcr_cd8,
                "offringa_TR_9_samples"=offringa_TR_9_samples,
                "Tirosh_Mel_Exh"=Tirosh_Mel_Exh,
                "Yost_CD8.Exh"=Yost_CD8.Exh,
                "db_30"=db_30,
                "db_100"=db_100),
           specificity_oliv,
           list(
             "REACTOME_CTLA4_INHIBITORY_SIGNALING"=all_pathways$REACTOME_CTLA4_INHIBITORY_SIGNALING,
             "WP_TCELL_RECEPTOR_TCR_SIGNALING_PATHWAY"=all_pathways$WP_TCELL_RECEPTOR_TCR_SIGNALING_PATHWAY,
             "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY"=all_pathways$REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY))


cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

set.seed(333)
auc_numbers=t(data.frame(cells_AUC@assays@data$AUC))
auc_numbers=data.frame(auc_numbers)
auc_numbers$cell_id=rownames(auc_numbers)
meta$cell_id=gsub("-",".",rownames(meta))
auc_numbers_merged=merge(auc_numbers, meta, by="cell_id")

saveRDS(auc_numbers_merged,"Results_AUC_object/auc_numbers_merged.rds")
saveRDS(geneSets,"Results_AUC_object/geneSets.rds")
