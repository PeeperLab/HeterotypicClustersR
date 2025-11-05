##################################################################
# Characterise the cluster signatures
##################################################################

# required input data in the working directory:
# # Tcells_Final.Rds output from Tcell_Analysis_Code.R
# # signatures created in s1_create_signature.R
# # signatures offringa_TR_9_samples.txt,rosenberg_neoTCR_cd8.txt and tumor_viral_specificity_oliveira.csv  from data/signatures_tcells
# # c5.go.bp.v2024.1.Hs.symbols.gmt from data/signatures_tcells

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)

fig_path = "Results_annotate_signature"
if(!dir.exists(fig_path)) dir.create(fig_path) 

# complete comparison
comp=fread("Results_create_signature/complete_comparison.txt")

# signatures of interest
top100=scan("Results_create_signature/top100_genes_DB_signature.txt", what="character")
top30=scan("Results_create_signature/top30_genes_DB_signature.txt", what="character")
offringa_9_samples=scan("extdata/signatures_tcells/offringa_TR_9_samples.txt", what = "character")
neotcr_cd8=scan("extdata/signatures_tcells/rosenberg_neoTCR_cd8.txt",what="character")
oliv_Tumor.specific=fread("extdata/signatures_tcells/tumor_viral_specificity_oliveira.csv")
oliv_Tumor.specific=oliv_Tumor.specific[Signature=="Tumor-specific",]$gene

comp[, top30:= ifelse(gene %in% top30, 1,0)]
comp[, top100:= ifelse(gene %in% top100, 1,0)]
comp[, offringa_9_samples:= ifelse(gene %in% offringa_9_samples, 1,0)]
comp[, neotcr_cd8:= ifelse(gene %in% neotcr_cd8, 1,0)]

# FORA
library(fgsea)
library(scales)
go_path=gmtPathways("extdata/signatures_tcells/c5.go.bp.v2024.1.Hs.symbols.gmt")

# # top30
go_res=fora(go_path,top30, universe=comp$gene, minSize = 5, maxSize = 300)

setDT(go_res)
setorder(go_res, padj)
go_res[, pathway:= paste0(pathway," (",size,")")]
go_res[, pathway:= gsub("GOBP_","", pathway)]
df <- do.call(rbind, go_res)
write.csv(df,"Results_annotate_signature/ORA_All.csv")

ggplot(go_res[1:10],
       aes(-log10(padj), reorder(gsub("_",' ',pathway),-log10(padj)), fill=overlap))+
  geom_col()+scale_fill_gradient(high="red4", low="#feb9b9")+ylab("")+
  scale_y_discrete(labels=label_wrap(35))+theme_bw()+
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6))
ggsave("Results_annotate_signature/ORA_top30.pdf", width=4, height=3.1)
df <- do.call(rbind, go_res[1:10])
write.csv(df,"Results_annotate_signature/ORA_top30.csv")

# # top100
go_res=fora(go_path,top100, universe=comp$gene, minSize = 5, maxSize = 300)

setDT(go_res)
setorder(go_res, padj)
go_res[, pathway:= paste0(pathway," (",size,")")]
go_res[, pathway:= gsub("GOBP_","", pathway)]
ggplot(go_res[1:10],
       aes(-log10(padj), reorder(gsub("_",' ',pathway),-log10(padj)), fill=overlap))+
  geom_col()+scale_fill_gradient(high="red4", low="#feb9b9")+ylab("")+
  scale_y_discrete(labels=label_wrap(35))+theme_bw()+
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6))
ggsave("Results_annotate_signature/ORA_top100.pdf", width=4, height=3.1)



# intersect heatmap
comp[, oliv_Tumor.specific:= ifelse(gene %in% oliv_Tumor.specific, 1,0)]

comp_30=comp[top30>0,]
m30=melt(comp_30, id.vars = c("gene", "avg_log2FC"),
         measure.vars = c("oliv_Tumor.specific","neotcr_cd8","offringa_9_samples","top30"))
setDT(m30)
m30[, value2:=ifelse(value==1, "Present", "Absent")]
ggplot(m30, aes(reorder(gene, -avg_log2FC), variable, fill=value2))+geom_tile()+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_fill_manual("",values =c("Present"="grey50", "Absent"="white") )+xlab("")+ylab("")
ggsave("Results_annotate_signature/top30_intersect_with_other_signatures.pdf", width=7.5, height=2.5)


comp_100=comp[top100>0,]
m100=melt(comp_100, id.vars = c("gene", "avg_log2FC"),
          measure.vars = c("oliv_Tumor.specific","neotcr_cd8","offringa_9_samples","top100"))
setDT(m100)
m100[, value2:=ifelse(value==1, "Present", "Absent")]
ggplot(m100, aes(reorder(gene, -avg_log2FC), variable, fill=value2))+geom_tile()+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_fill_manual("",values =c("Present"="grey50", "Absent"="white") )+xlab("")+ylab("")
ggsave("Results_annotate_signature/top100_intersect_with_other_signatures.pdf", width=16, height=2.5)
