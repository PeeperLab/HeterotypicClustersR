##################################################################
# Tumor cells - statistical enrichment analysis per pat
##################################################################

# required input data in the working directory:
# # Tumor_Annotated.Rds from 1.Main/Tumor_Analysis_Code.R


#setwd("/YOUR/PATH/")
library(data.table)

fig_path = "Results_stats"
if(!dir.exists(fig_path)) dir.create(fig_path)

per_pat_enrichment=function(OBJ){
  
  # clonotype enrichment
  md = OBJ@meta.data
  print(unique(md$patient_id))
  setDT(md)
  md= md[sample_id!="DB_APC_Tcell" & annotated_clusters!="Low Gene Count Cells",]
  
  md_split=split(md, by = c("patient_id", "sample_id"))
  md_split=lapply(md_split, function(x){
    setDT(x)
    x <- x[, .(pheno_tot = .N), by = c("patient_id", "sample_id", "annotated_clusters")]
    y=x[,lapply(.SD, mean), .SDcols="pheno_tot",
        by=c("annotated_clusters", "patient_id", "sample_id")]
    
    return(y)
  })
  
  md_split_dt=rbindlist(md_split)
  
  # fill in 
  sample_ids <- unique(md_split_dt$sample_id)
  annotated_clusters <- unique(md_split_dt$annotated_clusters)
  patient_ids <- unique(md_split_dt$patient_id)
  all_combinations <- expand.grid(sample_id = sample_ids,patient_id = patient_ids, annotated_clusters = annotated_clusters)
  complete_data <- merge(all_combinations, md_split_dt[,.(sample_id,patient_id,annotated_clusters,pheno_tot)],
                         by = c("sample_id","patient_id", "annotated_clusters"), all.x = TRUE)
  complete_data[is.na(complete_data)] <- 0
  setDT(complete_data)
  
  
  complete_data[, tot_sample:= sum(pheno_tot), by=c("sample_id","patient_id")]
  complete_data[, other_top_sample:= tot_sample-pheno_tot]
  
  results <- data.frame(
    patient = character(),
    group = character(),
    comparison = character(),
    p_value = numeric()
  )
  
  for (pat in unique(complete_data$patient_id)){
    for (sam in c("DB_Tumor_Tcell")){
      for (top in unique(complete_data$annotated_clusters)){
        y=complete_data[patient_id == pat,]
        y0=y[sample_id %in% c("SG_Singlets",sam),]
        y0[, tot_pat:= sum(pheno_tot)]
        y1=y0[annotated_clusters ==top,]
        
        print(pat)
        print(sam)
        print(top)
        
        A = ifelse(is.na(as.numeric(y1[sample_id == sam, ]$pheno_tot)), 0, as.numeric(y1[sample_id == sam, ]$pheno_tot)) # goi
        B = ifelse(is.na(as.numeric(y1[sample_id == sam, ]$other_top_sample)), 0, as.numeric(y1[sample_id == sam, ]$other_top_sample)) # other DB
        C = ifelse(is.na(as.numeric(y1[sample_id == "SG_Singlets", ]$pheno_tot)), 0, as.numeric(y1[sample_id == "SG_Singlets", ]$pheno_tot)) 
        D = ifelse(is.na(as.numeric(y1[sample_id == "SG_Singlets", ]$other_top_sample)), 0, as.numeric(y1[sample_id == "SG_Singlets", ]$other_top_sample))
        
        if (length(c(A, B, C, D)) == 4) {
          print(c(A+B+C+D))
          print(y1$tot_pat[1])
          
          table_fish=data.table(db=c(A,B), sg=c(C,D))
          table_fish = table_fish[, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]
          print(table_fish)
          
          pval.res=fisher.test(table_fish)$p.value
          results=rbind(results, c(pat, top, sam, pval.res))}else{
            print(paste0("skip: ",pat,", ",sam,", ",top))
          }
      }
    }
  }
  
  results=setDT(results)
  names(results)=c("patient_id","group","comparison","p_value")
  results[, p_value:= as.numeric(p_value)]
  results[, padj_1:= p.adjust(p_value, method = "fdr")] # across all
  return(results)
}

tum_res=per_pat_enrichment(readRDS("Results_tumor/Tumor_Annotated.Rds"))

tum_res[, sig_symbol := ifelse(padj_1 < 0.0001, "****",
                              ifelse(padj_1 < 0.001, "***",
                              ifelse(padj_1 < 0.01, "**",
                              ifelse(padj_1 < 0.05, "*", "ns"))))]

fwrite(tum_res, "Results_stats/tum_per_pat_stats.txt")


