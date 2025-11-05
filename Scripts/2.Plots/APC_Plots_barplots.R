##################################################################
# APC enrichment visualization
##################################################################

# required input data in the working directory:
# # apc_data_filtered_fvf_corr.rds output from APC_Analysis_Code.R
# # palette output from APC_Analysis_Code.R

#setwd("/YOUR/PATH/")
library(Seurat)
library(data.table)
library(ggplot2)

fig_path = "Results_APC_plots_barplots"
if(!dir.exists(fig_path)) dir.create(fig_path) 

# loadings
apc_all=readRDS("Results_APCs/apc_data_filtered_fvf_corr.rds")
palette_apc=readRDS("Results_APCs/palette_apc_fvf_corr.rds")

# format
meta_apc=apc_all@meta.data
meta_apc=setDT(data.frame(meta_apc))
meta_apc=meta_apc[sample_id!="DB_Tumor_Tcell",]
meta_apc$sample_id=factor(meta_apc$sample_id, levels=c("SG_Singlets","DB_APC_Tcell"))

ggplot(meta_apc[major_classes=="mono_mac"|
                (major_classes=="DCs"& patient_id %in% c("P2","P3","P5"))|
                (major_classes=="B_plasma_cells"& patient_id %in% c("P3","P5"))], aes(sample_id, fill=hi_res_clus))+
  geom_bar(position="fill")+
  scale_fill_manual("",values=palette_apc)+
  facet_grid(major_classes~patient_id)+
  theme_bw()+ylab("Frequency")+xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5), 
        legend.key.height = unit(4, "mm"), legend.key.width = unit(4, "mm"),
        legend.text = element_text(size=8))
ggsave("Results_APC_plots_barplots/all_pats_all_freqs_filtered.pdf", height=5, width=8)
meta_apc_filtered <- meta_apc %>%
  filter(
    major_classes == "mono_mac" |
      (major_classes == "DCs" & patient_id %in% c("P2", "P3", "P5")) |
      (major_classes == "B_plasma_cells" & patient_id %in% c("P3", "P5"))
  ) %>%
  group_by(major_classes, patient_id, sample_id, hi_res_clus) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(Frequency = count / sum(count))
write.csv(meta_apc_filtered,"Results_APC_plots_barplots/all_pats_all_freqs_filtered.csv")
# plot without filtering
meta_apc[, tot_sample_class:= .N, by=c("sample_id", "patient_id", "major_classes")]
meta_apc[, tot_sample_clus:= .N, by=c("sample_id", "patient_id", "hi_res_clus")]
meta_apc[, percent_clus_class:= 100*tot_sample_clus/tot_sample_class] # the percenatge that each cluster takes up in its class and interaction category

meta_apc_plot1=meta_apc[, lapply(.SD, mean), by=c("patient_id","sample_id","major_classes"), .SDcols = "tot_sample_class"]
ggplot(meta_apc_plot1,
       aes(patient_id,major_classes, fill=tot_sample_class))+
  geom_tile()+geom_text(aes(label=tot_sample_class), col="white", size=3)+
  scale_fill_viridis_b(breaks = c(0, 5, 10, 25, 100, 500, 1000), 
                       limits = c(0, 1000))+
  facet_wrap(~sample_id)+
  theme_bw()+ylab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggsave("Results_APC_plots_barplots/filter_out_criteria.pdf", height=3, width=6)

write.csv(meta_apc_plot1,"Results_APC_plots_barplots/filter_out_criteria.csv")

# # some groups are missing: we have to fill them with 0s
sample_ids <- unique(meta_apc$sample_id)
hi_res_clus <- unique(meta_apc$hi_res_clus)
patient_ids <- unique(meta_apc$patient_id)
major_classes <- unique(meta_apc$major_classes)
all_combinations <- expand.grid(sample_id = sample_ids,patient_id = patient_ids, hi_res_clus = hi_res_clus,major_classes=major_classes)
all_combinations

complete_data <- merge(all_combinations, meta_apc[,.(sample_id,patient_id,hi_res_clus,major_classes,percent_clus_class)],
                       by = c("sample_id","patient_id", "hi_res_clus","major_classes"), all.x = TRUE)
complete_data[is.na(complete_data)] <- 0
setDT(complete_data)

# # aggregate to %s by patient/class/cluster/sample_id
avg_pat_meta_1=complete_data[, lapply(.SD,mean), 
                             by=c("sample_id", "hi_res_clus", "major_classes","patient_id"),
                             .SDcols = "percent_clus_class"]

# If filtering - remove classes for some of the patients
avg_pat_meta_1_FILT=avg_pat_meta_1[ (patient_id %in% c("P3","P5"))|
                                      (patient_id=="P2" & major_classes!="B_plasma_cells")|
                                      (patient_id %in% c("P1","P4") & major_classes=="mono_mac")]
avg_pat_meta_2_FILT=avg_pat_meta_1_FILT[, lapply(.SD,mean), 
                              by=c("sample_id", "hi_res_clus", "major_classes"), 
                              .SDcols = "percent_clus_class"]
dcbc=avg_pat_meta_2_FILT[(grepl('DC',hi_res_clus)|grepl('B cells',hi_res_clus)|
                           grepl('plasma',hi_res_clus))&major_classes!="mono_mac",]
dcbc[, major_classes:=factor(major_classes, levels=c("B_plasma_cells", "DCs"))]
ggplot(dcbc, aes(sample_id,y=percent_clus_class, fill=hi_res_clus))+
  geom_bar(stat="identity")+theme_bw()+
  facet_wrap(~major_classes)+
  scale_fill_manual(values=palette_apc, name="")+
  ylab("Average frequency")+xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
ggsave("Results_APC_plots_barplots/bcell_dc_avg_freq_barplot.pdf", height=5, width=4.5)
write.csv(dcbc,"Results_APC_plots_barplots/bcell_dc_avg_freq_barplot.csv")
mm=avg_pat_meta_2_FILT[grepl('mac',hi_res_clus)|grepl('mono',hi_res_clus),]
ggplot(mm, aes(sample_id,y=percent_clus_class, fill=hi_res_clus))+
  geom_bar(stat="identity")+theme_bw()+
  scale_fill_manual(values=palette_apc, name="")+
  ylab("Average frequency")+xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
ggsave("Results_APC_plots_barplots/mm_avg_freq_barplot.pdf", height=5, width=4)
write.csv(mm,"Results_APC_plots_barplots/mm_avg_freq_barplot.csv")
