# HeterotypicClustersR (Latest Updates:2025-11-05)
This repository contains the scripts related to the analysis of scRNA seq and scTCR seq from the clusters datasets. 
![My Image](https://github.com/PeeperLab/HeterotypicClustersR/blob/main/extdata/workflow_fig.jpeg)
# Dependencies
R > 4.3.1  
Seurat == 4.4.0  
Harmony == 1.2.1  
AUCell == 1.24.0  
scRepertoire == 2.0.4  
Nichenet == 2.2.0  

# Overview
The Scripts contains following analysis:    
[1. Main](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/1.Main) - Workflow for processing cell ranger outputs of aligned sequencing samples and integrating samples for cell types identification.     
[2. Plots](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/2.Plots) - celltype specific analysis  
[3. Signatures](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/3.Signatures) - cluster signature analysis  
[4. Nichenet](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/4.Nichenet) - cell-cell interactions analysis among clusters groups  
[5. Reactivity](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/5.Reactivity) - Tcell reactivity analysis  
[6. Stats](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/6.Stats) - clusters enrichment statistics  
[7. TCR_activity](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/7.TCR_activity) - TCR activity comparison between clusters and singlets  
[8. Tumor_Analysis](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/8.Tumor_Analysis) - Tumor cluster enrichment analysis     
[9. Infer_CNV](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/9.Infer_CNV) - Inferred copy number analysis of tumor cells    
[10. Response](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/10.Response) - Cluster signature distribution analysis in Barras response data    
[11. CD39_Scripts](https://github.com/PeeperLab/HeterotypicClustersR/tree/main/Scripts/11.CD39_scripts) - CD39 sorted CD8 T cells comparisons with clusters data 
