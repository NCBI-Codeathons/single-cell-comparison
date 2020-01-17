# Results
For each project analyzed, we computed cell clusters using Seurat and asked whether the sample attributes as submitted by the authors corresponds to unique clusters. Shown below are scatter plots of clusters as computed by Seurat followed by a table with sample attributes that fall within those clusters. Ideally, each cluster should correspond to a unique sample attribute, though that does not always seems to be the case.  Sometimes clusters that appear separated in the plot but are marked as one cell type in reality have more than one cell type in the metadata.  Perhaps with more "resolution" in the cluster marking we could have more closely aligned results.

## Seurat Normalized
The following data were obtained using Seurat's default normalization procedure which is a log normalization.
### SRP011546
![SRP011546](images/SRP011546_counts.png)   

Seurat_Cluster|Sample_Attributes
:-------------|:----------------------------------------------------------
0             |hESC passage_0, hESC passage_10
1             |8-cell embryo, Morulae
2             |2-cell embryo, 4-cell embryo, 8-cell embryo, Oocyte, Zygote
3             |Late blastocyst

For this project, a heatmap with recognizable clusters of runs using spearman rank correlation between run pairs exists as well.  Five clusters are visible, with several outlier runs identifiable.  Quantitative elaboration tbd.
![SRP011546_ranks](images/SRP011546_rankcorr_counts.png)

### SRP057196
![SRP057196](images/SRP057196_counts.png)   

Seurat_Cluster|Tissue             |Cell_Type                                                   |Age
:-------------|:------------------|:-----------------------------------------------------------|:-------------------------------------------------------------------
0             |cortex, hippocampus|astrocytes, hybrid, neurons, oligodendrocytes               |21 years, 22 years, 37 years, 47 years, 50 years, 54 years, 63 years
1             |cortex             |fetal_quiescent, fetal_replicating, neurons                 |16-18 W, 22 years, 37 years
2             |cortex             |astrocytes, hybrid, neurons                                 |21 years, 37 years, 47 years, 50 years, 54 years, 63 years
3             |cortex, hippocampus|astrocytes, endothelial, hybrid, microglia, oligodendrocytes|21 years, 37 years, 50 years, 54 years, 63 years
4             |cortex, hippocampus|endothelial, fetal_replicating                              |16-18 W, 22 years, 37 years, 47 years, 63 years
5             |cortex, hippocampus|OPC, hybrid, microglia, neurons, oligodendrocytes           |21 years, 54 years
6             |cortex, hippocampus|OPC, hybrid, microglia, neurons                             |21 years, 37 years, 47 years, 50 years, 54 years, 63 years

### SRP066632
![SRP066632](images/SRP066632_counts.png)    

Seurat_Cluster|ptprc    |erbb2    |er                |pr                |her2
:-------------|:--------|:--------|:-----------------|:-----------------|:-------
0             |high, low|high, low|negative, positive|negative, positive|negative
1             |low      |high, low|negative, positive|negative, positive|negative

### SRP050499
![SRP050499](images/SRP050499_counts.png)    

Seurat_Cluster|Gender|Cell_Type|Dev_Stage
:-------------|:-----|:--------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0             |F, M  |PGC      |10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo, 4W-embryo, 7W-embryo, 8W-embryo, 8W-embryo1-ps1, 8W-embryo1-ps2, 8W-embryo1-ps4, 8W-embryo1-ps5, 8W-embryo1-ps6, 8W-embryo1-ps7, 8W-embryo1-ps8, 8W-embryo1-ps9
1             |F, M  |PGC, Soma|10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo, 4W-embryo, 7W-embryo, 8W-embryo
2             |F, M  |PGC      |10W-embryo, 11W-embryo, 19W-embryo, 7W-embryo, 8W-embryo
3             |M     |Soma     |10W-embryo, 11W-embryo, 7W-embryo
4             |F     |PGC      |17W-embryo

### SRP061549
![SRP061549](images/SRP061549_counts.png)  

Seurat_Cluster|Dev_Stage
:-------------|:---------------------
0             |GW19.5, GW20.5, GW23.5
1             |GW19.5, GW20.5, GW23.5
2             |GW19.5, GW20.5, GW23.5

![SRP061549_rank](images/SRP061549_rank_counts.png)

## TPM Nornmalization
The following data were obtained by normalizing TPM counts.

### SRP011546
![SRP011546](images/SRP011546_tpm.png)      

Seurat_Cluster|Sample_Attributes
:-------------|:-------------------------------------------
0             |8-cell embryo, Morulae
1             |hESC passage_0, hESC passage_10
2             |Late blastocyst
3             |2-cell embryo, 4-cell embryo, Oocyte, Zygote

### SRP057196
![SRP057196](images/SRP057196_tpm.png)      

Seurat_Cluster|Tissue             |Cell_Type                                                     |Age
:-------------|:------------------|:-------------------------------------------------------------|:----------------------------------------------------------------------------
0             |cortex, hippocampus|astrocytes, fetal_quiescent, hybrid, neurons, oligodendrocytes|16-18 W, 21 years, 22 years, 37 years, 47 years, 50 years, 54 years, 63 years
1             |cortex, hippocampus|fetal_quiescent, fetal_replicating, hybrid, microglia         |16-18 W, 22 years, 50 years, 54 years
2             |cortex             |astrocytes, hybrid, neurons                                   |21 years, 37 years, 47 years, 50 years, 54 years, 63 years
3             |cortex, hippocampus|astrocytes, endothelial, hybrid, microglia, oligodendrocytes  |21 years, 37 years, 50 years, 54 years, 63 years
4             |cortex, hippocampus|OPC, hybrid, microglia, oligodendrocytes                      |37 years, 47 years, 50 years, 54 years, 63 years
5             |cortex, hippocampus|endothelial, fetal_replicating                                |16-18 W, 22 years, 37 years, 47 years, 63 years

### SRP066632
![SRP066632](images/SRP066632_tpm.png)

Seurat_Cluster|ptprc    |erbb2    |er                |pr                |her2
:-------------|:--------|:--------|:-----------------|:-----------------|:-------
0             |high, low|high, low|negative, positive|negative, positive|negative
1             |high, low|high, low|negative, positive|negative, positive|negative
2             |high, low|high, low|negative, positive|negative, positive|negative


### SRP050499
![SRP050499](images/SRP050499_tpm.png)

Seurat_Cluster|Gender|Cell_Type|Dev_Stage
:-------------|:-----|:--------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0             |F, M  |PGC      |10W-embryo, 11W-embryo, 17W-embryo, 4W-embryo, 7W-embryo, 8W-embryo, 8W-embryo1-ps1, 8W-embryo1-ps2, 8W-embryo1-ps4, 8W-embryo1-ps5, 8W-embryo1-ps6, 8W-embryo1-ps7, 8W-embryo1-ps8, 8W-embryo1-ps9
1             |F, M  |PGC      |10W-embryo, 11W-embryo, 19W-embryo, 7W-embryo, 8W-embryo
2             |F, M  |Soma     |10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo, 4W-embryo, 7W-embryo, 8W-embryo
3             |F, M  |PGC      |10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo
4             |M     |Soma     |10W-embryo, 11W-embryo, 7W-embryo
5             |F     |PGC      |17W-embryo


### SRP061549
![SRP061549](images/SRP061549_tpm.png)   

Seurat_Cluster|Dev_Stage
:-------------|:---------------------
0             |GW19.5, GW20.5, GW23.5
1             |GW19.5, GW20.5, GW23.5
2             |GW19.5, GW20.5, GW23.5


## RPKM Normalization
The following data were obtained by normalizing RPKM counts.

### SRP011546
![SRP011546](images/SRP011546_rpkm.png)

Seurat_Cluster|Sample_Attributes
:-------------|:----------------------------------------------------------
0             |hESC passage_0, hESC passage_10
1             |8-cell embryo, Morulae
2             |2-cell embryo, 4-cell embryo, 8-cell embryo, Oocyte, Zygote
3             |Late blastocyst


### SRP057196
![SRP057196](images/SRP057196_rpkm.png)

Seurat_Cluster|Tissue             |Cell_Type                                                                                    |Age
:-------------|:------------------|:--------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------
0             |cortex, hippocampus|fetal_quiescent, hybrid, neurons                                                             |16-18 W, 21 years, 22 years, 37 years, 47 years, 50 years, 54 years, 63 years
1             |cortex             |astrocytes, fetal_quiescent, fetal_replicating, neurons, oligodendrocytes                    |16-18 W, 22 years, 50 years, 54 years, 63 years
2             |cortex, hippocampus|OPC, astrocytes, endothelial, fetal_replicating, hybrid, microglia, neurons, oligodendrocytes|16-18 W, 21 years, 22 years, 37 years, 47 years, 50 years, 54 years, 63 years
3             |cortex             |astrocytes                                                                                   |21 years, 37 years, 47 years, 50 years
4             |cortex, hippocampus|astrocytes, endothelial, hybrid, microglia, oligodendrocytes                                 |21 years, 37 years, 50 years, 54 years, 63 years
5             |cortex             |fetal_replicating                                                                            |16-18 W


### SRP066632
![SRP066632](images/SRP066632_rpkm.png)

Seurat_Cluster|ptprc    |erbb2    |er                |pr                |her2
:-------------|:--------|:--------|:-----------------|:-----------------|:-------
0             |high, low|high, low|negative, positive|negative, positive|negative
1             |high, low|high, low|negative, positive|negative, positive|negative

### SRP050499
![SRP050499](images/SRP050499_rpkm.png)

Seurat_Cluster|Gender|Cell_Type|Dev_Stage
:-------------|:-----|:--------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0             |F, M  |PGC      |10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo, 4W-embryo, 7W-embryo, 8W-embryo, 8W-embryo1-ps1, 8W-embryo1-ps2, 8W-embryo1-ps4, 8W-embryo1-ps5, 8W-embryo1-ps6, 8W-embryo1-ps7, 8W-embryo1-ps8, 8W-embryo1-ps9
1             |F, M  |Soma     |10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo, 4W-embryo, 7W-embryo, 8W-embryo
2             |F, M  |PGC, Soma|10W-embryo, 11W-embryo, 17W-embryo, 19W-embryo
3             |M     |Soma     |10W-embryo, 11W-embryo, 7W-embryo
4             |F     |PGC      |17W-embryo

### SRP061549
![SRP061549](images/SRP061549_rpkm.png)

Seurat_Cluster|Dev_Stage
:-------------|:---------------------
0             |GW19.5, GW20.5, GW23.5
1             |GW19.5, GW20.5, GW23.5
2             |GW19.5, GW20.5, GW23.5
