# The Impact of Ethnicity on the Single-cell Immune Response to Falciparum Malaria


This repository conatins the information related to **The Impact of Ethnicity on the Single-cell Immune Response to Falciparum Malaria**.

1. [Data repository](#data-information)
2. [Consortium analysis](#consortium-analysis)
3. [eQTL analysis](#eqtl-analysis)
4. [In house analysis](#independent-unsupervised-clustering)
5. [Comparison between 2 pipelines](#comparison-between-2-pipelines)

---

### Data information

The multiplexed 10x scRNA-seq using the 10x Genomics Chromium with Single Cell 3’ library was sequenced in 2 batches in Illumina NovaSeq 6000 instrument. The sequenced data is subsequently processed using cellranger (3.0.2) by aligning to the human transcriptome(GRCh38-1.2.0) for each pool

```
cellranger count \
--id=Pool007-G7_cellrangerCount \
--fastqs=. \
--sample=Pool007-G7 \
--transcriptome=refdata-cellranger-GRCh38-1.2.0 \
--jobmode=local \
--localcores=14 \
```

Cell ranger output per batch is available via Gene Expression Omnibus

**Batch1:** [GSE273781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273781)

**Batch2:** [GSE273785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273785)


---


### Consortium Analysis

**1) Genotype Quality Control and Phasing**
Genotype data was initially quality-controlled by comparison to the 1000 Genomes Project to validate sex and ancestry matches. Phasing of the genotype data was performed using Eagle, and imputation was carried out with Minimac4. The resulting phased and imputed genotype data was used for de-multiplexing samples and identification of doublet within each pool. Detailed information on the pipeline and scripts is defined [here](https://wg1-pipeline-qc.readthedocs.io/en/latest/Imputation/index.html#imputation).


![Comapare](plots/1Kg.png)

**2) Demultiplexing and Doublet Removal**
Demultiplexing and doublet identification were conducted using multiple tools, including Demuxlet and Scrublet. The initial dataset comprised 33,694 genes across ~196K cells, derived from 2 batches and 32 pools. Droplet-type assignment and associated confidence scores for each cell were computed and reported in output_all_Pools_Asignment.txt for all the different softwares used. Highly confident singlets of ~167K cells were selected for downstream analysis using an intersectional method as outlined in pipeline (results Final_all_Pools_Asignment.txt). Detailed information on the pipeline and scripts is defined [here](https://wg1-pipeline-qc.readthedocs.io/en/latest/Demultiplexing/index.html#demultiplexing).


**3) Cell Type Identification**
Cell type identification was performed using two independent methods: Azimuth (reference-based method within the Seurat framework) and scPred. The results from both methods were merged and compared. Detailed information on the pipeline and scripts is defined [here](https://powellgenomicslab.github.io/WG2-pipeline-classification-docs/).


![Comapare](plots/comparison_heatmap_counts.png)



---

### eQTL analysis

Pesuod bulk was created using the muscate tool (Muscat-DGEA.R)

---

### Independent unsupervised clustering 

The  Demultiplexed and Doublet removed singlecell data was used to perform unsupervised clustering (Louvain method) to identify cell clusters de novo and compare the cells to the Reference based method. We empolyed simmilar QC filtering in both the methods. The script `Harmony.R`  contains the method used for independent clusterings analysis. For comparative analysis we removed the cells with predicted score < 0.7 at predicted level2.


#### Clustering based on Resolution 0.8 and 40 dims

**UMAP :seurat Clustures**

![Clusture UMAP](plots/UmapPlot_Clusture_Filtered0.7.png)

**UMAP :Level1 annotation based on pbmc multimodal**

![Clusture Celltype1](plots/UmapPlot_predicted.celltype.l1_Filtered0.7.png)

**UMAP :Level2 annotation based on pbmc multimodal**

![Clusture Celltype2](plots/UmapPlot_predicted.celltype.l2_Filtered0.7.png)

**UMAP :Based on Ethnicity**

![Clusture Ethnicity](plots/UmapPlot_Ethnicity_Filtered0.7.png)


### Comparison between 2 pipelines

To measure the similarity between two clusterings methods while accounting for chance and differences we ustilized python module `adjusted_rand_score` from the `sklearn.metrics`. We obtained randscore of `0.834068931703432` while comparing Level2 annotation of louvain method clusturing with reference based method and `0.6690088808363425` randscore seurat clustures numbers with reference based method


**Find common cells between the 2 methods (Seurat Level 2 annotation)**

```{r FindCommonLevel2, echo=FALSE}
Idents(seurat_combined) <- "predicted.celltype.l2"

common_cells <- intersect(colnames(seurat_combined), colnames(Consortium))
harmony_clusters <- as.vector(Idents(seurat_combined)[common_cells])
Consortium_clusters <- as.vector(Idents(Consortium)[common_cells])
filtered_clusters <- tolower(Consortium_clusters)
df <- data.frame(cell = common_cells, harmony_clusters, filtered_clusters)
head(df)
write.csv(df, "Harmony_0.8_40Dims_common_cluster_assignments_Filtered_0.7.csv", row.names = FALSE)
```

```{python loadpylib, echo=FALSE}
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("Harmony_0.8_40Dims_common_cluster_assignments_Filtered_0.7.csv")
harmony_clusters = df['harmony_clusters']
filtered_clusters = df['filtered_clusters']

# Calculate Adjusted Rand Index
ari_score = adjusted_rand_score(harmony_clusters, filtered_clusters)
print("Level 2 annoatation", ari_score)
```

`Level 2 annoatation 0.834068931703432`


**Find common cells between the 2 methods (Seurat Clusture Id)**

```{r FindCommonID, echo=FALSE}
Idents(seurat_combined) <- "seurat_clusters"
common_cells <- intersect(colnames(seurat_combined), colnames(Consortium))
harmony_clusters <- as.vector(Idents(seurat_combined)[common_cells])
Consortium_clusters <- as.vector(Idents(Consortium)[common_cells])
filtered_clusters <- tolower(Consortium_clusters)
df <- data.frame(cell = common_cells, harmony_clusters, filtered_clusters)
head(df)
write.csv(df, "Harmony_0.8_40Dims_common_clusterID_Filtered_0.7.csv", row.names = FALSE)
```
Calculating Randscore based on Cluster ID

```{python randClustureID }
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import seaborn as sns
import matplotlib.pyplot as plt
df = pd.read_csv("Harmony_0.8_40Dims_common_clusterID_Filtered_0.7.csv")
harmony_clusters = df['harmony_clusters']
filtered_clusters = df['filtered_clusters']

# Calculate Adjusted Rand Index
ari_score = adjusted_rand_score(harmony_clusters, filtered_clusters)
print("Clusture Level", ari_score)
```

`Clusture Level: 0.6690088808363425` 
