library(Seurat)
library(dplyr)
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')

seurat_combined <- readRDS("Harmony.Rds")
Idents(seurat_combined) <- "seurat_clusters"

Consortium<-readRDS("cons_data_filtered_125.RDS")
Idents(Consortium) <- "predicted.celltype.l2"



consortium_new_cell_names <- paste0(
  "B", consortium_metadata$scRNA_Batch,
  "-Pool", sprintf("%03d", as.numeric(gsub("Pool", "", consortium_metadata$Pool))), # Force three-digit format
  "-", consortium_metadata$Assignment, 
  "_", sapply(strsplit(rownames(consortium_metadata), "_"), `[`, 1),
  "-1"
)

Consortium <- RenameCells(Consortium, new.names = consortium_new_cell_names)
common_cells <- intersect(colnames(seurat_combined), colnames(Consortium))
harmony_clusters <- as.vector(Idents(seurat_combined)[common_cells])
Consortium_clusters <- as.vector(Idents(Consortium)[common_cells])
filtered_clusters <- tolower(Consortium_clusters)
df <- data.frame(cell = common_cells, harmony_clusters, filtered_clusters)
write.csv(df, "common_cluster_assignments.csv", row.names = FALSE)
