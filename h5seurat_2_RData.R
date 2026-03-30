# We save our Seurat objects in an RData file for Theresa

library(Seurat)
library(SeuratDisk)

# Path where data is saved
data_path <- "/Users/juliaelrod/Desktop/cmu/research/scleroderma/jSSc_ASCT_monocytes_manuscript/code_for_manuscript/data/"

# Load in scler_pbmc
scler_pbmc <- LoadH5Seurat(paste0(data_path, "h5seurat/scler_pbmc_10.h5seurat"))

# Load in SSc_monocytes
SSc_monocytes <- LoadH5Seurat(paste0(data_path, "h5seurat/SSc_monocytes_2.h5seurat"))

# Load in healty_monocytes
healthy_monocytes <- LoadH5Seurat(paste0(data_path, "h5seurat/healthy_monocytes_2.h5seurat"))

# Save as RData
save.image(file = paste0(data_path, "RData/CITE_seurat.RData"))