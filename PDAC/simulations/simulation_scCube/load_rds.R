library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(devtools)

convert_seurat_to_h5ad <- function(file, output_prefix) {
  seurat_obj <- readRDS(file)
  
  if ("RNA" %in% names(seurat_obj@assays)) {
    seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
  }
  
  h5seurat_file <- paste0(output_prefix, ".h5seurat")
  SaveH5Seurat(seurat_obj, filename = h5seurat_file)
  Convert(h5seurat_file, dest = "h5ad")
}

convert_seurat_to_h5ad("~/projects/datashare/genref_hadaca3/00_peng_k_2019.rds", "peng")
convert_seurat_to_h5ad("~/projects/datashare/genref_hadaca3/00_raghavan_s_2021.rds", "raghavan")







