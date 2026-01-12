#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Matrix)
  library(data.table)
  library(SeuratData)
  library(igraph)
  library(RColorBrewer)
  library(SeuratDisk)
})

args <- commandArgs(trailingOnly = TRUE)

ref_path <- args[1]          # RDS ref scRNA-seq
spatial_data_dir <- args[2]  # Visium Folder
slide_output <- args[3]      # Out name

slide_name <- basename(spatial_data_dir)
if (slide_name == "outs") {
  slide_name <- basename(dirname(spatial_data_dir))
}
slice_name <- sub("^Visium_FFPE_", "", slide_name)

# Load ref
sc_ref <- readRDS(ref_path)
if (!"fine_consensus_annotation" %in% colnames(sc_ref@meta.data)) {
  stop("No fine_consensus_annotation column in reference")
}
sc_ref$cell_type <- sc_ref$fine_consensus_annotation
sc_ref$cell_type <- factor(sc_ref$cell_type)

# Ref preprocess (RNA assay)
DefaultAssay(sc_ref) <- "integrated"
sc_ref <- RunPCA(sc_ref, verbose = FALSE)
sc_ref <- RunUMAP(sc_ref, dims = 1:30)
sc_ref <- FindNeighbors(sc_ref, dims = 1:30)
sc_ref <- FindClusters(sc_ref, resolution = 0.5)

# Plot UMAP por cell_type (fine consensus)
p_ref <- DimPlot(sc_ref, group.by = "cell_type", label = TRUE) +
  ggtitle("scRNA reference (fine consensus)")
ggsave("scRNA_fine_consensus.png", p_ref, width = 7, height = 6)

# Load spatial data
seurat_obj <- Load10X_Spatial(
  data.dir = spatial_data_dir,
  slice = slice_name
)

# Process spatial (RNA assay)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Confirm intersection
common_genes <- intersect(rownames(sc_ref), rownames(seurat_obj))
sc_ref <- subset(sc_ref, features = common_genes)
seurat_obj <- subset(seurat_obj, features = common_genes)

# Assay fix
seurat_obj@assays$RNA <- seurat_obj@assays$Spatial
DefaultAssay(seurat_obj) <- "RNA"

# Annotation transfer
anchors <- FindTransferAnchors(
  reference = sc_ref,
  query = seurat_obj,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:30
)

# Transfer fine_consensus
seurat_obj <- TransferData(
  anchorset = anchors,
  refdata = sc_ref$cell_type,
  query = seurat_obj,
  dims = 1:30,
  k.weight = 5
)

# Visualisation
# Spatial plot transfered cell types
p_pred <- SpatialDimPlot(
  seurat_obj,
  group.by = "predicted.id",
  label = TRUE,
  label.size = 3
) + ggtitle("Predicted Cell Types (fine consensus)")
ggsave("spatial_predicted_fine_consensus.png", p_pred, width = 7, height = 6)

# UMAP  reference
p_umap <- DimPlot(sc_ref, group.by = "cell_type", label = TRUE) +
  ggtitle("scRNA reference UMAP (fine consensus)")
ggsave("scrna_umap_fine_consensus.png", p_umap, width = 7, height = 6)
saveRDS(sc_ref, "sc_ref_processed_fine.rds")
saveRDS(seurat_obj, paste0(slide_output, "_processed.rds"))
