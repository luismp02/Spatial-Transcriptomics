library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(Matrix)
library(data.table)
library(SeuratData)
library(gt)
library(igraph)
library(RColorBrewer)
library(SeuratDisk)

# Load annotated scRNA-seq reference
sc_ref <- readRDS("/home/martinpl/projects/datashare/genref_hadaca3/00_peng_k_2019.rds")
sc_ref <- NormalizeData(sc_ref)
sc_ref <- FindVariableFeatures(sc_ref, selection.method = "vst", nfeatures = 3000)
sc_ref <- ScaleData(sc_ref, verbose = FALSE)
sc_ref <- RunPCA(sc_ref, verbose = FALSE)
sc_ref <- RunUMAP(sc_ref, dims = 1:30)
sc_ref <- FindNeighbors(sc_ref, dims = 1:30)
sc_ref <- FindClusters(sc_ref, resolution = 0.5)

# Save plot of reference clustering
p1 <- DimPlot(sc_ref, group.by = "seurat_clusters", label = TRUE) + ggtitle("scRNA clusters")
ggsave("results/scRNA_clusters.png", p1, width = 6, height = 5)

# Load spatial data
seurat_obj <- Load10X_Spatial(
  data.dir = "/home/martinpl/projects/datashare/visium_PDAC/Visium_FFPE_V44L01-325_D/outs",
  slice = "V44L01-325_A"
)

# Prepare spatial object
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Match features
common_genes <- intersect(rownames(sc_ref), rownames(seurat_obj))
sc_ref <- subset(sc_ref, features = common_genes)
seurat_obj <- subset(seurat_obj, features = common_genes)

#Assay
seurat_obj@assays$RNA <- seurat_obj@assays$Spatial
DefaultAssay(seurat_obj) <- "RNA"

# Find anchors between scRNA-seq and spatial
anchors <- FindTransferAnchors(
  reference = sc_ref,
  query = seurat_obj,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:30
)

# Transfer cluster labels
seurat_obj <- TransferData(
  anchorset = anchors,
  refdata = sc_ref$seurat_clusters,
  query = seurat_obj,
  dims = 1:30
)

table(seurat_obj$predicted.id)
SpatialDimPlot(seurat_obj, group.by = "predicted.id", label = TRUE)


# CLUSTERING
p_predicted <- SpatialDimPlot(
  seurat_obj, 
  group.by = "predicted.id", 
  label = TRUE, 
  label.size = 3
) + ggtitle("Predicted Cell Types (TransferData)")
ggsave("results/spatial_predicted_labels.png", p_predicted, width = 7, height = 6)

p_umap <- DimPlot(sc_ref, group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("scRNA-seq UMAP (Peng et al.)")

ggsave("results/scrna_umap_clusters.png", p_umap, width = 6, height = 5)



# FIG 1 CLASSIFICATION WITH TISSUE IMAGE
sc_ref$cell_type <- factor(sc_ref$cell_type)
seurat_obj <- TransferData(
  anchorset = anchors, 
  refdata = sc_ref$cell_type,  
  weight.reduction = "pca",  
  query = seurat_obj,
  dims = 1:30
)

DimPlot(seurat_obj, group.by = "predicted.id") + ggtitle("Predicted Cell Types")
ggsave("predicted_cell_types.png", width = 6, height = 5)
SpatialDimPlot(seurat_obj, group.by = "predicted.id") + ggtitle("Predicted Cell Types in Tissue")
ggsave("spatial_predicted_cell_types.png", width = 8, height = 6)


# FIG 2 WITHOUT TISSUE IMAGE
coords <- GetTissueCoordinates(seurat_obj)
coords$cell_type <- seurat_obj$predicted.id[coords$cell]
coords <- GetTissueCoordinates(seurat_obj)
coords$cell_type <- seurat_obj$predicted.id[coords$cell]

# Transpose (x <-> y)
coords_rotated <- coords
coords_rotated$x_new <- coords$y      
coords_rotated$y_new <- max(coords$x) - coords$x  

ggplot(coords_rotated, aes(x = x_new, y = y_new, color = cell_type)) +
  geom_point(size = 1.5) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("Predicted Cell Types") +
  theme(legend.position = "right")
ggsave("spatial_predicted_cell_types_rotated.png", width = 8, height = 6)

# SAVE DATA
saveRDS(sc_ref, "sc_ref_processed.rds")
saveRDS(seurat_obj, "seurat_obj_processed.rds")
