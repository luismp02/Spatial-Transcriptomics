# --- SPOTlight deconvolution ---

library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(SPOTlight)
library(scatterpie)

# PREPROCESS SC REF
sc_ref = readRDS("sc_ref_processed.rds")
seurat_obj = readRDS("seurat_obj_processed.rds")
sce <- as.SingleCellExperiment(sc_ref)
sce <- logNormCounts(sce)


genes <- !grepl("^Rp[l|s]|^Mt", rownames(sce))
dec <- modelGeneVar(sce, subset.row = genes)

# Validation
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

hvg <- getTopHVGs(dec, n = 3000)
colnames(colData(sce)) 
colLabels(sce) <- sce$cell_type 

# Markers estimation

mgs <- scoreMarkers(sce, subset.row = genes)

mgs_df <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  x <- x[x$mean.AUC > 0.8, ]
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_df)

idx <- split(seq(ncol(sce)), sce$cell_type)  
n_cells <- 75  
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells) n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]


#PREPROCESS SPATIAL
spatial_counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")

spatial_coords <- GetTissueCoordinates(seurat_obj)
spatial_coords$x <- as.numeric(spatial_coords$x)
spatial_coords$y <- as.numeric(spatial_coords$y)

spatial_coords_matrix <- as.matrix(spatial_coords[, c("x", "y")])

spe <- SpatialExperiment(
  assays = list(counts = spatial_counts),
  spatialCoords = spatial_coords_matrix
)


# DECONVOLUTION
res <- SPOTlight(
  x = sce,
  y = spe,
  groups = as.character(sce$cell_type),  
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene"
)
                                                                                                                                                                                                                                                                                                                                                                                                                                          sce_cells_clean <- colnames(sce)
spatial_cells_clean <- rownames(spatial_coords_matrix)
sce_cells_clean <- sub("^.*?_","", sce_cells_clean)
spatial_cells_clean <- sub("-.*$", "", spatial_cells_clean)
colnames(sce) <- sce_cells_clean
rownames(spatial_coords_matrix) <- spatial_cells_clean

intersect_cells <- intersect(colnames(sce), rownames(spatial_coords_matrix))
spatial_df <- data.frame(
  x = spatial_coords_matrix[, 1],
  y = spatial_coords_matrix[, 2],
  cell_type = sample(as.character(sce$cell_type), size = nrow(spatial_coords_matrix), replace = TRUE)
)

#Verify
head(spatial_df)
head(res$cell_proportions)

# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
# Extract NMF model fit
mod <- res$NMF
sce$free_annotation <- sce$cell_type 



# FIGURE 1 TOPICS
plot1 <- plotTopicProfiles(
  x = mod,
  y = sce$free_annotation,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)
ggsave("topic_profiles_no_facet.png", plot = plot1, width = 6, height = 7)


# FIGURE 2 TOPICS
plot2 <- plotTopicProfiles(
  x = mod,
  y = sce$free_annotation,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)
ggsave("topic_profiles_facet.png", plot = plot2, width = 9, height = 6)


# FIGURE 3 INTERCATION NETWORK
plot4 <- plotInteractions(mat, which = "network")
ggsave("interactions_network.png", plot = plot4, width = 9, height = 9)

# FIGURE 4: DECONV
ct <- colnames(mat)
mat[mat < 0.1] <- 0

paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

library(png)
library(jpeg)
library(S4Vectors)


# SET SPATIAL AND IMG PATHS
img_path <- "/home/martinpl/projects/datashare/visium_PDAC/Visium_FFPE_V44L01-325_D/outs/spatial/tissue_hires_image.png"
spatial_dir <- "/home/martinpl/projects/datashare/visium_PDAC/Visium_FFPE_V44L01-325_D/outs/spatial"
img_list <- readImgData(
  path = spatial_dir,
  sample_id = "PDAC_sample"
)
imgData(spe) <- img_list

SpatialExperiment::imgData(spe)
SpatialExperiment::getImg(spe)

#scatter
plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))


# CT MAPPING DECONVOLUTION 
df <- as.data.frame(spatialCoords(spe))

df_rot <- data.frame(
  x = df$y,
  y = -df$x
)

df_rot$x <- df_rot$x - min(df_rot$x)
df_rot$y <- df_rot$y - min(df_rot$y)

predominant_type <- apply(mat, 1, function(x) colnames(mat)[which.max(x)])
df_rot$cell_type <- predominant_type

p <- ggplot(df_rot, aes(x = x, y = y, fill = cell_type)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_manual(values = pal) +
  coord_equal() +
  theme_minimal() +
  ggtitle("Predominant CT")

ggsave("predominant_cell_type_rotated.png", plot = p, width = 6, height = 5, dpi = 300)



# HEATMAP
for (ct in colnames(mat)) {
  df_ct <- as.data.frame(spatialCoords(spe))
  df_ct$score <- mat[, ct]
  
  p <- ggplot(df_ct, aes(x = x, y = y, fill = score)) +
    geom_point(shape = 21, size = 2) +
    scale_fill_gradientn(colours = c("white", pal[ct]), name = ct) +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste("Spatial distribution of", ct))
  
  ggsave(paste0("heatmap_", ct, ".png"), plot = p, width = 6, height = 5, dpi = 300)
}


# UMAP
library(umap)

umap_res <- umap(mat)
df_umap <- data.frame(UMAP1 = umap_res$layout[,1],
                      UMAP2 = umap_res$layout[,2],
                      cell_type = apply(mat, 1, function(x) colnames(mat)[which.max(x)]))

p_umap <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = pal) +
  theme_minimal() +
  ggtitle("UMAP of cell type composition per spot")

ggsave("umap_celltype.png", plot = p_umap, width = 6, height = 5, dpi = 300)

