suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
})

# SPOTlight deconvolution 
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(SPOTlight)
library(scatterpie)
library(Matrix)
library(jsonlite)
library(RColorBrewer)
library(SpatialExperiment)
library(pheatmap)
library(ineq)


# === base paths ===
results_base <- "/home/martinpl/projects/PDAC/deconv/test6"
visium_base <- "/home/martinpl/projects/datashare/PDAC/visium_PDAC"

slides <- c(
  "Visium_FFPE_V43T08-051_D", "Visium_FFPE_V44L22-378_D", "Visium_FFPE_V44L23-391_D",
  "Visium_FFPE_V43T08-041_A", "Visium_FFPE_V44L01-325_A", "Visium_FFPE_V44L23-362_A",
  "Visium_FFPE_V43T08-041_D", "Visium_FFPE_V44L01-325_D", "Visium_FFPE_V44L23-362_D",
  "Visium_FFPE_V43T08-051_A", "Visium_FFPE_V44L22-378_A", "Visium_FFPE_V44L23-391_A"
)

for (slide_name in slides) {
  message("Processing: ", slide_name)
  
  slide_dir <- file.path(results_base, slide_name)
  setwd(slide_dir)
  
  # === load data ===
  seurat_obj <- readRDS("seurat_obj_processed.rds")
  res <- readRDS(paste0("spotlight_result_", slide_name, ".rds"))
  mat <- res$mat
  mat <- mat[complete.cases(mat), ]
  
  # === spatial data ===
  spatial_dir <- file.path(visium_base, slide_name, "outs", "spatial")
  scalefactors_path <- file.path(spatial_dir, "scalefactors_json.json")
  scales <- jsonlite::fromJSON(readLines(scalefactors_path))
  spot_diameter <- scales$spot_diameter_lowres
  
  # === coords ===
  coords <- GetTissueCoordinates(seurat_obj)
  coords$x_new <- coords$y
  coords$y_new <- max(coords$x) - coords$x
  coords <- coords[rownames(mat), ]
  
  # === Scale ===
  scales <- jsonlite::fromJSON(readLines(scalefactors_path))
  spot_diameter <- scales$spot_diameter_lowres
  point_size <- spot_diameter * 1 
  coords <- GetTissueCoordinates(seurat_obj)
  coords$x_new <- coords$y
  coords$y_new <- max(coords$x) - coords$x
  coords <- coords[rownames(mat), ]
  df_rot <- GetTissueCoordinates(seurat_obj)


  # filter 
  abundance_mean <- colMeans(mat)
  mat_filtered <- mat[, abundance_mean > 0.01]
  
  #  Corelation matrix 
  cor_matrix <- cor(mat_filtered, method = "spearman")
  max_cor <- max(abs(cor_matrix), na.rm = TRUE)
  breaks <- seq(-max_cor, max_cor, length.out = 101)
  
  heatmap_plot <- pheatmap(
    cor_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks = breaks,
    main = "Cell type colinearity (SPOTlight)",
    border_color = NA
  )

  png("correlation_matrix_SPOTlight.png", width = 1800, height = 1500, res = 300)
  grid::grid.draw(heatmap_plot$gtable)
  dev.off()
  
  #  Codominance
  cell_types <- colnames(mat)
  mat_codom <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
                      dimnames = list(cell_types, cell_types))
  for (i in 1:nrow(mat)) {
    top2 <- sort(mat[i, ], decreasing = TRUE)[1:2]
    ct1 <- names(top2)[1]
    ct2 <- names(top2)[2]
    mat_codom[ct1, ct2] <- mat_codom[ct1, ct2] + 1
  }
  mat_codom <- mat_codom + t(mat_codom) - diag(diag(mat_codom))
  palette_codom <- colorRampPalette(c("white", "orange", "red3"))(100)
  png(paste0("codominance_", slide_name, ".png"), width = 2000, height = 1800, res = 300)
  pheatmap(mat_codom, cluster_rows = TRUE, cluster_cols = TRUE,
           color = palette_codom, main = "Cell type co-dominance per spot",
           border_color = NA)
  dev.off()
  
  #  Gini Index 
  gini_index <- apply(mat, 1, function(x) ineq::Gini(x))
  df_gini <- data.frame(spot = names(gini_index), gini = gini_index)
  df_gini$x <- coords$x_new
  df_gini$y <- coords$y_new
  plot_gini <- ggplot(df_gini, aes(x = x, y = y, color = gini)) +
    geom_point(size = 1.3) +
    scale_color_viridis_c(option = "plasma") +
    theme_void() +
    coord_fixed() +
    ggtitle("Gini index per spot")
  ggsave(paste0("gini_index_", slide_name, ".png"), plot = plot_gini, width = 6, height = 6, dpi = 300)
  
  # Spatial corelation 
  results <- lapply(colnames(mat), function(celltype) {
    abund <- mat[, celltype]
    data.frame(
      celltype = celltype,
      cor_x = cor(coords$x_new, abund, method = "spearman"),
      cor_y = cor(coords$y_new, abund, method = "spearman"),
      p_x = cor.test(coords$x_new, abund, method = "spearman")$p.value,
      p_y = cor.test(coords$y_new, abund, method = "spearman")$p.value
    )
  })
  
  
  df_cor <- do.call(rbind, results)
  df_melt <- reshape2::melt(df_cor[, c("celltype", "cor_x", "cor_y")], id.vars = "celltype")
  plot_cor_space <- ggplot(df_melt, aes(x= variable, y = reorder(celltype, value), fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Correlation with spatial coordinates",
         x = "Spatial axis", y = "Cell type") +
    theme_minimal(base_size = 12)
  ggsave(paste0("correlation_spatial_coordinates_", slide_name, ".png"),
         plot = plot_cor_space, width = 6, height = 5, dpi = 300)
}
