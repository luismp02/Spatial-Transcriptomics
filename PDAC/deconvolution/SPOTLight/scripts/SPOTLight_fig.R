suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(SPOTlight)
  library(scater)
  library(scran)
  library(dplyr)
  library(ggplot2)
  library(scatterpie)
  library(tidyr)
  library(scales)
  library(pheatmap)
  library(reshape2)
  library(RColorBrewer)
  library(png)
  library(gridExtra)
  library(tibble)
  library(jpeg)
  library(S4Vectors)
  library(jsonlite)
})

# config
root_dir <- "/home/martinpl/projects/PDAC/deconv/final"
base_dir <- "/home/martinpl/projects/datashare/PDAC/visium_PDAC"

# slides
slides <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)


for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  cat("\n Processing slide:", slide_name, "\n")
  
  setwd(slide_dir)
  cat("WD:", getwd(), "\n")
  
  # paths
  spatial_dir <- file.path(base_dir, slide_name, "outs", "spatial")
  scalefactors_path <- file.path(spatial_dir, "scalefactors_json.json")
  
  # verify files
  if (!file.exists("seurat_obj.rds") || !file.exists("spotlight_result.rds")) {
    cat("NO files", slide_dir, ", skip...\n")
    next
  }
  
  # load data
  seurat_obj <- readRDS("seurat_obj.rds")
  res <- readRDS("spotlight_result.rds")
  mat <- res$mat
  cat("Loaded:", nrow(mat), "spots x", ncol(mat), "cell types\n")
  
  # spatial
  spatial_counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")
  spatial_coords <- GetTissueCoordinates(seurat_obj)
  spatial_coords_matrix <- as.matrix(spatial_coords[, c("x", "y")])
  
  spe <- SpatialExperiment(
    assays = list(counts = spatial_counts),
    spatialCoords = spatial_coords_matrix
  )
  
  img_list <- readImgData(path = spatial_dir, sample_id = slide_name)
  imgData(spe) <- img_list
  
  
  # pal: set manually according to cell types
  pal_manual <- c(
    "basal_type_a"   = "#9E0000", "basal_type_b"   = "#E74C3C",
    "classic_type_a" = "#FFD700", "classic_type_b" = "#FFF59D",
    "ductal"         = "#FF8C00",
    "endocrine_alpha"= "#2E7D32", "endocrine_beta" = "#4CAF50", "endocrine_delta"= "#A5D6A7",
    "myCAF" = "#4575b4", "iCAF" = "#1976D2", "apCAF" = "#64B5F6", "PSCs_A" = "#AEDFF7",
    "b_cells" = "#7E57C2", "t_cells" = "#512DA8",
    "pi_TAM" = "#8D6E63", "m2_macrophages" = "#BCAAA4", "moMAC" = "#6D4C41",
    "monocytes" = "#8E8E8E", "trMAC" = "#455A64",
    "endothelium" = "#0097A7", "acinar" = "#4DD0E1"
  )
  pal_use <- pal_manual[colnames(mat)]
  
  # FIG 1: ScatterPie
  cat("Scatterpie...\n")
  p_scatter <- plotSpatialScatterpie(
    x = spe, y = mat, cell_types = colnames(mat),
    img = TRUE, scatterpie_alpha = 0.9, pie_scale = 0.4
  ) + scale_fill_manual(values = pal_use) +
    ggtitle(paste("SPOTlight proportions -", slide_name))
  ggsave("scatterpie_spatial.png", plot = p_scatter, width = 7, height = 6, dpi = 300)
  
  # Rotate coords
  cat("Coordenadas preprocess...\n")
  scales <- jsonlite::fromJSON(readLines(scalefactors_path))
  coords <- GetTissueCoordinates(seurat_obj)
  coords$x_new <- coords$y
  coords$y_new <- max(coords$x) - coords$x
  coords <- coords[rownames(mat), , drop = FALSE]
  
  # FIG 2: Heatmaps 
  cat("Heatmaps...\n")
  for (ct in colnames(mat)) {
    df_ct <- coords
    df_ct$score <- mat[, ct]
    safe_ct <- gsub("[/?]", "_", ct)
    p <- ggplot(df_ct, aes(x = x_new, y = y_new, fill = score)) +
      geom_point(shape = 21, size = 1.6) +
      scale_fill_gradientn(colours = c("white", pal_use[ct]), name = ct) +
      coord_fixed() + theme_minimal() +
      ggtitle(paste("Spatial distribution of", ct)) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title = element_blank())
    ggsave(paste0("heatmap_", safe_ct, ".png"), plot = p, width = 6, height = 5, dpi = 300)
  }
  
  # FIG 3: Predominant cell type 
  cat("Predominant cell type map...\n")
  pred_cell <- colnames(mat)[max.col(mat, ties.method = "first")]
  df_pred <- coords
  df_pred$cell_type <- factor(pred_cell, levels = names(pal_manual))
  p_pred <- ggplot(df_pred, aes(x = x_new, y = y_new, fill = cell_type)) +
    geom_point(shape = 21, size = 1.4) +
    scale_fill_manual(values = pal_use) +
    coord_fixed() + theme_minimal() +
    ggtitle("Predominant Cell Type per Spot") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank())
  ggsave("predominant_celltype_rotated.png", plot = p_pred, width = 7, height = 6, dpi = 300)
  
  cat("Slide", slide_name, "done")
  setwd(root_dir)
}
cat("\n Completed.\n")

#######################################################################################################
#######################################################################################################


# POPULATION ANALYSIS 

root_dir <- "/home/martinpl/projects/PDAC/deconv/final"
slides <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
df_all <- data.frame()

for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  res_path <- file.path(slide_dir, "spotlight_result.rds")
  
  if (!file.exists(res_path)) {
    cat(" spotlight_result.rds not found", slide_name, "\n")
    next
  }
  res <- readRDS(res_path)
  mat <- res$mat
  if (is.null(mat) || nrow(mat) == 0) next
  
  mat_norm <- sweep(mat, 1, rowSums(mat), FUN = "/")
  prop_means <- colMeans(mat_norm, na.rm = TRUE)
  
  df_slide <- data.frame(
    slide = slide_name,
    cell_type = names(prop_means),
    proportion = as.numeric(prop_means)
  )
  df_all <- rbind(df_all, df_slide)
}


# FIG: Stacked barplot 
p_bar <- ggplot(df_all, aes(x = slide, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_manual, na.value = "grey70") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "Cell type proportions per slide",
    x = "Slide",
    y = "Proportion"
  )

ggsave(file.path(root_dir, "fig_celltype_proportions.png"),
       plot = p_bar, width = 10, height = 6, dpi = 300)
cat("fig_celltype_proportions.png saved\n")

# FIGURES Basal vs Classic 
df_all <- df_all %>%
  mutate(group = case_when(
    grepl("basal_type", cell_type) ~ "Basal",
    grepl("classic_type", cell_type) ~ "Classic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

prop_df <- df_all %>%
  group_by(slide, group) %>%
  summarise(mean_prop = sum(proportion), .groups = "drop")

p_bar2 <- ggplot(prop_df, aes(x = slide, y = mean_prop, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Basal" = "#B71C1C", "Classic" = "#FFD54F")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Total Basal vs Classic proportions across slides",
    x = "Slide",
    y = "Mean proportion",
    fill = "Subtype"
  )

ggsave(file.path(root_dir, "fig_basal_classic_total_bar.png"),
       plot = p_bar2, width = 8, height = 5, dpi = 300)
cat(" fig_basal_classic_total_bar.png saved\n")

epsilon <- 1e-6

logratio_df <- prop_df %>%
  pivot_wider(names_from = group, values_from = mean_prop) %>%
  mutate(
    log2_basal_classic = log2((Basal + epsilon) / (Classic + epsilon))
  ) %>%
  arrange(log2_basal_classic)

p_logratio <- ggplot(logratio_df,
                     aes(x = log2_basal_classic,
                         y = reorder(slide, log2_basal_classic))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(fill = "#6A1B9A") +
  theme_bw() +
  labs(
    title = "Log2 ratio Basal / Classic per slide",
    x = "log2(Basal / Classic)",
    y = "Slide"
  )

df_spot_lr <- data.frame()
epsilon <- 1e-6

for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  res_path <- file.path(slide_dir, "spotlight_result.rds")
  
  if (!file.exists(res_path)) next
  
  res <- readRDS(res_path)
  mat <- res$mat
  if (is.null(mat) || nrow(mat) == 0) next

  mat_norm <- sweep(mat, 1, rowSums(mat), FUN = "/")
  basal_cols   <- grep("^basal_type", colnames(mat_norm), value = TRUE)
  classic_cols <- grep("^classic_type", colnames(mat_norm), value = TRUE)
  
  if (length(basal_cols) == 0 || length(classic_cols) == 0) next
  
  basal_prop   <- rowSums(mat_norm[, basal_cols, drop = FALSE])
  classic_prop <- rowSums(mat_norm[, classic_cols, drop = FALSE])
  
  log2_bc <- log2((basal_prop + epsilon) / (classic_prop + epsilon))
  
  df_slide_spot <- data.frame(
    slide = slide_name,
    barcode = rownames(mat_norm),
    log2_basal_classic = log2_bc
  )
  
  df_spot_lr <- rbind(df_spot_lr, df_slide_spot)
}

df_spot_props <- data.frame()

for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  res_path <- file.path(slide_dir, "spotlight_result.rds")
  if (!file.exists(res_path)) next
  
  res <- readRDS(res_path)
  mat <- res$mat
  if (is.null(mat) || nrow(mat) == 0) next
  
  mat_norm <- sweep(mat, 1, rowSums(mat), FUN = "/")
  
  basal_cols   <- grep("^basal_type", colnames(mat_norm), value = TRUE)
  classic_cols <- grep("^classic_type", colnames(mat_norm), value = TRUE)
  
  df_slide <- data.frame(
    slide   = slide_name,
    barcode = rownames(mat_norm),
    Basal   = rowSums(mat_norm[, basal_cols,   drop = FALSE]),
    Classic = rowSums(mat_norm[, classic_cols, drop = FALSE])
  )
  
  df_spot_props <- rbind(df_spot_props, df_slide)
}

ggplot(df_spot_props, aes(x = Basal, y = slide)) +
  geom_density_ridges(fill = "#B71C1C", alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Distribution of Basal program across spots",
    x = "Basal proportion",
    y = "Slide"
  )


ggplot(df_spot_props, aes(x = Classic, y = slide)) +
  geom_density_ridges(fill = "#FFD54F", alpha = 0.8) +
  theme_bw() +
  labs(
    title = "Distribution of Classic program across spots",
    x = "Classic proportion",
    y = "Slide"
  )

#######################################################################################################
#######################################################################################################


# FIG: Spot presence

norm_rows <- function(mat) {
  rs <- rowSums(mat)
  rs[rs == 0] <- NA
  sweep(mat, 1, rs, "/")
}

spot_presence_categories <- function(mat_norm, slide_name, th_pres = 0.10) {
  stopifnot(all(c("basal_type_a","basal_type_b","classic_type_a","classic_type_b") %in% colnames(mat_norm)))
  
  basal   <- mat_norm[, "basal_type_a"] + mat_norm[, "basal_type_b"]
  classic <- mat_norm[, "classic_type_a"] + mat_norm[, "classic_type_b"]
  
  basal_present   <- basal   >= th_pres
  classic_present <- classic >= th_pres
  
  catg <- dplyr::case_when(
    basal_present & !classic_present ~ "basal",
    !basal_present & classic_present ~ "classical",
    basal_present & classic_present ~ "classical and basal",
    TRUE ~ "others"
  )
  
  tibble(
    slide = slide_name,
    spot = rownames(mat_norm),
    basal = basal,
    classical = classic,
    category = factor(catg, levels = c("basal","classical","classical and basal","others"))
  )
}

plot_stacked_counts <- function(df_counts, title = NULL) {
  pal_cat <- c(
    "basal" = "#D32F2F",
    "classical" = "#FBC02D",
    "classical and basal" = "#F57C00",
    "others" = "#1F4E8C"
  )
  
  ggplot(df_counts, aes(x = slide, y = n_spots, fill = category)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = pal_cat) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    labs(title = title, x = NULL, y = "Number of spots")
}

root_dir <- "/home/martinpl/projects/PDAC/deconv/final"
slides <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)

th_pres <- 0.25

df_spots_all <- list()

for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  res_path <- file.path(slide_dir, "spotlight_result.rds")
  if (!file.exists(res_path)) next
  
  res <- readRDS(res_path)
  mat <- res$mat
  if (is.null(mat) || nrow(mat) == 0) next
  
  mat_norm <- norm_rows(mat)
  
  df_spots_all[[slide_name]] <- spot_presence_categories(
    mat_norm, slide_name, th_pres = th_pres
  )
}

df_spots_all <- bind_rows(df_spots_all)

df_counts <- df_spots_all %>%
  count(slide, category, name = "n_spots")

pA <- plot_stacked_counts(
  df_counts,
  title = paste0("A) SPOTlight – spot presence (th=", th_pres, ")")
)

ggsave(file.path(root_dir, paste0("A_spotlight_presence_stacked_counts_th", th_pres, ".png")),
       pA, width = 11, height = 5, dpi = 300)


df_counts_prop <- df_counts %>%
  group_by(slide) %>%
  mutate(prop = n_spots / sum(n_spots)) %>%
  ungroup()

pA_prop <- ggplot(df_counts_prop, aes(x = slide, y = prop, fill = category)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  labs(title = paste0("A) SPOTlight – proportion of spots (presence th=", th_pres, ")"),
       x = NULL, y = "Proportion of spots")


