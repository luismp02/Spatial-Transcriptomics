#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(SeuratDisk)
  library(igraph)
  library(RColorBrewer)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(scater)
  library(scran)
  library(SPOTlight)
  library(scatterpie)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Slides to process
slides <- c(
  "Visium_FFPE_V44L22-378_D", "Visium_FFPE_V44L23-391_D",
  "Visium_FFPE_V43T08-041_A", "Visium_FFPE_V44L23-362_A",
  "Visium_FFPE_V43T08-041_D", "Visium_FFPE_V44L23-362_D",
  "Visium_FFPE_V43T08-051_A", "Visium_FFPE_V44L22-378_A", "Visium_FFPE_V44L23-391_A"
)

# Paths
ref_base_path   <- "/home/martinpl/projects/datashare"
ref_file        <- "integrated_sce.rds"
visium_base     <- "/home/martinpl/projects/datashare/visium_PDAC"
slide_base      <- "/home/martinpl/projects/PDAC/deconv/final"
script_spotlight <- "/home/martinpl/projects/PDAC/deconv/SPOTlight_integrated.R"

slide_base <- path.expand(slide_base)

# Loop
for (slide in slides) {
  cat("\n Slide:", slide, "\n")

  slide_output     <- file.path(slide_base, slide)
  spatial_data_dir <- file.path(visium_base, slide, "outs")
  spatial_img_path <- file.path(spatial_data_dir, "spatial", "tissue_hires_image.png")
  spatial_dir      <- file.path(spatial_data_dir, "spatial")

  dir.create(slide_output, recursive = TRUE, showWarnings = FALSE)
  setwd(slide_output)

  args_spot <- c(file.path(ref_base_path, ref_file), slide_output, spatial_img_path, spatial_dir)
  cat("  CMD: Rscript", script_spotlight, paste(args_spot, collapse = " "), "\n")
  system2("Rscript", args = c(script_spotlight, args_spot))

  gc()
}
