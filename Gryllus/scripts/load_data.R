library(Matrix)
library(Seurat)
library(dplyr)
library(readr)
library(stringr)

# -------------------------------
# Step 1: Load raw Visium data
# -------------------------------

# Define paths
base_path <- "../../datashare/grillus"
matrix_path <- file.path(base_path, "outs", "filtered_feature_bc_matrix")
spatial_path <- file.path(base_path, "outs", "spatial")
data_path <- file.path(base_path, "files")

# Load matrix and metadata
expr_matrix <- readMM(file.path(matrix_path, "matrix.mtx.gz"))
barcodes <- read_tsv(file.path(matrix_path, "barcodes.tsv.gz"), col_names = "barcode")
features <- read_tsv(file.path(matrix_path, "features.tsv.gz"), 
                     col_names = c("gene_id", "gene_name", "feature_type"))

# Clean gene IDs
features <- as.data.frame(features)
features <- features %>%
  mutate(
    gene_id = str_replace_all(gene_id, "\t|Gene Expression", "") %>%
      str_trim() %>%
      str_to_lower() %>%
      str_sub(1, 9)
  )

# Apply consistent gene IDs to matrix
expr_matrix <- as(expr_matrix, "dgCMatrix")
rownames(expr_matrix) <- make.unique(features$gene_id)
colnames(expr_matrix) <- barcodes$barcode

cat("Matrix loaded with dimensions:", dim(expr_matrix)[1], "genes x", dim(expr_matrix)[2], "barcodes\n")

# -------------------------------
# Step 2: Create Seurat object
# -------------------------------
seurat_obj <- CreateSeuratObject(counts = expr_matrix, assay = "Spatial")

# -------------------------------
# Step 3: Load gene signature files
# -------------------------------
gene_list <- read_csv(file.path(data_path, "liste_gene_and_code.csv"))
genes_growth <- read_csv(file.path(data_path, "de_gene_names_growth.csv"))
genes_regen_up <- read_csv(file.path(data_path, "ur_regeneration.csv"))
genes_regen_down <- read_csv(file.path(data_path, "dr_regeneration.csv"))
genes_scarring <- read_csv(file.path(data_path, "gene_name_scaring.csv"))

# Clean gene names
clean_genes <- function(df) {
  df[[1]] <- df[[1]] %>% 
    str_replace_all("\t", "") %>% 
    str_trim() %>% 
    str_to_lower()
  df
}

genes_growth <- clean_genes(genes_growth)
genes_regen_up <- clean_genes(genes_regen_up)
genes_regen_down <- clean_genes(genes_regen_down)
genes_scarring <- clean_genes(genes_scarring)

# -------------------------------
# Step 4: Add metadata to features
# -------------------------------
features_fixed <- features
rownames(features_fixed) <- rownames(seurat_obj)  # Ensure match

features_fixed$growth_gene <- rownames(seurat_obj) %in% genes_growth[[1]]
features_fixed$regen_up_gene <- rownames(seurat_obj) %in% genes_regen_up[[1]]
features_fixed$regen_down_gene <- rownames(seurat_obj) %in% genes_regen_down[[1]]
features_fixed$scarring_gene <- rownames(seurat_obj) %in% genes_scarring[[1]]

seurat_obj@misc$gene_annotations <- features_fixed
gene_annot <- seurat_obj@misc$gene_annotations

# -------------------------------
# Step 5: Load spatial metadata
# -------------------------------
positions <- read_csv(file.path(spatial_path, "tissue_positions.csv"))
colnames(positions)[1] <- "barcode"

# Merge with Seurat object metadata
seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
seurat_obj@meta.data <- left_join(seurat_obj@meta.data, positions, by = "barcode")


# -------------------------------
# Step 6: Save object
# -------------------------------
saveRDS(seurat_obj, file = file.path(base_path, "grillus_seurat.rds"))

