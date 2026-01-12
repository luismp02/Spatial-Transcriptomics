
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
})

options(stringsAsFactors = FALSE)

# args
parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    if (!grepl("^--", x[i])) stop("Expected --key, got: ", x[i])
    key <- sub("^--", "", x[i])
    if (i == length(x)) stop("Missing value for --", key)
    val <- x[i + 1]
    if (grepl("^--", val)) stop("Missing value for --", key)
    out[[key]] <- val
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

req <- c("outs", "files", "out_rds")
missing <- req[!req %in% names(args)]
if (length(missing) > 0) {
  stop(
    "Missing args: ", paste(missing, collapse = ", "), "\n",
    "Usage:\n",
    "  Rscript load_data_nf.R --outs <.../outs> --files <dir_with_csvs> --out_rds <output.rds>\n"
  )
}

outs_dir  <- path.expand(args$outs)     # Space Ranger outs
files_dir <- path.expand(args$files)    # dir with CSVs
out_rds   <- path.expand(args$out_rds)

matrix_path  <- file.path(outs_dir, "filtered_feature_bc_matrix")
spatial_path <- file.path(outs_dir, "spatial")

if (!dir.exists(matrix_path)) stop("Missing filtered_feature_bc_matrix at: ", matrix_path)
if (!dir.exists(spatial_path)) stop("Missing spatial at: ", spatial_path)

# -------------------------------
# Step 1: Load raw Visium data
# -------------------------------
expr_matrix <- readMM(file.path(matrix_path, "matrix.mtx.gz"))
barcodes <- read_tsv(file.path(matrix_path, "barcodes.tsv.gz"), col_names = "barcode", show_col_types = FALSE)
features <- read_tsv(
  file.path(matrix_path, "features.tsv.gz"),
  col_names = c("gene_id", "gene_name", "feature_type"),
  show_col_types = FALSE
)

# Clean gene IDs 
features <- as.data.frame(features)
features <- features %>%
  mutate(
    gene_id = str_replace_all(gene_id, "\t|Gene Expression", "") %>%
      str_trim() %>%
      str_to_lower() %>%
      str_sub(1, 9)
  )

expr_matrix <- as(expr_matrix, "dgCMatrix")
rownames(expr_matrix) <- make.unique(features$gene_id)
colnames(expr_matrix) <- barcodes$barcode

cat("Matrix loaded with dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "barcodes\n")

# -------------------------------
# Step 2: Create Seurat object
# -------------------------------
seurat_obj <- CreateSeuratObject(counts = expr_matrix, assay = "Spatial")

# -------------------------------
# Step 3: Load gene signature files
# -------------------------------
gene_list_path <- file.path(files_dir, "liste_gene_and_code.csv")
if (file.exists(gene_list_path)) {
  gene_list <- read_csv(gene_list_path, show_col_types = FALSE)
} else {
  gene_list <- NULL
  warning("liste_gene_and_code.csv not found at: ", gene_list_path)
}

genes_growth_path    <- file.path(files_dir, "de_gene_names_growth.csv")
genes_regen_up_path  <- file.path(files_dir, "ur_regeneration.csv")
genes_regen_down_path<- file.path(files_dir, "dr_regeneration.csv")
genes_scarring_path  <- file.path(files_dir, "gene_name_scaring.csv")

stopifnot(file.exists(genes_growth_path))
stopifnot(file.exists(genes_regen_up_path))
stopifnot(file.exists(genes_regen_down_path))
stopifnot(file.exists(genes_scarring_path))

genes_growth     <- read_csv(genes_growth_path, show_col_types = FALSE)
genes_regen_up   <- read_csv(genes_regen_up_path, show_col_types = FALSE)
genes_regen_down <- read_csv(genes_regen_down_path, show_col_types = FALSE)
genes_scarring   <- read_csv(genes_scarring_path, show_col_types = FALSE)

clean_genes <- function(df) {
  df[[1]] <- df[[1]] %>%
    str_replace_all("\t", "") %>%
    str_trim() %>%
    str_to_lower()
  df
}

genes_growth     <- clean_genes(genes_growth)
genes_regen_up   <- clean_genes(genes_regen_up)
genes_regen_down <- clean_genes(genes_regen_down)
genes_scarring   <- clean_genes(genes_scarring)

# -------------------------------
# Step 4: Add metadata to features
# -------------------------------
features_fixed <- features
rownames(features_fixed) <- rownames(seurat_obj)

growth_vec     <- genes_growth[[1]]
regen_up_vec   <- genes_regen_up[[1]]
regen_down_vec <- genes_regen_down[[1]]
scarring_vec   <- genes_scarring[[1]]

features_fixed$growth_gene     <- rownames(seurat_obj) %in% growth_vec
features_fixed$regen_up_gene   <- rownames(seurat_obj) %in% regen_up_vec
features_fixed$regen_down_gene <- rownames(seurat_obj) %in% regen_down_vec
features_fixed$scarring_gene   <- rownames(seurat_obj) %in% scarring_vec

seurat_obj@misc$gene_annotations <- features_fixed

# Si querés también guardar gene_list en misc
if (!is.null(gene_list)) {
  seurat_obj@misc$gene_list <- gene_list
}

# -------------------------------
# Step 5: Load spatial metadata
# -------------------------------
pos1 <- file.path(spatial_path, "tissue_positions.csv")
pos2 <- file.path(spatial_path, "tissue_positions_list.csv")

if (file.exists(pos1)) {
  positions <- read_csv(pos1, show_col_types = FALSE)
} else if (file.exists(pos2)) {
  # formato viejo: sin header, 6 columnas
  positions <- read_csv(pos2, col_names = FALSE, show_col_types = FALSE)
  colnames(positions) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
} else {
  stop("No tissue_positions file found in: ", spatial_path)
}

if (!"barcode" %in% colnames(positions)) colnames(positions)[1] <- "barcode"

seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
seurat_obj@meta.data <- left_join(seurat_obj@meta.data, positions, by = "barcode")

# -------------------------------
# Step 6: Save object
# -------------------------------
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_obj, file = out_rds)
cat("Saved Seurat object to:", out_rds, "\n")
