suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(scater)
  library(scran)
  library(SPOTlight)
  library(Matrix)
  library(Seurat)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
gc()

# parse args
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

# Required:
# --ref_rds       : RDS with SCE ref (assay integrated/logcounts)
# --visium_outs   : .../<slide>/outs
# --sample_id     : slide name (Visium_FFPE_...)
# --out_base      : out dir
req <- c("ref_rds", "visium_outs", "sample_id", "out_base")
missing <- req[!req %in% names(args)]
if (length(missing) > 0) {
  stop(
    "Missing args: ", paste(missing, collapse = ", "), "\n",
    "Usage:\n",
    "  Rscript SPOTlight.R \\\n",
    "    --ref_rds <ref.rds> \\\n",
    "    --visium_outs <.../slide/outs> \\\n",
    "    --sample_id <slide_name> \\\n",
    "    --out_base <outdir>\n"
  )
}

ref_path   <- path.expand(args$ref_rds)
visium_outs <- path.expand(args$visium_outs)
sample_id  <- args$sample_id
out_base   <- path.expand(args$out_base)

# Output folder per slide
slide_output <- file.path(out_base, sample_id)
dir.create(slide_output, recursive = TRUE, showWarnings = FALSE)

cat(">> sample_id:", sample_id, "\n")
cat(">> ref_path:", ref_path, "\n")
cat(">> visium_outs:", visium_outs, "\n")
cat(">> slide_output:", slide_output, "\n")

# Load ref 
cat(">> Loading ref...\n")
sce_ref <- readRDS(ref_path)
stopifnot("fine_consensus_annotation" %in% colnames(colData(sce_ref)))

sce_ref$fine_consensus_annotation <- droplevels(factor(sce_ref$fine_consensus_annotation))
sce_ref <- sce_ref[, !is.na(sce_ref$fine_consensus_annotation)]
colLabels(sce_ref) <- sce_ref$fine_consensus_annotation

# Load Visium 

slice_name <- sub("^Visium_FFPE_", "", sample_id)

cat(">> Loading Visium...\n")
seurat_obj <- Load10X_Spatial(
  data.dir = visium_outs,
  slice = slice_name
)

# Fix assay RNA from Spatial and preprocess
seurat_obj[["RNA"]] <- seurat_obj[["Spatial"]]
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Intersección genes 
genes_ref     <- rownames(sce_ref)
genes_spatial <- rownames(seurat_obj)
common_genes  <- intersect(genes_ref, genes_spatial)

cat(">> Genes in ref:", length(genes_ref),
    " | genes in spatial:", length(genes_spatial),
    " | comunes:", length(common_genes), "\n")

if (length(common_genes) < 500) {
  stop("Not enough intersection genes with the reference: ", length(common_genes))
}

sce_ref    <- sce_ref[common_genes, , drop = FALSE]
seurat_obj <- subset(seurat_obj, features = common_genes)
gc()

# HVGs on SCE 
genes_keep <- !grepl("^Rp[l|s]|^Mt", rownames(sce_ref))
colLabels(sce_ref) <- sce_ref$fine_consensus_annotation

dec <- modelGeneVar(sce_ref, subset.row = genes_keep, assay.type = "logcounts")

n_hvg <- 2000L
max_hvg <- min(n_hvg, nrow(sce_ref))
hvg <- getTopHVGs(dec, n = max_hvg)

cat(">> HVGs selected:", length(hvg), "\n")

# Markers 
mgs_list <- scoreMarkers(sce_ref, subset.row = genes_keep)

mgs_fil <- lapply(names(mgs_list), function(i) {
  x <- mgs_list[[i]]
  x <- x[x$mean.AUC > 0.7, ]   # ajustable
  if (nrow(x) == 0) return(NULL)
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})

mgs_fil <- Filter(Negate(is.null), mgs_fil)
mgs_df <- do.call(rbind, mgs_fil)
mgs_df$cluster <- as.character(mgs_df$cluster)

cat(">> Final markers:", nrow(mgs_df),
    "genes in", length(unique(mgs_df$cluster)), "clusters\n")

# Downsampling 
cat(">> Downsampling of ref identity\n")
set.seed(123)
idx <- split(seq_len(ncol(sce_ref)), sce_ref$fine_consensus_annotation)
n_cells <- 80L
cs_keep <- lapply(idx, function(i) sample(i, min(length(i), n_cells)))
sce_ref <- sce_ref[, unlist(cs_keep)]

# SpatialExperiment 
spatial_counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")
spatial_coords <- GetTissueCoordinates(seurat_obj)
spatial_coords_matrix <- as.matrix(spatial_coords[, c("x", "y")])

spe <- SpatialExperiment(
  assays = list(counts = spatial_counts),
  spatialCoords = spatial_coords_matrix
)

#  logcounts ---
if (!"logcounts" %in% assayNames(sce_ref)) {
  message("logcounts not found, creating from counts.")
  assay(sce_ref, "logcounts") <- log1p(assay(sce_ref, "counts"))
}

# Normalize non-negative for NMF
logmat <- assay(sce_ref, "logcounts")
logmat[!is.finite(logmat)] <- 0
logmat[logmat < 0] <- 0
assay(sce_ref, "logcounts") <- logmat

if ("counts" %in% assayNames(sce_ref)) {
  counts_mat <- assay(sce_ref, "counts")
  counts_mat[!is.finite(counts_mat)] <- 0
  counts_mat[counts_mat < 0] <- 0
  assay(sce_ref, "counts") <- counts_mat
}

#  SPOTlight 
cat(">> Running SPOTlight...\n")
res <- SPOTlight(
  x = sce_ref,
  y = spe,
  groups    = as.character(sce_ref$fine_consensus_annotation),
  mgs       = mgs_df,
  hvg       = hvg,
  weight_id = "mean.AUC",
  group_id  = "cluster",
  gene_id   = "gene"
)

mat <- res$mat
mat_norm <- sweep(mat, 1, rowSums(mat), FUN = "/")
res$mat <- mat_norm

# Save outputs (in slide_output) 
out_rds   <- file.path(slide_output, "spotlight_result.rds")
out_mod   <- file.path(slide_output, "spotlight_model_nmf.rds")
out_seu   <- file.path(slide_output, "seurat_obj.rds")
out_ct    <- file.path(slide_output, "celltypes_ordered.rds")

saveRDS(res, out_rds);                cat(">> Saved:", out_rds, "\n")
saveRDS(res$mod, out_mod);            cat(">> Saved:", out_mod, "\n")
saveRDS(seurat_obj, out_seu);         cat(">> Saved:", out_seu, "\n")
saveRDS(colnames(res$mat), out_ct);   cat(">> Saved:", out_ct, "\n")

cat(">>> SPOTlight done for slide:", sample_id, "\n")