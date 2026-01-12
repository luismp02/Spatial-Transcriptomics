suppressPackageStartupMessages({
  library(SPADE)
  library(Matrix)
  library(SpatialExperiment)
  library(S4Vectors)
  library(zellkonverter)
  library(hdf5r)
})

options(stringsAsFactors = FALSE)

# utils: 
parse_args <- function(x) {
  if (length(x) == 0) return(list())
  if (any(grepl("^--", x)) == FALSE) stop("All arguments must be in '--key value' format.")
  keys <- gsub("^--", "", x[grepl("^--", x)])
  vals <- x[!grepl("^--", x)]

  out <- list()
  i <- 1
  while (i <= length(x)) {
    if (!grepl("^--", x[i])) stop("Unexpected token: ", x[i])
    key <- gsub("^--", "", x[i])
    if (i == length(x)) stop("Missing value for --", key)
    val <- x[i + 1]
    if (grepl("^--", val)) stop("Missing value for --", key)
    out[[key]] <- val
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

# required params 
req <- c("slides", "visium_base", "spagcn_out", "ref_rds", "markers_rds", "out_base")
missing <- req[!req %in% names(args)]
if (length(missing) > 0) {
  stop("Missing required args: ", paste(missing, collapse = ", "), "\n",
       "Usage:\n",
       "  Rscript SPADE.R --slides <slide1,slide2> --visium_base <dir> --spagcn_out <dir>\n",
       "                 --ref_rds <file.rds> --markers_rds <file.rds> --out_base <dir>\n")
}

slides <- strsplit(args$slides, ",")[[1]]
slides <- trimws(slides)
slides <- slides[nzchar(slides)]

visium_base  <- path.expand(args$visium_base)
spagcn_out   <- path.expand(args$spagcn_out)
path_ref_rds <- path.expand(args$ref_rds)
path_markers <- path.expand(args$markers_rds)
out_base     <- path.expand(args$out_base)

# load ref + markers ONCE
message(">>> Loading reference: ", path_ref_rds)
ref <- readRDS(path_ref_rds)

assay_int <- ref@assays[["integrated"]]
expr_ref  <- assay_int@data
meta_ref  <- ref@meta.data

message(">>> Loading fused markers: ", path_markers)
sign_list_fused <- readRDS(path_markers)

# Fusion endocrine/immune 
endocrine_types <- c("endocrine_alpha", "endocrine_beta", "endocrine_delta")
immune_types    <- c("b_cells", "t_cells", "moMAC", "monocytes", "m2_macrophages", "trMAC", "pi_TAM")

celltypes_fine  <- meta_ref$fine_consensus_annotation
celltypes_fused <- celltypes_fine
celltypes_fused[celltypes_fused %in% endocrine_types] <- "endocrine"
celltypes_fused[celltypes_fused %in% immune_types]    <- "immune"
meta_ref$celltype_fused <- celltypes_fused

# scref_fused: promedios por tipo fusionado
scref_fused <- sapply(names(sign_list_fused), function(ct) {
  idx <- which(meta_ref$celltype_fused == ct)
  if (length(idx) == 0) {
    rep(0, nrow(expr_ref))
  } else {
    Matrix::rowMeans(expr_ref[, idx, drop = FALSE])
  }
})
scref_fused <- as.matrix(scref_fused)
rownames(scref_fused) <- rownames(expr_ref)

# -------- loop slides --------
for (sample_id in slides) {
  message("\n>>> Processing slide: ", sample_id)

  path_spade <- file.path(out_base, sample_id)
  dir.create(path_spade, showWarnings = FALSE, recursive = TRUE)

  # Paths per slide
  visium_dir <- file.path(visium_base, sample_id, "outs")
  path_spagcn_h5ad <- file.path(spagcn_out, sample_id, paste0(sample_id, "_results.h5ad"))

  if (!file.exists(file.path(visium_dir, "filtered_feature_bc_matrix.h5"))) {
    stop("Missing Visium H5 for slide ", sample_id, " at: ", visium_dir)
  }
  if (!file.exists(path_spagcn_h5ad)) {
    stop("Missing SpaGCN h5ad for slide ", sample_id, " at: ", path_spagcn_h5ad)
  }

  ## Load Visium H5
  h5_path  <- file.path(visium_dir, "filtered_feature_bc_matrix.h5")
  pos_path <- file.path(visium_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(pos_path)) {
    stop("Missing tissue_positions.csv at: ", pos_path)
  }

  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5$close(), silent = TRUE), add = TRUE)

  data    <- h5[["matrix/data"]][]
  indices <- h5[["matrix/indices"]][]
  indptr  <- h5[["matrix/indptr"]][]
  shape   <- h5[["matrix/shape"]][]

  gene_names <- as.character(h5[["matrix/features/name"]][])
  gene_ids   <- as.character(h5[["matrix/features/id"]][])
  barcodes   <- as.character(h5[["matrix/barcodes"]][])
  h5$close()

  nr <- as.integer(shape[1]); nc <- as.integer(shape[2])
  stcount <- new("dgCMatrix",
    Dim      = c(nr, nc),
    Dimnames = list(gene_names, barcodes),
    x        = as.numeric(data),
    i        = as.integer(indices),
    p        = as.integer(indptr)
  )

  coords <- read.csv(pos_path, header = FALSE)
  colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
  coords <- coords[match(colnames(stcount), coords$barcode), , drop = FALSE]

  spatial_mat <- as.matrix(data.frame(
    pxl_col = as.numeric(coords$pxl_col),
    pxl_row = as.numeric(coords$pxl_row)
  ))

  spe <- SpatialExperiment(
    assays       = list(counts = stcount),
    rowData      = S4Vectors::DataFrame(gene_id = gene_ids, gene_name = gene_names),
    colData      = S4Vectors::DataFrame(
      barcode   = colnames(stcount),
      in_tissue = coords$in_tissue,
      sample_id = sample_id
    ),
    spatialCoords = spatial_mat
  )

  ## SpaGCN domains
  spagcn <- readH5AD(path_spagcn_h5ad)
  spagcn_barcodes <- colnames(spagcn)

  dom_col <- intersect(c("domain", "pred", "spatial_domain"), colnames(colData(spagcn)))
  if (length(dom_col) == 0) {
    stop("No domain column found in SpaGCN h5ad. Expected one of: domain/pred/spatial_domain")
  }

  spagcn_df <- data.frame(
    barcode = spagcn_barcodes,
    domain  = as.integer(colData(spagcn)[[dom_col[1]]])
  )

  # Match barcodes
  bar_common <- intersect(colnames(stcount), spagcn_df$barcode)
  if (length(bar_common) == 0) stop("No common barcodes between Visium and SpaGCN for ", sample_id)

  stcount_ready <- stcount[, bar_common, drop = FALSE]
  spagcn_df     <- spagcn_df[match(bar_common, spagcn_df$barcode), , drop = FALSE]

  # loc table (ACA arreglé tu bug grande)
  loc <- data.frame(
    barcode = colnames(stcount_ready),
    pxl_col = spatial_mat[match(colnames(stcount_ready), colnames(stcount)), 1],
    pxl_row = spatial_mat[match(colnames(stcount_ready), colnames(stcount)), 2]
  )

  loc_spade <- merge(loc, spagcn_df, by = "barcode")
  loc_spade <- loc_spade[match(colnames(stcount_ready), loc_spade$barcode), , drop = FALSE]
  loc_spade$location <- seq_len(nrow(loc_spade))
  loc_ready <- loc_spade

  ## Intersect genes
  common_genes <- intersect(rownames(scref_fused), rownames(stcount_ready))
  if (length(common_genes) < 50) {
    stop("Too few common genes (", length(common_genes), ") for slide ", sample_id)
  }

  scref_ready    <- scref_fused[common_genes, , drop = FALSE]
  stcount_ready  <- stcount_ready[common_genes, , drop = FALSE]
  sign_list_ready <- lapply(sign_list_fused, function(vec) intersect(vec, common_genes))

  stcount_for_spade <- as.matrix(stcount_ready)
  scref_for_spade   <- as.matrix(scref_ready)

  ## CTperDom
  loc_ready$domain <- loc_ready$domain - min(loc_ready$domain, na.rm = TRUE)
  nlay <- max(loc_ready$domain, na.rm = TRUE)
  if (!is.finite(nlay) || nlay < 0) stop("Invalid domain values for slide ", sample_id)

  CTperLayer_ready <- CTperDom(
    loc       = loc_ready,
    stcount   = stcount_for_spade,
    scref     = scref_for_spade,
    sign_list = unlist(sign_list_ready),
    lasso     = TRUE
  )

  saveRDS(CTperLayer_ready, file.path(path_spade, "CTperLayer_ready.rds"))
  saveRDS(loc_ready,         file.path(path_spade, "loc_ready.rds"))

  ## SPADE final
  CTest_ready <- SPADE(
    stcount   = stcount_for_spade,
    scref     = scref_for_spade,
    sign_list = sign_list_ready,
    loc       = loc_ready,
    ctData    = CTperLayer_ready,
    offset    = 0,
    yNorm     = "cpm",
    bNorm     = "cpm"
  )

  saveRDS(CTest_ready, file.path(path_spade, "CTest_ready.rds"))
  message(">>> SPADE done for slide: ", sample_id)
}