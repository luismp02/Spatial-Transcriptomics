suppressPackageStartupMessages({
  library(scatterpie)
  library(ggplot2)
  library(dplyr)
})
library(SPADE)


CTest <- readRDS("CTest_ready.rds")
loc   <- readRDS("loc_ready.rds")

all_spots <- unlist(lapply(CTest, rownames))
all_ct    <- unique(unlist(lapply(CTest, colnames)))

estCT <- matrix(
  0,
  nrow = length(all_spots),
  ncol = length(all_ct),
  dimnames = list(all_spots, all_ct)
)

for (i in seq_along(CTest)) {
  rr <- rownames(CTest[[i]])
  cc <- colnames(CTest[[i]])
  estCT[rr, cc] <- CTest[[i]]
}

## Align loc with estCT per barcode
loc_ord <- loc[match(rownames(estCT), loc$barcode), ]
stopifnot(all(loc_ord$barcode == rownames(estCT)))


## def coords
fuldat <- cbind(loc_ord, estCT)
fuldat$x <- fuldat$pxl_col
fuldat$y <- fuldat$pxl_row
fuldat$y <- max(fuldat$y, na.rm = TRUE) - fuldat$y


#################################################################################################
#################################################################################################

## Scatterpie
library(scatterpie)
library(ggplot2)
library(RColorBrewer)

ct_show <- colnames(estCT) #change according to ref CT

pal_spade <- c(
  "basal_type_a"   = "#9E0000",
  "basal_type_b"   = "#E74C3C",
  "classic_type_a" = "#FFD700",
  "classic_type_b" = "#FFF59D",
  "ductal"         = "#FF8C00",
  "endocrine"      = "#4CAF50",  
  "immune"         = "#6A4AA1",   
  "myCAF"          = "#4575b4",
  "iCAF"           = "#1976D2",
  "apCAF"          = "#64B5F6",
  "PSCs_A"         = "#AEDFF7",
  "endothelium"    = "#0097A7",
  "acinar"         = "#4DD0E1"
)

# spot size (adjust manually)
spot_r <- (max(fuldat$x) - min(fuldat$x)) / 110

ggplot() +
  geom_scatterpie(
    aes(x = x, y = y, r = spot_r),
    data = fuldat,
    cols = ct_show,
    color = NA,
    alpha = 0.90
  ) +
  coord_fixed() +
  scale_fill_manual(values = pal_spade) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

ggsave("SPADE_scatterpie.png", width = 10, height = 10, dpi = 300)


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(SPADE)
})


#################################################################################################
#################################################################################################

# Partition figure (basal/classic/other per slide)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

root_dir <- "/home/martinpl/projects/PDAC/deconv/SPADE"  
out_dir  <- root_dir
th_pres <- 0.25
CTest_candidates <- c("CTest_ready.rds", "CTest.rds", "CTest_spade.rds")
loc_candidates   <- c("loc_ready.rds", "loc.rds", "loc_spade.rds")


# helpers

norm_rows <- function(mat) {
  rs <- rowSums(mat)
  rs[rs == 0] <- NA
  sweep(mat, 1, rs, "/")
}

spot_presence_categories <- function(mat_norm, slide_name, th_pres = 0.25) {
  needed <- c("basal_type_a","basal_type_b","classic_type_a","classic_type_b")
  if (!all(needed %in% colnames(mat_norm))) {
    stop("Missing required columns: ",
         paste(setdiff(needed, colnames(mat_norm)), collapse = ", "))
  }
  
  basal   <- mat_norm[, "basal_type_a"] + mat_norm[, "basal_type_b"]
  classic <- mat_norm[, "classic_type_a"] + mat_norm[, "classic_type_b"]
  
  basal_present   <- basal   >= th_pres
  classic_present <- classic >= th_pres
  
  catg <- dplyr::case_when(
    basal_present & !classic_present ~ "basal",
    !basal_present & classic_present ~ "classical",
    basal_present & classic_present  ~ "classical and basal",
    TRUE ~ "others"
  )
  
  tibble(
    slide = slide_name,
    spot = rownames(mat_norm),
    category = factor(catg, levels = c("basal","classical","classical and basal","others"))
  )
}

find_first_existing <- function(dir_path, candidates) {
  for (fn in candidates) {
    p <- file.path(dir_path, fn)
    if (file.exists(p)) return(p)
  }
  return(NA_character_)
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

## loop slides

slides <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
cat("=== SPADE cohort loop:", length(slides), "slide dirs ===\n")

df_all <- list()

for (slide_dir in slides) {
  slide_name <- basename(slide_dir)
  cat("\nProcessing:", slide_name, "\n")
  
  ctest_path <- find_first_existing(slide_dir, CTest_candidates)
  loc_path   <- find_first_existing(slide_dir, loc_candidates)
  
  if (is.na(ctest_path) || is.na(loc_path)) {
    cat("  -> Missing SPADE files in", slide_name, "(skip)\n")
    next
  }
  
  CTest <- readRDS(ctest_path)
  loc   <- readRDS(loc_path)

  if (is.matrix(CTest)) {
    estCT <- CTest
    
  } else if (is.list(CTest) && all(vapply(CTest, is.matrix, logical(1)))) {
    all_spots <- unlist(lapply(CTest, rownames))
    all_ct    <- unique(unlist(lapply(CTest, colnames)))
    
    estCT <- matrix(
      0,
      nrow = length(all_spots),
      ncol = length(all_ct),
      dimnames = list(all_spots, all_ct)
    )
    
    for (i in seq_along(CTest)) {
      rr <- rownames(CTest[[i]])
      cc <- colnames(CTest[[i]])
      estCT[rr, cc] <- CTest[[i]]
    }
    
  } else {
    cat("  -> Unsupported CTest type in", slide_name, "(skip)\n")
    next
  }
  
  if (nrow(estCT) == 0) {
    cat("  -> Empty estCT (skip)\n")
    next
  }
  
  estCT_norm <- norm_rows(estCT)
  df_slide <- spot_presence_categories(estCT_norm, slide_name, th_pres = th_pres)
  
  df_all[[slide_name]] <- df_slide
  cat("  -> spots:", nrow(df_slide), "\n")
}

df_all <- bind_rows(df_all)

if (nrow(df_all) == 0) stop("No SPADE data found across slides. Check file names/paths.")

df_counts <- df_all %>%
  count(slide, category, name = "n_spots")

pB <- plot_stacked_counts(df_counts, title = paste0("B) SPADE – spot presence (th=", th_pres, ")"))

out_png <- file.path(out_dir, paste0("B_spade_presence_stacked_counts_th", th_pres, ".png"))
ggsave(out_png, pB, width = 11, height = 5, dpi = 300)

cat("\nSaved:", out_png, "\nDone.\n")
#################################################################################################
#################################################################################################