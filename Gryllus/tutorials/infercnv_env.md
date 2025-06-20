# Reproducible Environment Setup for `inferCNV` in R (conda-based)

Tested: **June 2025**, Linux `x86_64`, R `4.3.3`.

---

## 1. Create conda environment

```bash
conda create -n infercnv_env -c conda-forge -c bioconda \
    r-base=4.3.3 \
    r-devtools \
    r-remotes \
    r-optparse \
    r-reshape2 \
    r-igraph \
    r-cluster \
    r-seurat \
    bioconductor-biobase \
    bioconductor-s4vectors \
    bioconductor-singlecellexperiment \
    jags \
    r-rjags \
    git \
    -y

conda activate infercnv_env
```

---

## 2. R package installation (manual fix for known dependency issues)

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Core dependencies for inferCNV + SeuratObject
install.packages(c("Rcpp", "RcppArmadillo", "Matrix"), repos = "https://cloud.r-project.org")

# Downgrade Matrix for SeuratObject compatibility
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz",
                 repos = NULL, type = "source")
```

---

## 3. Solve compilation issues via conda (C libs)

```bash
conda install -c conda-forge libpng zlib hdf5 libxml2 glpk igraph r-igraph
```

```r
# These R packages fail to compile without system libs, so install after fixing
install.packages("png", type = "source")
install.packages("data.table", type = "source")
install.packages("scattermore")
install.packages("digest")
install.packages("cowplot")
install.packages("hdf5r", repos = "https://cloud.r-project.org")
```

---

## 4. Install correct versions of Seurat and SeuratDisk

```r
# Seurat v4 to avoid breaking changes in v5
remotes::install_version("Seurat", version = "4.3.0", repos = "https://cloud.r-project.org")

# Compatible SeuratDisk version (based on commit used in inferCNV Dockerfile)
remotes::install_github("mojaveazure/seurat-disk@877d4e1")
```

---

## 5. Additional packages for inferCNV scripts

```r
install.packages(c(
  "optparse",        
  "parallelDist",     
  "cluster",         
  "igraph",           
  "reshape2"
))
```

---

## Reference

This environment was based on the [official inferCNV Dockerfile](https://github.com/broadinstitute/infercnv/blob/master/Dockerfile)  
Adapted for conda + R without Docker.

---
