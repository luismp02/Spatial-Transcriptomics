# Deconvolution of 10X Visium samples, with pie plot visualisation and co-localisation matrix

These scripts form part of the ACACIA project and were created during a Master's internship with the Gene Expression and Oncogenesis team at IGDR in Rennes, France. 


### SpaGCN.py 
*SPADE requires clustered data. Therefore, SpaGCN was used for clustering based on the author's recommendation.*

[SpaGCN github repository](https://github.com/jianhuupenn/SpaGCN)
##### Input : 
- preprocessed Spatial expression matrix, .csv 
- spots locations, .csv
- histology image, .tiff
  
##### Output : 
- Spatial expression matrix with additional column = cluster annotation for each spot, results.h5ad
- Cluster visualisation, pred.png

#### Requirements
- Python: 3.7.9
- Python packages: pandas = 1.1.3, numpy = 1.20.2, python-igraph=0.8.3, torch=1.6.0,louvain=0.7.0, scipy = 1.5.2, scanpy = 1.6.0, anndata = 0.7.4, natsort = 7.0.1, sklearn = 0.23.3
-----
### ST_deconv_coloc.Rmd

[SPADE github repository](https://github.com/YyLu5/SPADE)
##### Input : 
- preprocessed + clustered spatial expression matrix, .rds
- spots locations, .csv
- single cell reference matrix, .rds
- Lists of marker genes (one list per cell type), .rds
##### output : 
- Deconvolution matrix (spots X cell types), .csv
- Pie plot visualisation, .pdf
- Co-localisation matrix, .pdf

#### Requirements
- R version >= 4.0.0.
- R packages: Biobase, rlist, caret, glmnet, tidyverse, NOISeq, EBImage

Exemple with SPADE's authors toy dataset is shown in ST_deconv_testdata.html, toy dataset can be downloaded [here](https://filesender.renater.fr/?s=download&token=d816f9b4-9613-49e7-8d6f-ab41efc4c32e)

----
#### Credits : 
SpaGCN : 	Hu, J. et al. SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. Nat. Methods 18, 1342–1351 (2021).

SPADE : Lu, Y., Chen, Q. M. & An, L. SPADE: spatial deconvolution for domain specific cell-type estimation. Commun. Biol. 7, 469 (2024).

