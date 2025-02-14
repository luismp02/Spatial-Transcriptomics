# load_data.py
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import json
import scipy.io
from scipy.io import mmread
from PIL import Image

# PATHS (change to your own paths)
res_path = "../results/Genes"
base_path = "../../datashare/grillus"
matrix_path = f"{base_path}/outs/filtered_feature_bc_matrix"
spatial_path = f"{base_path}/outs/spatial"
data_path = f"{base_path}/files"

# LOAD DATA
matrix = scipy.io.mmread(f"{matrix_path}/matrix.mtx.gz").tocsc()  
barcodes = pd.read_csv(f"{matrix_path}/barcodes.tsv.gz", header=None, names=["barcode"])  
features = pd.read_csv(f"{matrix_path}/features.tsv.gz", header=None, names=["gene_id", "gene_name", "feature_type"])  
gene_list = pd.read_csv(f"{data_path}/liste_gene_and_code.csv")  
genes_growth = pd.read_csv(f"{data_path}/de_gene_names_growth.csv") 
genes_regen_up = pd.read_csv(f"{data_path}/ur_regeneration.csv")  
genes_regen_down = pd.read_csv(f"{data_path}/dr_regeneration.csv")  
genes_scarring = pd.read_csv(f"{data_path}/gene_name_scaring.csv")  

# Clean id and names (indeixation may not work if you try to read every column directly)
features['gene_id'] = features['gene_id'].str.replace(r'\t', '', regex=True).str.replace('Gene Expression', '', regex=True).str.strip().str.lower()
features['gene_id'] = features['gene_id'].str.slice(0, 9)

for df in [genes_growth, genes_regen_up, genes_regen_down, genes_scarring]:
    df[df.columns[0]] = df[df.columns[0]].str.replace(r'\t', '', regex=True).str.strip().str.lower()

# Create AnnData object with index setting 
# IMPORTANT: verify if your data works with gene_name or gene_id
adata = anndata.AnnData(X=matrix.T, obs=barcodes, var=features.set_index('gene_id'))
adata.var.index = adata.var.index.str.slice(0, 9)
genes_regen_up = genes_regen_up.set_index('gene_name')
genes_regen_down = genes_regen_down.set_index('gene_name')
genes_scarring = genes_scarring.set_index('gene_name')

adata.var['growth_gene'] = adata.var.index.isin(genes_scarring.index)
adata.var['regen_up_gene'] = adata.var.index.isin(genes_regen_up.index)
adata.var['regen_down_gene'] = adata.var.index.isin(genes_regen_down.index)
adata.var['scarring_gene'] = adata.var.index.isin(genes_scarring.index)

# Load spatial data
tissue_positions = pd.read_csv(f"{spatial_path}/tissue_positions.csv")
adata.obs = adata.obs.merge(tissue_positions, left_on="barcode", right_on="barcode", how="left")

# Save AnnData object
adata.write_h5ad(f"{base_path}/grillus_anndata.h5ad")
