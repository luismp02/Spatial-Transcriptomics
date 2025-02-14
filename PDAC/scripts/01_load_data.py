import os
import scanpy as sc
import numpy as np
import pandas as pd
import json
import gzip
import matplotlib.pyplot as plt
from skimage import io
from pathlib import Path

################################ PARAMETERS #######################################
DATA_DIR = "../../datashare/visium_PDAC"
RESULTS_DIR = "../results"
REF_DIR = "../../datashare/sc_PDAC"
os.makedirs(RESULTS_DIR, exist_ok=True)
###################################################################################

# 1. Load data
def load_visium_data(sample_path):
    matrix_dir = os.path.join(sample_path, "outs", "filtered_feature_bc_matrix")
    if not os.path.exists(matrix_dir):
        matrix_dir = os.path.join(sample_path, "outs_old", "filtered_feature_bc_matrix")
    
    spatial_dir = os.path.join(sample_path, "outs", "spatial")
    if not os.path.exists(spatial_dir):
        spatial_dir = os.path.join(sample_path, "outs_old", "spatial")
    
    adata = sc.read_10x_mtx(
        matrix_dir,
        var_names='gene_symbols',
        cache=True
    )
    tissue_positions = pd.read_csv(
        os.path.join(spatial_dir, "tissue_positions.csv"),
        header=None, index_col=0
    )
    
    with open(os.path.join(spatial_dir, "scalefactors_json.json"), 'r') as f:
        scale_factors = json.load(f)
      
    adata.obs["x"] = tissue_positions[4]
    adata.obs["y"] = tissue_positions[5]
    adata.uns["scale_factors"] = scale_factors
    
    return adata

# 2. Adjust scale of each spot
def get_spot_diameter(sample_path):
    spatial_dir = os.path.join(sample_path, "outs", "spatial")
    if not os.path.exists(spatial_dir):
        spatial_dir = os.path.join(sample_path, "outs_old", "spatial")
    
    with open(os.path.join(spatial_dir, "scalefactors_json.json"), 'r') as f:
        scale_factors = json.load(f)
        spot_diameter = scale_factors.get("spot_diameter_fullres", None)
    if spot_diameter is None:
        raise ValueError(f"not found {sample_path}")
    
    return spot_diameter

# 3. Load samples
samples = [d for d in os.listdir(DATA_DIR) if d.startswith("Visium_FFPE_")]
for sample in samples:
    sample_path = os.path.join(DATA_DIR, sample)
    adata = load_visium_data(sample_path)
    adata.write(os.path.join(RESULTS_DIR, f"{sample}.h5ad"))


# 4. Load reference data
def load_sc_dataset(dataset_path):
    adatas = []
    for sample in os.listdir(dataset_path):
        sample_path = os.path.join(dataset_path, sample)
        
        barcodes_file = sorted([f for f in os.listdir(sample_path) if "barcodes.tsv.gz" in f])[0]
        features_file = sorted([f for f in os.listdir(sample_path) if "features.tsv.gz" in f])[0]
        matrix_file = sorted([f for f in os.listdir(sample_path) if "matrix.mtx.gz" in f])[0]

        barcodes = pd.read_csv(os.path.join(sample_path, barcodes_file), header=None, sep="\t")
        features = pd.read_csv(os.path.join(sample_path, features_file), header=None, sep="\t")
        matrix = sc.read_mtx(os.path.join(sample_path, matrix_file)).X

        adata = sc.AnnData(X=matrix)
        adata.obs_names = barcodes[0].values
        adata.var_names = features[1].values  #gene symbols are in the second column
        adata.var["gene_id"] = features[0].values  # Store original gene IDs
        
        adata.obs["sample"] = sample

        adatas.append(adata)

    return adatas

# subdirectories containing single-cell data
sc_datasets = [
    os.path.join(REF_DIR, "GSE194247_RAW"),
    os.path.join(REF_DIR, "GSE235449_RAW")
]
adata_list = []
for dataset in sc_datasets:
    adata_list.extend(load_sc_dataset(dataset))
adata_ref = sc.concat(adata_list, axis=0, join="outer")

# Save as h5ad file for cell2location input
output_file = os.path.join(REF_DIR, "sc_PDAC_reference.h5ad")
adata_ref.write(output_file)
print(f"Reference dataset saved as {output_file}")
