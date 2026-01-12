import scCube 
from scCube import scCube
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pandas as pd
import scanpy as sc
import numpy as np
import warnings
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.sparse import csr_matrix
import time
import os
import anndata
warnings.filterwarnings("ignore")

###################- Load data -#####################

sc_adata = sc.read("home/martinpl/projects/datashare/grillus/grillus_anndata.h5ad")
coordinates = pd.read_csv("../../datashare/grillus/outs/spatial/tissue_positions.csv")
sc_adata.obs['x'] = coordinates.set_index('barcode')['pxl_row_in_fullres']
sc_adata.obs['y'] = coordinates.set_index('barcode')['pxl_col_in_fullres']

print(sc_adata.obs[['x', 'y']].head())
output_dir = ".home/martinpl/projects/datashare/grillus/"
#####################################################

###################- Prep data -#####################
model = scCube()

#index cell type and barcode
sc_adata.obs['Cell'] = sc_adata.obs.index
sc_adata.obs['Cell_type'] = sc_adata.obs['sctype_classification']
sc_data = sc_adata.raw.X
sc_meta = sc_adata.obs

###################################################

###################- scCube -######################

# Train VAE and generate cells
generate_sc_meta, generate_sc_data = model.train_vae_and_generate_cell(
    sc_adata=sc_adata,
    celltype_key='Cell_type',
    cell_key='Cell',
    target_num=dict(sc_meta.Cell_type.value_counts() * 3),
    batch_size=512,
    epoch_num=500,
    lr=0.0001,
    hidden_size=128,
    save_model=False,
    used_device='cpu')

# Generate spatial patterns
generate_sc_data_new, generate_sc_meta_new = model.generate_pattern_random(
    generate_sc_data=generate_sc_data,
    generate_sc_meta=generate_sc_meta,
    set_seed=True,
    seed=12345,
    spatial_cell_type=None,
    ### PARAMETER DIMENSIONS ###
    spatial_dim=2,
    spatial_size=70,
    delta=100,
    lamda=1,
)

# Plot
fig, ax = plt.subplots(figsize=(8, 8))
unique_cell_types = generate_sc_meta_new['Cell_type'].unique()
color_map = plt.cm.get_cmap('Set2', len(unique_cell_types))  
cell_type_colors = {cell_type: color_map(i) for i, cell_type in enumerate(unique_cell_types)}
generate_sc_meta_new['color'] = generate_sc_meta_new['Cell_type'].map(cell_type_colors)
scatter = ax.scatter(generate_sc_meta_new['point_x'], generate_sc_meta_new['point_y'],
                     c=generate_sc_meta_new['color'], s=10, alpha=0.6)
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) 
           for color in cell_type_colors.values()]
labels = list(cell_type_colors.keys())

ax.legend(handles, labels, title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_title('Simulated Spatial Distribution of Cell Types')
plt.tight_layout()
plt.show()

###################################################

#Save sim
generate_sc_meta_new_fixed = generate_sc_meta_new.drop(columns=['color'], errors='ignore')
adata_simulated = anndata.AnnData(X=generate_sc_data_new.T, obs=generate_sc_meta_new_fixed)
adata_simulated.write(os.path.join(output_dir, "simulated_Visium_FFPE_V43T08-041_D_sctype.h5ad"))
generate_sc_meta_new.to_csv(os.path.join(output_dir, "simulated_Visium_FFPE_Visium_FFPE_V43T08-041_D_metadata.csv"), index=False)
