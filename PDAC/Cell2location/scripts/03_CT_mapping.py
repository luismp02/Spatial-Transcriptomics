import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from PIL import Image
from matplotlib import rcParams
from cell2location.plt import plot_spatial
rcParams['pdf.fonttype'] = 42  # Plotting of text for PDFs

#####################################################################################

################################-- LOAD DATA --######################################

#####################################################################################

results_folder = '/home/martinpl/projects/datashare/PDAC/results/'
spatial_data_dir = '/home/martinpl/projects/datashare/visium_PDAC/Visium_FFPE_V43T08-041_A/outs/spatial/'

adata_vis = sc.read(f'{results_folder}/cell2location_map/sp.h5ad')  
tissue_positions = pd.read_csv(f'{spatial_data_dir}tissue_positions.csv')

# Filter common barcodes
common_barcodes = np.intersect1d(tissue_positions['barcode'].values, adata_vis.obs.index.values)
tissue_positions_filtered = tissue_positions[tissue_positions['barcode'].isin(common_barcodes)]
adata_vis = adata_vis[adata_vis.obs.index.isin(common_barcodes)]

# Asign positions
adata_vis.obsm['spatial'] = tissue_positions_filtered[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
print(adata_vis.obsm['spatial'].shape)

cytassist_image = Image.open(f'{spatial_data_dir}cytassist_image.tiff')
tissue_hires_image = Image.open(f'{spatial_data_dir}tissue_hires_image.png')

with open(f'{spatial_data_dir}scalefactors_json.json') as f:
    scale_factors = json.load(f)

adata_vis.uns['spatial_image'] = {
    'cytassist': cytassist_image,
    'tissue_hires': tissue_hires_image
}

adata_vis.uns['scale_factors'] = scale_factors

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = np.log1p(
    adata_vis.obs[adata_vis.uns['mod']['factor_names']]
)

# plot in spatial coordinates
with plt.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_vis, cmap='magma',
                  color=['acinar', 'ductal', 'endothelial', 'epsilon',
                         'schwann', 't_cell', 'gamma', 'alpha'],
                  ncols=4, size=1.3,
                  img_key='tissue_hires',  
                  vmin=0, vmax='p99.2',
                  spot_size=14.907718000088739)
    plt.savefig(f"{results_folder}/spatial_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

# Labels
clust_labels = ['acinar', 'ductal', 't_cell']
clust_col = ['' + str(i) for i in clust_labels] 

adata_vis.uns["spatial"] = {
    "PDAC_sample": {  
        "images": {
             "hires": np.array(adata_vis.uns["spatial_image"]["tissue_hires"]),  
            "cytassist": np.array(adata_vis.uns["spatial_image"]["cytassist"])
        },
        "scalefactors": adata_vis.uns["scale_factors"]
    }
}

with plt.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,  
        style='fast',  
        max_color_quantile=0.992, 
        circle_diameter=6,  
        colorbar_position='right'
    )
    fig.savefig(f"{results_folder}/multiple_cell_types_spatial_plot.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
