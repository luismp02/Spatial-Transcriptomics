import sys
import subprocess
import os
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    subprocess.run(["pip", "install", "--quiet", "scvi-colab"], check=True)
    
    from scvi_colab import install
    install()
    
    subprocess.run(["pip", "install", "--quiet", "git+https://github.com/BayraktarLab/cell2location#egg=cell2location[tutorials]"], check=True)
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import pandas as pd
import json
from PIL import Image
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # Plotting of text for PDFs

#####################################################################################

#################################-- MAPPING --#######################################

#####################################################################################

#Paths
results_folder = '/home/martinpl/projects/datashare/PDAC/results/'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
os.makedirs(results_folder, exist_ok=True)  

###- Load ST -###
adata_vis = sc.read_h5ad("/home/martinpl/projects/datashare/PDAC/Visium_FFPE_V43T08-041_A.h5ad")
adata_vis.var.set_index("gene_ids", drop=True, inplace=True)
adata_vis.var["SYMBOL"] = adata_vis.var.index
adata_vis = sc.read_h5ad("/home/martinpl/projects/datashare/PDAC/Visium_FFPE_V43T08-041_A.h5ad")

if 'SYMBOL' not in adata_vis.var.columns:
    adata_vis.var['SYMBOL'] = adata_vis.var.index
if 'gene_ids' in adata_vis.var.columns:
    adata_vis.var.set_index("gene_ids", drop=True, inplace=True)

###- Load SC -###
adata_ref = sc.read_h5ad("/home/martinpl/projects/datashare/genref_hadaca3/baron.h5ad")
adata_ref.var['SYMBOL'] = adata_ref.var.index
print("adata_vis.obs.columns:", adata_vis.obs.columns)
#Filter SC
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

#####################################################################################

#################################-- ReMODEL --#######################################

#####################################################################################

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['ds_id']
                       )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
mod.train(max_epochs=1000)
fig = plt.figure()  # Create a new figure for plot_history
mod.plot_history(20)
plt.legend(labels=['full data training'])
plt.savefig(f"{results_folder}/plot_history.png")  # Save plot_history figure
plt.close()  # Close the figure to free memory

# Renombrar la columna '_index' a 'features' dentro de '_raw' para evitar conflicto
adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

#export
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 200}
)
# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

fig = plt.figure()  # Create a new figure for plot_QC
mod.plot_QC()
plt.savefig(f"{results_folder}/plot_QC.png")  # Save plot_QC figure
plt.close()  # Close the figure to free memory

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

#####################################################################################

#################################-- NB MODEL --######################################

#####################################################################################

# find shared genes and subset both anndata and reference signatures
adata_vis.var_names = adata_vis.var["SYMBOL"]
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)

adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=10,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=3000,
          # train using full data (batch_size=None)
          batch_size=200,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.savefig(f"{results_folder}/elbo_loss_history.png", dpi=300, bbox_inches='tight')
plt.close()

#In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 200}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file
mod.plot_QC()
plt.savefig(f"{results_folder}/qc_plot2.png", dpi=300, bbox_inches='tight')
plt.close()
