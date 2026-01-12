# gene_analysis.py
import scanpy as sc
import anndata
import numpy as np
import matplotlib.pyplot as plt
import json
from PIL import Image
import seaborn as sns
import pandas as pd
import scipy
import scipy.spatial
import sklearn
import sklearn.cluster
import networkx as nx
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu

# PATHS (change to your own paths)
res_path = "../results/Genes"
base_path = "../../datashare/grillus"
spatial_path = f"{base_path}/outs/spatial"
adata_path = f"{base_path}/grillus_anndata.h5ad"

# Load Data (AnnData object + spatial info)
adata = sc.read_h5ad(adata_path)
img = Image.open(f"{spatial_path}/tissue_hires_image.png")
with open(f"{spatial_path}/scalefactors_json.json") as f:
    scale_factors = json.load(f)
scale_factor = scale_factors['tissue_hires_scalef']
adata.obs['X_pixel'] = adata.obs['pxl_col_in_fullres'] * scale_factor  
adata.obs['Y_pixel'] = adata.obs['pxl_row_in_fullres'] * scale_factor  

# PARAMETERS #
#Set Thresold and minimum genes (for the moment it is arbitrary only for testing and debugging) 
def threshold_expression(adata, gene_category, quantile=0.75, min_genes=1):
    genes = adata.var[gene_category].values  
    if sum(genes) == 0:
        return np.zeros(adata.shape[0], dtype=bool) 
    
    expr = adata[:, genes].X  
    expr = expr.toarray() if hasattr(expr, "toarray") else expr 
    threshold = np.quantile(expr[expr > 0], quantile, axis=0)  
    expr_binary = (expr > threshold).sum(axis=1) >= min_genes 
    return expr_binary


# Filtering genes + get coordinates
growth_cells_expr = threshold_expression(adata, 'growth_gene', quantile=0.95, min_genes=2)
regen_up_cells_expr = threshold_expression(adata, 'regen_up_gene', quantile=0.80, min_genes=3)
regen_down_cells_expr = threshold_expression(adata, 'regen_down_gene', quantile=0.50, min_genes=3)
scarring_cells_expr = threshold_expression(adata, 'scarring_gene', quantile=0.95, min_genes=1)

growth_cells_positions = adata.obs.loc[growth_cells_expr, ['X_pixel', 'Y_pixel']]
regen_up_cells_positions = adata.obs.loc[regen_up_cells_expr, ['X_pixel', 'Y_pixel']]
regen_down_cells_positions = adata.obs.loc[regen_down_cells_expr, ['X_pixel', 'Y_pixel']]
scarring_cells_positions = adata.obs.loc[scarring_cells_expr, ['X_pixel', 'Y_pixel']]

# Set image orientation: it may change depending on the original rotation and transposition.
# HERE TRANSPOSITION FOR MIRRORED IMAGE + 90 DEGREES ROTATION
for positions in [growth_cells_positions, regen_up_cells_positions, regen_down_cells_positions, scarring_cells_positions]:
    positions[['X_pixel', 'Y_pixel']] = positions[['Y_pixel', 'X_pixel']].values
    positions['X_pixel'], positions['Y_pixel'] = positions['Y_pixel'], img.size[1] - positions['X_pixel']


# PLOT FUNCTIONAL REGIONS #
plt.figure(figsize=(8, 8))
plt.scatter(growth_cells_positions['X_pixel'], growth_cells_positions['Y_pixel'], color='green', label='Growth', alpha=0.5)
plt.scatter(regen_up_cells_positions['X_pixel'], regen_up_cells_positions['Y_pixel'], color='blue', label='Regeneration Up', alpha=0.5)
plt.scatter(regen_down_cells_positions['X_pixel'], regen_down_cells_positions['Y_pixel'], color='red', label='Regeneration Down', alpha=0.5)
plt.scatter(scarring_cells_positions['X_pixel'], scarring_cells_positions['Y_pixel'], color='purple', label='Scarring', alpha=0.5)
plt.imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
plt.legend()
plt.savefig(f"{res_path}/scatter_plot.png", dpi=300)
plt.show()

# 2. One image per category
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
img = Image.open(f"{spatial_path}/tissue_hires_image.png")
axes[0, 0].imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
axes[0, 0].scatter(growth_cells_positions['X_pixel'], growth_cells_positions['Y_pixel'], 
                   color='green', label='Growth', alpha=0.5)
axes[0, 0].set_title("Growth")
axes[0, 0].legend()
axes[0, 1].imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
axes[0, 1].scatter(regen_up_cells_positions['X_pixel'], regen_up_cells_positions['Y_pixel'], 
                   color='blue', label='Regeneration Up', alpha=0.5)
axes[0, 1].set_title("Regeneration Up")
axes[0, 1].legend()
axes[1, 0].imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
axes[1, 0].scatter(regen_down_cells_positions['X_pixel'], regen_down_cells_positions['Y_pixel'], 
                   color='red', label='Regeneration Down', alpha=0.5)
axes[1, 0].set_title("Regeneration Down")
axes[1, 0].legend()
axes[1, 1].imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
axes[1, 1].scatter(scarring_cells_positions['X_pixel'], scarring_cells_positions['Y_pixel'], 
                   color='purple', label='Scarring', alpha=0.5)
axes[1, 1].set_title("Scarring")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(f"{res_path}/scatter_plots_2x2.png", dpi=300)
plt.show()

'''
### EXPLORATORY ANLYSIS ###

# 1 Correlation (matrix, heatmap and mapping)

genes_to_analyze = adata.var.index[adata.var[['growth_gene', 'regen_up_gene', 'regen_down_gene', 'scarring_gene']].any(axis=1)]
expr_matrix = adata[:, genes_to_analyze].X.toarray() if hasattr(adata[:, genes_to_analyze].X, "toarray") else adata[:, genes_to_analyze].X
gene_correlation = np.corrcoef(expr_matrix.T)

plt.figure(figsize=(10, 8))
sns.heatmap(gene_correlation, xticklabels=genes_to_analyze, yticklabels=genes_to_analyze, cmap="coolwarm", center=0)
plt.title("Coexpression matrix")
plt.xticks(rotation=90)
plt.savefig(f"{res_path}/gene_correlation_heatmap.png", dpi=300, bbox_inches="tight")
plt.show()

# example for 2 pairs of genes
top_correlated_indices = np.unravel_index(np.argsort(gene_correlation, axis=None)[-10:], gene_correlation.shape)
top_gene_pairs = [(genes_to_analyze[i], genes_to_analyze[j]) for i, j in zip(*top_correlated_indices) if i != j][:4]
fig, axes = plt.subplots(1, 2, figsize=(12, 6)) 

for ax, (gene1, gene2) in zip(axes.flatten(), top_gene_pairs[:2]):  
    scatter1 = ax.scatter(adata.obs['X_pixel'], adata.obs['Y_pixel'], c=adata[:, gene1].X.toarray().flatten(), cmap="Reds", alpha=0.5, label=gene1)
    scatter2 = ax.scatter(adata.obs['X_pixel'], adata.obs['Y_pixel'], c=adata[:, gene2].X.toarray().flatten(), cmap="Blues", alpha=0.5, label=gene2)
    
    ax.imshow(img, extent=[0, img.size[0], 0, img.size[1]], alpha=0.2)
    ax.set_title(f"Coexp: {gene1} vs {gene2}")
    ax.legend()

plt.tight_layout()
plt.savefig(f"{res_path}/coexpression_top2.png", dpi=300)
plt.show()


# 2. Spatial distribution of categories
plt.figure(figsize=(10, 8))
category_colors = ['red', 'green', 'blue', 'purple']  
categories = ['growth_gene', 'regen_up_gene', 'regen_down_gene', 'scarring_gene']
cells_positions = [growth_cells_positions, regen_up_cells_positions, regen_down_cells_positions, scarring_cells_positions]

for idx, (category, positions) in enumerate(zip(categories, cells_positions)):
    rotated_X = positions['Y_pixel'] 
    rotated_Y = -positions['X_pixel'] 
    plt.scatter(rotated_X, rotated_Y, c=category_colors[idx], alpha=0.5, label=category)

plt.legend()
plt.title('Spatial distribution of spots by fonctional category')
plt.xlabel('Location X')
plt.ylabel('Location Y')
plt.savefig(f"{res_path}/spatial_distribution.png", dpi=300)
plt.show()

# 3. CLustering of expression 
expr_matrix = adata.X.toarray()
kmeans = KMeans(n_clusters=4)
clusters = kmeans.fit_predict(expr_matrix)
adata.obs['clusters'] = clusters
plt.figure(figsize=(10, 8))
rotated_X = adata.obs['X_pixel']
rotated_Y = -adata.obs['Y_pixel'] 

scatter = plt.scatter(rotated_X, rotated_Y, c=adata.obs['clusters'], cmap='tab10', alpha=0.7)
plt.title('Clustering of spots by gene expression')
plt.colorbar(scatter)
plt.xlabel('Location X')
plt.ylabel('Location Y')
plt.savefig(f"{res_path}/clustering.png", dpi=300)
plt.show()


# 4. PCA
expr_matrix = StandardScaler().fit_transform(adata.X.toarray())
pca = PCA(n_components=2)
pca_result = pca.fit_transform(expr_matrix)
plt.figure(figsize=(10, 8))
rotated_X = adata.obs['X_pixel']  
rotated_Y = -adata.obs['Y_pixel']
plt.scatter(rotated_X, rotated_Y, c=pca_result[:, 0], cmap='coolwarm', alpha=0.7)
plt.colorbar()
plt.title('PCA Gene expression')
plt.savefig(f"{res_path}/PCA.png", dpi=300)
plt.show()

# 5. Variabilty
var_expression = {}
for category, positions in zip(categories, cells_positions):
    gene_list = adata.var.index[adata.var[category]]
    expression_values = adata[:, gene_list].X.toarray()
    var_expression[category] = expression_values.std(axis=1)
plt.figure(figsize=(10, 8))
for idx, (category, values) in enumerate(var_expression.items()):
    rotated_X = adata.obs['X_pixel']
    rotated_Y = -adata.obs['Y_pixel']  
    plt.scatter(rotated_X, rotated_Y, c=values, cmap='plasma', alpha=0.5, label=category)
plt.legend()
plt.title('Variability of gene expression by category')
plt.savefig(f"{res_path}/variability.png", dpi=300)
plt.show()

# 6. Coexp network (example within the category)
categories = ['growth', 'regen_up', 'regen_down', 'scarring']
category_genes = {
    'growth': adata.var.index[adata.var['growth_gene']], 
    'regen_up': adata.var.index[adata.var['regen_up_gene']], 
    'regen_down': adata.var.index[adata.var['regen_down_gene']],  
    'scarring': adata.var.index[adata.var['scarring_gene']]  
}

expr_matrix = adata.X.toarray()
for category in categories:
    genes_in_category = category_genes[category]
    expr_subset = expr_matrix[:, [adata.var_names.get_loc(gene) for gene in genes_in_category]]
    correlation_matrix = np.corrcoef(expr_subset.T)
    G = nx.Graph()
    for i, gene_i in enumerate(genes_in_category):
        for j, gene_j in enumerate(genes_in_category):
            if i < j and correlation_matrix[i, j] > 0.7:  #CHANGE TRESHOLD
                G.add_edge(gene_i, gene_j, weight=correlation_matrix[i, j])

    plt.figure(figsize=(10, 10))
    nx.draw_networkx(G, with_labels=True, node_size=1000, node_color='skyblue')
    plt.title(f'Coexp network  - {category.capitalize()}')
    plt.savefig(f"{res_path}/networkot-{category.capitalize()}.png", dpi=300)
    plt.show()

'''
