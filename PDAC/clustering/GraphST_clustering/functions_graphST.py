from GraphST.GraphST.utils import clustering
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, silhouette_samples
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc

def graphST_clustering(adata, method, n_clusters_list, refinement):
    '''
    refinement (False or True) = label re-assignment as the same domain as the most common label of its surronding spots (based on radius)
    '''
    # radius specifies the number of neighbors considered during refinement
    if method == 'mclust':
      clustering(adata, n_clusters_list, radius=50, key='emb_pca', method=method, refinement=refinement) 
    elif method in ['leiden', 'louvain']:
      clustering(adata, n_clusters_list, radius=50, key='emb_pca', method=method, start=0.1, end=1.5, increment=0.01, refinement=refinement)    
    return adata


def clustering_plot(ax, adata, n_clusters):
    sc.pl.spatial(adata, img_key="hires", color=f"domain_{n_clusters}", show=False, ax=ax)


def umap(ax, adata, n_clusters):
    '''plotting predicted labels by UMAP'''
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=f'domain_{n_clusters}', title=['Predicted labels'], show=False, ax=ax)


def silhouette_metrics(ax, adata, n_clusters):
   '''
   metric to evaluate the clustering 
   '''

   ##### SILHOUETTE SCORE ######
   X = adata.obsm['emb_pca']
   cluster_labels = adata.obs[f'domain_{n_clusters}'].cat.codes.to_numpy()
   silhouette_avg = silhouette_score(X, cluster_labels)
   sample_silhouette_values = silhouette_samples(X, cluster_labels)
   
   ##### SILHOUETTE PLOT ######
   y_lower = 10
   for i in range(n_clusters):
      # Aggregate the silhouette scores for samples belonging to
      # cluster i, and sort them
      ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
      ith_cluster_silhouette_values.sort()

      size_cluster_i = ith_cluster_silhouette_values.shape[0]
      y_upper = y_lower + size_cluster_i

      color = cm.tab20(2 * i) # * permet de selectionner les bonnes couleurs (1 sur 2)
      #sc.pl.palettes.default_20[:n_clusters]
      ax.fill_betweenx(
         np.arange(y_lower, y_upper),
         0,
         ith_cluster_silhouette_values,
         facecolor=color,
         edgecolor=color,
         alpha=0.7,
      )

      # Label the silhouette plots with their cluster numbers at the middle
      ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i+1))

      # Compute the new y_lower for next plot
      y_lower = y_upper + 10  # 10 for the 0 samples

   ax.set_title("The silhouette plot for the various clusters.")
   ax.set_xlabel("The silhouette coefficient values")
   ax.set_ylabel("Cluster label")

   # The vertical line for average silhouette score of all the values
   ax.axvline(x=silhouette_avg, color="red", linestyle="--")

   ax.set_yticks([])  # Clear the yaxis labels / ticks
   ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

   ax.text(-0.77, 0.05, f"Average silhouette_score: {silhouette_avg:.2f}", fontsize=12, color='red')
   

def execute_all(adata, method, n_clusters_list, refinement):
    '''
    clustering and visualisation for a list of cluster numbers
    '''
    copy_list = n_clusters_list.copy()
    figures = []
    adata = graphST_clustering(adata, method, copy_list, refinement)
    for n_clusters in n_clusters_list:
        if f'domain_{n_clusters}' in adata.obs:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            clustering_plot(axes[0], adata, n_clusters)
            umap(axes[1], adata, n_clusters)
            silhouette_metrics(axes[2], adata, n_clusters)
            plt.suptitle(f"{method} : {n_clusters} clusters", fontsize=14, fontweight="bold")
            # plt.tight_layout()
            plt.show()
            figures.append(fig)
    return figures

def n_cluster_loop(adata, methods, n_clusters_list, refinement=True, filename=None):
    '''
    Execute the execute_all function for different methods at once
    Parameters:
        - methods can be a string or a list of string 
        - refinement: True/False
    '''
    if isinstance(methods, str):
        methods = [methods]

    if filename:
        with PdfPages(filename) as pdf:
            for method in methods:
                #removing of domain columns 
                domain_columns = [col for col in adata.obs.columns if col.startswith('domain_')]
                if domain_columns: #check if there is a column that begins with 'domain'
                    adata.obs.drop(columns=domain_columns, inplace=True, errors='ignore')
                
                figures = execute_all(adata, method, n_clusters_list, refinement)
                for fig in figures:
                    pdf.savefig(fig)
                    plt.close(fig)
    
    else:
        for method in methods:
            #removing of domain columns 
            domain_columns = [col for col in adata.obs.columns if col.startswith('domain_')]
            adata.obs.drop(columns=domain_columns, inplace=True)
            execute_all(adata, method, n_clusters_list)
            