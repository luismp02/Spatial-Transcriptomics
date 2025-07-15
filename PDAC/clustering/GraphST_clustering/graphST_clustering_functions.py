from GraphST.GraphST.utils import clustering
import anndata
import pandas as pd
import anndata

def data_preprocessing(count_data, spatial_data, typ='filtered'):
    """
    typ: '10X' 'SRTsim', 'filtered', 'scCube'
    seed: for SRTsim data (expect a number)
    col_label_true: ground true
    """

    if typ=='SRTsim' or typ=='filtered':
        adata = anndata.read_csv(count_data)
        spatial = pd.read_csv(spatial_data, index_col=0)
        #obs creation
        common_indices = adata.obs.index.intersection(spatial.index)
        filtered_spatial = spatial[spatial.index.isin(common_indices)]
        adata.obs = filtered_spatial
        col_label_true='label'
        #obsm creation
        coordinates = filtered_spatial.loc[:,['x','y']]
        adata.obsm['spatial'] = coordinates.astype(str).values
        
    elif typ=='scCube':
        adata = anndata.read_h5ad(count_data)
        adata.obs = adata.obs.iloc[:,1:]
        coord = adata.obs.iloc[:,[1,2]]
        adata.obsm['spatial'] = coord.astype(str).values
        col_label_true = 'Cell_type'
    
    else:
        raise ValueError(f"Type '{typ}' is not supported.")
        
    adata.obsm['spatial'] = adata.obsm['spatial'].astype('float64')
    return adata, col_label_true

def graphST_clustering(adata, method, n_clusters_list, start, end, refinement):
    '''
    refinement (False or True) = label re-assignment as the same domain as the most common label of its surronding spots (based on radius)
    if refinement = True, it generates a {method}_{n_clusters}_refined column in addition to the {method}_{n_clusters} column
    '''
    # radius specifies the number of neighbors considered during refinement
    if method == 'mclust':
      clustering(adata, n_clusters_list, radius=50, key='emb_pca', method=method, refinement=refinement) 
    elif method in ['leiden', 'louvain']:
      clustering(adata, n_clusters_list, radius=50, key='emb_pca', method=method, start=start, end=end, increment=0.01, refinement=refinement)
    adata.obs.drop(columns=[method], inplace=True) #remove the method column (redundant)
    return adata

def graphST_methods_loop(adata, methods, n_clusters_list, start, end, refinement):
    '''
    Execute the clustering and the visualization/metrics for the methods of interest 
    Parameters:
        - methods can be a string or a list of string 
        - refinement: True/False
    '''
    if refinement == 'true': #correct the nextflow parameter
        refinement = True
    elif refinement == 'false':
        refinement = False

    if isinstance(methods, str):
        methods = [methods]
        
    if isinstance(n_clusters_list, int):
        n_clusters_list = [n_clusters_list]

    for method in methods:
        copy_list = n_clusters_list.copy()
        # clustering
        adata = graphST_clustering(adata, method, copy_list, start, end, refinement)

    return adata