# GraphST

The GraphST tool (1) was implemeted using the https://github.com/JinmiaoChenLab/GraphST project. Nevertheless, in order to optimize the testing of Louvain and Leiden methods with different number of clusters, some modifications were done in GraphST/GraphST/utils.py. These changes avoid recalculating resolutions multiples times but determine the resolutions corresponding to a list of cluster numbers of interest. 
Different functions to visualize clustering (silhouette scores, UMAP) and automate the clustering with different methods can be found in functions_graphST.py.

The dependencies used are listed in the requirements.txt and can be installed using the following command:

```bash
pip install -r requirements.txt
```

## Preprocessing

GraphST generates a normalized gene expression matrix containing the top 3,000 highly variables gene expression in the different spots of interests using Scanpy (log tranformation + normalization).
Note: the Scanpy normalization divides the gene counts of each cell by the total count. Thus, it assumes that all regions have the same mRNA abundance. 


## References

1 - GraphST: Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST, 
Nature communication, Yahui Long et al. 2023
