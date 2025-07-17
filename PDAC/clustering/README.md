## Clustering Nextflow pipeline

![pipeline figure](figures/clustering_nf_pipeline.png)

### Options

```
--input_path (path to data directory)
--sample_name (expect a name that corresponds to an input_path subdirectory)
--thr_min  (minimum count cutoff to filter spots)
--thr_max  (maximum count cutoff to filter spots)
--basal_thr  (basal Uscore threshold to define 'basal' spots)
--classical_thr  (classical Uscore threshold to define 'classical' spots)
--method (GraphST clustering methods: mclust,louvain,leiden)
--refine (GraphST refinement step: True/False)
--min_clust (minimal number of clusters)
--max_clust (maximal number of clusters, such as min_clust to max_clust clusters will be generated)
--output (output directory)
--seeds (seeds used to generate simulated spatial transcriptomics  data)
--r_path (path to the directory where R is locally installed)
```

### Example

```
nextflow run filtering_clustering_pipeline.nf \
--input_path ../../../datashare/PDAC/visium_PDAC \
--sample_name Visium_FFPE_V43T08-051_D \
--thr_min 2500 \
--thr_max Inf \
--basal_thr 0.25 \
--classical_thr 0.25 \
--method mclust,louvain,leiden \
--refine True \
--min_clust 2 \
--max_clust 10 \
--output nf_output \
--seeds 1,2 \
--r_path /Library/Frameworks/R.framework/Resources
```
### Requirements

- anndata==0.8.0
- harmonypy==0.0.10
- igraph==0.11.8
- leidenalg==0.10.2
- louvain==0.8.2
- matplotlib==3.4.2
- numpy==1.22.3
- POT==0.9.5
- pandas==1.4.2
- python
- R
- rpy2==3.4.1
- scanpy=1.9.1
- scipy==1.8.1
- scikit-learn==1.1.1
- seaborn==0.13.2
- torch
- tqdm==4.64.0

### Modifications

#### Adding a new clustering method

To implement a new clustering method:
- Add a new **Nextflow process** in the `nf_modules` directory.
- Integrate the process in the `filtering_clustering_pipeline.nf` file.

#### Adding a new metric

To include a new metric for clustering assessment (in python):
- Add the method in `clustering_assessment/metrics_and_visualizations_functions.py` (see the implemented metrics, namely lisi_metrics, silhouette metrics, rand_index and adjusted_rand_index functions)
- The functions fron this file are exectued in `metrics_and_visualizations.ipynb`.
- If your new metric requires specific inputs specify these in the Nextflow pipeline and update the `metrics_and_visualization.nf` module accordingly.

### Methods

The different methods used are described in **M2-MEMOIRE-Marchand_Mehdi_2024-2025**.
> Note: GraphST proposes an optional method to reassign the spots to the same cluster as those of the surrounding spots within a defined radius (=refinement).

### Absence of clustering with Louvain and Leiden

Louvain and Leiden graph-based methods identify clusters of varying sizes depending
on the resolution parameter. A high resolution results in the identification of several
small clusters, while a low resolution leads to the identification of fewer but larger
clusters.
In the context of GraphST, the model computes the clustering by decreasing the
resolution from an initial value defined by the user until it identifies the number of
clusters specified.

Thus, if no clustering result is obtained try to increase the end parameter (i.e., the starting resolution) of the graphST_methods_loop function (`clustering.ipynb`).
