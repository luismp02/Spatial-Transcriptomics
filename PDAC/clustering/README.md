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

Nextflow process can be added in the nf_modules directory and implemented in the filtering_clustering_pipeline.nf to implement a new clustering method.  
Concerning metrics, python methods can be added in clustering_assesment/metrics_and_visualizations_functions.py (see the implemented ones, namely lisi_metrics, silhouette metrics, rand_index and adjusted_rand_index functions). The functions from metrics_and_visualizations_functions.py are executed in metrics_and_visualizaations.ipynb. If specific input are required to calculate the new metric, it must be specified in the nextflow pipeline and the metrics_and_visualization.nf module. 

### SCTransform 

The SCTransform function from Seurat (5.2.0) was used, as part of the preprocessing, to normalize and scale the 10X spatial transcriptomics expression data. 

It models the expression of each gene through a negative binomial distribution (i), using the counts of gene i from all spots and considering the sequencing depth as a covariate (ii).

(i) Y ~ NB(μ, θ)

(ii) log(expected count_j) = β₀ + β₁ × log₁₀(sequencing depth)

The parameters μ (mean) and θ (dispersion) were estimated independently for each gene, offering greater flexibility. Thus, gene expected counts for each spot were determined through the distribution, accounting for the sequencing depth of the spots.
Finally, Pearson residuals (iii), that correspond to the final normalized counts, were calculated, measuring the difference between observed and estimated counts.
  



