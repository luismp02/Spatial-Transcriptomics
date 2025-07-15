process graphst {
    tag "GraphST - ${filtered_count_data}"
    publishDir "${params.output}/data", mode: 'copy'

    input:
    tuple path(filtered_count_data), path(spatial_data), val(sample_name), val(path_to_dir), path(path_to_graphST_functions), val(r_path), val(min_clust), val(max_clust), val(method), val(refine)

    output:
    path "*_clustered.h5ad", emit : clustered_data

    script:
    """
    source ${path_to_dir}/clustering/envrt_clustering/bin/activate
    papermill ${path_to_dir}/clustering/GraphST_clustering/clustering.ipynb \
        -p filtered_count_data ${filtered_count_data} \
        -p spatial_data ${spatial_data} \
        -p sample_name ${sample_name} \
        -p r_path ${r_path} \
        -p min_clust ${min_clust} \
        -p max_clust ${max_clust} \
        -p method ${method} \
        -p refine ${refine}
    deactivate
    """
}