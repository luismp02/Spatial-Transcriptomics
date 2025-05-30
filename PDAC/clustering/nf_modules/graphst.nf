process graphst {
    tag "GraphST - ${filtered_count_data}"
    publishDir "${params.output}/data", mode: 'copy'

    input:
    tuple path(filtered_count_data), path(spatial_data), val(sample_name), path(path_to_project)

    output:
    path "*_clustered.h5ad", emit : clustered_data

    script:
    """
    source ~/Project/Spatial-Transcriptomics/PDAC/clustering/envrt_clustering/bin/activate
    papermill ~/Project/Spatial-Transcriptomics/PDAC/clustering/GraphST_clustering/clustering.ipynb \
        -p filtered_count_data ${filtered_count_data} \
        -p spatial_data ${spatial_data} \
        -p sample_name ${sample_name} \
        -p path_to_project ${path_to_project}
    deactivate
    """
}