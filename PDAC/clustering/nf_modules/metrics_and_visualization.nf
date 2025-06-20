process metrics_and_visualization {
    tag "metrics and visualization - ${clustered_data}"
    publishDir "${params.output}/results", mode: 'copy'

    input:
    path clustered_data
    val path_to_dir
    path path_to_visualization_functions

    output:
    path "*.pdf", emit : figures
    path "*.csv", emit : final_table

    script:
    """
    source ${path_to_dir}/clustering/envrt_clustering/bin/activate
    papermill ${path_to_dir}/clustering/clustering_assessment/metrics_and_visualization.ipynb \
        -p clustered_data ${clustered_data} \
    deactivate
    """
}