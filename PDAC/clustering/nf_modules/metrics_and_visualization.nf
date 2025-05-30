process metrics_and_visualization {
    tag "metrics and visualization - ${clustered_data}"
    publishDir "${params.output}/results", mode: 'copy'

    input:
    path clustered_data
    path path_to_project

    output:
    path "*.pdf", emit : figures
    path "*.csv", emit : final_table

    script:
    """
    source ~/Project/Spatial-Transcriptomics/PDAC/clustering/envrt_clustering/bin/activate
    papermill ~/Project/Spatial-Transcriptomics/PDAC/clustering/clustering_assessment/metrics_and_visualization.ipynb \
        -p clustered_data ${clustered_data} \
        -p path_to_project ${path_to_project}
    deactivate
    """
}