process final_figures{
    tag "final figures"
    publishDir "${params.output}/results", mode: 'copy'

    input:
    path(metrics_files)
    val(sample_name)
    val(refine)

    output:
    path "*final_graphs.pdf"

    script:
    """
    source ~/Project/Spatial-Transcriptomics/PDAC/clustering/envrt_clustering/bin/activate
    path_list=\$(echo ${metrics_files} | sed "s/ /','/g" | sed "s/^/['/" | sed "s/\$/']/")
    echo "metrics_files: \$path_list" 
    papermill ~/Project/Spatial-Transcriptomics/PDAC/clustering/clustering_assessment/final_graph.ipynb \
        -p sample_name ${sample_name} \
        -p table_list \$path_list \
        -p refine ${refine}
    deactivate
    """
}