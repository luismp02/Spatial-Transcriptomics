process spatial_preprocessing {
    tag "${params.sample_name}"
    publishDir "${params.output}/data", mode: 'copy'

    input:
    val path_to_dir
    val input_file
    val thr_min
    val thr_max
    val basal_thr
    val classical_thr

    output:
    path "${params.sample_name}_filtered.csv", emit: filtered_matrix
    path "${params.sample_name}_annotation.csv", emit: annotated_coordinates
    path "${params.sample_name}_table_sctype.csv", emit: contingency_table
    path "${params.sample_name}_preprocessing_figures.pdf", emit: preprocessing_figures

    script:
    """
    Rscript -e "rmarkdown::render('${path_to_dir}/analysis/spatial_preprocessing.Rmd',
        output_file = '${params.sample_name}_spatial_preprocessing_report.html',
        output_dir = '${path_to_dir}/analysis/report',
        params = list(
            input_file = '${input_file}',
            thr_min = ${thr_min},
            thr_max = ${thr_max},
            basal_thr = ${basal_thr},
            classical_thr = ${classical_thr},
            sample_name = '${params.sample_name}',
            working_dir = '\$(pwd)'
            )
        )"
    """
}