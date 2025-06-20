process srtsim {
    tag "Simulation seed: ${seed}"
    publishDir "${params.output}/data", mode: 'copy'


    input:
    tuple val(seed), path(filtered_matrix), path(annotated_coordinates), val(path_to_dir)

    output:
    path "${params.sample_name}_simulated_${seed}.csv", emit : simulated_data
    //val seed, emit : seed

    script:
    """
    Rscript -e "rmarkdown::render('${path_to_dir}/simulations/SRTsim/scripts/SRTsim_same_domains.Rmd',
        output_file = '${params.sample_name}_${seed}_SRTsim_same_domain_report.html',
        output_dir = '${path_to_dir}/simulations/SRTsim/report',
        params = list(
            seed = ${seed},
            filtered_matrix = '$filtered_matrix',
            annotated_coordinates = '$annotated_coordinates',
            sample_name = '${params.sample_name}',
            working_dir = '\$(pwd)'))"
    """
}