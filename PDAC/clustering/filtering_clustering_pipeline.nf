#!/usr/bin/env nextflow

log.info "Pipeline: Filtering - Annotation - Simulation - Clustering - Analysis"
log.info "======================================================================"
log.info "sample name : ${params.sample_name}"
log.info "count thresholds (min / max) : ${params.thr_min} / ${params.thr_max}"
log.info "Ucell thresholds (basal / classical) : ${params.basal_thr} / ${params.classical_thr}"
log.info "SRTsim seeds : ${params.seeds}"
log.info "ouput (path): ${params.output}"
log.info "\n"

def home_dir = System.getenv('HOME')
def path_to_project = "${home_dir}/Project"
def input_file = "${path_to_project}/datashare/PDAC/visium_PDAC"
log.info "Input file path (main.nf): ${input_file}"

include { spatial_preprocessing } from './nf_modules/spatial_preprocessing'
include { srtsim } from './nf_modules/srtsim'
include { graphst as graphst_ref } from './nf_modules/graphst'
include { graphst as graphst_sim } from './nf_modules/graphst'
include { metrics_and_visualization as mv_ref } from './nf_modules/metrics_and_visualization'
include { metrics_and_visualization as mv_sim } from './nf_modules/metrics_and_visualization'
include { final_figures} from './nf_modules/final_figures'

workflow {
    input = channel.fromPath(input_file, checkIfExists: true)
    spatial_preprocessing(input, params.thr_min, params.thr_max, params.basal_thr, params.classical_thr)

    if (params.seeds){
        if (params.seeds.contains(',')) {
            seed_ch = channel.from(params.seeds.split(',')*.trim()*.toInteger())
        }
        else{
            seed_ch = channel.of(params.seeds.trim.toInteger())
        }
        spatial_preprocessing.out.filtered_matrix
            .combine(spatial_preprocessing.out.annotated_coordinates)
            .combine(seed_ch)
            .map { filt, annot, seed -> tuple(seed, filt, annot) }
            | srtsim

        srtsim.out.simulated_data
            .combine(spatial_preprocessing.out.annotated_coordinates)
            .map { sim, annot -> tuple(sim, annot, params.sample_name, path_to_project) }
            | graphst_sim

        mv_sim(graphst_sim.out.clustered_data, path_to_project)
    }
    
    spatial_preprocessing.out.filtered_matrix
        .combine(spatial_preprocessing.out.annotated_coordinates)
        .map { filt, annot -> tuple(filt, annot, params.sample_name, path_to_project) }
        | graphst_ref

    mv_ref(graphst_ref.out.clustered_data, path_to_project)
    //graphst_ref.out.clustered_data
    //    | mv_ref

    path_list  = mv_ref.out.final_table
        .combine(mv_sim.out.final_table.collect())
    
    path_list.set { combinedChannel }
    
    final_figures(path_list, params.sample_name)
}
