#!/usr/bin/env nextflow

log.info "Pipeline: Filtering - Annotation - Simulation - Clustering - Analysis"
log.info "======================================================================"
log.info "path to data directory : ${params.input_path}"
log.info "sample name : ${params.sample_name}"
log.info "count thresholds (min / max) : ${params.thr_min} / ${params.thr_max}"
log.info "Ucell thresholds (basal / classical) : ${params.basal_thr} / ${params.classical_thr}"
log.info "Clustering methods (mclust, louvain, leiden) : ${params.method}"
log.info "Number of clusters (min / max) : ${params.min_clust} - ${params.max_clust}"
log.info "Refinement (True/False) : ${params.refine}"
log.info "SRTsim seeds : ${params.seeds}"
log.info "R installation path : ${params.r_path}"
log.info "ouput (path): ${params.output}"
log.info "\n"

def path_to_dir = file("..").toAbsolutePath().toString()
def path_to_graphST_functions = "${path_to_dir}/clustering/GraphST_clustering"
def path_to_visualisation_functions = "${path_to_dir}/clustering/clustering_assessment"

include { spatial_preprocessing } from './nf_modules/spatial_preprocessing'
include { srtsim } from './nf_modules/srtsim'
include { graphst as graphst_ref } from './nf_modules/graphst'
include { graphst as graphst_sim } from './nf_modules/graphst'
include { metrics_and_visualization as mv_ref } from './nf_modules/metrics_and_visualization'
include { metrics_and_visualization as mv_sim } from './nf_modules/metrics_and_visualization'
include { final_figures} from './nf_modules/final_figures'

workflow {
    input = channel.fromPath(params.input_path, checkIfExists: true)
    spatial_preprocessing(path_to_dir, input, params.thr_min, params.thr_max, params.basal_thr, params.classical_thr)

    if (params.seeds){
        if (params.seeds.toString().contains(',')) {
            seed_ch = channel.from(params.seeds.split(',')*.trim()*.toInteger())
        }
        else{
            seed_ch = channel.of(params.seeds.toInteger())
        }
        spatial_preprocessing.out.filtered_matrix
            .combine(spatial_preprocessing.out.annotated_coordinates)
            .combine(seed_ch)
            .map { filt, annot, seed -> tuple(seed, filt, annot, path_to_dir) }
            | srtsim

        srtsim.out.simulated_data
            .combine(spatial_preprocessing.out.annotated_coordinates)
            .map { sim, annot -> tuple(sim, annot, params.sample_name, path_to_dir, path_to_graphST_functions, params.r_path, params.min_clust, params.max_clust, params.method, params.refine) }
            | graphst_sim

        mv_sim(graphst_sim.out.clustered_data, path_to_dir, path_to_visualisation_functions)
    }
    
    spatial_preprocessing.out.filtered_matrix
        .combine(spatial_preprocessing.out.annotated_coordinates)
        .map { filt, annot -> tuple(filt, annot, params.sample_name, path_to_dir, path_to_graphST_functions, params.r_path, params.min_clust, params.max_clust, params.method, params.refine) }
        | graphst_ref

    mv_ref(graphst_ref.out.clustered_data, path_to_dir, path_to_visualisation_functions)

    if (params.seeds){
        path_list  = mv_ref.out.final_table
            .combine(mv_sim.out.final_table.collect())
        
        path_list.set { combinedChannel }
        
        final_figures(path_list, params.sample_name, params.refine)
    }
}
