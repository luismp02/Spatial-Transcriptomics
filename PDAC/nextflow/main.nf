nextflow.enable.dsl=2

include { PREPROCESS_INTEGRATED } from './modules/preprocess/main'
include { SPOTLIGHT_RUN }         from './modules/spotlight/run'
include { SPADE_RUN }             from './modules/spade/run'
include { SIMULATION_SCCUBE }     from './modules/simulation/run'

workflow {

  Channel
    .fromPath("${projectDir}/nextflow/assets/slides.txt")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .map { slide ->
      tuple(slide, file("${params.visium_base}/${slide}/outs"))
    }
    .set { CH_SLIDES }

  pre = PREPROCESS_INTEGRATED(CH_SLIDES)

  pre
    .map { sample, sc_ref, spatial_rds -> sample }
    .set { CH_SAMPLES }

  if (params.mode == 'spotlight' || params.mode == 'all') {
    SPOTLIGHT_RUN(CH_SAMPLES)
  }

  if (params.mode == 'spade' || params.mode == 'all') {
    SPADE_RUN(CH_SAMPLES)


  if (params.mode == 'simulation' || params.mode == 'all') {
    Channel.of(1).set { CH_SIM }
    SIMULATION_SCCUBE(CH_SIM)
  }
}
}
