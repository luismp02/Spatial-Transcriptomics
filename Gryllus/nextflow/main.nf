nextflow.enable.dsl=2

include { SPACERANGER_MKREF } from './modules/spaceranger_mkref/main'
include { SPACERANGER_COUNT } from './modules/spaceranger_count/main'
include { LOAD_GRyllus_PY }   from './modules/load_py/main'
include { LOAD_GRyllus_R }    from './modules/load_r/main'

workflow {

  /*
   * transcriptome source:
   * - if params.transcriptome_dir is set
   * - if not, if mode includes mkref => we generate it with mkref
   */
  Channel transcriptome_ch

  if (params.transcriptome_dir) {
    transcriptome_ch = Channel.value(file(params.transcriptome_dir))
  } else if (params.mode in ['mkref','count','all']) {
    Channel.of(1).set { CH_MK }
    mk = SPACERANGER_MKREF(CH_MK)
    transcriptome_ch = mk.out
  } else {
    transcriptome_ch = Channel.empty()
  }

  /*
   * Samplesheet  for count
   */
  Channel
    .fromPath("${projectDir}/Gryllus/nextflow/assets/samplesheet.tsv")
    .splitCsv(header:true, sep:'\t')
    .map { row ->
      tuple(
        row.sample_id as String,
        file(row.fastq_dir),
        row.sample_name as String,
        file(row.image_path),
        row.slide as String,
        row.area as String
      )
    }
    .set { CH_SAMPLES }

  /*
   * COUNT
   */
  if (params.mode in ['count','all']) {

    CH_SAMPLES
      .combine(transcriptome_ch)
      .map { t, transcriptome_dir ->
        tuple(t[0], t[1], t[2], t[3], t[4], t[5], transcriptome_dir)
      }
      .set { CH_COUNT_IN }

    cnt = SPACERANGER_COUNT(CH_COUNT_IN)

    // outs channel from count
    cnt.out
      .map { sample_id -> tuple(sample_id, file("${sample_id}/outs")) }
      .set { CH_OUTS_FROM_COUNT }

    if (params.mode == 'count') {
    } else {
      if (params.mode in ['all']) {
        if (params.run_load_py) LOAD_GRyllus_PY(CH_OUTS_FROM_COUNT)
        if (params.run_load_r)  LOAD_GRyllus_R(CH_OUTS_FROM_COUNT)
      }
    }
  }

  /*
   * LOAD ONLY (outs in base)
   * Existing ${params.outs_base}/${sample_id}/outs
   */
  if (params.mode in ['load_py','load_r','load_only']) {

    if (!params.outs_base) {
      error "For mode=${params.mode}, you must set --outs_base <dir containing sample_id/outs>"
    }

    CH_SAMPLES
      .map { t ->
        def sample_id = t[0]
        tuple(sample_id, file("${params.outs_base}/${sample_id}/outs"))
      }
      .set { CH_OUTS_ONLY }

    if (params.mode in ['load_py','load_only']) {
      LOAD_GRyllus_PY(CH_OUTS_ONLY)
    }
    if (params.mode in ['load_r','load_only']) {
      LOAD_GRyllus_R(CH_OUTS_ONLY)
    }
  }

  /*
   * MKREF only
   */
  if (params.mode == 'mkref') {
    // mkref run if mode==mkref
  }
}
