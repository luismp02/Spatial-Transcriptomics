process LOAD_GRyllus_R {

  tag "${sample_id}"
  publishDir "${params.outdir}/load_r", mode: 'copy'

  input:
    tuple val(sample_id), path(outs_dir)

  output:
    path("${sample_id}.rds")

  script:
  """
  Rscript ${projectDir}/Gryllus/scripts/load_data_nf.R \
    --outs ${outs_dir} \
    --files ${params.files_dir} \
    --out_rds ${sample_id}.rds
  """
}
