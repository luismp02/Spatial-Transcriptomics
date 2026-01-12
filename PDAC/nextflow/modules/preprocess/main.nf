process PREPROCESS_INTEGRATED {

  tag "$sample"
  publishDir "${params.outdir}/preprocess/${sample}", mode: 'copy'

  input:
    tuple val(sample), path(spatial_outs)

  output:
    tuple val(sample),
      path("sc_ref_processed_fine.rds"),
      path("${sample}_processed.rds")

  script:
  """
  Rscript ${projectDir}/Preprocessing/scripts/preprocess_integrated.R \
    ${params.ref_rds} \
    ${spatial_outs} \
    ${sample}
  """
}

