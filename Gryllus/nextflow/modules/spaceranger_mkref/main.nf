process SPACERANGER_MKREF {

  tag "mkref:${params.genome_name}"
  publishDir "${params.outdir}/mkref", mode: 'copy'

  input:
    val dummy

  output:
    path("${params.genome_name}")

  script:
  """
  bash ${projectDir}/Gryllus/scripts/spaceranger_mkref.sh \
    --spaceranger ${params.spaceranger_bin} \
    --genome ${params.genome_name} \
    --fasta ${params.fasta} \
    --genes ${params.gtf} \
    --outdir .
  """
}