process SPACERANGER_COUNT {

  tag "${sample_id}"
  publishDir "${params.outdir}/count", mode: 'copy'

  input:
    tuple val(sample_id),
          path(fastq_dir),
          val(sample_name),
          path(image_path),
          val(slide),
          val(area),
          path(transcriptome_dir)

  output:
    path("${sample_id}")

  script:
  """
  bash ${projectDir}/Gryllus/scripts/spaceranger_count.sh \
    --spaceranger ${params.spaceranger_bin} \
    --id ${sample_id} \
    --transcriptome ${transcriptome_dir} \
    --fastqs ${fastq_dir} \
    --sample ${sample_name} \
    --image ${image_path} \
    --slide ${slide} \
    --area ${area} \
    --create-bam ${params.create_bam} \
    ${ params.unknown_slide ? "--unknown-slide ${params.unknown_slide}" : "" } \
    --outdir .
  """
}