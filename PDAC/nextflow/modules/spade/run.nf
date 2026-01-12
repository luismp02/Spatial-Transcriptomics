process SPADE_RUN {

  tag "$sample"
  publishDir "${params.outdir}/spade", mode: 'copy'

  input:
    val sample

  output:
    path("${sample}/CTperLayer_ready.rds")
    path("${sample}/loc_ready.rds")
    path("${sample}/CTest_ready.rds")

  script:
  """
  Rscript ${projectDir}/deconvolution/SPADE/scripts/SPADE.R \
    --slides ${sample} \
    --visium_base ${params.visium_base} \
    --spagcn_out ${params.spagcn_out} \
    --ref_rds ${params.spade_ref_rds} \
    --markers_rds ${params.spade_markers_fused} \
    --out_base .
  """
}