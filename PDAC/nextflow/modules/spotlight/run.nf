process SPOTLIGHT_RUN {

  tag "$sample"
  publishDir "${params.outdir}/spotlight", mode: 'copy'

  input:
    tuple val(sample), path(sc_ref), path(spatial_rds)

  output:
    path("${sample}/spotlight_result.rds")
    path("${sample}/spotlight_model_nmf.rds")
    path("${sample}/seurat_obj.rds")
    path("${sample}/celltypes_ordered.rds")

  script:
  """
  Rscript ${projectDir}/deconvolution/SPOTlight/scripts/SPOTlight.R \
    --ref_rds ${params.ref_rds} \
    --visium_outs ${params.visium_base}/${sample}/outs \
    --sample_id ${sample} \
    --out_base .
  """
}