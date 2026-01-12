process LOAD_GRyllus_PY {

  tag "${sample_id}"
  publishDir "${params.outdir}/load_py", mode: 'copy'

  input:
    tuple val(sample_id), path(outs_dir)

  output:
    path("${sample_id}.h5ad")

  script:
  """
  python ${projectDir}/Gryllus/scripts/load_data_nf.py \
    --outs ${outs_dir} \
    --files ${params.files_dir} \
    --out_h5ad ${sample_id}.h5ad
  """
}
