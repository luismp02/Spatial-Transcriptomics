process SIMULATION_SCCUBE {

  tag "scCube"

  publishDir "${params.outdir}/simulation", mode: 'copy'

  input:
    val dummy

  output:
    path("*")

  script:
  """
  echo ">>> Running scCube simulation"

  python ${projectDir}/PDAC/simulation/run_scCube.py
  """
}