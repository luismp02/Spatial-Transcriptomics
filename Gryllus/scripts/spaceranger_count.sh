#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  spaceranger_count.sh --spaceranger <path/to/spaceranger> \
    --id <RUN_ID> \
    --transcriptome <mkref_output_dir/GENOME_NAME> \
    --fastqs <fastq_dir> \
    --sample <sample_name_in_fastqs> \
    --image <image.tif> \
    --slide <slide_id> \
    --area <A1|B1|C1|D1> \
    [--create-bam true|false] \
    [--unknown-slide visium-1] \
    --outdir <output_dir>

Example:
  spaceranger_count.sh --spaceranger ~/projects/spaceranger-3.1.2/spaceranger \
    --id SampleA_IAB --transcriptome /path/mkref_out/Gryllus_ref \
    --fastqs /path/fastqs/A --sample GMRS238 \
    --image /path/image.tif --slide V13B16-384 --area A1 \
    --create-bam false --outdir ./count_out
EOF
}

SPACERANGER=""
ID=""
TRANSCRIPTOME=""
FASTQS=""
SAMPLE=""
IMAGE=""
SLIDE=""
AREA=""
CREATE_BAM="false"
UNKNOWN_SLIDE=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --spaceranger) SPACERANGER="$2"; shift 2 ;;
    --id) ID="$2"; shift 2 ;;
    --transcriptome) TRANSCRIPTOME="$2"; shift 2 ;;
    --fastqs) FASTQS="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    --slide) SLIDE="$2"; shift 2 ;;
    --area) AREA="$2"; shift 2 ;;
    --create-bam) CREATE_BAM="$2"; shift 2 ;;
    --unknown-slide) UNKNOWN_SLIDE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$SPACERANGER" || -z "$ID" || -z "$TRANSCRIPTOME" || -z "$FASTQS" || -z "$SAMPLE" || -z "$IMAGE" || -z "$SLIDE" || -z "$AREA" || -z "$OUTDIR" ]]; then
  echo "Missing required arguments."
  usage
  exit 1
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

cmd=( "$SPACERANGER" count
  --id="$ID"
  --transcriptome="$TRANSCRIPTOME"
  --fastqs="$FASTQS"
  --sample="$SAMPLE"
  --image="$IMAGE"
  --slide="$SLIDE"
  --area="$AREA"
  --create-bam="$CREATE_BAM"
)

# opcional
if [[ -n "$UNKNOWN_SLIDE" ]]; then
  cmd+=( --unknown-slide="$UNKNOWN_SLIDE" )
fi

printf 'Running: %q ' "${cmd[@]}"; echo
"${cmd[@]}"
