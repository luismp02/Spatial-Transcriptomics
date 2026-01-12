#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  spaceranger_mkref.sh --spaceranger <path/to/spaceranger> \
    --genome <GENOME_NAME> \
    --fasta <genome.fasta> \
    --genes <genes.gtf> \
    --outdir <output_dir>

Example:
  spaceranger_mkref.sh --spaceranger ~/projects/spaceranger-3.1.2/spaceranger \
    --genome Gryllus_ref --fasta Gbi_genome.fa --genes Gbi_Genes_clean.gtf --outdir ./mkref_out
EOF
}

SPACERANGER=""
GENOME=""
FASTA=""
GENES=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --spaceranger) SPACERANGER="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --fasta) FASTA="$2"; shift 2 ;;
    --genes) GENES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$SPACERANGER" || -z "$GENOME" || -z "$FASTA" || -z "$GENES" || -z "$OUTDIR" ]]; then
  echo "Missing required arguments."
  usage
  exit 1
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

"$SPACERANGER" mkref \
  --genome="$GENOME" \
  --fasta="$FASTA" \
  --genes="$GENES"