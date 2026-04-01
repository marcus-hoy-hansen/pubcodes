#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bam_multicov_windows.sh sample.bam [windowsize]
#
# Example:
#   bam_multicov_windows.sh 00246-21_sup.bam 1000000

BAM="$1"
WINDOW="${2:-1000000}"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi

# Activate conda environment with samtools + bedtools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bedtools_env

# Derive names
BAM_BASE=$(basename "$BAM" .bam)
SIZES="${BAM_BASE}.genome.sizes"
BED="${BAM_BASE}_${WINDOW}.bed"
OUT="${BAM_BASE}_${WINDOW}.cov.tsv"

echo "▶ BAM        : $BAM"
echo "▶ Window size: $WINDOW bp"
echo "▶ Output     : $OUT"

# 1) Extract chromosome sizes from BAM header
samtools view -H "$BAM" \
| awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2"\t"$3}' \
> "$SIZES"

# 2) Generate windows
bedtools makewindows -g "$SIZES" -w "$WINDOW" > "$BED"

# 3) Multicov
bedtools multicov -bams "$BAM" -bed "$BED" > "$OUT"

echo "✔ Done"
