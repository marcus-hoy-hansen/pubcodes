#!/bin/bash
# # SBATCH --account <ACCOUNT_NAME> 
#SBATCH -c 1
#SBATCH --mem 6g
#SBATCH --time 24:00:00

# -------------------------------------------------------------------------
# Author: Marcus H. Hansen (marcus-hoy-hansen)
# -------------------------------------------------------------------------
# md5sum_catalogue.sh
# Generate a recursive, sorted MD5 checksum catalogue for a directory, 
# when TAR is not the obvious choice.
# Intended for cross-machine verification: run once per location.
# Check with: comm -3 <(sort md5sum_[file2].txt) <(sort md5sum_[file1].txt)
# -------------------------------------------------------------------------
# Usage:
#   ./md5sum_catalogue.sh /path/to/folder
#
# Output:
#   md5sum_<foldername>.txt
# -------------------------------------------------------------------------


set -euo pipefail

dir="$1"

if [[ -z "$dir" ]]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

if [[ ! -d "$dir" ]]; then
  echo "Error: '$dir' is not a directory."
  exit 1
fi

# Normalize directory path and output filename
absdir="$(realpath "$dir")"
basename="$(basename "$absdir")"
outfile="md5sum_${basename}.txt"

echo "[INFO] Generating recursive MD5 list for: $absdir"
echo "[INFO] Output file: $outfile"
echo "[INFO] This may take a while for large folders..."

# Use find+xargs (null-safe, sorted) for reproducibility
(
  cd "$absdir"
  find . -type f -print0 | LC_ALL=C sort -z | xargs -0 md5sum
) > "$outfile"

echo "[INFO] Done. Saved: $(realpath "$outfile")"
echo "[INFO] Lines: $(wc -l < "$outfile")"


