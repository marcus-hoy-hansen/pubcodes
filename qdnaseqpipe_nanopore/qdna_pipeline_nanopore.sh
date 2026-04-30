#!/usr/bin/env bash
#SBATCH --job-name=qdnapipe-nanopore
#SBATCH --account=nanopore_kga
#SBATCH -c 6
#SBATCH --mem 32g
#SBATCH --time 4:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

work_dir=${SLURM_SUBMIT_DIR:-$PWD}
script_dir=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
cd "$work_dir"

resolve_path() {
  local path=$1
  local base_dir=$2
  if [[ "$path" = /* ]]; then
    printf '%s\n' "$path"
  else
    printf '%s\n' "$base_dir/$path"
  fi
}

echo "PWD: $(pwd)"
echo "SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR:-}"
echo "SCRIPT_DIR: $script_dir"
echo "WORK_DIR: $work_dir"



usage() {
  cat <<EOF
Usage: qdna_pipeline.sh <input.cram> [sample_id]

Pipeline:
  1) accept BAM directly, or convert CRAM to BAM with samtools view -T
  2) sort BAM
  3) index BAM
  4) run QDNAseq CNV calling, plot, and VCF export

Defaults:
  reference: /faststorage/project/MomaDiagnosticSamples-KGA/BACKUP/reference/igv_genome_hg38/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.GRC_exclusions_masked.fasta
  threads:   6
  binsize:   100

Environment overrides:
  REFERENCE_FASTA   Reference FASTA for CRAM decoding
  THREADS           Samtools sort threads
  BINSIZE           QDNAseq bin size in kb
  QDNA_RSCRIPT      R script to run (default: qdna_stable5_cli.R)
  CONTINUE_ABORTED  Default 1; set to 0 for a clean QDNAseq rerun
  FORCE_BAM         Set to 1 to rebuild sorted BAM even if present
  FORCE_INDEX       Set to 1 to rebuild BAM index even if present
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 ]]; then
  usage >&2
  exit 1
fi

input_file=$(resolve_path "$1" "$work_dir")
input_basename=$(basename "$input_file")
sample_id=${2:-${input_basename%.cram}}
sample_id=${sample_id%.bam}

reference_fasta=${REFERENCE_FASTA:-/faststorage/project/MomaDiagnosticSamples-KGA/BACKUP/reference/igv_genome_hg38/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.GRC_exclusions_masked.fasta}
threads=${THREADS:-20}
binsize=${BINSIZE:-100}
qdna_rscript=$(resolve_path "${QDNA_RSCRIPT:-$script_dir/qdna_stable5_cli.R}" "$work_dir")
plot_script=$(resolve_path "${PLOT_CHR_SCRIPT:-$script_dir/plot_chr_boxwhisker_utest.py}" "$work_dir")
report_script=$(resolve_path "${REPORT_SCRIPT:-$script_dir/generate_sample_report.py}" "$work_dir")
control_sample=${CONTROL_SAMPLE:-02268-26_hg38_ASv2}
control_bed=$(resolve_path "${CONTROL_BED:-${control_sample}_log2.bed}" "$work_dir")

sorted_bam="${sample_id}.sorted.bam"
sorted_bai="${sorted_bam}.bai"
continue_aborted=${CONTINUE_ABORTED:-1}
force_bam=${FORCE_BAM:-0}
force_index=${FORCE_INDEX:-0}

if [[ ! -f "$input_file" ]]; then
  echo "Input file not found: $input_file" >&2
  exit 1
fi

if [[ ! -f "$reference_fasta" ]]; then
  echo "Reference FASTA not found: $reference_fasta" >&2
  exit 1
fi

if [[ ! -f "$qdna_rscript" ]]; then
  echo "QDNAseq R script not found: $qdna_rscript" >&2
  exit 1
fi

if [[ ! -f "$plot_script" ]]; then
  echo "Plot script not found: $plot_script" >&2
  exit 1
fi

if [[ ! -f "$report_script" ]]; then
  echo "Report script not found: $report_script" >&2
  exit 1
fi

input_ext=${input_basename##*.}

if [[ -f "$sorted_bam" && "$force_bam" != "1" ]]; then
  echo "[1/4] Reusing existing sorted BAM: $sorted_bam"
else
  case "$input_ext" in
    cram)
      echo "[1/4] Converting CRAM to sorted BAM..."
      samtools view -T "$reference_fasta" -b "$input_file" | \
        samtools sort -@ "$threads" -o "$sorted_bam" -
      ;;
    bam)
      echo "[1/4] Sorting BAM..."
      samtools sort -@ "$threads" -o "$sorted_bam" "$input_file"
      ;;
    *)
      echo "Unsupported input extension for $input_file (expected .bam or .cram)" >&2
      exit 1
      ;;
  esac
fi

if [[ -f "$sorted_bai" && "$force_index" != "1" ]]; then
  echo "[2/4] Reusing existing BAM index: $sorted_bai"
else
  echo "[2/4] Indexing BAM..."
  samtools index "$sorted_bam"
fi

echo "[3/4] Running QDNAseq..."
r_args=("$qdna_rscript" "$sorted_bam" "$sample_id" "$binsize")
if [[ "$continue_aborted" == "1" ]]; then
  r_args+=("--continue-aborted")
fi
Rscript "${r_args[@]}"

echo "[4/4] Generating chromosome comparison plot and report..."
python3 "$plot_script" "${sample_id}_log2.bed" "$control_bed" --plot-style dots
python3 "$report_script" "$sample_id"

echo "Done. Main outputs:"
echo "  - $sorted_bam"
echo "  - ${sorted_bam}.bai"
echo "  - ${sample_id}_QDNAseq_genome_plot.png"
echo "  - ${sample_id}_segments_log2.seg"
echo "  - ${sample_id}_segments_with_CN_class.seg"
echo "  - ${sample_id}_log2.vcf"
echo "  - ${sample_id}_chr_boxwhisker.png"
echo "  - ${sample_id}_chromosome_summary.tsv"
echo "  - ${sample_id}_report.md"
echo "  - ${sample_id}_report.pdf"
