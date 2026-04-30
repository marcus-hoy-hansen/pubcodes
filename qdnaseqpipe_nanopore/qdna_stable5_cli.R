#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(QDNAseq)
  library(QDNAseq.hg38)
})

args <- commandArgs(trailingOnly = TRUE)

continue_aborted <- "--continue-aborted" %in% args
args <- args[args != "--continue-aborted"]

if (length(args) < 2) {
  stop("Usage: qdna_stable5_cli.R <sorted.bam> <sampleID> [binSize] [--continue-aborted]", call. = FALSE)
}

bamfile <- args[[1]]
sampleID <- args[[2]]
binSize <- if (length(args) >= 3) as.numeric(args[[3]]) else 100

if (is.na(binSize) || binSize <= 0) {
  stop("binSize must be a positive number", call. = FALSE)
}

readcounts_rds <- paste0(sampleID, "_readCounts.rds")
copynumber_rds <- paste0(sampleID, "_copyNumbers.rds")

message("Loading hg38 bins...")
bins <- getBinAnnotations(binSize = binSize, genome = "hg38")

if (file.exists(readcounts_rds)) {
  message("Loading cached readCounts...")
  readCounts <- readRDS(readcounts_rds)
} else {
  message("Counting reads from BAM...")
  readCounts <- binReadCounts(bins, bamfiles = bamfile)
  saveRDS(readCounts, readcounts_rds)
}

if (continue_aborted && file.exists(copynumber_rds)) {
  message("Loading cached copyNumbers...")
  copyNumbers <- readRDS(copynumber_rds)
} else {

  message("Filtering + normalization (with X/Y)...")

  rc_auto <- applyFilters(readCounts)
  rc_auto <- estimateCorrection(rc_auto)
  rc_fit <- Biobase::assayDataElement(rc_auto, "fit")

  rc_all <- applyFilters(readCounts, chromosomes = NA)
  rc_all <- correctBins(rc_all, fit = rc_fit)
  rc_all <- normalizeBins(rc_all)
  rc_all <- smoothOutlierBins(rc_all)

  readCounts <- rc_all

  message("Segmenting...")
  copyNumbers <- segmentBins(
    readCounts,
    alpha = 0.01,
    undo.splits = "sdundo",
    undo.SD = 0.80
  )

  message("Calling CNVs...")
  copyNumbers <- callBins(copyNumbers)

  saveRDS(copyNumbers, copynumber_rds)
}

message("Plotting autosomes only...")
copyNumbers_plot <- copyNumbers[QDNAseq::chromosomes(copyNumbers) %in% as.character(1:22), ]
png(paste0(sampleID, "_QDNAseq_genome_plot.png"), width = 1400, height = 450)
plot(copyNumbers_plot, main = paste("QDNAseq CNV profile:", sampleID))
dev.off()

message("Exporting SEG (log2 segments)...")
seg_file <- paste0(sampleID, "_segments_log2.seg")

QDNAseq::exportBins(
  copyNumbers,
  file = seg_file,
  format = "seg",
  type = "segments",
  logTransform = TRUE
)

seg <- read.table(
  seg_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"LOG2_RATIO_MEAN" %in% colnames(seg)) {
  stop("LOG2_RATIO_MEAN column not found in SEG file.", call. = FALSE)
}

baseline_ploidy <- 2
seg$absolute_cn <- baseline_ploidy * 2^(seg$LOG2_RATIO_MEAN)

gain_threshold <- 0.07
loss_threshold <- -0.075

seg$CNV_class <- "Neutral"
seg$CNV_class[seg$LOG2_RATIO_MEAN >= gain_threshold] <- "Gain"
seg$CNV_class[seg$LOG2_RATIO_MEAN <= loss_threshold] <- "Loss"

write.table(
  seg,
  file = paste0(sampleID, "_segments_with_CN_class.seg"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Exporting BED (log2)...")
QDNAseq::exportBins(
  copyNumbers,
  file = paste0(sampleID, "_log2.bed"),
  format = "bed",
  logTransform = TRUE
)

message("Exporting IGV (log2)...")
QDNAseq::exportBins(
  copyNumbers,
  file = paste0(sampleID, "_log2.igv"),
  format = "igv",
  logTransform = TRUE
)

message("Exporting filtered IGV (log2)...")
QDNAseq::exportBins(
  copyNumbers,
  file = paste0(sampleID, "_log2_filtered.igv"),
  format = "igv",
  filter = TRUE,
  logTransform = TRUE
)

message("Exporting VCF (log2)...")
QDNAseq::exportBins(
  copyNumbers,
  file = paste0(sampleID, "_log2.vcf"),
  format = "vcf",
  logTransform = TRUE
)

message("Done.")

