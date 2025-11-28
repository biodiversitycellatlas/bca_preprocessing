#!/usr/bin/env Rscript

## Generate a QC report for a Salmon Alevin + alevin-fry run
## using the alevinQC package.
##
## Assumes the output directory contains three subfolders:
##   <sample_id>_run             -> alevin (map) output
##   <sample_id>_out_permit_knee -> permit list / knee step output
##   <sample_id>_counts          -> alevin-fry quantification output
##
## Usage:
##   Rscript alevin_QC.R /path/to/output_dir sample_id [report_dir]

## --- Load alevinQC -------------------------------------------------
library(alevinQC)

## --- Build alevin-fry directory paths ----------------------------------------

args <- commandArgs(trailingOnly = TRUE)
quant_dir  <- normalizePath(args[1])  # e.g. BCA..._counts
permit_dir <- normalizePath(args[2])  # e.g. BCA..._out_permit_knee
map_dir    <- normalizePath(args[3])  # e.g. BCA..._run
sample_id  <- args[4]
report_dir <- if (length(args) >= 5) normalizePath(args[5]) else map_dir

message("Using directories:")
message("  mapDir    : ", map_dir)
message("  permitDir : ", permit_dir)
message("  quantDir  : ", quant_dir)
message("  reportDir : ", report_dir)

## --- Check that required files are present -----------------------------------

alevinQC::checkAlevinFryInputFiles(
  mapDir    = map_dir,
  permitDir = permit_dir,
  quantDir  = quant_dir
)

## --- Remove existing Rmd file if present -------------------------------------
rmd_file <- file.path(report_dir, paste0(sample_id, "_alevinFry_QC.Rmd"))
if (file.exists(rmd_file)) {
  message("Existing Rmd found, removing: ", rmd_file)
  file.remove(rmd_file)
}

## --- Run QC report -----------------------------------------------------------

output_file <- paste0(sample_id, "_alevinFry_QC.html")

message("Generating QC report: ", file.path(report_dir, output_file))

alevinQC::alevinFryQCReport(
  mapDir       = map_dir,
  permitDir    = permit_dir,
  quantDir     = quant_dir,
  sampleId     = sample_id,
  outputFile   = output_file,
  outputDir    = report_dir,
  outputFormat = "html_document",
  forceOverwrite = TRUE,
  showCode     = FALSE,
  knitrProgress = FALSE,
  quiet        = FALSE
)

message("Done. Report written to: ", file.path(report_dir, output_file))
