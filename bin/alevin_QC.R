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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "Usage: Rscript run_alevin_fry_QC.R ",
      "/path/to/output_dir sample_id [report_dir]\n",
      "Example: Rscript run_alevin_fry_QC.R /data/alevin_out sample1"
    )
  )
}

base_dir   <- normalizePath(args[1])
sample_id  <- args[2]
report_dir <- if (length(args) >= 3) normalizePath(args[3]) else base_dir

## --- Load alevinQC -------------------------------------------------
library(alevinQC)

## --- Build alevin-fry directory paths ----------------------------------------

map_dir    <- file.path(base_dir, paste0(sample_id, "_run"))
permit_dir <- file.path(base_dir, paste0(sample_id, "_out_permit_knee"))
quant_dir  <- file.path(base_dir, paste0(sample_id, "_counts"))

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
