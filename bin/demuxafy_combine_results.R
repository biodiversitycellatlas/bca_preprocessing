#!/usr/bin/env Rscript
#
# Consensus doublet calls from the Scrublet and scDblFinder results, either or both.
#
# Ported from Demuxafy v3.0.0 (scripts/Combine_Results.R,
# https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs, commit
# 3c90f6a37323b2a03dbe826867dc077f596a84da), so the pipeline no longer needs the
# Demuxafy singularity image. Trimmed to the doublet-detecting path for the two
# tools the pipeline runs: same arguments for them (-r scrublet, -n scDblFinder,
# -m method, -o out), same left join on Barcode in Demuxafy's tool order (scDblFinder
# first), missing values set to "unassigned", the same four consensus rules, and the
# same output files:
#
#   <out>                               the joined per-barcode table
#   <base>_summary.tsv                  barcode counts per combination of tool calls
#   <base>_droplet_type_summary.tsv     barcode counts per consensus call (with -m)
#   <base>_w_combined_assignments.tsv   the joined table plus <method>_DropletType (with -m)
#
# With one tool given, every rule reduces to that tool's own call. Changes from the
# original: the demultiplexing tools, donor assignment and plots are gone (so no
# tidyverse, ComplexUpset or future.apply); barcodes are always read as text (fread
# would otherwise drop the leading zeros of all-digit barcodes); and errors exit
# with a non-zero status (Demuxafy's q() exited 0).
#
# Original licence:
#
# MIT License
#
# Copyright (c) 2024 drneavin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(argparse)))

METHODS <- c("MajoritySinglet", "AtLeastHalfSinglet", "AnySinglet", "AnyDoublet")

# create parser object
parser <- ArgumentParser()

parser$add_argument("-o", "--out", required = TRUE, type = "character", help = "The file where results will be saved")
parser$add_argument("-n", "--scDblFinder", required = FALSE, type = "character", default = NULL, help = "Path to scDblFinder results. Only use this option if you want to include the scDblFinder results.")
parser$add_argument("-r", "--scrublet", required = FALSE, type = "character", default = NULL, help = "Path to scrublet results. Only use this option if you want to include the scrublet results.")
parser$add_argument("-m", "--method", required = FALSE, type = "character", default = NULL, help = paste0("Combination method. Options are ", paste(METHODS, collapse = ", "), ". See https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html for detailed explanation of each intersectional method. Leave blank if you just want all the softwares to be merged into a single dataframe."))

args <- parser$parse_args()

fail <- function(...) {
    message(...)
    quit(save = "no", status = 1)
}

if (!is.null(args$method) && !(args$method %in% METHODS)) {
    fail("'", args$method, "' is not a combination method. Options are ", paste(METHODS, collapse = ", "), ".")
}

### Check to be sure that output directory already exists or create it for the user ###
out_dir <- dirname(args$out)
if (!dir.exists(out_dir)) {
    message(paste0("Creating directory as it doesn't exist. Complete path: ", out_dir))
    dir.create(out_dir, recursive = TRUE)
}
base <- tools::file_path_sans_ext(args$out)


# A tool's results: the file itself, or its default file name inside a directory
read_results <- function(path, tool, default_name) {
    message("Reading in ", tool, " results.")
    if (!file_test("-f", path)) {
        path <- file.path(path, default_name)
    }
    if (!file.exists(path)) {
        fail("Can't read in the ", tool, " report either from the file you provided or to find the '",
             default_name, "' file in the directory you provided. Please double check your path. Exiting.")
    }
    fread(path, sep = "\t", colClasses = list(character = "Barcode"))
}

### Make a list to store the results that will be combined in, in Demuxafy's tool order ###
results_list <- list()
if (!is.null(args$scDblFinder)) {
    results_list[["scDblFinder"]] <- read_results(args$scDblFinder, "scDblFinder", "scDblFinder_doublets_singlets.tsv")
}
if (!is.null(args$scrublet)) {
    results_list[["scrublet"]] <- read_results(args$scrublet, "scrublet", "scrublet_results.tsv")
}
if (length(results_list) < 1) {
    fail("You didn't provide any software results to combine. Please try again - this time providing at least one inputs.")
}
for (tool in names(results_list)) {
    if (!all(c("Barcode", paste0(tool, "_DropletType")) %in% colnames(results_list[[tool]]))) {
        fail("The ", tool, " results need 'Barcode' and '", tool, "_DropletType' columns; found: ",
             paste(colnames(results_list[[tool]]), collapse = ", "))
    }
}


message("\nCombining results.\n")
# A left join on Barcode in list order: the first tool's barcodes, in its order
left_join_dt <- function(x, y) {
    cols <- setdiff(colnames(y), "Barcode")
    joined <- copy(x)
    joined[y, (cols) := mget(paste0("i.", cols)), on = "Barcode"]
    joined
}
combined_results <- Reduce(left_join_dt, results_list)

### Update any missing data with "unassigned" ###
for (col in colnames(combined_results)) {
    missing <- is.na(combined_results[[col]])
    if (any(missing)) {
        set(combined_results, j = col, value = ifelse(missing, "unassigned", as.character(combined_results[[col]])))
    }
}

message("\nWriting output.\n")
fwrite(combined_results, args$out, sep = "\t", append = FALSE)

message("\nMaking and writing summary files.\n")
droplet_type_cols <- grep("DropletType", colnames(combined_results), value = TRUE)
combined_results_summary <- combined_results[, .(N = .N), by = droplet_type_cols]
combined_results_summary <- combined_results_summary[order(-N)]
fwrite(combined_results_summary, paste0(base, "_summary.tsv"), sep = "\t", append = FALSE)


## Call singlets and doublets using combined intersectional methods ##
if (is.null(args$method)) {
    message("\nNo method provided so the results will not be combined into a single assignment.\n")
} else {
    calls <- combined_results[, .SD, .SDcols = droplet_type_cols]
    n_singlet <- rowSums(calls == "singlet")
    n_tools <- length(droplet_type_cols)
    consensus <- switch(args$method,
        ### Only call a singlet if majority of softwares call a singlet
        MajoritySinglet    = ifelse(n_singlet > n_tools / 2, "singlet", "doublet"),
        ### Only call a singlet if at least half of softwares call a singlet
        AtLeastHalfSinglet = ifelse(n_singlet >= n_tools / 2, "singlet", "doublet"),
        ### Call a singlet if any softwares calls that droplet a singlet
        AnySinglet         = ifelse(n_singlet > 0, "singlet", "doublet"),
        ### Call a doublet if any software calls that droplet a doublet
        AnyDoublet         = ifelse(rowSums(calls == "doublet") > 0, "doublet", "singlet")
    )
    set(combined_results, j = paste0(args$method, "_DropletType"), value = consensus)

    droplet_type_summary <- as.data.frame(table(consensus))
    colnames(droplet_type_summary) <- c("Classification", "Droplet N")
    fwrite(droplet_type_summary, paste0(base, "_droplet_type_summary.tsv"), sep = "\t", append = FALSE)

    message("\nWriting output with combined calls.\n")
    fwrite(combined_results, paste0(base, "_w_combined_assignments.tsv"), sep = "\t", append = FALSE)
}
