#!/usr/bin/env Rscript
#
# scDblFinder doublet detection on a 10x matrix directory.
#
# Ported from Demuxafy v3.0.0 (scripts/scDblFinder.R,
# https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs, commit
# 3c90f6a37323b2a03dbe826867dc077f596a84da), so the pipeline no longer needs the
# Demuxafy singularity image. Same arguments, expected doublet rate, scDblFinder
# call and output files. Changes from the original:
#
# - the matrix is read with Matrix::readMM in place of Seurat's Read10X (gene IDs
#   from the first features column, made unique, as Read10X(gene.column = 1) does),
#   so the image needs no Seurat; a 10x .h5 file is no longer read;
# - the expected doublet rate per 1,000 cells is an argument (--doublet_rate_per_1000,
#   Demuxafy's 0.008 by default), capped at MAX_DOUBLET_RATE;
# - the random seed is fixed (--seed), so a rerun gives the same calls;
# - errors exit with a non-zero status (Demuxafy's q() exited 0), and the output
#   directory is created when missing;
# - the summary diagnostics' colapse= typo is fixed.
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

suppressMessages(suppressWarnings(library(argparse)))

# Above this the expected doublet rate stops meaning anything (Demuxafy's formula
# passes 1 at 125,000 cells at its default rate)
MAX_DOUBLET_RATE <- 0.5

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-o", "--out", required = TRUE, help = "The output directory where results will be saved")
parser$add_argument("-t", "--tenX_matrix", required = TRUE, type = "character", help = "Path to the 10x filtered matrix directory.")
parser$add_argument("-b", "--barcodes_filtered", required = FALSE, type = "character", help = "Path to a list of filtered barcodes to use for doublet detection.")
parser$add_argument("-d", "--doublet_rate_per_1000", required = FALSE, type = "double", default = 0.008, help = "Expected doublet rate per 1,000 cells; the expected rate is this times the number of cells / 1,000 (Demuxafy: 0.008).")
parser$add_argument("-s", "--seed", required = FALSE, type = "integer", default = 42, help = "Random seed.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

fail <- function(...) {
    message(...)
    quit(save = "no", status = 1)
}

suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(scDblFinder)))


# The first of `names` that exists in `dir`
find_file <- function(dir, names) {
    hits <- file.path(dir, names)
    hits <- hits[file.exists(hits)]
    if (length(hits) == 0) {
        fail("None of ", paste(names, collapse = ", "), " found in ", dir)
    }
    hits[1]
}

open_any <- function(path) {
    if (endsWith(path, ".gz")) gzfile(path) else file(path)
}

read_lines_any <- function(path) {
    con <- open_any(path)
    on.exit(close(con))
    readLines(con)
}

### Read in data as an sce object ###
if (!dir.exists(args$tenX_matrix)) {
    fail("Cannot find the counts matrix directory ", args$tenX_matrix)
}
print(paste0("Using the following counts: ", args$tenX_matrix))
matrix_file   <- find_file(args$tenX_matrix, c("matrix.mtx.gz", "matrix.mtx"))
barcodes_file <- find_file(args$tenX_matrix, c("barcodes.tsv.gz", "barcodes.tsv"))
features_file <- find_file(args$tenX_matrix, c("features.tsv.gz", "genes.tsv.gz", "features.tsv", "genes.tsv"))

counts <- as(readMM(matrix_file), "CsparseMatrix")
barcodes <- read_lines_any(barcodes_file)
# Base R rather than fread, which needs R.utils for a .gz file
features <- read.delim(open_any(features_file), header = FALSE, colClasses = "character", quote = "")
if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
    fail("Matrix is ", nrow(counts), " x ", ncol(counts), " but there are ", nrow(features),
         " features and ", length(barcodes), " barcodes")
}
rownames(counts) <- make.unique(as.character(features[[1]]))
colnames(counts) <- barcodes


if (!is.null(args$barcodes_filtered)) {
    if (!file.exists(args$barcodes_filtered)) {
        fail("Cannot find the provided filtered barcode list at ", args$barcodes_filtered, ".\n\nExiting.")
    }
    message("Reading in the filtered barcode list.")
    filtered_barcodes <- read_lines_any(args$barcodes_filtered)
    number_original_barcodes <- ncol(counts)
    barcodes_head <- head(colnames(counts))

    message("\nFiltering counts for barcodes included in the filtered barcode list.\n")
    ### Filter for the barcodes list of interest ###
    counts <- counts[, colnames(counts) %in% filtered_barcodes]

    ### Provide user informatin with the number of original barcodes, those in the filter list and those after filtering
    number_filtered_barcodes <- ncol(counts)

    if (number_filtered_barcodes == 0) {
        fail("\nThere were no barcodes left in counts matrix after filtering on the barcode list.\n",
             "Here is a list of the first few barcodes in the filter barcode list provided:\n",
             paste(head(filtered_barcodes), collapse = "\n"),
             "\nAnd here is a list of the first few barcodes in the original counts matrix:\n",
             paste(barcodes_head, collapse = "\n"),
             "\n\nPlease check for any inconsistencies between the barcode formatting.\n\nExiting\n")
    }
    message(paste0("\nThe original number of barcodes in the counts matrix: ", number_original_barcodes,
                   "\nThe number of barcodes in the user-provided barcode list :", length(filtered_barcodes),
                   "\nThe number of barcodes after filtering for user-provided barcodes: ", number_filtered_barcodes, "\n"))
}


sce <- SingleCellExperiment(list(counts = counts))

## Calculate doublet ratio ###
doublet_ratio <- ncol(sce) / 1000 * args$doublet_rate_per_1000
if (doublet_ratio > MAX_DOUBLET_RATE) {
    message(sprintf("WARNING: %d cells at %g doublets per 1,000 cells gives an expected doublet rate of %.3f; capped at %g.",
                    ncol(sce), args$doublet_rate_per_1000, doublet_ratio, MAX_DOUBLET_RATE))
    doublet_ratio <- MAX_DOUBLET_RATE
}
message(sprintf("Expected doublet rate: %.4f", doublet_ratio))


### Calculate Singlets and Doublets ###
set.seed(args$seed)
sce <- scDblFinder(sce, dbr = doublet_ratio)


### Make a dataframe of the results ###
results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)

dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
fwrite(results, file = file.path(args$out, "scDblFinder_doublets_singlets.tsv"), sep = "\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(results$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
fwrite(summary, file.path(args$out, "scDblFinder_doublet_summary.tsv"), sep = "\t")
