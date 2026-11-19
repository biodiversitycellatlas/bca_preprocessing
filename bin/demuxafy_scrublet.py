#!/usr/bin/env python3
"""
Scrublet doublet detection on a 10x matrix directory.

Ported from Demuxafy v3.0.0 (scripts/Scrublet.py and mods/read10x.py,
https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs, commit
3c90f6a37323b2a03dbe826867dc077f596a84da), so the pipeline no longer needs the
Demuxafy singularity image. Same arguments, defaults, expected doublet rate and
output files. Changes from the original:

- the UMAP embedding and its plot are gone (the pipeline never published them),
  and so is reading a 10x .h5 file; the image no longer needs umap, numba or scanpy;
- the expected doublet rate per 1,000 cells is an argument (--doublet_rate_per_1000,
  Demuxafy's 0.008 by default) and is capped at MAX_DOUBLET_RATE;
- -f subsets the barcodes along with the matrix (Demuxafy left them unsubset, so
  the calls were written against the wrong barcodes);
- when Scrublet finds no threshold it exits with NO_THRESHOLD_EXIT and a message,
  where Demuxafy crashed in plot_histogram();
- the summary is written with its Classification / Droplet N columns under every
  pandas version.

Original licence:

MIT License

Copyright (c) 2024 drneavin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse
import os
import sys

# Above this the expected doublet rate stops meaning anything (Demuxafy's formula
# passes 1 at 125,000 cells at its default rate)
MAX_DOUBLET_RATE = 0.5

# Exit status when Scrublet's automatic threshold finds no bimodal split, which
# run_scrublet_doublet.py answers with a re-run at a manual threshold
NO_THRESHOLD_EXIT = 3

parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required=True, help="cell ranger counts matrix directory containing matrix files or full path to matrix.mtx.")
parser.add_argument("-b", "--barcodes", required=False, help="barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("-f", "--filtered_barcodes", required=False, default=None, help="File containing a filtered list of droplet barcodes. This may be used if you want to use a filtered list of barcodes for doublet detection (ie need to remove droplets that are empty or high in ambient RNA).")
parser.add_argument("-r", "--sim_doublet_ratio", required=False, default=2, type=int, help="Number of doublets to simulate relative to the number of observed transcriptomes.")
parser.add_argument("-c", "--min_counts", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("-e", "--min_cells", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("-v", "--min_gene_variability_pctl", required=False, default=85, type=int, help="Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].")
parser.add_argument("-p", "--n_prin_comps", required=False, default=30, type=int, help="Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.")
parser.add_argument("-t", "--scrublet_doublet_threshold", required=False, default=None, type=float, help="Manually Set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time")
parser.add_argument("-d", "--doublet_rate_per_1000", required=False, default=0.008, type=float, help="Expected doublet rate per 1,000 cells; the expected rate is this times the number of cells / 1,000 (Demuxafy: 0.008).")
parser.add_argument("-o", "--outdir", required=False, default=os.getcwd(), help="The output directory")
args = parser.parse_args()

import scipy.io
import scrublet as scr
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pandas as pd


def import_cellranger_mtx(mtx_directory):
    """
    Reads in a sparse matrix (.mtx) exported by CellRanger in barcodes-by-genes format

    Recognizes output from CellRanger version 2 (files: matrix.mtx) and
    CellRanger v3 (files: matrix.mtx.gz)
    """
    if '.mtx' in mtx_directory:
        mtx_file = mtx_directory  ### If the file path was all the way to the .mtx file
        mtx_directory = os.path.abspath(os.path.join(mtx_file, os.pardir))
    else:
        mtx_file = os.path.join(mtx_directory, "matrix.mtx")

    if not os.path.exists(mtx_file):
        mtx_file = mtx_file + ".gz"
        if not os.path.exists(mtx_file):
            raise Exception("Directory {} does not contain a recognizable matrix file".format(mtx_directory))
    return scipy.io.mmread(mtx_file).T.tocsc()


def read_barcodes(barcodes):
    """Reads a barcodes.tsv(.gz) file into a one-column ('Barcode') data frame."""
    return pd.read_csv(barcodes, sep="\t", header=None, names=["Barcode"])


if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
if not os.path.exists(args.counts_matrix):
    sys.exit("Couldn't find the counts file " + args.counts_matrix)
counts_matrix = import_cellranger_mtx(args.counts_matrix)
matrix_dir = args.counts_matrix if os.path.isdir(args.counts_matrix) else os.path.dirname(args.counts_matrix)

if args.barcodes is None:
    if os.path.exists(os.path.join(matrix_dir, "barcodes.tsv.gz")):
        print("Reading in barcodes file")
        barcodes_df = read_barcodes(os.path.join(matrix_dir, "barcodes.tsv.gz"))
    elif os.path.exists(os.path.join(matrix_dir, "barcodes.tsv")):
        print("Reading in barcodes file")
        barcodes_df = read_barcodes(os.path.join(matrix_dir, "barcodes.tsv"))
    else:
        sys.exit("No barcode file in provided counts matrix directory")
else:
    print("Reading in barcodes file")
    barcodes_df = read_barcodes(args.barcodes)


if args.filtered_barcodes is None:
    print("Will not filter barcodes as no barcode filtering file was provided.")
else:
    if not os.path.exists(args.filtered_barcodes):
        sys.exit("Cannot read filtered barcode file, please check the directory path and try again.\n"
                 "Interpreted path for filtered barcodes file: " + args.filtered_barcodes)
    print("Filtered barcodes exist.")
    barcodes_filtered_df = read_barcodes(args.filtered_barcodes)
    original_head = barcodes_df['Barcode'].head()
    keep = barcodes_df['Barcode'].isin(barcodes_filtered_df['Barcode']).to_numpy()
    counts_matrix = counts_matrix[keep]
    barcodes_df = barcodes_df[keep].reset_index(drop=True)
    if counts_matrix.shape[0] == 0:
        sys.exit("There are no barcodes remaining in your dataframe after filtering on the provided --filtered_barcodes file.\n"
                 "This is what the top your original barcodes looks like:\n {}\n\n"
                 "And this is what the filtering barcodes look like:\n {}\n"
                 "Please check that the provided filter barcode file is accurate for this data and has the same format."
                 .format(original_head, barcodes_filtered_df['Barcode'].head()))
    print('\nThe number of barcodes in the user-provided barcode list: {}.\n'
          'The number of barcodes after filtering for user-provided barcodes: {}'
          .format(barcodes_filtered_df.shape[0], counts_matrix.shape[0]))


dbl_rate = counts_matrix.shape[0] / 1000 * args.doublet_rate_per_1000
if dbl_rate > MAX_DOUBLET_RATE:
    print("WARNING: {} cells at {} doublets per 1,000 cells gives an expected doublet rate of {:.3f}; "
          "capped at {}.".format(counts_matrix.shape[0], args.doublet_rate_per_1000, dbl_rate, MAX_DOUBLET_RATE))
    dbl_rate = MAX_DOUBLET_RATE
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Expected doublet rate: {:.4f}'.format(dbl_rate))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio=args.sim_doublet_ratio)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts,
                                                          min_cells=args.min_cells,
                                                          min_gene_variability_pctl=args.min_gene_variability_pctl,
                                                          n_prin_comps=args.n_prin_comps)


if args.scrublet_doublet_threshold is None:
    if scrub.predicted_doublets_ is None:
        # Scrublet found no minimum in the simulated score distribution to put its
        # threshold at, and so called nothing
        print("Scrublet could not set a doublet threshold automatically (no bimodal split in the "
              "doublet scores). Run again with -t to set one.", file=sys.stderr)
        sys.exit(NO_THRESHOLD_EXIT)
    histogram_name = 'doublet_score_histogram.png'
    results_name = 'scrublet_results.tsv'
else:
    print(args.scrublet_doublet_threshold)
    scrub.call_doublets(threshold=args.scrublet_doublet_threshold)
    histogram_name = 'doublet_score_histogram_manual_threshold.png'
    results_name = 'scrublet_doublets_singlets.tsv'

### Plotting and saving
scrub.plot_histogram()
plt.savefig(os.path.join(args.outdir, histogram_name))

results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([barcodes_df, results, scores], axis=1)
dataframe["scrublet_DropletType"] = dataframe["scrublet_DropletType"].map({True: "doublet", False: "singlet"})

print("Writing results to {}.".format(os.path.join(args.outdir, results_name)))
dataframe.to_csv(os.path.join(args.outdir, results_name), sep="\t", index=False)


### Make summary of singlets and doublets and write to file ###
summary = dataframe["scrublet_DropletType"].value_counts().rename_axis("Classification").reset_index(name="Droplet N")

print("Writing summary to {}.".format(os.path.join(args.outdir, 'scrublet_summary.tsv')))
summary.to_csv(os.path.join(args.outdir, 'scrublet_summary.tsv'), sep="\t", index=False)

print("Done!")
