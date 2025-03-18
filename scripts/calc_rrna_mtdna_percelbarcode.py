#!/usr/bin/env python3

import scipy.io
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict

# ############## SETTINGS ##############
# Define paths to the input files
workdir = "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/mapping_STARsolo_fastqsplitted_v6_merged"

# STARsolo matrix files for intronic calculation:
gene_full_mtx_path = workdir + "/TTTATCCTTG/TTTATCCTTG_Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx"
gene_mtx_path      = workdir + "/TTTATCCTTG/TTTATCCTTG_Solo.out/Gene/raw/matrix.mtx"
barcodes_path      = workdir + "/TTTATCCTTG/TTTATCCTTG_Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv"
features_path      = workdir + "/TTTATCCTTG/TTTATCCTTG_Solo.out/GeneFull_Ex50pAS/raw/features.tsv"

# BAM file for mtDNA calculation
bam_path           = workdir + "/TTTATCCTTG/TTTATCCTTG_Aligned.sortedByCoord.out.bam"

# Output CSV file.
output_csv_path    = "/users/asebe/bvanwaardenburg/git/data/250221_OAKseq_Nvec/mapping_STARsolo/metrics_per_CB/TTTATCCTTG_STARsolo_v6_merged_.csv"

# Name of the mitochondrial contig
mt_contig          = "OW052000.1"
# ######################################

# ---- Calculate Intronic Percentage from Matrices ----
print("Loading matrices...")

# Load matrices in Matrix Market format and convert to CSC for fast column operations
gene_full = scipy.io.mmread(gene_full_mtx_path).tocsc()
gene = scipy.io.mmread(gene_mtx_path).tocsc()

# Sum counts per cell (each column corresponds to a cell).
gene_full_counts = np.array(gene_full.sum(axis=0)).flatten()
gene_counts      = np.array(gene.sum(axis=0)).flatten()

# Calculate intronic percentage per cell as (gene_full - gene) / gene_full,
# while handling cells with zero counts
intronic_percent = np.divide(
    gene_full_counts - gene_counts, 
    gene_full_counts, 
    out=np.zeros_like(gene_full_counts, dtype=float), 
    where=gene_full_counts != 0
)

# Replace any negative values with 0
intronic_percent[intronic_percent < 0] = 0
print("Intronic percentage per cell (first few values):", intronic_percent[:5])

# ---- Load Cell Barcodes ----
print("Loading barcodes...")
# Load the barcode file and force it to a 1D series
barcodes = pd.read_csv(barcodes_path, header=None)
if barcodes.shape[1] > 1:
    barcodes = barcodes.iloc[:, 0]
else:
    barcodes = barcodes.squeeze()

# ---- Calculate Mitochondrial Percentage from BAM ----
print("Processing BAM file for mitochondrial read counts per cell...")

# Dictionaries to store total read counts and mitochondrial read counts per cell
total_counts_bam = defaultdict(int)
mt_counts_bam = defaultdict(int)

# Open the BAM file.
with pysam.AlignmentFile(bam_path, "rb") as bam:
    # Iterate over all reads in the BAM
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        try:
            cell_barcode = read.get_tag("CB")
        except KeyError:
            continue  # Skip reads without a cell barcode
        total_counts_bam[cell_barcode] += 1
        if read.reference_name == mt_contig:
            mt_counts_bam[cell_barcode] += 1

# Compute mitochondrial percentage per cell
mt_percent_dict = {}
for cell, tot in total_counts_bam.items():
    mt = mt_counts_bam.get(cell, 0)
    mt_percent_dict[cell] = mt / tot if tot > 0 else 0

# Build an array of mitochondrial percentages aligned with the barcodes from the matrix
mt_percent = np.array([mt_percent_dict.get(cell, 0) for cell in barcodes])
print("Mitochondrial percentage per cell (first few values):", mt_percent[:5])

# ---- Combine the Results and Save to CSV ----
result = pd.DataFrame({
    "Cell": barcodes,
    "IntronicPercent": intronic_percent,
    "MTPercent": mt_percent
})

print("Per-cell metrics (first few rows):")
print(result.head())

result.to_csv(output_csv_path, index=False)
print(f"Results saved to {output_csv_path}")
