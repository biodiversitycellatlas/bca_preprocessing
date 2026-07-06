#!/usr/bin/env python3
"""
STARsolo matrix filter and QC Statistics Generator
"""

import argparse
import os
import json
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.sparse import csc_matrix

def main() -> None:
    parser = argparse.ArgumentParser(description="Filter STARsolo raw matrices and report QC statistics.")
    parser.add_argument('-d', '--dir', required=True, type=str, help="Path to raw STARsolo directory (containing matrix.mtx, barcodes.tsv, features.tsv)")
    parser.add_argument('-c', '--cutoff', required=True, type=int, help="UMI cutoff threshold for filtering cells")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory for filtered matrices")
    args = parser.parse_args()

    # Define paths
    matrix_path = os.path.join(args.dir, 'matrix.mtx')
    barcodes_path = os.path.join(args.dir, 'barcodes.tsv')

    # Handle variations in STARsolo feature file naming
    features_path = os.path.join(args.dir, 'features.tsv')
    if not os.path.exists(features_path):
        features_path = os.path.join(args.dir, 'genes.tsv')

    try:
        # ------ Load Data ------

        # Load as CSC (Compressed Sparse Column) for efficient column/cell slicing
        mat = sio.mmread(matrix_path).tocsc()
        barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')
        features = pd.read_csv(features_path, header=None, sep='\t')

        # Filter Cells based on UMI cutoff
        # Sum across rows (axis=0) to get UMIs per cell. .A1 flattens the matrix output to a 1D numpy array.
        umis_per_cell = mat.sum(axis=0).A1

        valid_cells_idx = np.where(umis_per_cell >= args.cutoff)[0]

        if len(valid_cells_idx) == 0:
            print(f"Error: No cells found meeting the threshold of {args.cutoff} UMIs.")
            return

        # Slice matrix and barcodes
        filtered_mat = mat[:, valid_cells_idx]
        filtered_barcodes = barcodes.iloc[valid_cells_idx]
        filtered_umis_per_cell = umis_per_cell[valid_cells_idx]

        # ------ Calculate Statistics ------

        # Mean and median UMIs
        mean_umi = np.mean(filtered_umis_per_cell)
        median_umi = np.median(filtered_umis_per_cell)

        # Genes per cell (count elements > 0 per column)
        genes_per_cell = (filtered_mat > 0).sum(axis=0).A1
        median_genes = np.median(genes_per_cell)

        # Total genes detected (count rows where sum across filtered cells > 0)
        total_genes_detected = np.sum(filtered_mat.sum(axis=1).A1 > 0)

        # ------ Save Filtered Output ------

        os.makedirs(args.outdir, exist_ok=True)

        # Write matrix
        sio.mmwrite(os.path.join(args.outdir, 'matrix.mtx'), filtered_mat)

        # Write barcodes and features
        filtered_barcodes.to_csv(os.path.join(args.outdir, 'barcodes.tsv'), sep='\t', header=False, index=False)
        features.to_csv(os.path.join(args.outdir, 'features.tsv'), sep='\t', header=False, index=False)

        # Report statistics saved as json
        json_out_path = os.path.join("secondderiv_statistics.json")
        json_data = {
            "estimated_cells": int(len(valid_cells_idx)),
            "umi_threshold_applied": int(args.cutoff),
            "mean_umis_per_cell": float(mean_umi),
            "median_umis_per_cell": float(median_umi),
            "median_genes_per_cell": float(median_genes),
            "total_genes_detected": int(total_genes_detected)
        }
        with open(json_out_path, "w") as f:
            json.dump(json_data, f, indent=4)

        print(f"Stats saved to: {json_out_path}")

    except Exception as e:
        print(f"Critical Failure: {e}")

if __name__ == "__main__":
    main()
