#!/usr/bin/env python3
"""
STARsolo matrix filter and QC Statistics Generator

Filters a raw STARsolo count matrix on a UMI cutoff (typically the one derived
by ``secondderiv_cellcalling.py``) and recomputes the cell-level summary
statistics on the resulting matrix, so that the dashboard reports numbers that
describe the cells that were actually retained.
"""

import argparse
import json
import os
import sys
from typing import Dict, Optional, Set

import numpy as np
import pandas as pd
import scipy.io as sio


# Column of CellReads.stats holding the number of unique-gene reads per barcode, in preference order.
_UNIQUE_READ_COLUMNS = ["countedU", "featureU"]

# Column of CellReads.stats holding the number of deduplicated (unique) UMIs per barcode.
_UNIQUE_UMI_COLUMN = "nUMIunique"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Filter STARsolo raw matrices and report QC statistics.")
    parser.add_argument("-d", "--dir", required=True, help="Path to raw STARsolo directory (containing matrix.mtx, barcodes.tsv, features.tsv)")
    parser.add_argument("-c", "--cutoff", required=True, type=int, help="UMI cutoff threshold for filtering cells")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for filtered matrices")
    parser.add_argument("-s", "--stats", default="secondderiv_statistics.json", help="Output JSON file with the recomputed statistics")
    parser.add_argument("--cellreads", default=None, help="STARsolo CellReads.stats, used to recompute read-level statistics for the retained cells")
    return parser.parse_args()


def read_level_stats(cellreads_path: Optional[str], cell_barcodes: Set[str]) -> Dict[str, float]:
    """Recompute read-level statistics for *cell_barcodes* from ``CellReads.stats``.

    Returns a dict with ``fraction_unique_reads_in_cells``, ``reads_in_cells``
    and ``sequencing_saturation``, or ``{}`` when the statistics cannot be
    derived.  STARsolo reports all three in ``Summary.csv``, but against its own
    cell set; recomputing them here keeps them consistent with the filtered
    matrix.
    """
    if not cellreads_path or not os.path.exists(cellreads_path):
        return {}

    try:
        stats = pd.read_csv(cellreads_path, sep="\t", index_col=0)
    except Exception as exc:
        print(f"Warning: could not read {cellreads_path}: {exc}", file=sys.stderr)
        return {}

    column = next((c for c in _UNIQUE_READ_COLUMNS if c in stats.columns), None)
    if column is None:
        print(
            f"Warning: none of {_UNIQUE_READ_COLUMNS} present in {cellreads_path}; "
            "skipping read-level statistics.",
            file=sys.stderr,
        )
        return {}

    total_unique = float(stats[column].sum())
    if total_unique <= 0:
        return {}

    in_cells = stats.loc[stats.index.intersection(cell_barcodes)]
    reads_in_cells = float(in_cells[column].sum())

    out: Dict[str, float] = {
        "fraction_unique_reads_in_cells": reads_in_cells / total_unique,
        "unique_reads_in_cells": int(reads_in_cells),
    }

    if "cbMatch" in stats.columns:
        out["reads_in_cells"] = int(in_cells["cbMatch"].sum())

    # Sequencing saturation is the duplication rate of the reads counted into the
    # matrix, so it is a property of the cell set: 1 - unique UMIs / counted reads
    # over the retained cells. 
    if _UNIQUE_UMI_COLUMN in stats.columns and reads_in_cells > 0:
        umis_in_cells = float(in_cells[_UNIQUE_UMI_COLUMN].sum())
        out["sequencing_saturation"] = 1.0 - (umis_in_cells / reads_in_cells)

    return out


def main() -> None:
    args = parse_args()

    matrix_path = os.path.join(args.dir, "matrix.mtx")
    barcodes_path = os.path.join(args.dir, "barcodes.tsv")

    # Handle variations in STARsolo feature file naming
    features_path = os.path.join(args.dir, "features.tsv")
    if not os.path.exists(features_path):
        features_path = os.path.join(args.dir, "genes.tsv")

    # ------ Load Data ------

    # Load as CSC (Compressed Sparse Column) for efficient column/cell slicing
    mat = sio.mmread(matrix_path).tocsc()
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
    features = pd.read_csv(features_path, header=None, sep="\t")

    # Filter Cells based on UMI cutoff
    # Sum across rows (axis=0) to get UMIs per cell. .A1 flattens the matrix output to a 1D numpy array.
    umis_per_cell = mat.sum(axis=0).A1

    valid_cells_idx = np.where(umis_per_cell >= args.cutoff)[0]

    if len(valid_cells_idx) == 0:
        raise SystemExit(
            f"Error: no cells found meeting the threshold of {args.cutoff} UMIs."
        )

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
    sio.mmwrite(os.path.join(args.outdir, "matrix.mtx"), filtered_mat)

    # Write barcodes and features
    filtered_barcodes.to_csv(os.path.join(args.outdir, "barcodes.tsv"), sep="\t", header=False, index=False)
    features.to_csv(os.path.join(args.outdir, "features.tsv"), sep="\t", header=False, index=False)

    # Report statistics saved as json
    json_data = {
        "estimated_cells": int(len(valid_cells_idx)),
        "umi_threshold_applied": int(args.cutoff),
        "mean_umis_per_cell": float(mean_umi),
        "median_umis_per_cell": float(median_umi),
        "median_genes_per_cell": float(median_genes),
        "total_genes_detected": int(total_genes_detected),
    }
    json_data.update(read_level_stats(args.cellreads, set(filtered_barcodes.iloc[:, 0].astype(str))))

    with open(args.stats, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"Stats saved to: {args.stats}")


if __name__ == "__main__":
    main()
