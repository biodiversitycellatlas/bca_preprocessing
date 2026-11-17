#!/usr/bin/env python3
"""
STARsolo matrix filter and QC Statistics Generator

Three subcommands, the first two mirroring ``secondderiv_alevin.py``:

``umis``
    Write the per-cell UMI totals of a raw STARsolo matrix sorted descending, one
    integer per line, in the format of STARsolo's ``UMIperCellSorted.txt``. Needed
    when the analysed matrix is a multi-gene one (``UniqueAndMult-EM.mtx``):
    STAR's own curve counts unique-gene UMIs only, and the cutoff has to be read
    off totals on the same basis as the matrix it filters.

``filter``
    Filter the raw matrix on a UMI cutoff (typically the one derived by
    ``secondderiv_cellcalling.py``) and recompute the cell-level summary
    statistics on the result, so that the dashboard reports numbers that describe
    the cells that were actually retained.

``stats``
    Compute the matrix-level statistics of an already filtered matrix, without
    any cutoff or read-level statistics. Used for STARsolo's own cell call applied
    to a multi-gene matrix (SUBSET_STARSOLO_CELLS): STAR's ``Summary.csv``
    describes the unique-gene counts of those cells, not the matrix analysed.

``--matrix`` names the matrix inside the raw directory (``matrix.mtx`` or
``UniqueAndMult-<type>.mtx``). A multi-gene matrix holds fractional counts.
Totals are truncated to whole UMIs for the curve; the filter compares the exact
totals with the integer cutoff, which keeps the same cells (floor(x) >= c exactly
when x >= c). A gene counts as detected in a cell from one UMI on, so the small
fractions EM spreads over a gene's paralogues do not inflate the gene counts; for
integer matrices this is the same as any non-zero count.

The filtered matrix is always written as ``matrix.mtx``, whichever raw matrix it
was taken from.
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

# Minimum count at which a gene is detected in a cell
_DETECTED = 1


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Second-derivative cell calling helpers for STARsolo raw matrices.")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_matrix(p: argparse.ArgumentParser) -> None:
        p.add_argument("-d", "--dir", required=True, help="Path to raw STARsolo directory (containing the matrix, barcodes.tsv, features.tsv)")
        p.add_argument("-m", "--matrix", default="matrix.mtx", help="Matrix file inside --dir (default: matrix.mtx)")

    umis = sub.add_parser("umis", help="Write per-cell UMI totals, sorted descending, one per line.")
    add_matrix(umis)
    umis.add_argument("-o", "--output", required=True, help="Output UMIperCellSorted-style text file")

    filt = sub.add_parser("filter", help="Filter cells on a UMI cutoff and recompute cell-level statistics.")
    add_matrix(filt)
    filt.add_argument("-c", "--cutoff", required=True, type=int, help="UMI cutoff threshold for filtering cells")
    filt.add_argument("-o", "--outdir", required=True, help="Output directory for filtered matrices")
    filt.add_argument("-s", "--stats", default="secondderiv_statistics.json", help="Output JSON file with the recomputed statistics")
    filt.add_argument("--cellreads", default=None, help="STARsolo CellReads.stats, used to recompute read-level statistics for the retained cells")

    stats = sub.add_parser("stats", help="Compute the matrix-level statistics of a filtered matrix.")
    add_matrix(stats)
    stats.add_argument("-s", "--stats", required=True, help="Output JSON file with the statistics")
    stats.add_argument("--label", default=None, help="Raw matrix the filtered one was taken from, recorded as 'matrix' (default: --matrix)")
    return parser.parse_args()


def load_matrix(dirpath: str, matrix_name: str):
    """Load a raw STARsolo matrix as CSC (genes x cells), plus its barcodes and features."""
    matrix_path = os.path.join(dirpath, matrix_name)
    if not os.path.exists(matrix_path):
        raise SystemExit(
            f"Error: {matrix_path} not found. A multi-gene matrix (UniqueAndMult-<type>.mtx) is only "
            "written when star_soloMultiMappers lists that type."
        )
    barcodes_path = os.path.join(dirpath, "barcodes.tsv")

    # Handle variations in STARsolo feature file naming
    features_path = os.path.join(dirpath, "features.tsv")
    if not os.path.exists(features_path):
        features_path = os.path.join(dirpath, "genes.tsv")

    # CSC (Compressed Sparse Column) for efficient column/cell slicing
    mat = sio.mmread(matrix_path).tocsc()
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
    features = pd.read_csv(features_path, header=None, sep="\t")
    return mat, barcodes, features


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


def matrix_stats(mat) -> Dict[str, float]:
    """Cell-level summary statistics of a filtered matrix (genes x cells)."""
    umis_per_cell = mat.sum(axis=0).A1

    # Genes per cell (count elements >= 1 per column)
    genes_per_cell = (mat >= _DETECTED).sum(axis=0).A1

    return {
        "estimated_cells": int(mat.shape[1]),
        "mean_umis_per_cell": float(np.mean(umis_per_cell)),
        "median_umis_per_cell": float(np.median(umis_per_cell)),
        "median_genes_per_cell": float(np.median(genes_per_cell)),
        # Total genes detected (count rows where sum across cells >= 1)
        "total_genes_detected": int(np.sum(mat.sum(axis=1).A1 >= _DETECTED)),
    }


def cmd_umis(args: argparse.Namespace) -> None:
    """Write the descending per-cell UMI totals, truncated to whole UMIs."""
    mat, barcodes, _features = load_matrix(args.dir, args.matrix)
    totals = np.sort(np.floor(mat.sum(axis=0).A1).astype(np.int64))[::-1]

    with open(args.output, "w") as fh:
        for value in totals:
            fh.write(f"{int(value)}\n")

    print(f"Wrote {args.matrix} UMI totals for {len(barcodes)} barcodes to {args.output}", file=sys.stderr)


def cmd_filter(args: argparse.Namespace) -> None:
    """Filter cells on the UMI cutoff and write the matrix plus statistics."""
    mat, barcodes, features = load_matrix(args.dir, args.matrix)

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

    # ------ Save Filtered Output ------

    os.makedirs(args.outdir, exist_ok=True)

    # Write matrix
    sio.mmwrite(os.path.join(args.outdir, "matrix.mtx"), filtered_mat)

    # Write barcodes and features
    filtered_barcodes.to_csv(os.path.join(args.outdir, "barcodes.tsv"), sep="\t", header=False, index=False)
    features.to_csv(os.path.join(args.outdir, "features.tsv"), sep="\t", header=False, index=False)

    # Report statistics saved as json
    cell_stats = matrix_stats(filtered_mat)
    json_data = {
        "estimated_cells": cell_stats["estimated_cells"],
        "umi_threshold_applied": int(args.cutoff),
        **cell_stats,
    }
    json_data.update(read_level_stats(args.cellreads, set(filtered_barcodes.iloc[:, 0].astype(str))))

    with open(args.stats, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"Kept {len(valid_cells_idx)} of {len(barcodes)} barcodes at >= {args.cutoff} UMIs in {args.matrix}")
    print(f"Stats saved to: {args.stats}")


def cmd_stats(args: argparse.Namespace) -> None:
    """Write the matrix-level statistics of a filtered matrix, labelled with its source."""
    mat, _barcodes, _features = load_matrix(args.dir, args.matrix)
    json_data = {"matrix": args.label or args.matrix, **matrix_stats(mat)}

    with open(args.stats, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"Wrote statistics of {json_data['estimated_cells']} cells ({json_data['matrix']}) to {args.stats}")


def main() -> None:
    args = parse_args()
    if args.command == "umis":
        cmd_umis(args)
    elif args.command == "filter":
        cmd_filter(args)
    else:
        cmd_stats(args)


if __name__ == "__main__":
    main()
