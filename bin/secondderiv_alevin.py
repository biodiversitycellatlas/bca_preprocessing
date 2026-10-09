#!/usr/bin/env python3
"""
alevin-fry counterpart of the second-derivative cell calling and matrix filtering.

alevin-fry writes its counts as a *cells x genes* matrix in ``quants_mat.mtx``
(the transpose of STARsolo's layout), with ``quants_mat_rows.txt`` naming the
cell barcodes and ``quants_mat_cols.txt`` naming the gene columns.  In USA mode
-- which this pipeline always runs, since ``alevin-fry quant`` is given a
3-column ``t2g_3col.tsv`` -- every gene occupies three columns: its spliced,
unspliced and ambiguous counts, in three equal blocks named ``<gene>``,
``<gene>-U`` and ``<gene>-A`` (see ``alevin_usa.py``).

Two subcommands, so the USA-aware matrix loader has a single definition:

``umis``
    Write the per-cell UMI totals sorted descending, one integer per line: the
    same format as STARsolo's ``UMIperCellSorted.txt``, so that
    ``secondderiv_cellcalling.py`` derives the cutoff for both mappers.

``filter``
    Keep the cells at or above a UMI cutoff and recompute the cell-level summary
    statistics, in the same JSON schema ``secondderiv_filter_matrices.py`` writes
    for STARsolo, so the dashboard reads both through one code path.

Which blocks count towards a cell's UMI total is set by ``--counts``, and must be
the same choice that is made for the downstream matrix (``collapse_alevin_usa.py``,
``params.alevin_usa_counts``): the cutoff, the statistics reported against it and
the matrix that is analysed all have to be on one basis.  The default ``SUA``
matches STARsolo's GeneFull_Ex50pAS, which counts reads over exons *and* introns,
so a cutoff means the same thing on either mapper.

The filtered matrix itself keeps every USA column -- only its *cells* are
selected.  That leaves it a faithful subset of alevin-fry's own output, usable
anywhere the unfiltered matrix is, with the block selection applied once further
downstream.

A matrix that is not a USA column set is an error: alevin-fry always runs in
USA mode here, and counting every column instead would count each gene up to
three times in the gene statistics.
"""

import argparse
import json
import os
import sys
from typing import List, Tuple

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

from alevin_usa import sum_blocks, usa_genes

_MATRIX_FILE = "quants_mat.mtx"
_ROWS_FILE = "quants_mat_rows.txt"
_COLS_FILE = "quants_mat_cols.txt"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Second-derivative cell calling helpers for alevin-fry USA matrices."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_counts(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--counts", default="SUA", choices=["SUA", "SA", "S", "UA", "U"],
            help="USA blocks that count towards a cell's UMI total; must match "
                 "params.alevin_usa_counts (default: SUA)",
        )

    umis = sub.add_parser(
        "umis", help="Write per-cell UMI totals, sorted descending, one per line."
    )
    umis.add_argument("-d", "--dir", required=True, help="alevin-fry quant matrix directory")
    umis.add_argument("-o", "--output", required=True, help="Output UMIperCellSorted-style text file")
    add_counts(umis)

    filt = sub.add_parser(
        "filter", help="Filter cells on a UMI cutoff and recompute cell-level statistics."
    )
    filt.add_argument("-d", "--dir", required=True, help="alevin-fry quant matrix directory")
    filt.add_argument("-c", "--cutoff", required=True, type=int, help="UMI cutoff threshold for filtering cells")
    filt.add_argument("-o", "--outdir", required=True, help="Output directory for the filtered matrix")
    filt.add_argument("-s", "--stats", default="secondderiv_statistics.json", help="Output JSON file with the recomputed statistics")
    add_counts(filt)

    return parser.parse_args()


def load_matrix(dirpath: str) -> Tuple[sp.csr_matrix, List[str], List[str]]:
    """Load an alevin-fry quant matrix as ``(cells x genes, barcodes, columns)``.

    CSR is used because every operation here slices or sums over cells, which are
    the matrix rows.
    """
    matrix_path = os.path.join(dirpath, _MATRIX_FILE)
    rows_path = os.path.join(dirpath, _ROWS_FILE)
    cols_path = os.path.join(dirpath, _COLS_FILE)

    for path in (matrix_path, rows_path, cols_path):
        if not os.path.exists(path):
            raise SystemExit(f"Error: {path} not found; is {dirpath} an alevin-fry quant directory?")

    mat = sio.mmread(matrix_path).tocsr()
    barcodes = _read_lines(rows_path)
    columns = _read_lines(cols_path)

    if mat.shape[0] != len(barcodes) or mat.shape[1] != len(columns):
        raise SystemExit(
            f"Error: {_MATRIX_FILE} is {mat.shape[0]}x{mat.shape[1]} but "
            f"{_ROWS_FILE} has {len(barcodes)} barcodes and {_COLS_FILE} has "
            f"{len(columns)} columns."
        )

    return mat, barcodes, columns


def _read_lines(path: str) -> List[str]:
    """Read a one-name-per-line text file, dropping blank lines."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def gene_level(
    mat: sp.csr_matrix, columns: List[str], blocks: str
) -> Tuple[sp.csr_matrix, List[str]]:
    """The ``cells x genes`` matrix the cutoff and the statistics are derived from.

    Both need the blocks summed per gene: the UMI totals so that they are on the
    same basis as the analysed matrix (``collapse_alevin_usa.py`` sums the same
    blocks), and the gene counts so that a gene detected as spliced *and*
    unspliced is not counted twice.
    """
    try:
        gene_names = usa_genes(columns)
    except ValueError as err:
        raise SystemExit(f"Error: {_COLS_FILE} is not a USA column set: {err}")

    return sum_blocks(mat, len(gene_names), blocks).tocsr(), gene_names


def umis_per_cell(mat: sp.csr_matrix) -> np.ndarray:
    """Total UMIs per cell (row sums)."""
    return np.asarray(mat.sum(axis=1)).ravel().astype(np.int64)


def cmd_umis(args: argparse.Namespace) -> None:
    """Write the descending per-cell UMI totals."""
    mat, barcodes, columns = load_matrix(args.dir)
    by_gene, _gene_names = gene_level(mat, columns, args.counts)
    totals = np.sort(umis_per_cell(by_gene))[::-1]

    with open(args.output, "w") as fh:
        for value in totals:
            fh.write(f"{int(value)}\n")

    print(
        f"Wrote {args.counts} UMI totals for {len(barcodes)} barcodes to {args.output}",
        file=sys.stderr,
    )


def cmd_filter(args: argparse.Namespace) -> None:
    """Filter cells on the UMI cutoff and write the matrix plus statistics."""
    mat, barcodes, columns = load_matrix(args.dir)
    by_gene, _gene_names = gene_level(mat, columns, args.counts)
    totals = umis_per_cell(by_gene)

    keep = np.where(totals >= args.cutoff)[0]
    if len(keep) == 0:
        raise SystemExit(
            f"Error: no cells found meeting the threshold of {args.cutoff} UMIs."
        )

    # Only the cells are selected: every USA column is kept, so the result stays a faithful subset of alevin-fry's own output
    filtered = mat[keep, :]
    filtered_barcodes = [barcodes[i] for i in keep]
    filtered_totals = totals[keep]

    kept_by_gene = by_gene[keep, :]
    genes_per_cell = np.asarray((kept_by_gene > 0).sum(axis=1)).ravel()
    total_genes_detected = int(np.sum(np.asarray(kept_by_gene.sum(axis=0)).ravel() > 0))

    os.makedirs(args.outdir, exist_ok=True)

    # Written under alevin-fry's own file names so the filtered directory can be consumed anywhere the unfiltered one is.
    sio.mmwrite(os.path.join(args.outdir, _MATRIX_FILE), filtered)
    _write_lines(os.path.join(args.outdir, _ROWS_FILE), filtered_barcodes)
    _write_lines(os.path.join(args.outdir, _COLS_FILE), columns)

    json_data = {
        "estimated_cells": int(len(keep)),
        "umi_threshold_applied": int(args.cutoff),
        "mean_umis_per_cell": float(np.mean(filtered_totals)),
        "median_umis_per_cell": float(np.median(filtered_totals)),
        "median_genes_per_cell": float(np.median(genes_per_cell)),
        "total_genes_detected": total_genes_detected,
    }

    with open(args.stats, "w") as fh:
        json.dump(json_data, fh, indent=4)

    print(f"Kept {len(keep)} of {len(barcodes)} barcodes at >= {args.cutoff} {args.counts} UMIs")
    print(f"Stats saved to: {args.stats}")


def _write_lines(path: str, names: List[str]) -> None:
    """Write one name per line."""
    with open(path, "w") as fh:
        for name in names:
            fh.write(f"{name}\n")


def main() -> None:
    args = parse_args()
    if args.command == "umis":
        cmd_umis(args)
    else:
        cmd_filter(args)


if __name__ == "__main__":
    main()
