#!/usr/bin/env python3
"""
alevin-fry counterpart of the second-derivative cell calling and matrix filtering.

alevin-fry writes its counts as a *cells x genes* matrix in ``quants_mat.mtx``
(the transpose of STARsolo's layout), with ``quants_mat_rows.txt`` naming the
cell barcodes and ``quants_mat_cols.txt`` naming the gene columns.  In USA mode
-- which this pipeline always runs, since ``alevin-fry quant`` is given a
3-column ``t2g_3col.tsv`` -- every gene occupies three columns, suffixed ``-S``,
``-U`` and ``-A`` for its spliced, unspliced and ambiguous counts.

Two subcommands, so the USA-aware matrix loader has a single definition:

``umis``
    Write the per-cell UMI totals sorted descending, one integer per line: the
    same format as STARsolo's ``UMIperCellSorted.txt``, so that
    ``secondderiv_cellcalling.py`` derives the cutoff for both mappers.

``filter``
    Keep the cells at or above a UMI cutoff and recompute the cell-level summary
    statistics, in the same JSON schema ``secondderiv_filter_matrices.py`` writes
    for STARsolo, so the dashboard reads both through one code path.

UMI totals sum all three USA blocks.  The STARsolo side of this pipeline calls
cells on GeneFull_Ex50pAS, which counts reads over exons *and* introns, so
spliced+unspliced+ambiguous is the comparable quantity and a cutoff means the
same thing on either mapper.
"""

import argparse
import json
import os
import re
import sys
from typing import List, Tuple

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

# Suffixes alevin-fry appends to the spliced / unspliced / ambiguous column
# blocks of a USA-mode count matrix.
_USA_SUFFIX_RE = re.compile(r"-[SUA]$")

_MATRIX_FILE = "quants_mat.mtx"
_ROWS_FILE = "quants_mat_rows.txt"
_COLS_FILE = "quants_mat_cols.txt"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Second-derivative cell calling helpers for alevin-fry USA matrices."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    umis = sub.add_parser(
        "umis", help="Write per-cell UMI totals, sorted descending, one per line."
    )
    umis.add_argument("-d", "--dir", required=True, help="alevin-fry quant matrix directory")
    umis.add_argument("-o", "--output", required=True, help="Output UMIperCellSorted-style text file")

    filt = sub.add_parser(
        "filter", help="Filter cells on a UMI cutoff and recompute cell-level statistics."
    )
    filt.add_argument("-d", "--dir", required=True, help="alevin-fry quant matrix directory")
    filt.add_argument("-c", "--cutoff", required=True, type=int, help="UMI cutoff threshold for filtering cells")
    filt.add_argument("-o", "--outdir", required=True, help="Output directory for the filtered matrix")
    filt.add_argument("-s", "--stats", default="secondderiv_statistics.json", help="Output JSON file with the recomputed statistics")

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


def collapse_usa(mat: sp.csr_matrix, columns: List[str]) -> sp.csr_matrix:
    """Sum the spliced / unspliced / ambiguous blocks of each gene into one column.

    Returns a ``cells x genes`` matrix.  The collapse is applied only when
    stripping the USA suffixes yields exactly one third as many distinct names,
    so a non-USA reference -- or a gene name legitimately ending in ``-S`` --
    leaves the matrix untouched.  Needed for the gene-level statistics: without
    it a gene detected as spliced *and* unspliced would be counted twice.
    """
    if len(columns) % 3 != 0:
        return mat

    stripped = [_USA_SUFFIX_RE.sub("", name) for name in columns]
    gene_names = sorted(set(stripped))
    if len(gene_names) != len(columns) // 3:
        return mat

    gene_index = {name: i for i, name in enumerate(gene_names)}
    col_to_gene = np.fromiter((gene_index[name] for name in stripped), dtype=np.int64, count=len(stripped))

    # (n_columns x n_genes) 0/1 aggregation matrix: one row per matrix column,
    # marking the gene that column belongs to.
    aggregator = sp.csr_matrix(
        (np.ones(len(stripped), dtype=mat.dtype), (np.arange(len(stripped)), col_to_gene)),
        shape=(len(stripped), len(gene_names)),
    )
    return (mat @ aggregator).tocsr()


def umis_per_cell(mat: sp.csr_matrix) -> np.ndarray:
    """Total UMIs per cell (row sums)."""
    return np.asarray(mat.sum(axis=1)).ravel().astype(np.int64)


def cmd_umis(args: argparse.Namespace) -> None:
    """Write the descending per-cell UMI totals."""
    mat, barcodes, _columns = load_matrix(args.dir)
    totals = np.sort(umis_per_cell(mat))[::-1]

    with open(args.output, "w") as fh:
        for value in totals:
            fh.write(f"{int(value)}\n")

    print(
        f"Wrote UMI totals for {len(barcodes)} barcodes to {args.output}",
        file=sys.stderr,
    )


def cmd_filter(args: argparse.Namespace) -> None:
    """Filter cells on the UMI cutoff and write the matrix plus statistics."""
    mat, barcodes, columns = load_matrix(args.dir)
    totals = umis_per_cell(mat)

    keep = np.where(totals >= args.cutoff)[0]
    if len(keep) == 0:
        raise SystemExit(
            f"Error: no cells found meeting the threshold of {args.cutoff} UMIs."
        )

    filtered = mat[keep, :]
    filtered_barcodes = [barcodes[i] for i in keep]
    filtered_totals = totals[keep]

    # Gene-level statistics need the USA blocks collapsed, so that a gene is
    # counted once however its reads were assigned.
    by_gene = collapse_usa(filtered, columns)
    genes_per_cell = np.asarray((by_gene > 0).sum(axis=1)).ravel()
    total_genes_detected = int(np.sum(np.asarray(by_gene.sum(axis=0)).ravel() > 0))

    os.makedirs(args.outdir, exist_ok=True)

    # Written under alevin-fry's own file names so the filtered directory can be
    # consumed anywhere the unfiltered one is.
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

    print(f"Kept {len(keep)} of {len(barcodes)} barcodes at >= {args.cutoff} UMIs")
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
