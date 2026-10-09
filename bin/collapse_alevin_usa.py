#!/usr/bin/env python3
"""
Collapse an alevin-fry USA count matrix to one column per gene.

In USA mode -- which this pipeline always runs, since ``alevin-fry quant`` is
given a 3-column ``t2g_3col.tsv`` -- every gene occupies three columns of
``quants_mat.mtx``: its spliced, unspliced and ambiguous counts, in three
equal blocks named ``<gene>``, ``<gene>-U`` and ``<gene>-A`` (see
``alevin_usa.py``). Handed to a downstream tool unchanged, that matrix presents
each gene three times as three correlated features, which distorts the
highly-variable-gene selection and PCA that Scrublet, scDblFinder and the
ambient-RNA callers all rely on.

This sums the requested blocks per gene, so the result is a gene-level matrix
whose feature names are plain gene IDs -- directly comparable to STARsolo's
``features.tsv``.

Which blocks to sum is the choice of what counts as expression:

  SUA  spliced + unspliced + ambiguous. The counterpart of STARsolo's
       GeneFull_Ex50pAS, which counts reads over exons *and* introns, and
       therefore the only setting under which the two mappers' matrices are
       comparable.
  SA   spliced + ambiguous, the conventional single-cell count, which discards
       intronic signal. The counterpart of STARsolo's Gene.
  S    spliced only.
  UA   unspliced + ambiguous.
  U    unspliced only -- the intronic matrix, for RNA velocity. The counterpart
       of the unspliced matrix in STARsolo's Velocyto feature.

Barcodes are never touched: cells stay in their original order, and no cell or
gene is dropped, so the matrix keeps its full feature axis whichever blocks
were summed.

A matrix that is not a USA column set is an error, not something to copy
through: alevin-fry always runs in USA mode here, so such a matrix can only
mean an upstream bug, and passing it on would hand every downstream step each
gene three times without anything failing.
"""

import argparse
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
        description="Collapse an alevin-fry USA count matrix to one column per gene."
    )
    parser.add_argument("-d", "--dir", required=True, help="alevin-fry quant matrix directory")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for the gene-level matrix")
    parser.add_argument(
        "-c", "--counts", default="SUA", choices=["SUA", "SA", "S", "UA", "U"],
        help="USA blocks to sum into each gene's count (default: SUA)",
    )
    return parser.parse_args()


def read_lines(path: str) -> List[str]:
    """Read a one-name-per-line text file, dropping blank lines."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def write_lines(path: str, names: List[str]) -> None:
    """Write one name per line."""
    with open(path, "w") as fh:
        for name in names:
            fh.write(f"{name}\n")


def load_matrix(dirpath: str) -> Tuple[sp.csr_matrix, List[str], List[str]]:
    """Load an alevin-fry quant matrix as ``(cells x columns, barcodes, columns)``."""
    matrix_path = os.path.join(dirpath, _MATRIX_FILE)
    rows_path = os.path.join(dirpath, _ROWS_FILE)
    cols_path = os.path.join(dirpath, _COLS_FILE)

    for path in (matrix_path, rows_path, cols_path):
        if not os.path.exists(path):
            raise SystemExit(f"Error: {path} not found; is {dirpath} an alevin-fry quant directory?")

    mat = sio.mmread(matrix_path).tocsr()
    barcodes = read_lines(rows_path)
    columns = read_lines(cols_path)

    if mat.shape[0] != len(barcodes) or mat.shape[1] != len(columns):
        raise SystemExit(
            f"Error: {_MATRIX_FILE} is {mat.shape[0]}x{mat.shape[1]} but "
            f"{_ROWS_FILE} has {len(barcodes)} barcodes and {_COLS_FILE} has "
            f"{len(columns)} columns."
        )

    return mat, barcodes, columns


def main() -> None:
    args = parse_args()
    mat, barcodes, columns = load_matrix(args.dir)

    try:
        gene_names = usa_genes(columns)
    except ValueError as err:
        raise SystemExit(f"Error: {os.path.join(args.dir, _COLS_FILE)} is not a USA column set: {err}")

    gene_mat = sum_blocks(mat, len(gene_names), args.counts).tocsr()

    os.makedirs(args.outdir, exist_ok=True)
    sio.mmwrite(os.path.join(args.outdir, _MATRIX_FILE), gene_mat)
    write_lines(os.path.join(args.outdir, _ROWS_FILE), barcodes)
    write_lines(os.path.join(args.outdir, _COLS_FILE), gene_names)

    empty_cells = int(np.sum(np.asarray(gene_mat.sum(axis=1)).ravel() == 0))
    if empty_cells:
        print(
            f"Warning: {empty_cells} of {len(barcodes)} cells have no counts left after "
            f"keeping only {args.counts}",
            file=sys.stderr,
        )

    print(
        f"Collapsed {len(columns)} USA columns to {len(gene_names)} genes "
        f"({args.counts}) for {len(barcodes)} cells -> {args.outdir}"
    )


if __name__ == "__main__":
    main()
