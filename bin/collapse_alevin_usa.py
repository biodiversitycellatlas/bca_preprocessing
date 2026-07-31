#!/usr/bin/env python3
"""
Collapse an alevin-fry USA count matrix to one column per gene.

In USA mode -- which this pipeline always runs, since ``alevin-fry quant`` is
given a 3-column ``t2g_3col.tsv`` -- every gene occupies three columns of
``quants_mat.mtx``, suffixed ``-S``, ``-U`` and ``-A`` for its spliced,
unspliced and ambiguous counts. Handed to a downstream tool unchanged, that
matrix presents each gene three times as three correlated features, which
distorts the highly-variable-gene selection and PCA that Scrublet, scDblFinder
and the ambient-RNA callers all rely on.

This sums the requested blocks per gene and strips the suffixes, so the result
is a gene-level matrix whose feature names are plain gene IDs -- directly
comparable to STARsolo's ``features.tsv``.

Which blocks to sum is the choice of what counts as expression:

  SUA  spliced + unspliced + ambiguous. The counterpart of STARsolo's
       GeneFull_Ex50pAS, which counts reads over exons *and* introns, and
       therefore the only setting under which the two mappers' matrices are
       comparable.
  SA   spliced + ambiguous, the conventional single-cell count, which discards
       intronic signal.
  S    spliced only. The counterpart of STARsolo's Gene.

Barcodes are never touched: cells stay in their original order, and no cell or
gene is dropped, so the matrix keeps its full feature axis whichever blocks
were summed.
"""

import argparse
import os
import re
import shutil
import sys
from typing import List, Optional, Tuple

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

# Suffixes alevin-fry appends to the spliced / unspliced / ambiguous column
# blocks of a USA-mode count matrix.
_USA_SUFFIX_RE = re.compile(r"-([SUA])$")

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
        "-c", "--counts", default="SUA", choices=["SUA", "SA", "S"],
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


def collapse_usa(
    mat: sp.csr_matrix, columns: List[str], blocks: str
) -> Optional[Tuple[sp.csr_matrix, List[str]]]:
    """Sum the *blocks* of each gene into a single column.

    Returns ``(cells x genes, gene_names)``, or ``None`` when *columns* is not a
    USA column set -- detected by requiring that stripping the suffixes yields
    exactly one third as many distinct names, so a non-USA reference, or a gene
    name legitimately ending in ``-S``, is left alone rather than mangled.

    Genes keep the order in which they first appear, which for alevin-fry's
    ``-S``-first layout is the reference's own gene order.
    """
    if len(columns) % 3 != 0:
        return None

    matches = [_USA_SUFFIX_RE.search(name) for name in columns]
    if not all(matches):
        return None

    stripped = [_USA_SUFFIX_RE.sub("", name) for name in columns]
    gene_names = list(dict.fromkeys(stripped))
    if len(gene_names) != len(columns) // 3:
        return None

    gene_index = {name: i for i, name in enumerate(gene_names)}
    keep = [i for i, match in enumerate(matches) if match.group(1) in blocks]
    if not keep:
        raise SystemExit(f"Error: no {blocks} columns found among {len(columns)} USA columns.")

    # (n_columns x n_genes) 0/1 aggregation matrix, carrying only the selected
    # blocks: one row per kept column, marking the gene that column belongs to.
    aggregator = sp.csr_matrix(
        (
            np.ones(len(keep), dtype=mat.dtype),
            (keep, [gene_index[stripped[i]] for i in keep]),
        ),
        shape=(len(columns), len(gene_names)),
    )
    return (mat @ aggregator).tocsr(), gene_names


def main() -> None:
    args = parse_args()
    mat, barcodes, columns = load_matrix(args.dir)

    os.makedirs(args.outdir, exist_ok=True)
    collapsed = collapse_usa(mat, columns, args.counts)

    if collapsed is None:
        # Not a USA matrix, so there is nothing to collapse and no block to
        # select: copy it through unchanged rather than guess at its layout.
        print(
            f"Warning: {len(columns)} columns in {args.dir} are not a USA column set; "
            f"copying the matrix through unchanged and ignoring --counts {args.counts}",
            file=sys.stderr,
        )
        for name in (_MATRIX_FILE, _ROWS_FILE, _COLS_FILE):
            shutil.copyfile(os.path.join(args.dir, name), os.path.join(args.outdir, name))
        return

    gene_mat, gene_names = collapsed

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
