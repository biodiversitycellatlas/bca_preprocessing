#!/usr/bin/env python3
"""
Build a velocity-ready AnnData object from either mapper's splicing-aware output.

RNA velocity tools (scVelo, velocyto, CellRank) expect one object carrying
``layers["spliced"]`` and ``layers["unspliced"]``, not three loose matrices. Both
mappers in this pipeline can produce the underlying counts, in different shapes:

  --starsolo-dir  a STARsolo ``Velocyto`` directory, holding ``spliced.mtx``,
                  ``unspliced.mtx`` and ``ambiguous.mtx`` over a shared
                  ``barcodes.tsv`` / ``features.tsv``.
  --alevin-dir    an alevin-fry USA-mode quant directory, where the three blocks
                  are column ranges of one ``quants_mat.mtx``, suffixed ``-S``,
                  ``-U`` and ``-A`` in ``quants_mat_cols.txt``.

``X`` is the sum of the three layers -- the total count per gene per cell under
the velocity model -- so the object is usable directly, and the layers carry the
splicing breakdown. Note that ``X`` is *not* expected to equal STARsolo's
``GeneFull_Ex50pAS`` matrix: the two features assign reads to genes by different
rules.
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# The three layers, in the order they are reported. Keys are the AnnData layer names.
_LAYERS = ("spliced", "unspliced", "ambiguous")

# Suffixes alevin-fry appends to the spliced / unspliced / ambiguous column blocks
# of a USA-mode count matrix.
_USA_SUFFIX_RE = re.compile(r"-([SUA])$")

# alevin-fry USA block letter -> layer name
_USA_BLOCK_TO_LAYER = {"S": "spliced", "U": "unspliced", "A": "ambiguous"}

_ALEVIN_MATRIX_FILE = "quants_mat.mtx"
_ALEVIN_ROWS_FILE = "quants_mat_rows.txt"
_ALEVIN_COLS_FILE = "quants_mat_cols.txt"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a velocity-ready AnnData object with spliced/unspliced/ambiguous layers."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--starsolo-dir", help="STARsolo Velocyto directory (spliced.mtx, unspliced.mtx, ambiguous.mtx)")
    source.add_argument("--alevin-dir", help="alevin-fry USA-mode quant directory (quants_mat.mtx)")
    parser.add_argument("-o", "--out", required=True, help="Output .h5ad file")
    parser.add_argument("--sample-id", default=None, help="Value for adata.obs['sample_id']")
    return parser.parse_args()


def read_lines(path: str) -> List[str]:
    """Read a one-name-per-line text file, dropping blank lines."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_starsolo(dirpath: str) -> Tuple[Dict[str, sp.csr_matrix], List[str], List[str]]:
    """Load a STARsolo ``Velocyto`` directory as ``(layers, barcodes, features)``.

    Layers come back oriented cells x genes. STARsolo writes genes x cells, but the
    orientation is decided from the axis lengths rather than assumed, so a matrix
    that has already been transposed upstream is handled too.
    """
    barcodes_path = os.path.join(dirpath, "barcodes.tsv")
    features_path = os.path.join(dirpath, "features.tsv")
    if not os.path.exists(features_path):
        features_path = os.path.join(dirpath, "genes.tsv")

    for path in (barcodes_path, features_path):
        if not os.path.exists(path):
            raise SystemExit(f"Error: {path} not found; is {dirpath} a STARsolo Velocyto directory?")

    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t").iloc[:, 0].astype(str).tolist()
    features = pd.read_csv(features_path, header=None, sep="\t").iloc[:, 0].astype(str).tolist()

    layers: Dict[str, sp.csr_matrix] = {}
    for layer in _LAYERS:
        matrix_path = os.path.join(dirpath, f"{layer}.mtx")
        if not os.path.exists(matrix_path):
            raise SystemExit(f"Error: {matrix_path} not found; is {dirpath} a STARsolo Velocyto directory?")

        mat = sio.mmread(matrix_path).tocsr()

        if mat.shape == (len(features), len(barcodes)):
            mat = mat.T.tocsr()
        elif mat.shape != (len(barcodes), len(features)):
            raise SystemExit(
                f"Error: {layer}.mtx is {mat.shape[0]}x{mat.shape[1]}, which matches neither "
                f"{len(features)} genes x {len(barcodes)} cells nor its transpose."
            )

        layers[layer] = mat

    return layers, barcodes, features


def load_alevin(dirpath: str) -> Tuple[Dict[str, sp.csr_matrix], List[str], List[str]]:
    """Load an alevin-fry USA quant directory as ``(layers, barcodes, features)``.

    alevin-fry writes cells x columns, with three suffixed columns per gene. The
    blocks are split back apart by suffix; the gene order is the order in which
    each gene first appears, which for the ``-S``-first layout is the reference's
    own gene order.
    """
    matrix_path = os.path.join(dirpath, _ALEVIN_MATRIX_FILE)
    rows_path = os.path.join(dirpath, _ALEVIN_ROWS_FILE)
    cols_path = os.path.join(dirpath, _ALEVIN_COLS_FILE)

    for path in (matrix_path, rows_path, cols_path):
        if not os.path.exists(path):
            raise SystemExit(f"Error: {path} not found; is {dirpath} an alevin-fry quant directory?")

    mat = sio.mmread(matrix_path).tocsr()
    barcodes = read_lines(rows_path)
    columns = read_lines(cols_path)

    if mat.shape != (len(barcodes), len(columns)):
        raise SystemExit(
            f"Error: {_ALEVIN_MATRIX_FILE} is {mat.shape[0]}x{mat.shape[1]} but "
            f"{_ALEVIN_ROWS_FILE} has {len(barcodes)} barcodes and {_ALEVIN_COLS_FILE} "
            f"has {len(columns)} columns."
        )

    matches = [_USA_SUFFIX_RE.search(name) for name in columns]
    if len(columns) % 3 != 0 or not all(matches):
        raise SystemExit(
            f"Error: the {len(columns)} columns in {cols_path} are not a USA column set; "
            "velocity layers need alevin-fry run in USA mode against a splici reference."
        )

    stripped = [_USA_SUFFIX_RE.sub("", name) for name in columns]
    features = list(dict.fromkeys(stripped))
    if len(features) != len(columns) // 3:
        raise SystemExit(
            f"Error: stripping the USA suffixes from {len(columns)} columns yields "
            f"{len(features)} genes, expected {len(columns) // 3}."
        )

    gene_index = {name: i for i, name in enumerate(features)}

    layers: Dict[str, sp.csr_matrix] = {}
    for block, layer in _USA_BLOCK_TO_LAYER.items():
        keep = [i for i, match in enumerate(matches) if match.group(1) == block]
        if len(keep) != len(features):
            raise SystemExit(f"Error: found {len(keep)} '-{block}' columns, expected {len(features)}.")

        # (n_columns x n_genes) 0/1 selection matrix, carrying only this block: one row
        # per kept column, marking the gene that column belongs to. Reordering through it
        # keeps every layer on the same gene axis regardless of the column layout.
        selector = sp.csr_matrix(
            (
                np.ones(len(keep), dtype=mat.dtype),
                (keep, [gene_index[stripped[i]] for i in keep]),
            ),
            shape=(len(columns), len(features)),
        )
        layers[layer] = (mat @ selector).tocsr()

    return layers, barcodes, features


def main() -> None:
    args = parse_args()

    if args.starsolo_dir:
        layers, barcodes, features = load_starsolo(args.starsolo_dir)
    else:
        layers, barcodes, features = load_alevin(args.alevin_dir)

    # X is the total under the velocity model, so the object works without unpacking layers
    X = layers["spliced"] + layers["unspliced"] + layers["ambiguous"]

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(barcodes, name=None)),
        var=pd.DataFrame(index=pd.Index(features, name=None)),
        layers={name: layers[name] for name in _LAYERS},
    )
    adata.var_names_make_unique()
    if args.sample_id:
        adata.obs["sample_id"] = args.sample_id

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")

    totals = {name: int(layers[name].sum()) for name in _LAYERS}
    grand_total = sum(totals.values())
    if grand_total == 0:
        print("Warning: every layer is empty; the velocity object carries no counts.", file=sys.stderr)
    else:
        breakdown = ", ".join(f"{name} {100.0 * totals[name] / grand_total:.1f}%" for name in _LAYERS)
        print(f"Layer breakdown: {breakdown}")

    print(f"Wrote {adata.n_obs} cells x {adata.n_vars} genes with layers {list(_LAYERS)} -> {args.out}")


if __name__ == "__main__":
    main()
