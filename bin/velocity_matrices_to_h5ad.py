#!/usr/bin/env python3
"""
Build a velocity-ready AnnData object from either mapper's splicing-aware output.

RNA velocity tools (scVelo, velocyto, CellRank) expect one object carrying
``layers["spliced"]`` and ``layers["unspliced"]``, not three loose matrices. Both
mappers in this pipeline can produce the underlying counts, in different shapes:

  --starsolo-dir  a STARsolo ``Velocyto`` directory, holding ``spliced.mtx``,
                  ``unspliced.mtx`` and ``ambiguous.mtx`` over a shared
                  ``barcodes.tsv`` / ``features.tsv``.
  --alevin-dir    an alevin-fry USA-mode quant directory, where the spliced (S),
                  unspliced (U) and ambiguous (A) counts are three column blocks
                  of one ``quants_mat.mtx`` (see ``alevin_usa.py``).

Each mapper's layers follow that mapper's own convention, which is also what
makes the two comparable as velocity input:

  STARsolo     ``spliced``, ``unspliced`` and ``ambiguous`` are Velocyto's three
               disjoint matrices, which follow velocyto.py's rules. scVelo reads
               ``spliced`` and ``unspliced`` and ignores ``ambiguous``.
  alevin-fry   ``spliced`` is S + A and ``unspliced`` is U, as alevin-fry's
               velocity tutorial and pyroe's ``velocity`` output format define
               them: with a splici reference, A holds reads that fit an exon and
               its retained-intron flank equally well, and the alevin-fry paper
               found counting them as spliced changes the velocity graph only
               slightly. ``ambiguous`` is A on its own, for reference -- it is
               already inside ``spliced``.

``X`` is the total count per gene per cell (S + U + A for both mappers), so the
object is usable directly. It is *not* expected to equal STARsolo's
``GeneFull_Ex50pAS`` matrix: the two features assign reads to genes by different
rules. ``uns['velocity_layers']`` records what each layer holds.
"""

import argparse
import os
import sys
from typing import Dict, List, Tuple

import anndata as ad
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from alevin_usa import sum_blocks, usa_genes

# The three layers, in the order they are reported. Keys are the AnnData layer names.
_LAYERS = ("spliced", "unspliced", "ambiguous")

# What each layer holds, per source; written to uns['velocity_layers']
_STARSOLO_LAYERS = {
    "spliced": "STARsolo Velocyto spliced",
    "unspliced": "STARsolo Velocyto unspliced",
    "ambiguous": "STARsolo Velocyto ambiguous",
}

# alevin-fry USA blocks summed into each layer
_ALEVIN_LAYER_BLOCKS = {"spliced": "SA", "unspliced": "U", "ambiguous": "A"}

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


def load_alevin(
    dirpath: str,
) -> Tuple[Dict[str, sp.csr_matrix], sp.csr_matrix, Dict[str, int], List[str], List[str]]:
    """Load an alevin-fry USA quant directory as ``(layers, X, block_totals, barcodes, features)``.

    alevin-fry writes cells x columns, three column blocks per gene. The layers
    are built from the blocks as ``_ALEVIN_LAYER_BLOCKS`` defines, ``X`` is all
    three blocks, and *block_totals* are the disjoint S / U / A totals for the
    breakdown printed at the end.
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

    try:
        features = usa_genes(columns)
    except ValueError as err:
        raise SystemExit(
            f"Error: {cols_path} is not a USA column set ({err}); velocity layers need "
            "alevin-fry run in USA mode against a splici reference."
        )

    n_genes = len(features)
    layers = {
        layer: sum_blocks(mat, n_genes, blocks).tocsr()
        for layer, blocks in _ALEVIN_LAYER_BLOCKS.items()
    }
    X = sum_blocks(mat, n_genes, "SUA").tocsr()
    block_totals = {block: float(sum_blocks(mat, n_genes, block).sum()) for block in "SUA"}

    return layers, X, block_totals, barcodes, features


def main() -> None:
    args = parse_args()

    if args.starsolo_dir:
        layers, barcodes, features = load_starsolo(args.starsolo_dir)
        # Velocyto's three matrices are disjoint, so their sum is the total
        X = layers["spliced"] + layers["unspliced"] + layers["ambiguous"]
        layer_definitions = _STARSOLO_LAYERS
        breakdown_totals = {name: float(layers[name].sum()) for name in _LAYERS}
    else:
        layers, X, block_totals, barcodes, features = load_alevin(args.alevin_dir)
        layer_definitions = {
            layer: "alevin-fry " + "+".join(blocks) for layer, blocks in _ALEVIN_LAYER_BLOCKS.items()
        }
        breakdown_totals = {"S": block_totals["S"], "U": block_totals["U"], "A": block_totals["A"]}

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(barcodes, name=None)),
        var=pd.DataFrame(index=pd.Index(features, name=None)),
        layers={name: layers[name] for name in _LAYERS},
    )
    adata.var_names_make_unique()
    adata.uns["velocity_layers"] = dict(layer_definitions)
    if args.sample_id:
        adata.obs["sample_id"] = args.sample_id

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")

    grand_total = sum(breakdown_totals.values())
    if grand_total == 0:
        print("Warning: every layer is empty; the velocity object carries no counts.", file=sys.stderr)
    else:
        breakdown = ", ".join(
            f"{name} {100.0 * total / grand_total:.1f}%" for name, total in breakdown_totals.items()
        )
        print(f"Count breakdown: {breakdown}")

    layer_summary = ", ".join(f"{name} = {definition}" for name, definition in layer_definitions.items())
    print(f"Wrote {adata.n_obs} cells x {adata.n_vars} genes ({layer_summary}) -> {args.out}")


if __name__ == "__main__":
    main()
