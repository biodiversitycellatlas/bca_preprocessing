#!/usr/bin/env python3
"""
Build a synthetic STARsolo ``Solo.out`` tree carrying velocity matrices, plus an
alevin-fry USA quant directory, for tests/checks/velocity_matrix.sh.

Every count is known by construction, and -- crucially -- the three velocity
layers hold *different* values for the same gene/cell. A subset or a layer
assignment that silently swaps spliced for unspliced would otherwise pass every
shape check, which is exactly the failure this fixture exists to catch.

The counts follow one rule:

    spliced[g, c]   = 100 * (c + 1) + (g + 1)
    unspliced[g, c] =  10 * (c + 1) + (g + 1)
    ambiguous[g, c] =            (g + 1)

so a value identifies its layer, its gene and its cell on sight, and the three
layers can never be confused with one another.

The GeneFull_Ex50pAS matrix is written with per-cell totals that fall away
steeply across the barcodes, so a UMI cutoff selects a known prefix of cells.

Usage: make_velocyto_fixture.py OUTDIR
"""

import os
import sys
from typing import List

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

# Deliberately small: the assertions name individual cells and genes.
N_GENES = 6
N_CELLS = 8

# The barcodes the GeneFull cell call keeps, as a prefix of the full list. The
# velocity subset is asserted against exactly this set.
N_CALLED_CELLS = 3


def barcode(i: int) -> str:
    """A distinguishable, fixed-width barcode."""
    return f"CELL{i:04d}"


def gene(i: int) -> str:
    """A distinguishable gene ID."""
    return f"GENE{i:04d}"


def layer_value(layer: str, g: int, c: int) -> int:
    """The count for *layer* at gene *g*, cell *c*; see the module docstring."""
    if layer == "spliced":
        return 100 * (c + 1) + (g + 1)
    if layer == "unspliced":
        return 10 * (c + 1) + (g + 1)
    if layer == "ambiguous":
        return g + 1
    raise ValueError(f"unknown layer: {layer}")


def write_lines(path: str, names: List[str]) -> None:
    """Write one name per line."""
    with open(path, "w") as fh:
        for name in names:
            fh.write(f"{name}\n")


def write_velocyto(dirpath: str, barcodes: List[str], features: List[str]) -> None:
    """Write a STARsolo Velocyto ``raw`` directory: three matrices, genes x cells."""
    os.makedirs(dirpath, exist_ok=True)

    for layer in ("spliced", "unspliced", "ambiguous"):
        dense = np.array(
            [[layer_value(layer, g, c) for c in range(len(barcodes))] for g in range(len(features))],
            dtype=np.int64,
        )
        # STARsolo writes genes x cells
        sio.mmwrite(os.path.join(dirpath, f"{layer}.mtx"), sp.csr_matrix(dense))

    write_lines(os.path.join(dirpath, "barcodes.tsv"), barcodes)
    write_lines(os.path.join(dirpath, "features.tsv"), features)


def write_genefull(dirpath: str, barcodes: List[str], features: List[str]) -> None:
    """Write a GeneFull_Ex50pAS ``raw`` directory whose per-cell totals fall away steeply.

    Cell 0 gets the most counts, cell 1 fewer, and so on, so that a UMI cutoff
    keeps a known prefix of the barcode list.
    """
    os.makedirs(dirpath, exist_ok=True)

    dense = np.array(
        [[(len(barcodes) - c) * 1000 + (g + 1) for c in range(len(barcodes))] for g in range(len(features))],
        dtype=np.int64,
    )
    sio.mmwrite(os.path.join(dirpath, "matrix.mtx"), sp.csr_matrix(dense))
    write_lines(os.path.join(dirpath, "barcodes.tsv"), barcodes)
    write_lines(os.path.join(dirpath, "features.tsv"), features)


def write_alevin_usa(dirpath: str, barcodes: List[str], features: List[str]) -> None:
    """Write an alevin-fry USA quant directory holding the same counts.

    alevin-fry writes cells x columns with three suffixed columns per gene, all
    ``-S`` first, then all ``-U``, then all ``-A`` -- the layout the collapse and
    velocity scripts both have to unpick.
    """
    os.makedirs(dirpath, exist_ok=True)

    columns = (
        [f"{g}-S" for g in features]
        + [f"{g}-U" for g in features]
        + [f"{g}-A" for g in features]
    )

    blocks = []
    for layer in ("spliced", "unspliced", "ambiguous"):
        blocks.append(
            np.array(
                [[layer_value(layer, g, c) for g in range(len(features))] for c in range(len(barcodes))],
                dtype=np.int64,
            )
        )

    sio.mmwrite(os.path.join(dirpath, "quants_mat.mtx"), sp.csr_matrix(np.hstack(blocks)))
    write_lines(os.path.join(dirpath, "quants_mat_rows.txt"), barcodes)
    write_lines(os.path.join(dirpath, "quants_mat_cols.txt"), columns)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {os.path.basename(sys.argv[0])} OUTDIR")

    outdir = sys.argv[1]
    barcodes = [barcode(i) for i in range(N_CELLS)]
    features = [gene(i) for i in range(N_GENES)]

    solo = os.path.join(outdir, "sampleA_starsolo_Solo.out")
    write_velocyto(os.path.join(solo, "Velocyto", "raw"), barcodes, features)
    write_genefull(os.path.join(solo, "GeneFull_Ex50pAS", "raw"), barcodes, features)
    write_alevin_usa(os.path.join(outdir, "alevin"), barcodes, features)

    # The cell call the velocity subset has to reproduce
    write_lines(os.path.join(outdir, "called_barcodes.tsv"), barcodes[:N_CALLED_CELLS])

    # A Solo.out with no Velocyto directory at all, for the "absent" case
    novelo = os.path.join(outdir, "sampleB_starsolo_Solo.out")
    write_genefull(os.path.join(novelo, "GeneFull_Ex50pAS", "raw"), barcodes, features)

    print(f"Wrote fixture with {N_CELLS} cells, {N_GENES} genes to {outdir}")
    print(f"Called cells: {', '.join(barcodes[:N_CALLED_CELLS])}")


if __name__ == "__main__":
    main()
