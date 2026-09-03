#!/usr/bin/env python3
"""
Build a synthetic count-matrix triplet in both mappers' orientations, plus the
doublet and CellSweep annotations that get merged into the published .h5ad, for
tests/checks/annotated_h5ad.sh.

Every count is known by construction and identifies its own gene and cell:

    raw[g, c]      = 100 * (c + 1) + (g + 1)
    denoised[g, c] =  10 * (c + 1) + (g + 1)

so a raw count can never be mistaken for a denoised one, and a transposed matrix
or an annotation joined by position instead of by name shows up as a wrong value
rather than as a wrong shape -- which is the failure this fixture exists to catch,
since every step here would pass a shape check either way.

The two triplets hold the *same* matrix in opposite orientations:

    starsolo/   matrix.mtx        genes x cells, features.tsv with three columns
    alevin/     quants_mat.mtx    cells x genes, quants_mat_cols.txt with one

The CellSweep object is deliberately awkward: its barcodes are shuffled and its
genes are a shuffled subset, so a merge that trusts positions cannot pass.

Usage: make_annotation_fixture.py OUTDIR
"""

import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# Deliberately small: the assertions name individual cells and genes.
N_GENES = 6
N_CELLS = 8

# The barcodes a cell call kept, as a prefix of the full list: doublets are only ever
# called on those, exactly as in the pipeline.
N_EVALUATED = 5

# Positions in the barcode list that the consensus calls doublets. Both are inside the
# evaluated prefix; DOUBLET_ONE_TOOL is called by a single tool, so AnySinglet keeps it.
DOUBLET_CELLS = [1, 3]
DOUBLET_ONE_TOOL = 4

# The gene CellSweep's output does not carry, so the merge has to reindex rather than
# assume the two axes line up.
CS_DROPPED_GENE = 4


def barcode(i):
    return f"CELL{i:03d}"


def gene(i):
    return f"GENE{i:03d}"


def raw_value(g, c):
    return 100 * (c + 1) + (g + 1)


def denoised_value(g, c):
    return 10 * (c + 1) + (g + 1)


def alpha_hat(c):
    return round(0.01 * (c + 1), 4)


def ambient_hat(g):
    return round(0.5 * (g + 1), 4)


def is_empty(c):
    return c >= N_EVALUATED


def celltype(c):
    return "empty" if is_empty(c) else str(c % 2)


def raw_matrix():
    """The genes x cells count matrix, dense: the fixture is tiny by design."""
    return np.array(
        [[raw_value(g, c) for c in range(N_CELLS)] for g in range(N_GENES)],
        dtype=np.int64,
    )


def barcodes():
    return [barcode(i) for i in range(N_CELLS)]


def genes():
    return [gene(i) for i in range(N_GENES)]


def doublet_barcodes():
    """The barcodes the AnySinglet consensus calls doublets."""
    return [barcode(i) for i in DOUBLET_CELLS]


def kept_barcodes():
    """The barcodes a doublet filter leaves behind, in source order."""
    return [bc for bc in barcodes() if bc not in doublet_barcodes()]


def write_starsolo(outdir):
    """STARsolo's genes x cells triplet."""
    os.makedirs(outdir, exist_ok=True)

    sio.mmwrite(os.path.join(outdir, "matrix.mtx"), sp.csr_matrix(raw_matrix()))

    with open(os.path.join(outdir, "barcodes.tsv"), "w") as f:
        f.write("\n".join(barcodes()) + "\n")

    with open(os.path.join(outdir, "features.tsv"), "w") as f:
        for name in genes():
            f.write(f"{name}\t{name}_sym\tGene Expression\n")


def write_alevin(outdir):
    """Alevin-fry's cells x genes triplet, holding the same counts transposed."""
    os.makedirs(outdir, exist_ok=True)

    sio.mmwrite(os.path.join(outdir, "quants_mat.mtx"), sp.csr_matrix(raw_matrix().T))

    with open(os.path.join(outdir, "quants_mat_rows.txt"), "w") as f:
        f.write("\n".join(barcodes()) + "\n")

    with open(os.path.join(outdir, "quants_mat_cols.txt"), "w") as f:
        f.write("\n".join(genes()) + "\n")


def write_combined_results(path):
    """
    Demuxafy's Combine_Results.R output: one row per evaluated barcode, one
    DropletType column per tool plus the AnySinglet consensus.
    """
    rows = []
    for c in range(N_EVALUATED):
        scrublet = "doublet" if c in DOUBLET_CELLS + [DOUBLET_ONE_TOOL] else "singlet"
        scdblfinder = "doublet" if c in DOUBLET_CELLS else "singlet"
        rows.append({
            "Barcode": barcode(c),
            "scrublet_DropletType": scrublet,
            "scDblFinder_DropletType": scdblfinder,
            # AnySinglet: a doublet only where every tool agrees
            "AnySinglet_DropletType": "doublet" if c in DOUBLET_CELLS else "singlet",
        })

    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_cellsweep_h5ad(path, aligned=False):
    """
    A stand-in for CellSweep's ``_cs_full.h5ad``: the denoised counts plus the obs
    and var columns MTX_TO_H5AD lifts out of it.

    Two variants, because the merge has two paths. ``aligned=True`` gives the axes
    CellSweep actually produces now that it reads the same triplet, which the merge
    takes as it stands; the default gives shuffled barcodes and a missing gene, so
    the general path cannot get away with merging by position.
    """
    import anndata as ad

    if aligned:
        obs_order = barcodes()
        var_order = genes()
    else:
        kept_genes = [gene(g) for g in range(N_GENES) if g != CS_DROPPED_GENE]
        # Shuffled deterministically: the same fixture every run, but never in axis order
        rng = np.random.default_rng(0)
        obs_order = list(rng.permutation(barcodes()))
        var_order = list(rng.permutation(kept_genes))

    cell_of = {barcode(c): c for c in range(N_CELLS)}
    gene_of = {gene(g): g for g in range(N_GENES)}

    X = np.array(
        [[denoised_value(gene_of[gn], cell_of[bc]) for gn in var_order] for bc in obs_order],
        dtype=np.int64,
    )

    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(
            {
                "is_empty": [is_empty(cell_of[bc]) for bc in obs_order],
                "celltype": [celltype(cell_of[bc]) for bc in obs_order],
                "alpha_hat": [alpha_hat(cell_of[bc]) for bc in obs_order],
            },
            index=pd.Index(obs_order),
        ),
        var=pd.DataFrame(
            {"ambient_hat": [ambient_hat(gene_of[gn]) for gn in var_order]},
            index=pd.Index(var_order),
        ),
    )
    adata.write_h5ad(path)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {os.path.basename(sys.argv[0])} OUTDIR")

    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    write_starsolo(os.path.join(outdir, "starsolo"))
    write_alevin(os.path.join(outdir, "alevin"))
    write_combined_results(os.path.join(outdir, "combined_results.tsv"))

    # The h5ads need anndata, which not every environment that runs the rest has.
    try:
        write_cellsweep_h5ad(os.path.join(outdir, "cellsweep_full.h5ad"))
        write_cellsweep_h5ad(os.path.join(outdir, "cellsweep_aligned.h5ad"), aligned=True)
    except ImportError:
        print("anndata not available; skipping the CellSweep fixtures", file=sys.stderr)

    print(f"Fixture written to {outdir}")


if __name__ == "__main__":
    main()
