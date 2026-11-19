#!/usr/bin/env python3
"""
Fixtures for tests/checks/doublet_detection.sh: small 10x matrix directories
(matrix.mtx.gz genes x cells, barcodes.tsv.gz, features.tsv.gz), as MTX_TO_10X
writes them for the doublet callers.

    <outdir>/mixed/          N_TYPES cell types with their own marker genes, plus
                             N_DOUBLETS heterotypic doublets (two cells of different
                             types summed); the doublets' barcodes are listed in
                             <outdir>/mixed_doublets.txt
    <outdir>/homogeneous/    one cell type and no doublets, so Scrublet's simulated
                             doublets look like the cells and its score distribution
                             has no second mode to put a threshold at

Seeded, so every run writes the same matrices.
"""

import gzip
import os
import sys

import numpy as np
import scipy.io
import scipy.sparse as sp

N_GENES = 800
N_TYPES = 4
MARKERS_PER_TYPE = 60
N_SINGLETS = 600
N_DOUBLETS = 60
N_HOMOGENEOUS = 400
SEED = 7


def barcode(i):
    """A 16-nt barcode, unique per index."""
    bases = "ACGT"
    return "".join(bases[(i >> (2 * k)) & 3] for k in range(16)) + "-1"


def profiles(rng, n_types):
    """Per-type expression profiles: a shared background plus each type's markers."""
    base = rng.gamma(0.3, 1.0, N_GENES)
    out = []
    for t in range(n_types):
        p = base.copy()
        p[t * MARKERS_PER_TYPE:(t + 1) * MARKERS_PER_TYPE] *= 25
        out.append(p / p.sum())
    return out


def sample_cells(rng, profile, n, depth=(1500, 3000)):
    sizes = rng.integers(*depth, n)
    return np.stack([rng.poisson(profile * s) for s in sizes])


def write_10x(outdir, counts, barcodes):
    """counts: cells x genes."""
    os.makedirs(outdir, exist_ok=True)
    with gzip.open(os.path.join(outdir, "matrix.mtx.gz"), "wb") as f:
        scipy.io.mmwrite(f, sp.coo_matrix(counts.T))
    with gzip.open(os.path.join(outdir, "barcodes.tsv.gz"), "wt") as f:
        f.write("\n".join(barcodes) + "\n")
    with gzip.open(os.path.join(outdir, "features.tsv.gz"), "wt") as f:
        for g in range(N_GENES):
            f.write(f"GENE{g:05d}\tgene{g}\tGene Expression\n")


def main(outdir):
    rng = np.random.default_rng(SEED)

    # Mixed: singlets of every type, then heterotypic doublets, shuffled together
    types = profiles(rng, N_TYPES)
    per_type = N_SINGLETS // N_TYPES
    singlets = np.concatenate([sample_cells(rng, types[t], per_type) for t in range(N_TYPES)])
    labels = np.repeat(np.arange(N_TYPES), per_type)
    doublets = []
    for _ in range(N_DOUBLETS):
        a = rng.integers(len(singlets))
        b = rng.choice(np.flatnonzero(labels != labels[a]))
        doublets.append(singlets[a] + singlets[b])
    counts = np.concatenate([singlets, np.stack(doublets)])
    is_doublet = np.r_[np.zeros(len(singlets), bool), np.ones(N_DOUBLETS, bool)]
    order = rng.permutation(len(counts))
    counts, is_doublet = counts[order], is_doublet[order]
    barcodes = [barcode(i) for i in range(len(counts))]
    write_10x(os.path.join(outdir, "mixed"), counts, barcodes)
    with open(os.path.join(outdir, "mixed_doublets.txt"), "w") as f:
        f.write("\n".join(b for b, d in zip(barcodes, is_doublet) if d) + "\n")

    # Homogeneous: one type, no doublets
    one = profiles(rng, 1)[0]
    counts = sample_cells(rng, one, N_HOMOGENEOUS)
    write_10x(os.path.join(outdir, "homogeneous"), counts, [barcode(i) for i in range(N_HOMOGENEOUS)])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} OUTDIR")
    main(sys.argv[1])
