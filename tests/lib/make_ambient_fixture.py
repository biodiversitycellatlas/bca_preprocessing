#!/usr/bin/env python3
"""
Build synthetic STARsolo raw directories with a known cell/empty boundary, for the
matrix-aware selection in ``secondderiv_cellcalling_v5.py``.

``make_knee_fixture.py`` plants boundaries in a curve; v5 also reads what the
barcodes on either side of a bend look like, so these fixtures plant expression too.
Cells are drawn from distinct cell-type profiles, empty droplets and debris from a
Dirichlet-multinomial around an ambient profile that mixes the cell types (as lysed
cells do) with an extra mitochondrial-like share.  Two archetypes:

``two_tier``
    The 2026_Q3 ``spolac01_lib1`` failure: 300 high-RNA cells of one type above 2,700
    cells of two others, then the empty tail.  The curve bends at both 300 and 3,000.
    v4 at an expected count of 450 (0.15x) takes the bend at 300, which is a bend
    among cells; the cell/empty boundary is 3,000 at any expected count.

``shoulder``
    2,500 cells above 12,000 debris droplets of about 250 UMIs drawn from the ambient
    profile, then the empty tail: a bend whose two sides are both ambient, which an
    expected count set too high would hand the whole shoulder.

``gradual``
    The 2026_Q3 ``spolac01_lib2`` shape: 150 high-RNA cells, then 3,850 cells whose
    depths spread so widely (sigma 1.0 around 400 UMIs) that they thin out into a deep
    empty tail without a qualifying bend at the boundary.  The only bends that qualify
    are among cells, so the boundary has to be read off the ambient test.  Cells and
    empties overlap in depth, so no cutoff keeps every cell and only cells; the checks
    score it against the planted cells in ``truth/cells.tsv`` instead of a rank.

``split_pool``
    The Parse Biosciences shape: no plateau of empty droplets, only a thin tail of
    debris at a few UMIs, and a low-RNA cell type centred near 90 UMIs, so cells tail
    off to about 20 UMIs and many sit inside the 100-UMI ambient pool.  The cliff below
    the cells never recovers into a plateau, so no bend passes v4's recovery test.  v5
    then scans for the boundary, finds the ambient test's "empties" level among the
    shallow cells, and cuts far above them; scored against ``truth/cells.tsv``.

Importable, so the checks name the same constants that generated the data.
"""

import argparse
import os
from typing import Dict, Tuple

import numpy as np


SEED = 20261102
N_GENES = 300
MARKERS_PER_TYPE = 15
MARKER_BOOST = 30.0
CELL_CONCENTRATION = 200.0
AMBIENT_CONCENTRATION = 500.0
TYPE_SHARES_IN_AMBIENT = (0.2, 0.4, 0.4)

# Populations as (count, median UMIs, lognormal sigma, source); the source is a cell
# type, or "ambient" for droplets drawn from the ambient profile.
ARCHETYPES: Dict[str, Tuple[Tuple[int, float, float, object], ...]] = {
    "two_tier": ((300, 8000, 0.25, 0), (1350, 1200, 0.3, 1), (1350, 1200, 0.3, 2),
                 (30000, 25, 0.6, "ambient")),
    "shoulder": ((800, 1500, 0.3, 0), (850, 1500, 0.3, 1), (850, 1500, 0.3, 2),
                 (12000, 250, 0.15, "ambient"), (30000, 25, 0.6, "ambient")),
    "gradual": ((150, 6000, 0.2, 0), (1900, 400, 1.0, 1), (1950, 400, 1.0, 2),
                (30000, 50, 0.8, "ambient")),
    "split_pool": ((1500, 600, 0.5, 0), (1500, 350, 0.6, 1), (1000, 90, 0.6, 2),
                   (40000, 4, 0.7, "ambient")),
    "shallow_split_pool": ((1500, 200, 0.4, 0), (1500, 120, 0.4, 1), (1000, 45, 0.4, 2),
                           (40000, 4, 0.7, "ambient")),
}
PLANTED_BOUNDARY = {"two_tier": 3000, "shoulder": 2500, "gradual": 4000, "split_pool": 4000,
                    "shallow_split_pool": 4000}
DECOY_BEND = {"two_tier": 300, "shoulder": 14500, "gradual": 150, "split_pool": None, "shallow_split_pool": None}
EXPECTED_CELLS = {"two_tier": 3000, "shoulder": 2500, "gradual": 4000, "split_pool": 4000,
                  "shallow_split_pool": 4000}


def profiles(rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Three cell-type profiles with their own marker genes, and the ambient mixture of them."""
    base = rng.dirichlet(np.full(N_GENES, 0.5))
    types = []
    for k in range(3):
        prof = base.copy()
        markers = np.arange(k * MARKERS_PER_TYPE, (k + 1) * MARKERS_PER_TYPE) + 10
        prof[markers] *= MARKER_BOOST
        types.append(prof / prof.sum())
    types = np.vstack(types)
    ambient = np.asarray(TYPE_SHARES_IN_AMBIENT) @ types
    ambient[:5] *= 5.0  # the share lysis adds: mitochondrial transcripts and the like
    return types, ambient / ambient.sum()


def build(kind: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(counts as barcodes x genes, totals, is_cell), barcodes in random order as STAR writes them."""
    rng = np.random.default_rng(SEED)
    types, ambient = profiles(rng)
    rows, labels = [], []
    for n, median, sigma, source in ARCHETYPES[kind]:
        depth = np.maximum(np.round(np.exp(np.log(median) + sigma * rng.standard_normal(n))), 1).astype(np.int64)
        if source == "ambient":
            alpha = AMBIENT_CONCENTRATION * ambient
        else:
            alpha = CELL_CONCENTRATION * types[source]
        for d in depth:
            rows.append(rng.multinomial(d, rng.dirichlet(alpha)))
        labels.extend([source != "ambient"] * n)
    counts = np.vstack(rows)
    order = rng.permutation(len(counts))
    counts = counts[order]
    return counts, counts.sum(axis=1), np.asarray(labels)[order]


def precision_recall(totals: np.ndarray, is_cell: np.ndarray, cutoff: int) -> Tuple[float, float]:
    """Of the barcodes a cutoff keeps, the share that are planted cells; of the cells, the share kept."""
    kept = totals >= cutoff
    return float(is_cell[kept].mean()) if kept.any() else 0.0, float(kept[is_cell].mean())


def write(outdir: str, kind: str) -> str:
    """
    A ``GeneFull_Ex50pAS``-shaped directory: raw/{matrix.mtx,barcodes.tsv,features.tsv},
    the curve, and truth/cells.tsv listing the planted cells.
    """
    counts, totals, is_cell = build(kind)
    solo = os.path.join(outdir, f"{kind}_Solo.out", "GeneFull_Ex50pAS")
    raw = os.path.join(solo, "raw")
    os.makedirs(raw, exist_ok=True)
    with open(os.path.join(raw, "barcodes.tsv"), "w") as fh:
        fh.write("\n".join(f"BC{i:06d}" for i in range(len(counts))) + "\n")
    with open(os.path.join(raw, "features.tsv"), "w") as fh:
        fh.write("\n".join(f"g{g}\tgene{g}\tGene Expression" for g in range(N_GENES)) + "\n")
    cell, gene = np.nonzero(counts)
    with open(os.path.join(raw, "matrix.mtx"), "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        fh.write(f"{N_GENES} {len(counts)} {len(cell)}\n")
        fh.write("\n".join(f"{g + 1} {c + 1} {v}" for g, c, v in zip(gene, cell, counts[cell, gene])) + "\n")
    with open(os.path.join(solo, "UMIperCellSorted.txt"), "w") as fh:
        fh.write("\n".join(str(int(t)) for t in np.sort(totals[totals > 0])[::-1]) + "\n")
    os.makedirs(os.path.join(solo, "truth"), exist_ok=True)
    with open(os.path.join(solo, "truth", "cells.tsv"), "w") as fh:
        fh.write("\n".join(f"BC{i:06d}" for i in np.flatnonzero(is_cell)) + "\n")
    return solo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write synthetic STARsolo raw directories for the v5 checks.")
    parser.add_argument("-o", "--outdir", required=True, help="Directory the fixtures are written to.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for kind in ARCHETYPES:
        solo = write(args.outdir, kind)
        print(f"{kind}: {sum(n for n, *_ in ARCHETYPES[kind])} barcodes, boundary {PLANTED_BOUNDARY[kind]} -> {solo}")


if __name__ == "__main__":
    main()
