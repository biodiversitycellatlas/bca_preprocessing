#!/usr/bin/env python3
"""
The column layout of an alevin-fry USA-mode count matrix, defined once.

``alevin-fry quant`` runs in USA mode whenever it is given a 3-column
``t2g_3col.tsv`` -- which this pipeline always does -- and then writes three
columns per gene into ``quants_mat.mtx``. ``quants_mat_cols.txt`` names them in
three equal blocks, in a fixed order:

  1. the spliced counts (S), under the bare gene IDs;
  2. the unspliced counts (U), the same IDs in the same order plus ``-U``;
  3. the ambiguous counts (A), the same IDs in the same order plus ``-A``.

There is no ``-S`` suffix: alevin-fry's ``quant.rs`` has written this layout
since at least v0.8.0, and its own loaders (pyroe's ``load_fry``, simpleaf's
``af-anndata``) split the columns by position into thirds rather than by name.
This module does the same, and uses the names only to confirm the layout. A
gene whose own ID ends in ``-U`` or ``-A`` (``HLA-A``) is therefore handled like
any other, where stripping suffixes would have mangled it.

Only the standard library is used, and matrices are only sliced and added, so
the module imports into any of the environments its callers run in --
``collapse_alevin_usa.py``, ``secondderiv_alevin.py``,
``velocity_matrices_to_h5ad.py``, ``generate_dashboard.py`` and
``dashboard_mappingstats.py``. Nextflow puts ``bin/`` on PATH and Python puts a
script's own directory on ``sys.path``, so a plain ``import alevin_usa`` works.
"""

from typing import List, Sequence

# The blocks in the order alevin-fry writes them, with the suffix each block's names carry
BLOCKS = "SUA"
_SUFFIXES = {"S": "", "U": "-U", "A": "-A"}


def usa_genes(columns: Sequence[str]) -> List[str]:
    """The gene IDs of a USA column set: its first third.

    Raises ``ValueError`` when *columns* is not laid out as alevin-fry writes a
    USA matrix, naming the first column that breaks the layout -- a matrix that
    is not a USA set has no blocks to select, and guessing would hand every
    consumer each gene three times.
    """
    n_columns = len(columns)
    if n_columns == 0 or n_columns % 3 != 0:
        raise ValueError(
            f"{n_columns} columns cannot be a USA column set, which holds three columns per gene"
        )

    n_genes = n_columns // 3
    genes = list(columns[:n_genes])

    for block_index, block in enumerate(BLOCKS[1:], start=1):
        suffix = _SUFFIXES[block]
        offset = block_index * n_genes
        for i, gene in enumerate(genes):
            if columns[offset + i] != gene + suffix:
                raise ValueError(
                    f"column {offset + i + 1} is '{columns[offset + i]}' where a USA column set "
                    f"has '{gene}{suffix}': expected {n_genes} gene IDs, then the same IDs with "
                    f"'-U', then with '-A'"
                )

    return genes


def sum_blocks(mat, n_genes: int, blocks: str):
    """Sum the *blocks* (any of ``S``, ``U``, ``A``) of a cells x (3 x genes) matrix.

    Returns a cells x genes matrix of the same kind as *mat*. The blocks are
    column ranges, so this only slices and adds and works on any matrix that
    supports both (scipy sparse, numpy).
    """
    unknown = set(blocks) - set(BLOCKS)
    if not blocks or unknown:
        raise ValueError(f"'{blocks}' is not a selection of the USA blocks {BLOCKS}")

    total = None
    for block in BLOCKS:
        if block not in blocks:
            continue
        start = BLOCKS.index(block) * n_genes
        part = mat[:, start:start + n_genes]
        total = part if total is None else total + part
    return total
