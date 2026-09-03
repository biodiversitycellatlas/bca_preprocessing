#!/usr/bin/env python3
"""
Reading and writing the (mtx, barcodes, features) count-matrix triplet that every
step of the filtering workflow now works on.

The two mappers hand over the same three files in opposite orientations:
STARsolo's ``matrix.mtx`` is genes (rows) x cells (columns), the 10x convention,
while alevin-fry's ``quants_mat.mtx`` is cells (rows) x genes (columns) --
``quants_mat_rows.txt`` holds the barcodes. Resolving that is the only
mapper-specific decision left downstream of ``COLLAPSE_ALEVIN_USA``, so it lives
here once rather than in each consumer: an orientation resolved differently in
two places would give an object of exactly the right shape carrying the
transposed numbers, and nothing would raise.

Both readers build their result straight from the coordinate arrays
MatrixMarket is parsed into, in whichever orientation the caller wants. A raw
matrix is the largest thing this pipeline holds in memory, and the obvious
spelling -- read, convert to CSR, transpose, convert back -- keeps three copies
of its nonzeros alive to produce one.

Counts are handed back at the precision they were written at. ``downcast_counts``
would halve that again, but it is opt-in per call (``downcast=True``) rather than
the default, so what the pipeline publishes is bit for bit what the mapper
quantified.

Imported by ``mtx_to_10x.py``, ``filter_doublets.py``, ``run_cellsweep.py`` and
``mtx_to_h5ad.py``. Nextflow puts ``bin/`` on PATH and Python puts a script's own
directory on ``sys.path``, so a plain ``import mtx_io`` works from any of them
(as ``scirocket_demux_rocket.py`` already relies on).

``anndata`` is imported lazily: not every module that needs the reader carries it.
"""

import logging
import os

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

logger = logging.getLogger("MtxIO")

_INT32_MAX = np.iinfo(np.int32).max


def downcast_counts(data):
    """
    Halve the bytes per nonzero, halving the peak memory of a step that reads a
    matrix and speeding up the h5ad write with it. Not applied by default: the
    readers take ``downcast=True`` for a caller that wants it.

    Counts do not need 64 bits. int32 holds 2.1e9, against a per-cell-per-gene
    count that tops out in the tens of thousands. alevin-fry's counts are
    fractional -- ``quant -r cr-like-em`` splits multi-mapping UMIs by
    expectation-maximisation -- and float32 represents every integer below
    16,777,216 exactly and any value to ~7 significant digits, i.e. a relative
    error under 6e-8. On an EM-split count of 3.43 that is an absolute error
    around 2e-7, some seven orders of magnitude below the Poisson noise on the
    count itself, and float32 is the precision the downstream scanpy steps
    normalise and decompose in anyway.

    An input that would not survive the narrowing is left as it is rather than
    wrapped around.
    """
    if np.issubdtype(data.dtype, np.integer):
        if data.size and np.abs(data).max() > _INT32_MAX:
            logger.warning("Counts exceed int32; keeping %s", data.dtype)
            return data
        return data.astype(np.int32, copy=False)

    if np.issubdtype(data.dtype, np.floating):
        return data.astype(np.float32, copy=False)

    return data


def downcast_matrix(matrix):
    """
    ``downcast_counts`` applied to a matrix that has already been assembled --
    a layer read back out of an .h5ad, typically. Sparse matrices are narrowed in
    place, so nothing is copied.

    Nothing in the pipeline calls this at present; it is the assembled-matrix
    counterpart of ``downcast_counts``, for a caller that opts in.
    """
    if sp.issparse(matrix):
        matrix.data = downcast_counts(matrix.data)
        return matrix

    return downcast_counts(np.asarray(matrix))


def read_barcodes(path):
    """The barcode axis, as strings, in file order."""
    return pd.read_csv(path, header=None, sep="\t")[0].astype(str).tolist()


def read_features(path):
    """
    The feature axis as a DataFrame, keeping every column the file carries.

    STARsolo's ``features.tsv`` is gene_id / gene_name / feature_type; alevin-fry's
    ``quants_mat_cols.txt`` is a single column of gene names. Column 0 is the
    identifier in both cases.
    """
    return pd.read_csv(path, header=None, sep="\t")


def _read_coordinates(mtx, n_features, n_barcodes, downcast=False):
    """
    Parse the matrix once and hand back its coordinate arrays as
    ``(gene_index, cell_index, data)``, whichever way round the file was written.

    Raises when the matrix fits neither orientation, rather than silently
    producing a transposed object.
    """
    logger.info(f"Loading {mtx}")
    coo = sio.mmread(str(mtx))

    if coo.shape == (n_features, n_barcodes):
        logger.info(f"Input matrix is genes x cells ({n_features} x {n_barcodes})")
        gene_index, cell_index = coo.row, coo.col
    elif coo.shape == (n_barcodes, n_features):
        logger.info(f"Input matrix is cells x genes ({n_barcodes} x {n_features})")
        gene_index, cell_index = coo.col, coo.row
    else:
        raise ValueError(
            f"Matrix shape {coo.shape} does not match {n_features} features / "
            f"{n_barcodes} barcodes in either orientation"
        )

    data = downcast_counts(coo.data) if downcast else coo.data
    return gene_index, cell_index, data


def read_triplet(mtx, barcodes, features, downcast=False):
    """
    Read a triplet and return ``(matrix, barcodes, features)`` with the matrix
    normalised to genes (rows) x cells (columns) -- the orientation the 10x export
    and the doublet filter work in.
    """
    barcode_list = read_barcodes(barcodes)
    feature_df = read_features(features)

    gene_index, cell_index, data = _read_coordinates(
        mtx, len(feature_df), len(barcode_list), downcast=downcast
    )

    matrix = sp.csr_matrix(
        (data, (gene_index, cell_index)), shape=(len(feature_df), len(barcode_list))
    )
    return matrix, barcode_list, feature_df


def read_triplet_anndata(mtx, barcodes, features, sample_id=None, downcast=False):
    """
    The same triplet as an AnnData object: cells x genes, unique var names, and
    ``obs['sample_id']`` when a sample id is given.

    Assembled in AnnData's orientation directly, so a matrix that arrived as
    genes x cells is never transposed through a second full copy.
    """
    import anndata as ad

    barcode_list = read_barcodes(barcodes)
    feature_df = read_features(features)

    gene_index, cell_index, data = _read_coordinates(
        mtx, len(feature_df), len(barcode_list), downcast=downcast
    )

    adata = ad.AnnData(
        X=sp.csr_matrix(
            (data, (cell_index, gene_index)), shape=(len(barcode_list), len(feature_df))
        ),
        obs=pd.DataFrame(index=pd.Index(barcode_list, name=None)),
        var=pd.DataFrame(index=pd.Index(feature_df[0].astype(str).values, name=None)),
    )
    adata.var_names_make_unique()

    if sample_id is not None:
        adata.obs["sample_id"] = sample_id

    return adata


def write_triplet(outdir, matrix, barcodes, features,
                  mtx_name="matrix.mtx", barcodes_name="barcodes.tsv", features_name="features.tsv"):
    """
    Write a genes x cells matrix back out as a triplet.

    Always written in the 10x orientation, whichever mapper the input came from, so
    a rewritten triplet reads back the same way for both. The file names default to
    STARsolo's, which is what the workflow resolves a matrix directory into.
    """
    os.makedirs(outdir, exist_ok=True)

    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"Refusing to write a {matrix.shape} matrix as {len(features)} features x "
            f"{len(barcodes)} barcodes"
        )

    sio.mmwrite(os.path.join(outdir, mtx_name), matrix)

    pd.Series(barcodes).to_csv(
        os.path.join(outdir, barcodes_name), sep="\t", header=False, index=False
    )
    features.to_csv(
        os.path.join(outdir, features_name), sep="\t", header=False, index=False
    )

    logger.info(f"Triplet written to {outdir}: {matrix.shape[0]} genes x {matrix.shape[1]} cells")
