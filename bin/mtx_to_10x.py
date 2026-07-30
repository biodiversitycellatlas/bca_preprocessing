#!/usr/bin/env python3
"""
Repacks a raw mtx/barcodes/features triplet (as produced by STARsolo or
alevin-fry) into a 10x Genomics-style matrix directory (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz) for the Demuxafy CLI wrappers
(Scrublet.py, scDblFinder.R). This works directly on the pre-h5ad triplet
to avoid an unnecessary mtx -> h5ad -> mtx round-trip through MTX_TO_H5AD.

STARsolo's raw matrix is already genes (rows) x cells (columns), matching
the 10x convention. alevin-fry's quants_mat.mtx is cells (rows) x genes
(columns), so it is transposed -- the same orientation check used in
MTX_TO_H5AD (modules/local/tools/scanpy/main.nf).
"""

import argparse
import gzip
import logging
import os
import shutil
import sys

import pandas as pd
from scipy.io import mmread, mmwrite

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("MtxTo10x")


def main():
    parser = argparse.ArgumentParser(description="Repack an mtx/barcodes/features triplet into 10x format.")
    parser.add_argument("--mtx", type=str, required=True)
    parser.add_argument("--barcodes", type=str, required=True)
    parser.add_argument("--features", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    logger.info(f"Loading {args.mtx}")
    matrix = mmread(args.mtx).tocsr()

    barcodes = pd.read_csv(args.barcodes, header=None)[0].astype(str).values
    features = pd.read_csv(args.features, header=None, sep="\t")

    n_features, n_barcodes = len(features), len(barcodes)

    # Normalize orientation to genes (rows) x cells (columns), the 10x convention.
    if matrix.shape[0] == n_barcodes and matrix.shape[1] == n_features:
        logger.info("Input matrix is cells x genes; transposing to genes x cells")
        matrix = matrix.T.tocsr()
    elif matrix.shape[0] != n_features or matrix.shape[1] != n_barcodes:
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match {n_features} features / "
            f"{n_barcodes} barcodes in either orientation"
        )

    logger.info(f"Writing {n_barcodes} barcodes and {n_features} features")

    barcodes_path = os.path.join(args.outdir, "barcodes.tsv.gz")
    with gzip.open(barcodes_path, "wt") as f:
        f.write("\n".join(barcodes) + "\n")

    features_path = os.path.join(args.outdir, "features.tsv.gz")
    with gzip.open(features_path, "wt") as f:
        for row in features.itertuples(index=False):
            gene_id = str(row[0])
            gene_name = str(row[1]) if len(row) > 1 else gene_id
            feature_type = str(row[2]) if len(row) > 2 else "Gene Expression"
            f.write(f"{gene_id}\t{gene_name}\t{feature_type}\n")

    matrix_path = os.path.join(args.outdir, "matrix.mtx")
    mmwrite(matrix_path, matrix)
    with open(matrix_path, "rb") as f_in, gzip.open(f"{matrix_path}.gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(matrix_path)

    logger.info(f"10x-format export written to {args.outdir}")


if __name__ == "__main__":
    main()
