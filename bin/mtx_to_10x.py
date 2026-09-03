#!/usr/bin/env python3
"""
Repacks a raw mtx/barcodes/features triplet (as produced by STARsolo or
alevin-fry) into a 10x Genomics-style matrix directory (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz) for the Demuxafy CLI wrappers
(Scrublet.py, scDblFinder.R). This works directly on the pre-h5ad triplet,
like every other step of the filtering workflow.

The genes x cells orientation the 10x convention expects is resolved by
``mtx_io.read_triplet``, which handles alevin-fry's transposed
``quants_mat.mtx`` for every consumer alike.
"""

import argparse
import gzip
import logging
import os
import sys

from scipy.io import mmwrite

import mtx_io

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

    matrix, barcodes, features = mtx_io.read_triplet(args.mtx, args.barcodes, args.features)

    logger.info(f"Writing {len(barcodes)} barcodes and {len(features)} features")

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

    # Straight into the gzip stream: mmwrite takes a file object, so the matrix does
    # not need writing out in full and reading back to be compressed.
    with gzip.open(os.path.join(args.outdir, "matrix.mtx.gz"), "wb") as f:
        mmwrite(f, matrix)

    logger.info(f"10x-format export written to {args.outdir}")


if __name__ == "__main__":
    main()
