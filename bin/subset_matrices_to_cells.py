#!/usr/bin/env python3
"""
Subset one or more matrices in a directory to a given list of cell barcodes.

Written for STARsolo's ``Velocyto/raw`` directory, which holds ``spliced.mtx``,
``unspliced.mtx`` and ``ambiguous.mtx`` over a single shared ``barcodes.tsv`` /
``features.tsv`` pair. The cell set is taken from a barcode list rather than
recomputed from a UMI cutoff, so that the velocity matrices describe exactly the
cells that were called on the main ``GeneFull_Ex50pAS`` matrix -- see
``secondderiv_filter_matrices.py``, which produces that call and that barcode
list. Recomputing a cutoff here instead would give a different cell set for
every feature, and the published matrices would no longer be comparable.

Barcodes are emitted in the order they appear in the barcode list, so every
output matrix shares one column order.
"""

import argparse
import glob
import os
import sys
from typing import List

import pandas as pd
import scipy.io as sio


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Subset the matrices in a directory to a list of cell barcodes."
    )
    parser.add_argument("-d", "--dir", required=True, help="Directory holding the matrices, barcodes.tsv and features.tsv")
    parser.add_argument("-b", "--barcodes", required=True, help="One barcode per line; the cells to keep")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for the subset matrices")
    parser.add_argument(
        "-m", "--matrices", nargs="+", default=None,
        help="Matrix file names to subset (default: every *.mtx in --dir)",
    )
    return parser.parse_args()


def resolve_features(dirpath: str) -> str:
    """Return the features file, tolerating STARsolo's older ``genes.tsv`` name."""
    features_path = os.path.join(dirpath, "features.tsv")
    if not os.path.exists(features_path):
        features_path = os.path.join(dirpath, "genes.tsv")
    return features_path


def resolve_matrices(dirpath: str, requested: List[str]) -> List[str]:
    """Return the matrix file names to subset, erroring when none can be found."""
    if requested:
        missing = [name for name in requested if not os.path.exists(os.path.join(dirpath, name))]
        if missing:
            raise SystemExit(f"Error: {', '.join(missing)} not found in {dirpath}")
        return requested

    found = sorted(os.path.basename(p) for p in glob.glob(os.path.join(dirpath, "*.mtx")))
    if not found:
        raise SystemExit(f"Error: no .mtx files found in {dirpath}")
    return found


def main() -> None:
    args = parse_args()

    barcodes_path = os.path.join(args.dir, "barcodes.tsv")
    features_path = resolve_features(args.dir)

    for path in (barcodes_path, features_path):
        if not os.path.exists(path):
            raise SystemExit(f"Error: {path} not found; is {args.dir} a matrix directory?")

    # ------ Resolve the cells to keep ------

    source_barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
    keep_barcodes = pd.read_csv(args.barcodes, header=None, sep="\t")

    # Position of each barcode in the source matrix, so the requested order is the output order
    position = {bc: i for i, bc in enumerate(source_barcodes.iloc[:, 0].astype(str))}

    wanted = keep_barcodes.iloc[:, 0].astype(str).tolist()
    keep_idx = [position[bc] for bc in wanted if bc in position]

    # A cell called on another feature's matrix need not appear here: the barcode axis is the
    # whitelist in both cases, but a barcode with no counts at all can be absent. Report it
    # rather than fail, so one odd barcode does not sink the run.
    absent = len(wanted) - len(keep_idx)
    if absent:
        print(
            f"Warning: {absent} of {len(wanted)} requested barcodes are not present in "
            f"{barcodes_path}; they are dropped from the subset.",
            file=sys.stderr,
        )

    if not keep_idx:
        raise SystemExit(
            f"Error: none of the {len(wanted)} barcodes in {args.barcodes} are present in {barcodes_path}."
        )

    kept_barcodes = source_barcodes.iloc[keep_idx]
    features = pd.read_csv(features_path, header=None, sep="\t")

    # ------ Subset and write ------

    os.makedirs(args.outdir, exist_ok=True)

    for name in resolve_matrices(args.dir, args.matrices):
        # CSC (Compressed Sparse Column) for efficient column/cell slicing
        mat = sio.mmread(os.path.join(args.dir, name)).tocsc()

        if mat.shape[1] != len(source_barcodes):
            raise SystemExit(
                f"Error: {name} has {mat.shape[1]} columns but {barcodes_path} has "
                f"{len(source_barcodes)} barcodes."
            )

        sio.mmwrite(os.path.join(args.outdir, name), mat[:, keep_idx])
        print(f"{name}: {mat.shape[0]} x {mat.shape[1]} -> {mat.shape[0]} x {len(keep_idx)}")

    kept_barcodes.to_csv(os.path.join(args.outdir, "barcodes.tsv"), sep="\t", header=False, index=False)
    features.to_csv(os.path.join(args.outdir, "features.tsv"), sep="\t", header=False, index=False)

    print(f"Subset {len(keep_idx)} cells to: {args.outdir}")


if __name__ == "__main__":
    main()
