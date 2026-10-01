#!/usr/bin/env python3
"""
Packs a count-matrix triplet into the published .h5ad, merging in whatever the
detections upstream produced.

This is the last step of the filtering chain, not the first: doublet detection,
optional doublet filtering and ambient-RNA denoising all work on the triplet, and
their results converge here so that one object per matrix carries all of them
instead of a succession of partial h5ads.

    X                                    counts from the triplet, doublet-filtered
                                         when params.perform_doublet_filtering ran
    obs['sample_id']                     always
    obs['doublet_status']                --doublet_results, the Combine_Results.R
                                         consensus (singlet / doublet)
    layers['cellsweep']                  --cellsweep_h5ad, CellSweep's denoised counts
    obs['is_empty' / 'celltype' /        --cellsweep_h5ad
        'alpha_hat']
    var['ambient_hat']                   --cellsweep_h5ad

Both annotation inputs are optional; without either, this is the plain triplet ->
h5ad conversion it has always been.

CellSweep reads the same triplet, so its axes normally match this object's exactly
and its matrix is taken as it stands. Nothing in its API promises they will, so the
general case is merged by name -- never by position, since a positional merge that
silently slipped would produce an object of exactly the right shape carrying the
wrong numbers.
"""

import argparse
import logging
import sys

import scipy.sparse as sp

import doublet_calls
import mtx_io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("MtxToH5ad")

# CellSweep's per-cell and per-gene results, as they are named in its own output
CS_OBS_COLS = ["is_empty", "celltype", "alpha_hat"]
CS_VAR_COLS = ["ambient_hat"]


def denoised_layer(adata, cs):
    """
    CellSweep's denoised counts, on `adata`'s axes.

    CellSweep reads the same triplet this object was built from, so its axes
    normally match exactly and its matrix can be taken as it stands. Nothing in
    CellSweep's API promises that, though, so the general case is still handled --
    by name, never by position, since a positional merge that slipped would give a
    layer of exactly the right shape carrying the wrong cells' counts.
    """
    if cs.obs_names.equals(adata.obs_names) and cs.var_names.equals(adata.var_names):
        logger.info("CellSweep's axes match the count matrix; taking its counts as they are")
        return cs.X

    # Scatter onto this object's axes. Genes CellSweep dropped stay zero in the
    # layer, which a sparse matrix stores for free.
    sub = cs[:, cs.var_names.isin(adata.var_names)]
    coo = sp.coo_matrix(sub.X)
    rows = adata.obs_names.get_indexer(sub.obs_names)[coo.row]
    cols = adata.var_names.get_indexer(sub.var_names)[coo.col]

    logger.info(
        f"Denoised counts placed for {sub.n_obs} / {adata.n_obs} barcodes and "
        f"{sub.n_vars} / {adata.n_vars} genes"
    )

    return sp.coo_matrix(
        (coo.data, (rows, cols)), shape=adata.shape
    ).tocsr()


def merge_cellsweep(adata, cellsweep_h5ad):
    """
    Lift CellSweep's denoised counts into a layer and its annotations into obs/var,
    aligned to `adata`'s axes by name.
    """
    import anndata as ad

    logger.info(f"Loading CellSweep results from {cellsweep_h5ad}")
    cs = ad.read_h5ad(cellsweep_h5ad)
    cs.var_names_make_unique()

    unknown = cs.obs_names.difference(adata.obs_names)
    if len(unknown):
        raise ValueError(
            f"{len(unknown)} of CellSweep's {cs.n_obs} barcodes are absent from the count "
            f"matrix (e.g. {list(unknown[:3])}); the two are not describing the same sample."
        )

    adata.layers["cellsweep"] = denoised_layer(adata, cs)

    for col in CS_OBS_COLS:
        if col in cs.obs:
            adata.obs[col] = cs.obs[col].reindex(adata.obs_names)
        else:
            logger.warning(f"CellSweep output carries no obs['{col}']; skipped")

    for col in CS_VAR_COLS:
        if col in cs.var:
            adata.var[col] = cs.var[col].reindex(adata.var_names)
        else:
            logger.warning(f"CellSweep output carries no var['{col}']; skipped")

    return adata


def write_versions(path, process):
    import platform

    import anndata as ad
    import scanpy as sc

    with open(path, "w") as f:
        f.write(f'"{process}":\n')
        f.write(f"    python: {platform.python_version()}\n")
        f.write(f"    anndata: {ad.__version__}\n")
        f.write(f"    scanpy: {sc.__version__}\n")


def main():
    parser = argparse.ArgumentParser(description="Pack a count-matrix triplet into an annotated h5ad.")
    parser.add_argument("--mtx", type=str, required=True, help="Count matrix (mtx)")
    parser.add_argument("--barcodes", type=str, required=True, help="Barcode axis of --mtx")
    parser.add_argument("--features", type=str, required=True, help="Feature axis of --mtx")
    parser.add_argument("--sample_id", type=str, required=True, help="Written to .obs['sample_id']")
    parser.add_argument("--output_h5ad", type=str, required=True, help="Output path for the h5ad")
    parser.add_argument("--doublet_results", type=str, default=None,
                        help="Optional Combine_Results.R '_w_combined_assignments.tsv' output; "
                             "its consensus calls become obs['doublet_status']")
    parser.add_argument("--doublet_method", type=str, default=None,
                        help="Consensus method name used in Combine_Results.R (e.g. AnySinglet). "
                             "Required with --doublet_results")
    parser.add_argument("--cellsweep_h5ad", type=str, default=None,
                        help="Optional CellSweep '_cs_full.h5ad'; its denoised counts and "
                             "ambient-RNA annotations are merged in")
    parser.add_argument("--compression", type=str, default="gzip",
                        help="HDF5 compression for the output ('gzip', 'lzf', or 'none'). "
                             "gzip is the interoperable default; lzf writes several times "
                             "faster for ~15%% more disk but is only readable through h5py")
    parser.add_argument("--versions_yml", type=str, default=None, help="Optional versions.yml path")
    parser.add_argument("--process_name", type=str, default="MTX_TO_H5AD", help="Name for versions.yml")
    args = parser.parse_args()

    if args.doublet_results and not args.doublet_method:
        parser.error("--doublet_method is required when --doublet_results is given")

    adata = mtx_io.read_triplet_anndata(
        args.mtx, args.barcodes, args.features, sample_id=args.sample_id
    )
    logger.info(f"Read {adata.n_obs} barcodes x {adata.n_vars} genes")

    if args.doublet_results:
        logger.info(f"Loading consensus doublet calls from {args.doublet_results}")
        barcodes = doublet_calls.read_doublet_barcodes(args.doublet_results, args.doublet_method)
        adata = doublet_calls.annotate(adata, barcodes)

    if args.cellsweep_h5ad:
        adata = merge_cellsweep(adata, args.cellsweep_h5ad)

    compression = None if args.compression.lower() in ("", "none", "null") else args.compression
    adata.write_h5ad(args.output_h5ad, compression=compression)
    logger.info(f"Wrote {args.output_h5ad} (compression: {compression or 'none'})")

    if args.versions_yml:
        write_versions(args.versions_yml, args.process_name)


if __name__ == "__main__":
    main()
