#!/usr/bin/env python3
"""
Group one sample's called cells into Metacell2 metacells, and summarise every group
for the metacell filtering report.

Input is the cell-called h5ad MTX_TO_H5AD published: raw counts in X and, where
detection ran, the consensus calls in ``obs['doublet_status']``. Doublets are
annotated, not removed, so that they gather into metacells of their own and can be
excluded there, by the user, with the rest of the metacell in view. CellSweep only
ever runs on the unfiltered matrix, so its per-cell ``alpha_hat``, per-gene
``ambient_hat`` and denoised layer are lifted from the raw (or alevin-fry 'full')
h5ad for the called barcodes.

Metacell2 groups on a cleaned copy (metacells 0.9):

    exclude_genes       mitochondrial genes (from the gene table) plus
                        --excluded_gene_patterns, so they cannot drive the grouping
    exclude_cells       cells below --min_cell_umis or above --max_cell_umis. The
                        excluded-gene fraction cap is off by default: with it, the
                        cells with most mitochondrial reads would be dropped before
                        grouping, and the report's per-metacell mito percentage --
                        the thing it filters on -- would be censored.
    mark_lateral_genes  --lateral_gene_patterns (e.g. cell cycle), optional
    divide_and_conquer_pipeline, collect_metacells, compute_umap_by_markers

All QC numbers are then computed on the full matrix by ``metacell_qc_summary.py``.

Outputs, under ``--prefix``:

    _mc2_cells.h5ad       every called cell: X raw counts, obs['metacell'] plus the
                          per-cell QC, the doublet and CellSweep annotations,
                          layers['cellsweep'] when present, uns['mc_fingerprint']
    _mc2_metacells.h5ad   one row per group: X summed raw UMIs over all genes, obs
                          the per-group QC and UMAP coordinates
    _mc_summary.json      the report's payload for this sample

A sample with too few cells to group writes only the summary, with status 'skipped'
and the reason, and exits 0: an expected outcome must not fail the task, or Nextflow
would never cache it and every -resume would try again.

``--metacell_obs_col`` takes an existing assignment from obs instead of running
Metacell2, for metacells computed elsewhere and for the tests.
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

import metacell_qc_summary as qcs

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("RunMetacells")

CS_LAYER = "cellsweep"
CS_OBS_COLS = [qcs.ALPHA_COL, "is_empty", "celltype"]

# Below this many metacells a UMAP has fewer points than neighbours to place them by
MIN_METACELLS_FOR_UMAP = 15


class TooFewCells(Exception):
    pass


# --------------------------------------------------------------------------
# Ambient-RNA annotations
# --------------------------------------------------------------------------

def _anndata_io():
    try:
        from anndata.io import read_elem, sparse_dataset
    except ImportError:
        from anndata.experimental import read_elem, sparse_dataset
    return read_elem, sparse_dataset


def lift_ambient(adata, ambient_h5ad):
    """
    Copy CellSweep's annotations for ``adata``'s cells and genes out of the raw h5ad.

    The raw object holds every droplet, so it is read through h5py: obs and var in
    full, and only the called barcodes' rows of the denoised layer. Alignment is by
    name on both axes.
    """
    import h5py

    read_elem, sparse_dataset = _anndata_io()

    with h5py.File(ambient_h5ad, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])

        rows = obs.index.get_indexer(adata.obs_names)
        missing = int((rows < 0).sum())
        if missing:
            logger.warning(f"{missing} / {adata.n_obs} called barcodes are absent from {ambient_h5ad}")

        for col in CS_OBS_COLS:
            if col in obs.columns:
                adata.obs[col] = obs[col].reindex(adata.obs_names).to_numpy()
        if qcs.AMBIENT_VAR_COL in var.columns:
            adata.var[qcs.AMBIENT_VAR_COL] = var[qcs.AMBIENT_VAR_COL].reindex(adata.var_names).to_numpy()

        if "layers" in f and CS_LAYER in f["layers"]:
            adata.layers[CS_LAYER] = _rows_of_layer(
                sparse_dataset(f["layers"][CS_LAYER]), rows, var.index, adata.var_names
            )
            logger.info(f"Lifted CellSweep's denoised counts for {adata.n_obs - missing} cells")

    present = [c for c in CS_OBS_COLS if c in adata.obs.columns]
    logger.info(f"Lifted CellSweep annotations: obs {present}, "
                f"var {[qcs.AMBIENT_VAR_COL] if qcs.AMBIENT_VAR_COL in adata.var.columns else []}")
    return adata


def _rows_of_layer(dataset, rows, source_var, target_var):
    """``rows`` of a backed sparse layer, on ``target_var``'s gene axis; -1 rows stay zero."""
    present = np.flatnonzero(rows >= 0)
    src = rows[present]
    order = np.argsort(src)
    block = sp.csr_matrix(dataset[src[order]])

    # back to the called cells' order
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    block = block[inverse]

    cols = source_var.get_indexer(target_var)
    if not (cols >= 0).all() or not np.array_equal(cols, np.arange(len(target_var))):
        keep = np.flatnonzero(cols >= 0)
        block = sp.csr_matrix(block[:, cols[keep]])
        coo = block.tocoo()
        block = sp.csr_matrix((coo.data, (coo.row, keep[coo.col])), shape=(block.shape[0], len(target_var)))

    coo = block.tocoo()
    return sp.csr_matrix((coo.data, (present[coo.row], coo.col)), shape=(len(rows), len(target_var)))


# --------------------------------------------------------------------------
# Metacell2
# --------------------------------------------------------------------------

def _patterns(value):
    return [p for p in (value or "").replace(",", " ").split() if p]


def run_mc2(adata, info, args):
    """
    Per-cell group names (see metacell_qc_summary), the metacells' UMAP or None, and
    a few numbers worth recording.
    """
    import anndata as ad
    import metacells as mc

    if hasattr(mc.ut, "set_processors_count"):
        mc.ut.set_processors_count(args.cpus)

    X = qcs.as_csr(adata.X).astype(np.float32)
    if X.nnz and not np.array_equal(X.data, np.rint(X.data)):
        # alevin-fry's EM splits UMIs into fractions; Metacell2 downsamples whole UMIs.
        # Rounded for the grouping only -- every number reported is on the original.
        logger.info("Counts are fractional; rounding them for Metacell2's grouping only")
        X.data = np.rint(X.data)
        X.eliminate_zeros()

    full = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )
    mc.ut.set_name(full, args.sample_id)

    excluded_names = list(info.index[info["is_mito"].to_numpy()])
    mc.pl.exclude_genes(
        full,
        excluded_gene_names=excluded_names or None,
        excluded_gene_patterns=_patterns(args.excluded_gene_patterns) or None,
        random_seed=args.random_seed,
    )
    mc.tl.compute_excluded_gene_umis(full)
    # None disables a bound (metacells 0.9 takes all three as Optional)
    mc.pl.exclude_cells(
        full,
        properly_sampled_min_cell_total=int(args.min_cell_umis) if args.min_cell_umis else None,
        properly_sampled_max_cell_total=int(args.max_cell_umis) if args.max_cell_umis else None,
        properly_sampled_max_excluded_genes_fraction=args.max_excluded_genes_fraction,
    )

    clean = mc.pl.extract_clean_data(full, name=f"{args.sample_id}.clean")
    n_clean = 0 if clean is None else clean.n_obs
    min_cells = max(2 * args.target_metacell_size, args.min_cells)
    if n_clean < min_cells:
        raise TooFewCells(
            f"{n_clean} cells passed Metacell2's cell filters, fewer than the {min_cells} "
            f"needed (max(2 x target_metacell_size, min_cells))"
        )

    # Always set, even empty: the pile selection reads the lateral_gene mask unconditionally
    lateral = _patterns(args.lateral_gene_patterns)
    mc.pl.mark_lateral_genes(clean, lateral_gene_names=[], lateral_gene_patterns=lateral)

    if hasattr(mc.pl, "guess_max_parallel_piles"):
        mc.pl.set_max_parallel_piles(mc.pl.guess_max_parallel_piles(clean))

    mc.pl.divide_and_conquer_pipeline(
        clean, target_metacell_size=args.target_metacell_size, random_seed=args.random_seed
    )
    metacells = mc.pl.collect_metacells(clean, name=f"{args.sample_id}.metacells", random_seed=args.random_seed)

    index = clean.obs["metacell"].to_numpy().astype(np.int64)
    names = np.asarray(metacells.obs_names, dtype=object)
    clean_groups = np.where(index >= 0, names[np.clip(index, 0, None)], qcs.PSEUDO_OUTLIERS)

    groups = pd.Series(qcs.PSEUDO_EXCLUDED, index=adata.obs_names, dtype=object)
    groups.loc[clean.obs_names] = clean_groups

    umap = None
    if metacells.n_obs >= MIN_METACELLS_FOR_UMAP:
        try:
            # The UMAP is built on the metacells' marker genes, which nothing upstream marks
            # (compute_for_mcview would, along with much the report does not need)
            mc.tl.find_metacells_marker_genes(metacells)
            mc.pl.compute_umap_by_markers(metacells, random_seed=args.random_seed)
            umap = pd.DataFrame(
                {"x": metacells.obs["x"].to_numpy(), "y": metacells.obs["y"].to_numpy()},
                index=pd.Index(names.astype(str)),
            )
        except Exception as exc:  # the report falls back to a QC scatter
            logger.warning(f"Metacell2's UMAP failed ({exc}); the report will plot QC metrics instead")
    else:
        logger.info(f"{metacells.n_obs} metacells; too few for a UMAP")

    stats = {
        "n_clean_cells": int(n_clean),
        "n_excluded_cells": int(adata.n_obs - n_clean),
        "n_outlier_cells": int((index < 0).sum()),
        "n_excluded_genes": int(np.asarray(full.var["excluded_gene"]).sum()) if "excluded_gene" in full.var else None,
        "metacells_version": package_version("metacells"),
    }
    return groups.to_numpy(), umap, stats


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

def _uns_safe(d):
    """anndata cannot write None into uns."""
    return {k: ("" if v is None else v) for k, v in d.items()}


def annotate_cells(adata, info, groups, umap, fp, mode, params):
    qc = qcs.per_cell_qc(adata.X, info["is_mito"].to_numpy(), info["is_rrna"].to_numpy())
    adata.obs[qcs.METACELL_COL] = pd.Categorical(groups, categories=qcs.group_order(groups))
    for col in qc.columns:
        adata.obs[col] = qc[col].to_numpy()

    adata.var["gene_name"] = info["gene_name"].to_numpy()
    adata.var["is_mito"] = info["is_mito"].to_numpy()
    adata.var["is_rrna"] = info["is_rrna"].to_numpy()
    adata.var["pfam"] = [",".join(d) for d in info["pfam"]]

    adata.uns["mc_fingerprint"] = fp
    adata.uns["mc_doublet_mode"] = mode
    adata.uns["mc_params"] = _uns_safe(params)
    if umap is not None:
        adata.uns["mc_umap"] = {
            "id": umap.index.to_numpy().astype(str),
            "x": umap["x"].to_numpy(dtype=float),
            "y": umap["y"].to_numpy(dtype=float),
        }
    return adata


def metacell_object(adata, records, groups):
    """One row per group: summed raw UMIs over every gene, the group QC in obs."""
    import anndata as ad

    order = qcs.group_order(groups)
    sums = qcs.group_indicator(groups, order) @ qcs.as_csr(adata.X)

    obs = pd.DataFrame(records).set_index("id").reindex(order)
    obs["markers"] = obs["markers"].map(lambda m: ",".join(m))
    for col in obs.columns:
        if obs[col].dtype == object:
            obs[col] = obs[col].map(lambda v: "" if v is None else v)
            obs[col] = obs[col].astype(str) if col == "markers" else pd.to_numeric(obs[col], errors="coerce")

    return ad.AnnData(X=sp.csr_matrix(sums), obs=obs, var=adata.var.copy())


def package_version(name):
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def write_versions(path, process):
    import platform

    with open(path, "w") as f:
        f.write(f'"{process}":\n')
        f.write(f"    python: {platform.python_version()}\n")
        f.write(f"    anndata: {package_version('anndata')}\n")
        f.write(f"    metacells: {package_version('metacells')}\n")


def main():
    parser = argparse.ArgumentParser(description="Group one sample's called cells into Metacell2 metacells.")
    parser.add_argument("--cells_h5ad", required=True, help="Cell-called h5ad (MTX_TO_H5AD)")
    parser.add_argument("--ambient_h5ad", default=None,
                        help="Optional raw/full h5ad carrying CellSweep's annotations for these barcodes")
    parser.add_argument("--gene_table", required=True, help="build_gene_table.py output")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--mapping_method", default="")
    parser.add_argument("--prefix", required=True, help="Output prefix, e.g. the sample id")
    parser.add_argument("--target_metacell_size", type=int, default=96)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--min_cell_umis", type=float, default=200)
    parser.add_argument("--max_cell_umis", type=float, default=None)
    parser.add_argument("--max_excluded_genes_fraction", type=float, default=None,
                        help="Metacell2's cap on excluded-gene UMIs per cell; off by default (see above)")
    parser.add_argument("--excluded_gene_patterns", default="", help="Extra genes kept out of the grouping")
    parser.add_argument("--lateral_gene_patterns", default="", help="Genes marked lateral (e.g. cell cycle)")
    parser.add_argument("--min_cells", type=int, default=300,
                        help="Fewest clean cells to group; with fewer the sample is skipped")
    parser.add_argument("--doublets_removed", action="store_true",
                        help="The matrix was doublet-filtered (params.perform_doublet_filtering)")
    parser.add_argument("--metacell_obs_col", default=None,
                        help="Take the assignment from this obs column instead of running Metacell2")
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--compression", default="gzip")
    parser.add_argument("--versions_yml", default=None)
    parser.add_argument("--process_name", default="METACELL2")
    args = parser.parse_args()

    os.environ.setdefault("METACELLS_PROCESSORS_COUNT", str(args.cpus))

    import anndata as ad

    params = {
        "target_metacell_size": args.target_metacell_size,
        "random_seed": args.random_seed,
        "min_cell_umis": args.min_cell_umis,
        "max_cell_umis": args.max_cell_umis,
        "max_excluded_genes_fraction": args.max_excluded_genes_fraction,
        "excluded_gene_patterns": args.excluded_gene_patterns,
        "lateral_gene_patterns": args.lateral_gene_patterns,
        "source": f"obs['{args.metacell_obs_col}']" if args.metacell_obs_col else "metacell2",
    }
    summary_path = f"{args.prefix}_mc_summary.json"

    adata = ad.read_h5ad(args.cells_h5ad)
    adata.var_names_make_unique()
    adata.X = qcs.as_csr(adata.X)
    logger.info(f"Read {adata.n_obs} called cells x {adata.n_vars} genes from {args.cells_h5ad}")

    if args.ambient_h5ad:
        lift_ambient(adata, args.ambient_h5ad)

    info, how = qcs.gene_info(adata.var_names, qcs.load_gene_table(args.gene_table))
    logger.info(
        f"Gene table matched {how['match_rate']:.1%} of features "
        f"(gene_id {how['gene_id']}, gene_name {how['gene_name']}, cut id {how['cut_id']}, "
        f"unmatched {how['unmatched']}); {int(info['is_mito'].sum())} mitochondrial"
    )
    if how["match_rate"] < 0.5:
        logger.warning("Fewer than half the features matched the gene table: mito flags and domains "
                       "will be missing for the rest. Is this the GTF the sample was counted against?")

    umap = None
    try:
        if args.metacell_obs_col:
            if args.metacell_obs_col not in adata.obs.columns:
                raise SystemExit(f"--metacell_obs_col '{args.metacell_obs_col}' is not in obs")
            groups = adata.obs[args.metacell_obs_col].astype(str).to_numpy()
            if "mc_umap" in adata.uns:
                umap = pd.DataFrame(adata.uns["mc_umap"]).set_index("id")
        else:
            if adata.n_obs < args.min_cells:
                raise TooFewCells(f"only {adata.n_obs} cells were called, fewer than min_cells = {args.min_cells}")
            groups, umap, stats = run_mc2(adata, info, args)
            params.update(stats)
    except TooFewCells as exc:
        reason = str(exc)
        logger.warning(f"No metacells for {args.sample_id}: {reason}")
        qcs.write_json(summary_path, qcs.skipped_summary(
            args.sample_id, args.mapping_method, reason, n_cells=int(adata.n_obs), params=params))
        if args.versions_yml:
            write_versions(args.versions_yml, args.process_name)
        return

    mode = qcs.doublet_mode(adata.obs, removed=args.doublets_removed)
    fp = qcs.fingerprint(adata.obs_names, groups, adata.var_names)
    annotate_cells(adata, info, groups, umap, fp, mode, params)

    summary = qcs.build_summary(adata, info, how, args.sample_id, args.mapping_method,
                                removed=args.doublets_removed, umap=umap, params=params)
    if summary["sample"]["fingerprint"] != fp:
        raise RuntimeError("The fingerprint changed while annotating the cells; refusing to write a summary "
                           "a selection could not be applied against")

    compression = None if args.compression.lower() in ("", "none", "null") else args.compression
    adata.write_h5ad(f"{args.prefix}_mc2_cells.h5ad", compression=compression)
    metacell_object(adata, summary["metacells"], groups).write_h5ad(
        f"{args.prefix}_mc2_metacells.h5ad", compression=compression)
    qcs.write_json(summary_path, summary)

    logger.info(
        f"{summary['sample']['n_metacells']} metacells for {args.sample_id}; "
        f"{int((groups == qcs.PSEUDO_OUTLIERS).sum())} outlier and "
        f"{int((groups == qcs.PSEUDO_EXCLUDED).sum())} excluded cells"
    )

    if args.versions_yml:
        write_versions(args.versions_yml, args.process_name)


if __name__ == "__main__":
    main()
