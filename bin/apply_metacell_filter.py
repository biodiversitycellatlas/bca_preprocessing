#!/usr/bin/env python3
"""
Apply a selection made in ``filtering_report.html`` to one sample's metacells, and
write the filtered UMI matrices for downstream analysis.

The selection is the JSON the report's export button downloads. Its per-sample
``keep_metacells`` and ``excluded_genes`` are authoritative: they were resolved in
the report from the thresholds and the blacklist the user set, and
are applied here as they stand. The rules are carried along as a record of how
the lists were made; of them, only the sample's ``gene_rules.min_total_umis`` is
read, to drop the genes with no UMI that the report could not list.

A selection is only meaningful for the assignment it was made on: metacell names
are reused between runs, so the same name applied to a regrouped sample would keep
a different set of cells and nothing would look wrong. Each sample in the selection
therefore carries the fingerprint of the assignment it was made against
(``metacell_qc_summary.fingerprint``), and it is checked against the one recomputed
from the cells h5ad here before anything is written.

Outputs, under ``--outdir``:

    matrix.mtx / barcodes.tsv / features.tsv   kept cells x genes, raw UMIs, in the
                                               10x orientation (genes x cells)
    <sample>_final.h5ad                        the same cells, with obs['metacell'],
                                               the per-cell QC and annotations,
                                               layers['cellsweep'] when present
    <sample>_final_metacells.h5ad              kept metacells x genes, summed raw UMIs
    <sample>_filter_summary.json               what was kept, and from what

``--gene_mode drop`` (default) removes the excluded genes from every output;
``--gene_mode flag`` keeps the full gene axis -- convenient when samples will be
integrated later -- and records the decision in ``var['pass_filter']`` and
``var['exclusion_reasons']`` instead.

Runs inside the pipeline (APPLY_METACELL_FILTER, with -resume and
--metacell_selection) and standalone on the published ``metacells/<id>/`` outputs:

    apply_metacell_filter.py --cells_h5ad S1_starsolo_mc2_cells.h5ad \\
        --selection selection.json --sample_id S1_starsolo --outdir S1_starsolo_final
"""

import argparse
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

import metacell_qc_summary as qcs
import mtx_io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("ApplyMetacellFilter")

SELECTION_SCHEMA = "bca_metacell_selection"
# Version 2 gives each sample its own gene_rules; version 1 had one set for the run
SELECTION_VERSIONS = (1, 2)

EXIT_INVALID_SELECTION = 2


class SelectionError(Exception):
    pass


def load_selection(path):
    with open(path) as fh:
        try:
            selection = json.load(fh)
        except json.JSONDecodeError as exc:
            raise SelectionError(f"{path} is not valid JSON: {exc}")

    if selection.get("schema") != SELECTION_SCHEMA:
        raise SelectionError(
            f"{path} is not a metacell selection (schema '{selection.get('schema')}', "
            f"expected '{SELECTION_SCHEMA}'). Export it from filtering_report.html."
        )
    if selection.get("schema_version") not in SELECTION_VERSIONS:
        raise SelectionError(
            f"{path} has schema_version {selection.get('schema_version')}; this pipeline "
            f"reads {list(SELECTION_VERSIONS)}."
        )
    if not isinstance(selection.get("samples"), dict):
        raise SelectionError(f"{path} carries no 'samples' object")
    return selection


def sample_entry(selection, sample_id):
    samples = selection["samples"]
    if sample_id not in samples:
        raise SelectionError(
            f"Sample '{sample_id}' is not in the selection, which covers: {sorted(samples)}. "
            "Samples left out of a selection are not filtered."
        )
    entry = samples[sample_id]
    for key in ("fingerprint", "keep_metacells", "excluded_genes"):
        if key not in entry:
            raise SelectionError(f"The selection for '{sample_id}' has no '{key}'")
    return entry


def check_fingerprint(adata, entry, sample_id):
    """
    The selection must have been made on this exact assignment, and the h5ad must
    still hold the assignment it says it does.
    """
    recomputed = qcs.fingerprint(adata.obs_names, adata.obs[qcs.METACELL_COL], adata.var_names)
    stored = adata.uns.get("mc_fingerprint")

    if stored and stored != recomputed:
        raise SelectionError(
            f"The cells h5ad for '{sample_id}' no longer matches the fingerprint it was written "
            "with; it has been modified since Metacell2 ran."
        )
    if entry["fingerprint"] != recomputed:
        raise SelectionError(
            f"The selection for '{sample_id}' was made on a different metacell assignment "
            f"(selection {entry['fingerprint'][:19]}..., this run {recomputed[:19]}...). Metacell "
            "names are reused between runs, so applying it would keep the wrong cells. Make a new "
            "selection in this run's filtering_report.html."
        )
    return recomputed


def gene_rules_for(selection, entry):
    """The sample's gene rules: its own (version 2), else the run's (version 1)."""
    return entry.get("gene_rules") or selection.get("gene_rules")


def resolve(adata, entry, sample_id, gene_rules=None):
    """
    Boolean masks of the cells and genes to keep, after checking every name exists.

    Genes with no UMI in any called cell are not in the report's payload, so they
    cannot be listed; a minimum-UMI rule excludes them by definition, and they are
    added here with the same reason.
    """
    groups = adata.obs[qcs.METACELL_COL].astype(str)
    known = set(groups)

    keep = set(map(str, entry["keep_metacells"]))
    unknown = sorted(keep - known)
    if unknown:
        raise SelectionError(
            f"The selection for '{sample_id}' keeps {len(unknown)} metacell(s) this sample does not "
            f"have, e.g. {unknown[:5]}"
        )

    excluded = {str(g): reasons for g, reasons in entry["excluded_genes"].items()}
    unknown = sorted(set(excluded) - set(adata.var_names))
    if unknown:
        raise SelectionError(
            f"The selection for '{sample_id}' excludes {len(unknown)} gene(s) this sample does not "
            f"have, e.g. {unknown[:5]}"
        )

    min_umis = (gene_rules or {}).get("min_total_umis") or 0
    if min_umis > 0:
        totals = np.asarray(qcs.as_csr(adata.X).sum(axis=0)).ravel()
        for g in adata.var_names[totals == 0]:
            excluded.setdefault(str(g), ["min_umi"])

    cell_mask = groups.isin(keep).to_numpy()
    gene_mask = ~adata.var_names.isin(list(excluded))
    return cell_mask, gene_mask, excluded


def features_frame(var):
    """A 10x-style features.tsv: id, name (the id where there is none), feature type."""
    names = var["gene_name"].astype(str) if "gene_name" in var.columns else pd.Series("", index=var.index)
    names = names.where(names != "", var.index.to_series())
    return pd.DataFrame({0: var.index.astype(str), 1: names.to_numpy(), 2: "Gene Expression"})


def metacell_sums(adata):
    import anndata as ad

    groups = adata.obs[qcs.METACELL_COL].astype(str).to_numpy()
    order = qcs.group_order(groups)
    sums = qcs.group_indicator(groups, order) @ qcs.as_csr(adata.X)
    obs = pd.DataFrame(index=pd.Index(order))
    obs["n_cells"] = pd.Series(groups).value_counts().reindex(order).to_numpy()
    return ad.AnnData(X=sp.csr_matrix(sums), obs=obs, var=adata.var.copy())


def write_versions(path, process):
    import platform
    from importlib.metadata import version

    with open(path, "w") as f:
        f.write(f'"{process}":\n')
        f.write(f"    python: {platform.python_version()}\n")
        f.write(f"    anndata: {version('anndata')}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Apply a filtering_report.html selection and write the filtered UMI matrices."
    )
    parser.add_argument("--cells_h5ad", required=True, help="The sample's _mc2_cells.h5ad")
    parser.add_argument("--selection", required=True, help="selection.json exported from the report")
    parser.add_argument("--sample_id", required=True, help="The sample's id in the selection")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--gene_mode", choices=["drop", "flag"], default="drop",
                        help="Remove excluded genes (drop) or keep them flagged in var (flag)")
    parser.add_argument("--compression", default="gzip")
    parser.add_argument("--versions_yml", default=None)
    parser.add_argument("--process_name", default="APPLY_METACELL_FILTER")
    args = parser.parse_args()

    import anndata as ad

    try:
        selection = load_selection(args.selection)
        entry = sample_entry(selection, args.sample_id)

        adata = ad.read_h5ad(args.cells_h5ad)
        if qcs.METACELL_COL not in adata.obs.columns:
            raise SelectionError(f"{args.cells_h5ad} has no obs['{qcs.METACELL_COL}']")

        fp = check_fingerprint(adata, entry, args.sample_id)
        cell_mask, gene_mask, excluded = resolve(adata, entry, args.sample_id, gene_rules_for(selection, entry))
    except SelectionError as exc:
        logger.error(str(exc))
        sys.exit(EXIT_INVALID_SELECTION)

    if not cell_mask.any():
        logger.error(f"The selection for '{args.sample_id}' keeps no cells")
        sys.exit(EXIT_INVALID_SELECTION)

    total_umis = float(qcs.as_csr(adata.X).sum())
    kept = adata[cell_mask].copy()

    reasons = pd.Series([";".join(excluded.get(g, [])) for g in kept.var_names], index=kept.var_names)
    if args.gene_mode == "drop":
        kept = kept[:, gene_mask].copy()
    else:
        kept.var["pass_filter"] = gene_mask
        kept.var["exclusion_reasons"] = reasons.to_numpy()

    if "pfam" in kept.var.columns:
        kept.var["pfam"] = kept.var["pfam"].astype(str)
    kept.obs[qcs.METACELL_COL] = kept.obs[qcs.METACELL_COL].astype(str).astype("category")
    kept.X = qcs.as_csr(kept.X)
    kept.uns["mc_selection"] = {
        "selection_created": str(selection.get("created", "")),
        "fingerprint": fp,
        "gene_mode": args.gene_mode,
    }

    os.makedirs(args.outdir, exist_ok=True)
    mtx_io.write_triplet(args.outdir, kept.X.T.tocsr(), list(kept.obs_names), features_frame(kept.var))

    compression = None if args.compression.lower() in ("", "none", "null") else args.compression
    kept.write_h5ad(os.path.join(args.outdir, f"{args.sample_id}_final.h5ad"), compression=compression)
    metacell_sums(kept).write_h5ad(
        os.path.join(args.outdir, f"{args.sample_id}_final_metacells.h5ad"), compression=compression)

    groups = adata.obs[qcs.METACELL_COL].astype(str)
    kept_umis = float(qcs.as_csr(kept.X).sum())
    summary = {
        "sample_id": args.sample_id,
        "fingerprint": fp,
        "selection_created": selection.get("created"),
        "gene_mode": args.gene_mode,
        "cells_in": int(adata.n_obs),
        "cells_kept": int(cell_mask.sum()),
        "groups_in": int(groups.nunique()),
        "groups_kept": int(groups[cell_mask].nunique()),
        "kept_groups": sorted(set(groups[cell_mask])),
        "genes_in": int(adata.n_vars),
        "genes_kept": int(gene_mask.sum()),
        "genes_excluded": int((~gene_mask).sum()),
        "umis_in": round(total_umis, 2),
        "umis_kept": round(kept_umis, 2),
        "umis_kept_pct": round(100.0 * kept_umis / total_umis, 3) if total_umis else 0.0,
        "cell_rules": entry.get("cell_rules"),
        "gene_rules": gene_rules_for(selection, entry),
    }
    with open(os.path.join(args.outdir, f"{args.sample_id}_filter_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    logger.info(
        f"{args.sample_id}: kept {summary['cells_kept']} / {summary['cells_in']} cells in "
        f"{summary['groups_kept']} / {summary['groups_in']} groups and {summary['genes_kept']} / "
        f"{summary['genes_in']} genes ({summary['umis_kept_pct']}% of UMIs)"
    )

    if args.versions_yml:
        write_versions(args.versions_yml, args.process_name)


if __name__ == "__main__":
    main()
