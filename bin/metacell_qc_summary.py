#!/usr/bin/env python3
"""
The per-metacell and per-gene numbers the metacell filtering report filters on,
computed from a cells h5ad whose ``obs['metacell']`` assigns every cell to a group.

Shared by ``run_metacells.py``, which calls it right after grouping so the h5ad is
never staged twice, and by ``apply_metacell_filter.py``, which has to recompute the
fingerprint the report's selection was made against. It also runs on its own, on any
cells h5ad that already carries an assignment, which is what the tests do.

Every number is computed on the full count matrix of the called cells, before the
exclusions Metacell2 applies to choose its grouping genes: a mitochondrial fraction
computed after the mitochondrial genes were excluded would read 0 everywhere.

Groups:

    <name>          a metacell, named as Metacell2 names it (e.g. "M12.34")
    __outliers__    cells Metacell2 grouped with no metacell (metacell -1)
    __excluded__    cells Metacell2 never grouped (too few UMIs, see exclude_cells)

The two pseudo-groups are reported like metacells, so the report can show what is in
them and the user can choose to keep them.

The summary JSON is the report's only input per sample:

    sample        id, mapping method, status, doublet mode, fingerprint, counts
    metacells     one record per group (see summarize_groups)
    genes         per-gene arrays aligned to one another (see summarize_genes)
"""

import argparse
import hashlib
import json
import logging
import math
import re
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger("MetacellQC")

PSEUDO_OUTLIERS = "__outliers__"
PSEUDO_EXCLUDED = "__excluded__"
PSEUDO_GROUPS = (PSEUDO_OUTLIERS, PSEUDO_EXCLUDED)

METACELL_COL = "metacell"
DOUBLET_COL = "doublet_status"
ALPHA_COL = "alpha_hat"
AMBIENT_VAR_COL = "ambient_hat"

SUMMARY_SCHEMA = "bca_metacell_summary"
SUMMARY_VERSION = 1

N_MARKERS = 5
MARKER_MIN_UMIS = 5
MARKER_EPS = 1e-5


# --------------------------------------------------------------------------
# Genes
# --------------------------------------------------------------------------

def load_gene_table(path):
    """``build_gene_table.py``'s output, with flags as booleans and domains as lists."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    table["is_mito"] = table["is_mito"].astype(int).astype(bool)
    table["is_rrna"] = table["is_rrna"].astype(int).astype(bool)
    table["pfam"] = table["pfam"].map(lambda s: [d for d in s.split(",") if d])
    return table


def _unique_index(keys):
    """key -> row for the keys that occur exactly once; an ambiguous key matches nothing."""
    counts = pd.Series(keys).value_counts()
    unique = set(counts.index[counts == 1])
    return {k: i for i, k in enumerate(keys) if k in unique and k != ""}


def match_genes(var_names, table):
    """
    The gene-table row each feature corresponds to (-1 when none), and how it matched.

    Tried in order, so a feature takes the most specific match available:

        gene_id     STARsolo's features.tsv column 0
        gene_name   a feature axis of names
        cut_id      the part of a GTF gene_id before its first '-', which is what
                    salmon_create_splici_ref.R reduces alevin-fry's gene ids to
    """
    by_id = {gid: i for i, gid in enumerate(table["gene_id"])}
    by_name = _unique_index(list(table["gene_name"]))
    by_cut = _unique_index([gid.split("-", 1)[0] for gid in table["gene_id"]])

    rows = np.full(len(var_names), -1, dtype=np.int64)
    how = {"gene_id": 0, "gene_name": 0, "cut_id": 0, "unmatched": 0}

    for j, name in enumerate(map(str, var_names)):
        for kind, index in (("gene_id", by_id), ("gene_name", by_name), ("cut_id", by_cut)):
            if name in index:
                rows[j] = index[name]
                how[kind] += 1
                break
        else:
            how["unmatched"] += 1

    return rows, how


def gene_info(var_names, table):
    """
    Per-feature annotation aligned to ``var_names``, plus the match statistics.

    Features with no gene-table row -- GeneExt's new genes, or a feature axis from a
    different annotation -- are kept, unflagged and without domains.
    """
    rows, how = match_genes(var_names, table)
    matched = rows >= 0
    take = np.where(matched, rows, 0)

    info = pd.DataFrame(index=pd.Index(list(map(str, var_names))))
    info["gene_name"] = np.where(matched, table["gene_name"].to_numpy()[take], "")
    info["is_mito"] = np.where(matched, table["is_mito"].to_numpy()[take], False).astype(bool)
    info["is_rrna"] = np.where(matched, table["is_rrna"].to_numpy()[take], False).astype(bool)
    pfam = table["pfam"].to_numpy()
    info["pfam"] = [pfam[r] if r >= 0 else [] for r in rows]

    n = len(var_names)
    how["match_rate"] = (n - how["unmatched"]) / n if n else 0.0
    return info, how


# --------------------------------------------------------------------------
# Cells
# --------------------------------------------------------------------------

def as_csr(X):
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))


def per_cell_qc(X, is_mito, is_rrna):
    """Total, mitochondrial and rRNA UMIs and detected genes per cell, on the full matrix."""
    X = as_csr(X)
    total = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    mito = np.asarray(X[:, np.flatnonzero(is_mito)].sum(axis=1)).ravel() if is_mito.any() else np.zeros_like(total)
    rrna = np.asarray(X[:, np.flatnonzero(is_rrna)].sum(axis=1)).ravel() if is_rrna.any() else np.zeros_like(total)
    n_genes = np.diff(X.indptr)

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mito = np.where(total > 0, 100.0 * mito / total, 0.0)
        pct_rrna = np.where(total > 0, 100.0 * rrna / total, 0.0)

    return pd.DataFrame({
        "total_umis": total,
        "mito_umis": mito,
        "rrna_umis": rrna,
        "n_genes": n_genes,
        "pct_mito": pct_mito,
        "pct_rrna": pct_rrna,
    })


def doublet_mode(obs, removed=False):
    """
    ``annotated`` when the calls are in obs, ``removed`` when the matrix was already
    doublet-filtered, ``absent`` when detection was off or a caller failed. The
    report shows a doublet percentage only in the first case: 0% would be a claim.
    """
    if removed:
        return "removed"
    return "annotated" if DOUBLET_COL in obs.columns else "absent"


def fingerprint(barcodes, groups, var_names):
    """
    A hash of the assignment and the gene axis, which a selection is bound to.

    Built from sorted ``barcode<TAB>group`` lines, so it does not depend on the
    order cells were written in, only on which cell went where.
    """
    h = hashlib.sha256()
    for bc, grp in sorted(zip(map(str, barcodes), map(str, groups))):
        h.update(f"{bc}\t{grp}\n".encode())
    h.update(b"#genes\n")
    for g in map(str, var_names):
        h.update(f"{g}\n".encode())
    return "sha256:" + h.hexdigest()


# --------------------------------------------------------------------------
# Groups
# --------------------------------------------------------------------------

def group_order(groups):
    """Real metacells in Metacell2's own order (natural sort on the name), pseudo-groups last."""
    def natural(name):
        return [(0, int(c), "") if c.isdigit() else (1, 0, c) for c in re.split(r"(\d+)", name)]

    names = pd.unique(pd.Series(groups, dtype=str))
    real = sorted((g for g in names if g not in PSEUDO_GROUPS), key=natural)
    return real + [g for g in PSEUDO_GROUPS if g in set(names)]


def group_indicator(groups, order):
    """A groups x cells 0/1 matrix, so group sums are one sparse product."""
    position = {g: i for i, g in enumerate(order)}
    rows = np.array([position[g] for g in map(str, groups)], dtype=np.int64)
    cols = np.arange(len(rows))
    return sp.csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(len(order), len(rows)))


def marker_genes(sums, names, usable, n_markers=N_MARKERS):
    """
    The genes most enriched in each metacell over the mean of all metacells.

    Enrichment is log2 of the gene's fraction in the metacell over its mean fraction,
    on genes with at least MARKER_MIN_UMIS in that metacell so a single UMI in a small
    metacell cannot top the list. Rows are densified in chunks, never all at once.
    """
    sums = as_csr(sums)
    totals = np.asarray(sums.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    fractions = sp.diags(1.0 / totals) @ sums
    mean_frac = np.asarray(fractions.mean(axis=0)).ravel()
    names = np.asarray(names)

    markers = []
    for start in range(0, sums.shape[0], 256):
        block_frac = fractions[start:start + 256].toarray()
        block_sums = sums[start:start + 256].toarray()
        score = np.log2((block_frac + MARKER_EPS) / (mean_frac + MARKER_EPS))
        score[:, ~usable] = -np.inf
        score[block_sums < MARKER_MIN_UMIS] = -np.inf
        for row in score:
            top = np.argsort(-row)[:n_markers]
            markers.append([str(names[j]) for j in top if np.isfinite(row[j])])
    return markers


def _r(value, digits=4):
    """A float for the JSON, NaN and inf as null."""
    if value is None:
        return None
    value = float(value)
    return round(value, digits) if math.isfinite(value) else None


def summarize_groups(adata, info, qc, mode, umap=None):
    """
    One record per group:

        id, pseudo, n_cells, umis, median_umis, mito_pct (pooled: mito UMIs over all
        UMIs of the group), median_cell_mito_pct, rrna_pct, n_doublets, doublet_pct
        (null unless the calls are annotated), alpha_hat_mean (null without CellSweep),
        x, y (Metacell2's UMAP, null when not computed), markers
    """
    groups = adata.obs[METACELL_COL].astype(str).to_numpy()
    order = group_order(groups)
    ind = group_indicator(groups, order)

    n_cells = np.asarray(ind.sum(axis=1)).ravel()
    umis = ind @ qc["total_umis"].to_numpy()
    mito = ind @ qc["mito_umis"].to_numpy()
    rrna = ind @ qc["rrna_umis"].to_numpy()

    is_doublet = None
    if mode == "annotated":
        is_doublet = (adata.obs[DOUBLET_COL].astype(str).str.lower() == "doublet").to_numpy().astype(float)
        n_doublets = ind @ is_doublet

    alpha = None
    if ALPHA_COL in adata.obs.columns:
        alpha = pd.to_numeric(adata.obs[ALPHA_COL], errors="coerce").to_numpy(dtype=float)

    sums = ind @ as_csr(adata.X)
    real = [i for i, g in enumerate(order) if g not in PSEUDO_GROUPS]
    usable = ~info["is_mito"].to_numpy()
    labels = np.where(info["gene_name"].to_numpy() != "", info["gene_name"].to_numpy(), info.index.to_numpy())
    markers = dict(zip(real, marker_genes(sums[real], labels, usable))) if real else {}

    cell_mito = qc["pct_mito"].to_numpy()
    cell_umis = qc["total_umis"].to_numpy()
    members = pd.Series(np.arange(len(groups))).groupby(groups).apply(np.asarray)

    records = []
    for i, g in enumerate(order):
        idx = members[g]
        rec = {
            "id": g,
            "pseudo": g in PSEUDO_GROUPS,
            "n_cells": int(n_cells[i]),
            "umis": _r(umis[i], 2),
            "median_umis": _r(np.median(cell_umis[idx]), 2),
            "mito_pct": _r(100.0 * mito[i] / umis[i]) if umis[i] > 0 else 0.0,
            "median_cell_mito_pct": _r(np.median(cell_mito[idx])),
            "rrna_pct": _r(100.0 * rrna[i] / umis[i]) if umis[i] > 0 else 0.0,
            "n_doublets": int(n_doublets[i]) if is_doublet is not None else None,
            "doublet_pct": _r(100.0 * n_doublets[i] / n_cells[i]) if is_doublet is not None else None,
            "alpha_hat_mean": _r(np.nanmean(alpha[idx])) if alpha is not None and np.isfinite(alpha[idx]).any() else None,
            "x": None,
            "y": None,
            "markers": markers.get(i, []),
        }
        if umap is not None and g in umap.index:
            rec["x"], rec["y"] = _r(umap.loc[g, "x"]), _r(umap.loc[g, "y"])
        records.append(rec)

    return records, order, sums


def summarize_genes(adata, info):
    """
    Per-gene arrays, all aligned to ``ids``:

        ids, names, mito, rrna, pfam (domain lists), total_umi, n_cells,
        ambient_hat (null where CellSweep gave none, or everywhere without it)
    """
    X = as_csr(adata.X)
    total = np.asarray(X.sum(axis=0)).ravel()
    n_cells = np.bincount(X.indices, minlength=X.shape[1])

    ambient = None
    if AMBIENT_VAR_COL in adata.var.columns:
        ambient = [_r(v, 6) for v in pd.to_numeric(adata.var[AMBIENT_VAR_COL], errors="coerce")]

    return {
        "ids": list(info.index),
        "names": list(info["gene_name"]),
        "mito": [int(v) for v in info["is_mito"]],
        "rrna": [int(v) for v in info["is_rrna"]],
        "pfam": list(info["pfam"]),
        "total_umi": [_r(v, 2) for v in total],
        "n_cells": [int(v) for v in n_cells],
        "ambient_hat": ambient,
    }


def build_summary(adata, info, how, sample_id, mapping_method, removed=False, umap=None,
                  status="ok", reason=None, params=None):
    """The report's payload for one sample (see the module docstring)."""
    qc = per_cell_qc(adata.X, info["is_mito"].to_numpy(), info["is_rrna"].to_numpy())
    mode = doublet_mode(adata.obs, removed=removed)
    groups, order, _ = summarize_groups(adata, info, qc, mode, umap=umap)
    fp = fingerprint(adata.obs_names, adata.obs[METACELL_COL], adata.var_names)

    return {
        "schema": SUMMARY_SCHEMA,
        "schema_version": SUMMARY_VERSION,
        "sample": {
            "id": sample_id,
            "mapping_method": mapping_method,
            "status": status,
            "reason": reason,
            "doublet_mode": mode,
            "ambient_available": ALPHA_COL in adata.obs.columns,
            "fingerprint": fp,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_metacells": sum(1 for g in order if g not in PSEUDO_GROUPS),
            "total_umis": _r(qc["total_umis"].sum(), 2),
            "mito_pct": _r(100.0 * qc["mito_umis"].sum() / qc["total_umis"].sum()) if qc["total_umis"].sum() > 0 else 0.0,
            "gene_match": how,
            "params": params or {},
        },
        "metacells": groups,
        "genes": summarize_genes(adata, info),
    }


def skipped_summary(sample_id, mapping_method, reason, n_cells=None, params=None):
    """A sample that produced no metacells, so the report can say why instead of leaving it out."""
    return {
        "schema": SUMMARY_SCHEMA,
        "schema_version": SUMMARY_VERSION,
        "sample": {
            "id": sample_id, "mapping_method": mapping_method, "status": "skipped",
            "reason": reason, "n_cells": n_cells, "params": params or {},
        },
        "metacells": [],
        "genes": None,
    }


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, separators=(",", ":"), allow_nan=False)


# --------------------------------------------------------------------------
# Standalone
# --------------------------------------------------------------------------

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])

    parser = argparse.ArgumentParser(
        description="Summarise a cells h5ad with obs['metacell'] for the metacell filtering report."
    )
    parser.add_argument("--cells_h5ad", required=True, help="Cells h5ad carrying obs['metacell']")
    parser.add_argument("--gene_table", required=True, help="build_gene_table.py output")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--mapping_method", default="")
    parser.add_argument("--doublets_removed", action="store_true",
                        help="The matrix was doublet-filtered (params.perform_doublet_filtering)")
    parser.add_argument("--output", required=True, help="Output summary JSON")
    args = parser.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.cells_h5ad)
    if METACELL_COL not in adata.obs.columns:
        raise SystemExit(f"{args.cells_h5ad} has no obs['{METACELL_COL}'] to summarise")

    umap = None
    if "mc_umap" in adata.uns:
        umap = pd.DataFrame(adata.uns["mc_umap"]).set_index("id")

    info, how = gene_info(adata.var_names, load_gene_table(args.gene_table))
    summary = build_summary(adata, info, how, args.sample_id, args.mapping_method,
                            removed=args.doublets_removed, umap=umap)
    write_json(args.output, summary)
    logger.info(f"Wrote {args.output}: {summary['sample']['n_metacells']} metacells")


if __name__ == "__main__":
    main()
