#!/usr/bin/env python3
"""
Pool the cell-called objects of a samplesheet 'filtering_report' group into one h5ad,
the input of the group's Metacell2 run.

Every member has been called, denoised and doublet-annotated on its own, since the
ambient profile and the doublets belong to one capture. Pooling only follows that:

    obs_names       '<sample>:<barcode>', so cells of different captures that share a
                    barcode stay apart. A barcode never holds ':' (STARsolo joins the
                    segments of a complex barcode with '_'), so rpartition(':') undoes it.
    obs['sample']   the member's samplesheet name; obs['sample_id'] keeps its analytical
                    run id, as MTX_TO_H5AD wrote it
    CellSweep       lifted per member from that member's own raw (or full) h5ad, as
                    run_metacells.py would (its lift_ambient): alpha_hat, is_empty and
                    celltype in obs, the denoised counts in layers['cellsweep']
    ambient_hat     per gene, the members' mean weighted by their cells, with each
                    member's own value kept as var['ambient_hat:<sample>']

The members must be counted against the same genes; a group mixing annotations is
refused rather than reconciled.

A pair of members sharing most of their cell barcodes is warned about: that is one
library sequenced twice, whose cells would enter Metacell2 twice. Its runs belong under
one sample id, which merges the reads before mapping.
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd

import metacell_qc_summary as qcs
from run_metacells import CS_LAYER, lift_ambient, package_version

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("MergeGroup")

MEMBER_COL = "sample"
SEPARATOR = ":"
# Two members sharing more than this share of the smaller one's barcodes are one library.
# Unrelated captures overlap far less, even on a split-pool protocol restricted to a few
# RT wells (a few percent).
SHARED_BARCODES_WARN = 0.5


def plain_strings(adata):
    """
    Names and text columns as plain object strings. An h5ad written under pandas >= 3
    reads back with pandas' nullable 'string' dtype, and concat can produce it too; the
    anndata of the metacells environment refuses to write it, and metacells rejects it
    (as in run_metacells.py). Same values, so nothing downstream changes.
    """
    adata.obs_names = pd.Index(np.asarray(adata.obs_names, dtype=object))
    adata.var_names = pd.Index(np.asarray(adata.var_names, dtype=object))
    for frame in (adata.obs, adata.var):
        for col in frame.columns:
            if isinstance(frame[col].dtype, pd.StringDtype):
                frame[col] = frame[col].astype(object)
    return adata


def read_member(member_id, sample, cells_h5ad, ambient_h5ad):
    import anndata as ad

    adata = plain_strings(ad.read_h5ad(cells_h5ad))
    adata.var_names_make_unique()
    adata.X = qcs.as_csr(adata.X)
    logger.info(f"{sample} ({member_id}): {adata.n_obs} cells x {adata.n_vars} genes from {cells_h5ad}")
    if ambient_h5ad:
        lift_ambient(adata, ambient_h5ad)
    return adata


def check_genes(objs, samples):
    """Every member on the first one's gene axis, in its order; refuse different genes."""
    ref = objs[0].var_names
    for i, (adata, sample) in enumerate(zip(objs, samples)):
        if adata.var_names.equals(ref):
            continue
        if set(adata.var_names) != set(ref):
            only_here = adata.var_names.difference(ref)
            only_ref = ref.difference(adata.var_names)
            raise SystemExit(
                f"{sample} is counted against other genes than {samples[0]} "
                f"({len(only_here)} genes only in {sample}, e.g. {list(only_here[:3])}; "
                f"{len(only_ref)} only in {samples[0]}, e.g. {list(only_ref[:3])}). "
                "The samples of a filtering_report group must share one reference annotation."
            )
        objs[i] = adata[:, ref].copy()
    return objs


def shared_barcodes(objs, samples):
    """``[(a, b, shared, fraction of the smaller)]`` for every pair of members."""
    out = []
    for i in range(len(objs)):
        for j in range(i + 1, len(objs)):
            a, b = set(objs[i].obs_names), set(objs[j].obs_names)
            shared = len(a & b)
            smaller = max(1, min(len(a), len(b)))
            out.append((samples[i], samples[j], shared, shared / smaller))
    return out


def clean_obs(obs):
    """
    The columns one member lacks come back as object arrays with gaps, which h5ad cannot
    write: flags become nullable booleans, numbers floats, the rest categorical strings.
    """
    for col in obs.columns:
        s = obs[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            s = s.astype(object)
        if s.dtype != object:
            continue
        present = s.dropna()
        if len(present) and present.map(lambda v: isinstance(v, (bool, np.bool_))).all():
            obs[col] = s.astype("boolean")
            continue
        numeric = pd.to_numeric(present, errors="coerce")
        if len(present) and numeric.notna().all() and not present.map(lambda v: isinstance(v, str)).any():
            obs[col] = pd.to_numeric(s, errors="coerce").astype(float)
            continue
        obs[col] = pd.Categorical(s.map(lambda v: "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)))
    return obs


def pool_ambient(objs, samples):
    """
    ``(pooled, per_member)``: the cell-weighted mean of the members' ambient_hat and
    each member's own column, or ``(None, {})`` when no member carries one.
    """
    col = qcs.AMBIENT_VAR_COL
    per_member = {}
    for adata, sample in zip(objs, samples):
        if col in adata.var.columns:
            per_member[sample] = pd.to_numeric(adata.var[col], errors="coerce").to_numpy(dtype=float)
    if not per_member:
        return None, {}

    weights = {s: o.n_obs for o, s in zip(objs, samples) if s in per_member}
    values = np.vstack([per_member[s] for s in weights])
    w = np.array([weights[s] for s in weights], dtype=float)[:, None]
    present = np.isfinite(values)
    total_w = (w * present).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = np.where(total_w > 0, (np.where(present, values, 0.0) * w).sum(axis=0) / total_w, np.nan)
    if len(per_member) < len(objs):
        logger.warning(f"Only {len(per_member)} of {len(objs)} samples carry CellSweep's ambient_hat; "
                       "the pooled value averages those")
    return pooled, per_member


def merge(objs, member_ids, samples):
    import anndata as ad

    objs = check_genes(objs, samples)

    for a, b, shared, frac in shared_barcodes(objs, samples):
        if frac > SHARED_BARCODES_WARN:
            logger.warning(
                f"{a} and {b} share {shared} cell barcodes ({frac:.0%} of the smaller): they look like "
                "one library sequenced twice, whose cells would enter Metacell2 twice. Give their runs "
                "one sample id instead, which merges the reads before mapping."
            )
        else:
            logger.info(f"{a} and {b} share {shared} cell barcodes ({frac:.1%} of the smaller)")

    with_layer = [CS_LAYER in o.layers for o in objs]
    if any(with_layer) and not all(with_layer):
        logger.warning(
            f"CellSweep's denoised counts are missing for "
            f"{[s for s, h in zip(samples, with_layer) if not h]}; the pooled object leaves the layer out"
        )
        for o in objs:
            if CS_LAYER in o.layers:
                del o.layers[CS_LAYER]

    pooled_ambient, per_member_ambient = pool_ambient(objs, samples)

    for adata, member_id, sample in zip(objs, member_ids, samples):
        adata.obs_names = [f"{sample}{SEPARATOR}{b}" for b in adata.obs_names]
        adata.obs[MEMBER_COL] = sample
        adata.obs["sample_id"] = member_id
        if qcs.AMBIENT_VAR_COL in adata.var.columns:
            del adata.var[qcs.AMBIENT_VAR_COL]

    pooled = plain_strings(ad.concat(objs, axis=0, join="outer", merge="same", uns_merge="same"))
    pooled.X = qcs.as_csr(pooled.X)
    pooled.obs = clean_obs(pooled.obs)
    pooled.obs[MEMBER_COL] = pd.Categorical(pooled.obs[MEMBER_COL].astype(str), categories=list(samples))

    if pooled_ambient is not None:
        pooled.var[qcs.AMBIENT_VAR_COL] = pooled_ambient
        for sample, values in per_member_ambient.items():
            pooled.var[f"{qcs.AMBIENT_VAR_COL}:{sample}"] = values

    pooled.uns["group_members"] = {
        "sample": np.asarray(samples, dtype=str),
        "id": np.asarray(member_ids, dtype=str),
        "n_cells": np.asarray([o.n_obs for o in objs], dtype=np.int64),
    }
    return pooled


def write_versions(path, process):
    import platform

    with open(path, "w") as f:
        f.write(f'"{process}":\n')
        f.write(f"    python: {platform.python_version()}\n")
        f.write(f"    anndata: {package_version('anndata')}\n")


def main():
    parser = argparse.ArgumentParser(description="Pool the cell-called objects of a filtering_report group.")
    parser.add_argument("--group_id", required=True)
    parser.add_argument("--member", action="append", nargs="+", required=True,
                        metavar=("ID SAMPLE CELLS_H5AD", "AMBIENT_H5AD"),
                        help="One member: its analytical id, its samplesheet name, its cell-called h5ad "
                             "and optionally its raw/full h5ad carrying CellSweep's annotations")
    parser.add_argument("--output", required=True)
    parser.add_argument("--members_tsv", default=None)
    parser.add_argument("--compression", default="gzip")
    parser.add_argument("--versions_yml", default=None)
    parser.add_argument("--process_name", default="MERGE_GROUP_H5AD")
    args = parser.parse_args()

    for m in args.member:
        if len(m) not in (3, 4):
            parser.error(f"--member takes ID SAMPLE CELLS_H5AD [AMBIENT_H5AD], got {m}")
    member_ids = [m[0] for m in args.member]
    samples = [m[1] for m in args.member]
    if len(set(samples)) != len(samples):
        parser.error(f"--member names a sample twice: {samples}")
    bad = [s for s in samples if SEPARATOR in s]
    if bad:
        parser.error(f"sample names cannot hold '{SEPARATOR}', which separates them from the barcode: {bad}")

    objs = [read_member(m[0], m[1], m[2], m[3] if len(m) == 4 and m[3] else None) for m in args.member]
    pooled = merge(objs, member_ids, samples)

    compression = None if args.compression.lower() in ("", "none", "null") else args.compression
    pooled.write_h5ad(args.output, compression=compression)
    logger.info(f"{args.group_id}: {pooled.n_obs} cells from {len(samples)} samples x {pooled.n_vars} genes")

    if args.members_tsv:
        counts = pooled.obs[MEMBER_COL].value_counts()
        with open(args.members_tsv, "w") as fh:
            fh.write("id\tsample\tn_cells\n")
            for member_id, sample in zip(member_ids, samples):
                fh.write(f"{member_id}\t{sample}\t{int(counts.get(sample, 0))}\n")

    if args.versions_yml:
        write_versions(args.versions_yml, args.process_name)


if __name__ == "__main__":
    main()
